// spida-worker
//
// Runs exactly one SPIDA simulation job from a JSON config and reports its
// outcome the way the job-service backend expects to find it (see
// docs/api/openapi.yaml, docs/api/events.schema.json): a status.json written
// before and after the run, a manifest.json summarizing whatever report
// files SPIDA itself wrote into the output directory, and an events.ndjson
// (one SimulationEvent per line) appended as the run progresses.
//
// Usage: spida-worker <config.json> <output-dir>
//
// This is the in-tree successor to a worker originally authored and
// numerically verified in spida-console's services/worker — see
// docs/adr/0003-worker-relocation-and-cooperative-cancellation.md for why it
// moved here and exactly what changed. In short: construction now goes
// through spida::config::SimulationRun instead of hand-duplicating each
// model's setup; progress comes from the real
// BasePropagator::setProgressObserver() (exact stepsTaken/currentStepSize)
// instead of an estimate derived from report cadence; stopReason comes from
// the real spida::StopReason instead of a manifest-scanning heuristic; and
// SIGTERM now triggers real cooperative cancellation via requestCancel()
// instead of always requiring a hard SIGKILL from the caller. See
// ADR-0003 for a real consequence this has on api-server's own exit-handling
// logic that still needs to be picked up there (out of scope for this repo).

#include <spida/config/simulationbuilder.h>
#include <spida/config/simulationconfig.h>

#include <nlohmann/json.hpp>

#include <chrono>
#include <csignal>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <regex>
#include <string>

namespace fs = std::filesystem;
using json = nlohmann::json;

namespace {

std::string nowIso8601()
{
    auto now = std::chrono::system_clock::now();
    auto t = std::chrono::system_clock::to_time_t(now);
    std::tm tm{};
    gmtime_r(&t, &tm);
    char buf[32];
    std::strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%SZ", &tm);
    return buf;
}

void writeStatus(const fs::path& dir, const json& status)
{
    std::ofstream os(dir / "status.json");
    os << status.dump(2);
}

/// String values chosen to stay compatible with the original worker's
/// "tf_reached"/"max_reports_reached" — api-server already passes
/// status.json's "stopReason" through opaquely (see jobs.ts), so the two new
/// values below (real cancellation/divergence, neither reachable before this
/// rewrite) don't need any change on that side to show up correctly.
std::string stopReasonToString(spida::StopReason reason)
{
    switch (reason) {
    case spida::StopReason::None:
        return "tf_reached";
    case spida::StopReason::MaxReportsReached:
        return "max_reports_reached";
    case spida::StopReason::CancelRequested:
        return "cancel_requested";
    case spida::StopReason::Diverged:
        return "diverged";
    }
    return "tf_reached";
}

// ---- events.ndjson: SimulationEvent, one JSON object per line ------------
// Appended (not rewritten) so api-server can tail it with a byte offset
// while the worker is still running (see events.ts's poll()). Unchanged
// from the original worker — see docs/api/events.schema.json for the shape.
class EventLog {
public:
    EventLog(fs::path dir, std::string simulationId)
        : m_path(std::move(dir) / "events.ndjson"), m_simId(std::move(simulationId))
    {
    }

    void status(const std::string& status_)
    {
        append({{"kind", "status"}, {"status", status_}});
    }

    void log(const std::string& level, const std::string& message)
    {
        append({{"kind", "log"}, {"level", level}, {"message", message}});
    }

    void progress(const spida::ProgressSnapshot& s)
    {
        json p = {{"t", s.t}, {"stepsTaken", s.stepsTaken}};
        if (s.tf.has_value())
            p["tf"] = *s.tf;
        if (s.currentStepSize.has_value())
            p["currentStepSize"] = *s.currentStepSize;
        append({{"kind", "progress"}, {"progress", p}});
    }

private:
    void append(json evt)
    {
        evt["id"] = std::to_string(m_seq++);
        evt["simulationId"] = m_simId;
        evt["at"] = nowIso8601();
        std::ofstream os(m_path, std::ios::app);
        os << evt.dump() << "\n";
    }

    fs::path m_path;
    std::string m_simId;
    std::size_t m_seq{0};
};

// ---- Manifest: summarize whatever report files SPIDA itself wrote --------
// Unchanged from the original worker. Convention verified against
// src/propagator/propagator.cpp and src/utils/report.hpp:
// BasePropagator::report1D()/report2D() write a numbered "<name>_<i>.json"
// family (a "field1d"/"field2d" series, one frame per report);
// BasePropagator::reportTrack() overwrites a single un-numbered
// "<name>.json" every call (a "track" series) — not distinguishable from a
// one-shot Report1D by JSON shape alone, only by this naming convention.
json buildManifest(const fs::path& dir)
{
    struct Series {
        std::string kind;
        std::string valueType;
        int maxIndex = -1;
        bool single = false;
    };
    std::map<std::string, Series> series;
    static const std::regex numbered(R"(^(.+)_(\d+)\.json$)");

    for (auto const& entry : fs::directory_iterator(dir)) {
        if (!entry.is_regular_file() || entry.path().extension() != ".json")
            continue;
        auto fname = entry.path().filename().string();
        if (fname == "status.json" || fname == "manifest.json" || fname == "config.json")
            continue;

        std::smatch m;
        bool isNumbered = std::regex_match(fname, m, numbered);
        std::string name = isNumbered ? m[1].str() : entry.path().stem().string();

        std::ifstream is(entry.path());
        json j;
        is >> j;
        std::string type = j.value("type", "xy");
        std::string valueType = (type == "xy_complex" || type == "xyz_complex") ? "complex" : "real";
        std::string kind = (type == "xyz" || type == "xyz_complex") ? "field2d" : "field1d";

        auto& s = series[name];
        s.valueType = valueType;
        if (isNumbered) {
            int idx = std::stoi(m[2].str());
            s.maxIndex = std::max(s.maxIndex, idx);
            s.kind = kind;
        }
        else {
            s.single = true;
            s.kind = "track";
        }
    }

    json out = json::array();
    for (auto const& [name, s] : series) {
        out.push_back({
            {"name", name},
            {"kind", s.kind},
            {"valueType", s.valueType},
            {"frameCount", s.single ? 1 : s.maxIndex + 1},
        });
    }
    return out;
}

// ---- SIGTERM -> cooperative cancel ---------------------------------------
// api-server sends SIGTERM to cancel a running job (see jobs.ts's
// hardKill()) — this now takes effect cooperatively via requestCancel()
// (see propagator.h) instead of unconditionally requiring a hard kill.
// requestCancel() only stores to an atomic bool, which is async-signal-safe;
// nothing else happens inside the handler itself.
spida::BasePropagator* g_propagator = nullptr;

extern "C" void handleSigterm(int)
{
    if (g_propagator != nullptr)
        g_propagator->requestCancel();
}

} // namespace

int main(int argc, char** argv)
{
    if (argc != 3) {
        std::cerr << "usage: spida-worker <config.json> <output-dir>\n";
        return 2;
    }
    const fs::path configPath = argv[1];
    const fs::path outDir = argv[2];
    fs::create_directories(outDir);

    json configJson;
    {
        std::ifstream is(configPath);
        if (!is) {
            std::cerr << "cannot open config: " << configPath << "\n";
            return 2;
        }
        is >> configJson;
    }

    EventLog events(outDir, outDir.filename().string());

    try {
        const auto cfg = configJson.get<spida::config::SimulationConfig>();

        writeStatus(outDir, {{"status", "running"}, {"startedAt", nowIso8601()}});
        events.status("running");

        spida::config::SimulationRun run(cfg, outDir);
        g_propagator = &run.propagator();
        std::signal(SIGTERM, handleSigterm);

        // Progress events at report cadence (stepsPerOutput1D), matching the
        // original worker's rate — setProgressObserver() itself fires every
        // accepted solver step, which is much finer-grained than report
        // cadence and would otherwise multiply events.ndjson's write volume
        // several-fold for no benefit api-server currently uses.
        const auto stepsPerOutput1D = cfg.reporting.stepsPerOutput1D;
        run.propagator().setProgressObserver([&events, stepsPerOutput1D](const spida::ProgressSnapshot& s) {
            if (stepsPerOutput1D == 0 || s.stepsTaken % stepsPerOutput1D == 0)
                events.progress(s);
        });

        const bool stepOk = run.run();

        auto manifest = buildManifest(outDir);
        std::ofstream mf(outDir / "manifest.json");
        mf << json{{"simulationId", outDir.filename().string()}, {"series", manifest}}.dump(2);

        if (!stepOk) {
            const std::string detail = "solver.evolve() returned false (a solver step failed)";
            writeStatus(outDir,
                        {{"status", "failed"},
                         {"failureReason", "runtime_exception"},
                         {"failureDetail", detail},
                         {"finishedAt", nowIso8601()}});
            events.log("error", detail);
            events.status("failed");
            return 1;
        }

        const auto reason = run.propagator().stopReason();
        if (reason == spida::StopReason::Diverged) {
            const std::string detail = "solution diverged";
            writeStatus(outDir,
                        {{"status", "failed"},
                         {"failureReason", "divergence"},
                         {"stopReason", stopReasonToString(reason)},
                         {"failureDetail", detail},
                         {"finishedAt", nowIso8601()}});
            events.log("error", detail);
            events.status("failed");
            return 1;
        }

        // MaxReportsReached / CancelRequested / None (tf reached) all land
        // on "completed" from the worker's own point of view — a cancelled
        // run now finishes through this same path instead of being killed
        // mid-flight. See ADR-0003 for what that means for api-server's
        // exit handler, which still assumes a cancelling job never gets to
        // write its own terminal status.json.
        writeStatus(outDir,
                    {{"status", "completed"},
                     {"stopReason", stopReasonToString(reason)},
                     {"finishedAt", nowIso8601()}});
        events.status("completed");
        return 0;
    }
    catch (const nlohmann::json::exception& e) {
        writeStatus(outDir,
                    {{"status", "failed"},
                     {"failureReason", "config_validation"},
                     {"failureDetail", e.what()},
                     {"finishedAt", nowIso8601()}});
        events.log("error", e.what());
        events.status("failed");
        return 1;
    }
    catch (const std::invalid_argument& e) {
        writeStatus(outDir,
                    {{"status", "failed"},
                     {"failureReason", "config_validation"},
                     {"failureDetail", e.what()},
                     {"finishedAt", nowIso8601()}});
        events.log("error", e.what());
        events.status("failed");
        return 1;
    }
    catch (const std::exception& e) {
        writeStatus(outDir,
                    {{"status", "failed"},
                     {"failureReason", "runtime_exception"},
                     {"failureDetail", e.what()},
                     {"finishedAt", nowIso8601()}});
        events.log("error", e.what());
        events.status("failed");
        return 1;
    }
}
