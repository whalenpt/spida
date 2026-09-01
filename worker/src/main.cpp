// spida-worker
//
// Runs exactly one SPIDA simulation job from a JSON config and reports its
// outcome the way the job-service backend expects to find it (see
// docs/api/openapi.yaml, docs/api/events.schema.json): a status.json written
// before and after the run, a manifest.json summarizing whatever report
// files SPIDA itself wrote into the output directory -- rewritten LIVE on
// every new frame, not just once at the end, so GET /results can answer
// mid-run instead of only after the process exits -- and a SimulationEvent
// stream (one JSON object per line, EventLog below, now including a
// "report" kind fired alongside each manifest rewrite) delivered through
// TWO independent EventSinks as the run progresses:
//   - FileEventSink -> events.ndjson, appended -- durable; replayable by a
//     caller after this process exits, or after an api-server restart.
//   - StdoutEventSink -> this process's stdout, flushed per line -- the
//     LIVE channel. A caller that spawns spida-worker already holds a pipe
//     to its stdout the instant it does, so this needs no polling and no
//     filesystem watching to get real-time delivery: read a line from the
//     child's stdout, forward it (e.g. as one SSE frame per line to a
//     browser), done. See EventSink's own header comment below for why
//     these are two sinks fed by one producer, not two separate code paths.
//
// Usage: spida-worker <config.json> <output-dir> [timeout-seconds]
//        spida-worker --describe   (prints capabilities.h's JSON, exits 0)
//
// timeout-seconds is an optional OPERATOR wall-clock cap (proposal's error
// taxonomy: "timeout" — "optional operator wall-clock cap exceeded"), not
// part of SimulationConfig itself — a deployment concern, not something the
// submitted config carries. Enforced the same way SIGTERM-driven
// cancellation already is: reuses requestCancel(), checked at the next
// stepUpdate() checkpoint, not an immediate kill. Omitted or <= 0 disables
// it entirely (no wall-clock cap, the only behavior before this).
//
// This is the in-tree successor to a worker originally authored and
// numerically verified in spida-console's services/worker — see
// docs/adr/0002-worker-model-coverage-and-config-registry.md for why it
// moved here and exactly what changed. In short: construction now goes
// through spida::config::SimulationRun instead of hand-duplicating each
// model's setup; progress comes from the real
// BasePropagator::setProgressObserver() (exact stepsTaken/currentStepSize)
// instead of an estimate derived from report cadence; stopReason comes from
// the real spida::StopReason instead of a manifest-scanning heuristic; and
// SIGTERM now triggers real cooperative cancellation via requestCancel()
// instead of always requiring a hard SIGKILL from the caller. See
// docs/adr/0003-transport-and-live-events.md for a real consequence this
// has on api-server's own exit-handling logic that still needs to be
// picked up there (out of scope for this repo).
//
// config.json is now checked with spida::config::validate() before
// SimulationRun is ever constructed, producing a structured
// "validationErrors" array (field + message per entry) in status.json
// instead of a single free-text exception message — see
// docs/adr/0002-worker-model-coverage-and-config-registry.md.
//
// Two more fields since: status.json's "config" now echoes the fully-
// resolved SimulationConfig (every field defaulted, not just what the
// caller's config.json actually set) on every write from "running" onward
// — a client no longer has to re-derive spida's own default table (see
// modelregistry.h) to know what a run actually used. manifest.json's
// per-series entries now also carry "gridCoords" (and "gridCoordsY" for a
// 2D series) — the x/y axis values are static across a series' frames and
// captured once from the report sink instead of requiring a caller to
// open and re-parse a report file to get them.

#include <spida/config/capabilities.h>
#include <spida/config/modelregistry.h>
#include <spida/config/simulationbuilder.h>
#include <spida/config/simulationconfig.h>
#include <spida/config/validation.h>

#include <nlohmann/json.hpp>

#include <chrono>
#include <csignal>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

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

// Rewrites manifest.json in full each time -- called from the report sink
// on every new frame (see main()'s setReportSink() call), not just once
// after the run finishes. Without this, GET /simulations/:id/results (the
// endpoint that tells a frontend what series exist and their axes) can't
// return anything until the ENTIRE run is over, which defeats the point of
// the live status/progress channel: a caller could see "73% done" but not
// what the field actually looks like at 73%, only at 100%. Cheap to do
// unconditionally -- a few KB, called at most a few hundred times over a
// run, trivial next to actual solve cost -- and each frame's own file is
// already on disk by the time this runs (report1D()/report2D() write the
// file and call the sink synchronously, in that order, at the same
// checkpoint), so nothing this writes ever points at a file that isn't
// there yet.
void writeManifest(const fs::path& dir, std::string_view simulationId, const json& series)
{
    std::ofstream mf(dir / "manifest.json");
    mf << json{{"simulationId", std::string(simulationId)}, {"series", series}}.dump(2);
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

// ---- EventSink: where a built SimulationEvent goes -----------------------
// An implementation's only job is delivery -- it never builds event content
// (id/simulationId/at/kind) or knows what "progress" vs "log" means. That's
// EventLog's job (below), which fans one already-built event out to every
// registered sink. Splitting these apart (rather than one class that both
// assembles and appends-to-a-file, as before) is what makes adding a second
// destination -- stdout, for a future api-server to read as a live pipe,
// alongside the existing durable file -- purely additive: a new EventSink
// implementation and one more entry in main()'s sink list, no change to
// EventLog's status()/log()/progress() or their callers at all.
class EventSink {
public:
    virtual ~EventSink() = default;
    virtual void write(const json& event) = 0;
};

// events.ndjson, appended (not rewritten) — the durable copy: still there
// for a caller to replay after the process exits (GET .../events?since=)
// or to recover from after an api-server restart, independent of whether
// anything was reading the live stdout feed at the time.
class FileEventSink : public EventSink {
public:
    explicit FileEventSink(fs::path path) : m_path(std::move(path))
    {
    }

    void write(const json& event) override
    {
        std::ofstream os(m_path, std::ios::app);
        os << event.dump() << "\n";
    }

private:
    fs::path m_path;
};

// One NDJSON line per event on stdout, explicitly flushed. This is the
// live channel: a parent process that spawns spida-worker already holds a
// pipe to its stdout the instant it does, with no polling, no inotify, and
// none of the partial-write races that come with tailing a file mid-append
// -- see this file's own header comment for the fuller argument. The
// explicit flush() matters: stdout is fully-buffered (not line-buffered)
// the moment it's redirected to a pipe rather than a tty, so without it
// events would sit in libstdc++'s internal buffer and only reach the
// parent in a burst at process exit -- silently defeating the entire
// point of using stdout as a live channel instead of a batch one.
class StdoutEventSink : public EventSink {
public:
    void write(const json& event) override
    {
        std::cout << event.dump() << "\n";
        std::cout.flush();
    }
};

// ---- EventLog: builds a SimulationEvent, fans it out to every sink -------
// Public API (status()/log()/progress()) unchanged from before this
// refactor — only how an assembled event is delivered changed. See
// docs/api/events.schema.json for the shape.
class EventLog {
public:
    EventLog(std::string simulationId, std::vector<std::unique_ptr<EventSink>> sinks)
        : m_simId(std::move(simulationId)), m_sinks(std::move(sinks))
    {
    }

    void status(const std::string& status_)
    {
        emit({{"kind", "status"}, {"status", status_}});
    }

    void log(const std::string& level, const std::string& message)
    {
        emit({{"kind", "log"}, {"level", level}, {"message", message}});
    }

    void progress(const spida::ProgressSnapshot& s)
    {
        json p = {{"t", s.t}, {"stepsTaken", s.stepsTaken}};
        if (s.tf.has_value())
            p["tf"] = *s.tf;
        if (s.currentStepSize.has_value())
            p["currentStepSize"] = *s.currentStepSize;
        emit({{"kind", "progress"}, {"progress", p}});
    }

    // A new frame of `seriesName` is available -- both as its own file
    // (<seriesName>_<frameIndex>.json) and reflected into the just-rewritten
    // manifest.json (see writeManifest()) -- by the time this event is
    // observed. Lets a live subscriber fetch/render the new frame the
    // instant it exists instead of polling GET /results to notice
    // frameCount grew, the same reasoning that motivated pushing
    // status/progress live rather than leaving them poll-only.
    void report(const std::string& seriesName, std::size_t frameIndex)
    {
        emit({{"kind", "report"}, {"seriesName", seriesName}, {"frameIndex", frameIndex}});
    }

private:
    void emit(json evt)
    {
        evt["id"] = std::to_string(m_seq++);
        evt["simulationId"] = m_simId;
        evt["at"] = nowIso8601();
        for (auto const& sink : m_sinks)
            sink->write(evt);
    }

    std::string m_simId;
    std::vector<std::unique_ptr<EventSink>> m_sinks;
    std::size_t m_seq{0};
};

// ---- Manifest: summarized from ReportHandler's sink, not by re-reading ---
// the files SPIDA itself just wrote. The original worker re-parsed each
// "<name>_<i>.json" off disk and regex-matched filenames to recover what
// this JSON already told us the moment it was produced — exactly the
// file-tailing/re-parsing pattern
// docs/adr/0003-transport-and-live-events.md rejects for the live events
// channel too ("added latency, races on partial writes, and re-parsing
// JSON the same process just serialized"). Same argument applies here:
// setReportSink() (docs/adr/0001-library-extension-seams.md) delivers each
// report's name and JSON in-process, at the moment BasePropagator::report1D
// ()/report2D() call it, with no filesystem round-trip.
//
// "type" in that JSON ("xy"/"xy_complex"/"xyz"/"xyz_complex") distinguishes
// real vs complex and 1D vs 2D unambiguously. It does NOT distinguish a
// Track series from a one-shot Report1D — both serialize "xy". None of the
// three models wired today (burgers/kdv_rv/ks) register a Track report (see
// each model's own header), so that ambiguity is moot for now; it becomes
// live the day a model adds one, at which point this needs either a
// separate sink per report kind or a richer Sink signature that says which
// kind called it — flagged here rather than silently assumed away.
class ManifestBuilder {
public:
    // Returns the series' new frameCount (i.e. 1-based count including the
    // just-observed frame) so a caller can derive that frame's 0-based file
    // index (frameCount - 1) without duplicating any counting of its own --
    // matches BasePropagator's own m_report_count1D/m_report_count2D
    // numbering exactly, since both this series and every other 1D (or
    // both every other 2D) series registered on the same propagator share
    // one counter, incremented once per checkpoint (see propagator.cpp).
    std::size_t observe(std::string_view name, const json& j)
    {
        auto& s = m_series[std::string(name)];
        const std::string type = j.value("type", "xy");
        s.valueType = (type == "xy_complex" || type == "xyz_complex") ? "complex" : "real";
        s.kind = (type == "xyz" || type == "xyz_complex") ? "field2d" : "field1d";
        // Grid coordinates ("x", and "y" for a 2D series) are static across
        // every frame of a series -- report.hpp's buildJson() re-embeds them
        // redundantly in each frame's JSON, so capturing them once here (on
        // first observation) is enough. Done specifically so a future
        // api-server's ResultSeriesDescriptor.gridCoords/gridCoordsY (see
        // docs/api/openapi.yaml, docs/api/binary-frame-spec.md's "x is not
        // repeated per frame -- served once") doesn't have to re-open and
        // re-parse a report file to reconstruct them -- the same
        // sink-not-file-scan argument this class's own header comment
        // already makes for the manifest itself.
        if (s.frameCount == 0) {
            if (j.contains("x"))
                s.gridCoords = j.at("x");
            if (s.kind == "field2d" && j.contains("y"))
                s.gridCoordsY = j.at("y");
        }
        s.frameCount++;
        return s.frameCount;
    }

    [[nodiscard]] json build() const
    {
        json out = json::array();
        for (auto const& [name, s] : m_series) {
            json entry = {
                {"name", name},
                {"kind", s.kind},
                {"valueType", s.valueType},
                {"frameCount", s.frameCount},
                {"gridCoords", s.gridCoords},
            };
            if (s.kind == "field2d")
                entry["gridCoordsY"] = s.gridCoordsY;
            out.push_back(std::move(entry));
        }
        return out;
    }

private:
    struct Series {
        std::string kind;
        std::string valueType;
        std::size_t frameCount = 0;
        json gridCoords = json::array();
        json gridCoordsY = json::array();
    };
    std::map<std::string, Series> m_series;
};

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
    if (argc == 2 && std::string_view(argv[1]) == "--describe") {
        std::cout << spida::config::describeCapabilities().dump(2) << "\n";
        return 0;
    }
    if (argc != 3 && argc != 4) {
        std::cerr << "usage: spida-worker <config.json> <output-dir> [timeout-seconds]\n"
                     "       spida-worker --describe\n";
        return 2;
    }
    const fs::path configPath = argv[1];
    const fs::path outDir = argv[2];
    double timeoutSeconds = 0.0; // <= 0 means disabled
    if (argc == 4) {
        try {
            timeoutSeconds = std::stod(argv[3]);
        }
        catch (const std::exception&) {
            std::cerr << "invalid timeout-seconds: " << argv[3] << "\n";
            return 2;
        }
    }
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

    const std::string simId = outDir.filename().string();

    // Two sinks: events.ndjson (durable, replayable after this process
    // exits) and stdout (live -- the pipe an api-server already holds the
    // instant it spawns this process). See EventSink's own header comment
    // above for why this is a fan-out, not a choice between the two.
    std::vector<std::unique_ptr<EventSink>> sinks;
    sinks.push_back(std::make_unique<FileEventSink>(outDir / "events.ndjson"));
    sinks.push_back(std::make_unique<StdoutEventSink>());
    EventLog events(simId, std::move(sinks));

    try {
        auto cfg = configJson.get<spida::config::SimulationConfig>();

        // Backfill modelParams with this model's own registry defaults for
        // any key the caller's config.json omitted, BEFORE anything reads
        // cfg.modelParams (including SimulationRun's construction below) or
        // it's echoed into status.json's "config" field. Without this,
        // "config" would show modelParams as {} whenever a caller relied on
        // a default -- the struct-level fields (grid/solver/reporting)
        // really are resolved by NLOHMANN_DEFINE_TYPE_..._WITH_DEFAULT the
        // moment configJson.get<SimulationConfig>() ran above, but
        // modelParams is an opaque JSON blob whose defaults were previously
        // only ever applied ad hoc, inside simulationbuilder.h's
        // cfg.modelParams.value(name, ...) calls -- invisible to anything
        // that only reads cfg.modelParams itself. Caught by actually running
        // this end to end (a burgers config with no modelParams echoed back
        // "modelParams": {}), not by the unit tests, which never inspected
        // status.json's "config" field this closely.
        if (const auto* desc = spida::config::describe(cfg.model); desc != nullptr) {
            for (auto const& p : desc->modelParams) {
                if (!cfg.modelParams.contains(p.name))
                    cfg.modelParams[p.name] = p.defaultValue;
            }
        }

        // Checked before anything is constructed (or "running" is even
        // written) so a bad config never claims to have started. Structured
        // per-field errors, not just one free-text exception message — see
        // this file's header comment and validation.h's own header comment
        // for why (the proposal's own UX goal: inline field errors on the
        // config form, no run created).
        if (auto errors = spida::config::validate(cfg); !errors.empty()) {
            const std::string detail =
                "config validation failed (" + std::to_string(errors.size()) + " error(s))";
            writeStatus(outDir,
                        {{"status", "failed"},
                         {"failureReason", "config_validation"},
                         {"failureDetail", detail},
                         {"validationErrors", errors},
                         {"config", cfg},
                         {"finishedAt", nowIso8601()}});
            events.log("error", detail);
            events.status("failed");
            return 1;
        }

        // "config" below is the fully-resolved SimulationConfig (every
        // omitted field already filled by NLOHMANN_DEFINE_TYPE_..._WITH_
        // DEFAULT the moment configJson.get<SimulationConfig>() ran above),
        // not just an echo of the caller's raw config.json -- written into
        // every status.json from here on so a client can see exactly what
        // ran (e.g. the real "mu" used when the caller never set one)
        // without re-deriving spida's own default table itself.
        writeStatus(outDir,
                    {{"status", "running"}, {"config", cfg}, {"startedAt", nowIso8601()}});
        events.status("running");

        spida::config::SimulationRun run(cfg, outDir);
        g_propagator = &run.propagator();
        // Known, accepted gap: SIGTERM sent to this process BEFORE this
        // point (construction above, or the config-parse/validate steps
        // before it) gets the default disposition -- immediate termination,
        // no status.json written -- rather than cooperative cancellation.
        // Narrow in practice (construction is milliseconds at these problem
        // sizes) but real; not closed here because doing so needs a
        // pre-construction cancel flag checked at a couple of points, which
        // is meaningfully more complexity than the size of the actual risk
        // justifies right now. Documented in worker/README.md's
        // "Cancellation" section too -- flagged rather than silently
        // assumed away.
        std::signal(SIGTERM, handleSigterm);

        // Progress events at report cadence (stepsPerOutput1D), matching the
        // original worker's rate — setProgressObserver() itself fires every
        // accepted solver step, which is much finer-grained than report
        // cadence and would otherwise multiply events.ndjson's write volume
        // several-fold for no benefit api-server currently uses.
        //
        // The same callback also enforces the wall-clock timeout, checked on
        // EVERY invocation (not throttled by stepsPerOutput1D like the event
        // forwarding below it) — a coarse report cadence shouldn't delay
        // noticing a timeout. requestCancel() takes effect at the run's next
        // stepUpdate() checkpoint, same as SIGTERM — see handleSigterm above.
        const auto stepsPerOutput1D = cfg.reporting.stepsPerOutput1D;
        const auto startTime = std::chrono::steady_clock::now();
        bool timedOut = false;
        run.propagator().setProgressObserver(
            [&events, stepsPerOutput1D, timeoutSeconds, startTime, &timedOut](
                const spida::ProgressSnapshot& s) {
                if (timeoutSeconds > 0.0 && !timedOut) {
                    const double elapsed =
                        std::chrono::duration<double>(std::chrono::steady_clock::now() - startTime)
                            .count();
                    if (elapsed >= timeoutSeconds) {
                        timedOut = true;
                        g_propagator->requestCancel();
                    }
                }
                if (stepsPerOutput1D == 0 || s.stepsTaken % stepsPerOutput1D == 0)
                    events.progress(s);
            });

        // Report files are still written to outDir as before (setWriteReportFiles
        // defaults to true) — the sink only adds an in-process observation of
        // the same data, replacing the manifest's old post-run file scan.
        // manifest.json is rewritten on every frame, not just once at the
        // end (see writeManifest()'s own comment for why), and a "report"
        // event fires alongside it -- both derived from the same observe()
        // call, so a live subscriber and a GET /results poller can never
        // see mutually inconsistent states.
        ManifestBuilder manifest;
        run.propagator().setReportSink([&manifest, &outDir, &simId, &events](std::string_view name,
                                                                             const json& j) {
            const std::size_t frameCount = manifest.observe(name, j);
            writeManifest(outDir, simId, manifest.build());
            events.report(std::string(name), frameCount - 1); // 0-based, matches "<name>_<i>.json"
        });

        const bool stepOk = run.run();

        // Final write: a safety net for a hypothetical model that registers
        // zero report series (none of the six wired ones do -- every model
        // has at least one), so manifest.json always exists once a run
        // reaches this point either way. For every real model today this is
        // a no-op re-write of exactly what the sink above already wrote on
        // the last frame.
        writeManifest(outDir, simId, manifest.build());

        if (!stepOk) {
            const std::string detail = "solver.evolve() returned false (a solver step failed)";
            writeStatus(outDir,
                        {{"status", "failed"},
                         {"failureReason", "runtime_exception"},
                         {"failureDetail", detail},
                         {"config", cfg},
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
                         {"config", cfg},
                         {"finishedAt", nowIso8601()}});
            events.log("error", detail);
            events.status("failed");
            return 1;
        }

        // Checked before the generic "completed" path below: timedOut is
        // only ever set by this worker's own wall-clock check above, so
        // it's unambiguous — a run that also happened to receive SIGTERM
        // around the same time still reports "timeout" here, which is a
        // defensible call either way given both were true. reason is
        // reported alongside for the underlying mechanism (almost always
        // CancelRequested, since that's what requestCancel() produces).
        if (timedOut) {
            const std::string detail =
                "exceeded operator wall-clock time budget of " + std::to_string(timeoutSeconds) +
                "s";
            writeStatus(outDir,
                        {{"status", "failed"},
                         {"failureReason", "timeout"},
                         {"stopReason", stopReasonToString(reason)},
                         {"failureDetail", detail},
                         {"config", cfg},
                         {"finishedAt", nowIso8601()}});
            events.log("error", detail);
            events.status("failed");
            return 1;
        }

        // MaxReportsReached / CancelRequested / None (tf reached) all land
        // on "completed" from the worker's own point of view — a cancelled
        // run now finishes through this same path instead of being killed
        // mid-flight. See docs/adr/0003-transport-and-live-events.md for
        // what that means for api-server's exit handler, which still
        // assumes a cancelling job never gets to write its own terminal
        // status.json.
        writeStatus(outDir,
                    {{"status", "completed"},
                     {"stopReason", stopReasonToString(reason)},
                     {"config", cfg},
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
