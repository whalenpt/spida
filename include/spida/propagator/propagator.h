#pragma once

#include "spida/helper/constants.h"

#include <atomic>
#include <cstddef>
#include <filesystem>
#include <functional>
#include <memory>
#include <optional>
#include <string_view>
#include <vector>

#include <nlohmann/json.hpp>

namespace detail {
class StatCenter;
} // namespace detail

namespace spida {

class ReportHandler;
class ReportData1D;
class ReportData2D;
class TrackData;

/// Why a BasePropagator-driven evolve() loop stopped early, from the
/// propagator's point of view. `None` covers both "still running" and
/// "reached tf normally" — the solver's own evolve() return value is what
/// distinguishes a clean finish from a step failure; this enum only
/// classifies the propagator-level checkpoint that stepUpdate() reports
/// through (see stepUpdate()'s bool return, which every SolverCV::evolve()
/// loop already checks each step).
enum class StopReason { None, MaxReportsReached, CancelRequested, Diverged };

/// Snapshot of run progress, reported at the same checkpoint as console
/// progress logging (see setLogProgress()/setLogFrequency()). currentStepSize
/// is unset unless the solver reports it via the stepUpdate(t, dt) overload.
struct ProgressSnapshot {
    double t{0.0};
    std::optional<double> tf;
    std::size_t stepsTaken{0};
    std::optional<double> currentStepSize;
};

class BasePropagator {
public:
    explicit BasePropagator(const std::filesystem::path& dir_path);
    virtual ~BasePropagator();
    virtual void updateFields(double t) = 0;

    /// Optional hook for a subclass to detect a diverging solution (e.g. a
    /// non-finite value in its own field data) right after updateFields(t)
    /// runs. Default: never diverged. Checked from stepUpdate(), which is
    /// the same checkpoint SolverCV_AS/SolverCV_CS::evolve() already call
    /// once per accepted step — no solver-loop changes needed.
    virtual bool checkDiverged(double t)
    {
        (void)t;
        return false;
    }

    [[nodiscard]] bool hasData1D() const;

    void setDirPath(const std::filesystem::path& dir_path)
    {
        this->m_dir_path = dir_path;
    }

    void setLogProgress(bool val)
    {
        this->m_log_progress = val;
    }

    void setLogFrequency(std::size_t val);

    [[nodiscard]] bool logProgress() const
    {
        return this->m_log_progress;
    }

    /// Final time of the current run, used only to fill ProgressSnapshot::tf.
    /// Purely informational — does not affect stepUpdate()/evolve() behavior.
    void setFinalTime(double tf)
    {
        this->m_tf = tf;
    }

    /// Called from another thread (or a signal handler) to cooperatively stop
    /// the run. Takes effect at the next stepUpdate() checkpoint, not
    /// immediately — evolve()'s in-flight step still completes.
    void requestCancel()
    {
        this->m_cancel_requested.store(true, std::memory_order_relaxed);
    }

    [[nodiscard]] bool cancelRequested() const
    {
        return this->m_cancel_requested.load(std::memory_order_relaxed);
    }

    /// Reason stepUpdate() last returned false, or StopReason::None if it
    /// hasn't (yet). Only meaningful after evolve() returns.
    [[nodiscard]] StopReason stopReason() const
    {
        return this->m_stop_reason;
    }

    /// Invoked once per stepUpdate() checkpoint with the current
    /// ProgressSnapshot — an alternative (or complement) to the spdlog
    /// console output gated by setLogProgress()/setLogFrequency(), for a
    /// caller that wants structured progress instead of/alongside log lines.
    void setProgressObserver(std::function<void(const ProgressSnapshot&)> observer)
    {
        this->m_progress_observer = std::move(observer);
    }

    void setStepsPerOutput(std::size_t val);
    void setStepsPerOutput1D(std::size_t val);
    void setStepsPerOutput2D(std::size_t val);
    void setStepsPerOutputTrack(std::size_t val);
    void setMaxReports(std::size_t val);
    void setMaxReports1D(std::size_t val);
    void setMaxReports2D(std::size_t val);

    void addReport(std::unique_ptr<ReportData1D> def);
    void addReport(std::unique_ptr<ReportData2D> def);
    void addReport(std::unique_ptr<TrackData> def);

    /// Forwards to the internal ReportHandler's sink (see reporthandler.h) —
    /// called with each report definition's name and JSON payload, in
    /// addition to (or, with setWriteReportFiles(false), instead of) the
    /// default filesystem write.
    void setReportSink(std::function<void(std::string_view, const nlohmann::json&)> sink);
    void setWriteReportFiles(bool val);

    [[nodiscard]] const std::filesystem::path& dirPath() const
    {
        return this->m_dir_path;
    }

    void report1D(double t);
    void report2D(double t);
    void reportTrack(double t);
    void report(double t);
    void reportStats() const;
    [[nodiscard]] bool stepUpdate(double t);

    /// Same checkpoint as stepUpdate(double), plus the step size the solver
    /// just took — carried through to ProgressSnapshot::currentStepSize.
    /// SolverCV_AS/SolverCV_CS::evolve() call this overload; the single-arg
    /// overload above delegates here with currentStepSize left unset.
    [[nodiscard]] bool stepUpdate(double t, double dt);

private:
    static constexpr std::size_t DEFAULT_MAX_REPORTS_1D = 500;
    static constexpr std::size_t DEFAULT_MAX_REPORTS_2D = 200;

    std::unique_ptr<ReportHandler> m_report_handler;

    [[nodiscard]] bool ready1D(std::size_t step) const;
    [[nodiscard]] bool ready2D(std::size_t step) const;
    [[nodiscard]] bool readyTrack(std::size_t step) const;
    [[nodiscard]] bool readyForReport(std::size_t step) const;
    [[nodiscard]] bool maxReportReached() const;

    std::filesystem::path m_dir_path;
    std::size_t m_steps_taken{0};
    std::size_t m_report_count1D{0};
    std::size_t m_steps_per_out1D{1};
    std::size_t m_max_reports1D{DEFAULT_MAX_REPORTS_1D};

    std::size_t m_report_count2D{0};
    std::size_t m_steps_per_out2D{1};
    std::size_t m_max_reports2D{DEFAULT_MAX_REPORTS_2D};

    std::size_t m_steps_per_track{1};

    bool m_log_progress{false};
    std::size_t m_log_freq{1};
    std::unique_ptr<detail::StatCenter> m_stat;

    std::optional<double> m_tf;
    std::atomic<bool> m_cancel_requested{false};
    StopReason m_stop_reason{StopReason::None};
    std::function<void(const ProgressSnapshot&)> m_progress_observer;
};

class PropagatorCV : public BasePropagator {
public:
    using BasePropagator::BasePropagator;
    [[nodiscard]] virtual std::vector<dcmplx>& propagator() = 0;
};

} // namespace spida
