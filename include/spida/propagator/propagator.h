#pragma once

#include "spida/helper/constants.h"
#include "spida/propagator/reporthandler.h"

#include <cstddef>
#include <filesystem>
#include <memory>
#include <vector>

namespace pw {
class StatCenter;
class ReportData1D;
class ReportData2D;
class TrackData;
} // namespace pw

namespace spida {

class BasePropagator {
public:
    explicit BasePropagator(const std::filesystem::path& dir_path);
    virtual ~BasePropagator();
    virtual void updateFields(double t) = 0;

    [[nodiscard]] ReportHandler& reportHandler()
    {
        return this->m_report_handler;
    }

    [[nodiscard]] const ReportHandler& reportHandler() const
    {
        return this->m_report_handler;
    }

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

    void setStepsPerOutput(std::size_t val);
    void setStepsPerOutput1D(std::size_t val);
    void setStepsPerOutput2D(std::size_t val);
    void setStepsPerOutputTrack(std::size_t val);
    void setMaxReports(std::size_t val);
    void setMaxReports1D(std::size_t val);
    void setMaxReports2D(std::size_t val);

    void addReport(std::unique_ptr<pw::ReportData1D> def);
    void addReport(std::unique_ptr<pw::ReportData2D> def);
    void addReport(std::unique_ptr<pw::TrackData> def);

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

private:
    static constexpr std::size_t DEFAULT_MAX_REPORTS_1D = 500;
    static constexpr std::size_t DEFAULT_MAX_REPORTS_2D = 200;

    ReportHandler m_report_handler;

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
    std::unique_ptr<pw::StatCenter> m_stat;
};

class PropagatorCV : public BasePropagator {
public:
    using BasePropagator::BasePropagator;
    [[nodiscard]] virtual std::vector<dcmplx>& propagator() = 0;
};

} // namespace spida
