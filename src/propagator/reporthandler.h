#pragma once

#include <cstddef>
#include <filesystem>
#include <functional>
#include <memory>
#include <string_view>
#include <vector>

#include <nlohmann/json.hpp>

namespace detail {
class StatCenter;
} // namespace detail

namespace spida {

class ReportData1D;
class ReportData2D;
class TrackData;

class ReportHandler {
public:
    /// Called once per report definition, per report1D/report2D/reportTrack
    /// call, with that definition's name and the same nlohmann::json object
    /// the default filesystem write serializes (ReportBase::toJson() in
    /// utils/report.hpp) — a seam for a non-filesystem sink (e.g. an
    /// artifact store / event stream) to observe every report without
    /// touching BasePropagator. Runs in addition to the file write unless
    /// setWriteFiles(false) turns the latter off.
    using Sink = std::function<void(std::string_view name, const nlohmann::json&)>;

    ReportHandler() = default;
    ~ReportHandler();
    void report1D(const std::filesystem::path& dir_path, std::size_t rep_num);
    void report2D(const std::filesystem::path& dir_path, std::size_t rep_num);
    void reportData(const std::filesystem::path& dir_path, std::size_t rep_num);
    void reportTrack(const std::filesystem::path& dir_path, double t);
    void addReport(std::unique_ptr<ReportData1D> def);
    void addReport(std::unique_ptr<ReportData2D> def);
    void addReport(std::unique_ptr<TrackData> def);
    void setItem(std::string_view key, double val);

    void setSink(Sink sink)
    {
        this->m_sink = std::move(sink);
    }

    /// Default true (matches existing behavior). Set false to skip the
    /// std::ofstream file write entirely and rely solely on the sink —
    /// e.g. a worker that doesn't want redundant local disk I/O.
    void setWriteFiles(bool val)
    {
        this->m_write_files = val;
    }

    [[nodiscard]] bool hasData1D() const
    {
        return !m_defs_1D.empty();
    }

    [[nodiscard]] bool hasData2D() const
    {
        return !m_defs_2D.empty();
    }

    [[nodiscard]] bool hasDataTrack() const
    {
        return !m_tracker_defs.empty();
    }

private:
    using vec1D = std::vector<std::unique_ptr<ReportData1D>>;
    using vec2D = std::vector<std::unique_ptr<ReportData2D>>;
    using vecTrack = std::vector<std::unique_ptr<TrackData>>;

    vec1D m_defs_1D;
    vec2D m_defs_2D;
    vecTrack m_tracker_defs;
    Sink m_sink;
    bool m_write_files{true};
};

} // namespace spida
