#pragma once

#include <filesystem>
#include <memory>
#include <string_view>
#include <vector>

namespace pw{
    class StatCenter;
    class ReportData1D;
    class ReportData2D;
    class TrackData;
}

namespace spida{

class ReportHandler{
    using vec1D = std::vector<std::unique_ptr<pw::ReportData1D>>;
    using vec2D = std::vector<std::unique_ptr<pw::ReportData2D>>;
    using vecTrack = std::vector<std::unique_ptr<pw::TrackData>>;
    public:
        ReportHandler() = default;
        ~ReportHandler();
        void report1D(const std::filesystem::path& dir_path,unsigned repNum) const;
        void report2D(const std::filesystem::path& dir_path,unsigned repNum) const;
        void reportData(const std::filesystem::path& dir_path,unsigned repNum) const;
        void reportTrack(const std::filesystem::path& dir_path,double t) const;
        void addReport(std::unique_ptr<pw::ReportData1D> def);
        void addReport(std::unique_ptr<pw::ReportData2D> def);
        void addReport(std::unique_ptr<pw::TrackData> def);
        void setItem(std::string_view key,double val);
        [[nodiscard]] bool hasData1D() const {return !m_defs_1D.empty();}
        [[nodiscard]] bool hasData2D() const {return !m_defs_2D.empty();}
        [[nodiscard]] bool hasDataTrack() const {return !m_tracker_defs.empty();}
    private:
        vec1D m_defs_1D;
        vec2D m_defs_2D;
        vecTrack m_tracker_defs;
};

}
