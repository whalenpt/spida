#pragma once

#include <filesystem>
#include <memory>
#include <vector> 
#include <string>
#include <pwutils/report/basereport.h>
// basereport has defintions for ReportData1D, ReportData2D, and TrackData

namespace spida{

class ReportHandler{
    using vec1D = std::vector<std::unique_ptr<pw::ReportData1D>>;
    using vec2D = std::vector<std::unique_ptr<pw::ReportData2D>>;
    using vecTrack = std::vector<std::unique_ptr<pw::TrackData>>;
    public:
        ReportHandler() = default;
        void report1D(const std::filesystem::path& dir_path,int repNum) const;
        void report2D(const std::filesystem::path& dir_path,int repNum) const;
        void reportData(const std::filesystem::path& dir_path,int repNum) const;
        void reportTrack(const std::filesystem::path& dir_path,double t) const;
		void addReport(std::unique_ptr<pw::ReportData1D> def);
		void addReport(std::unique_ptr<pw::ReportData2D> def);
		void addReport(std::unique_ptr<pw::TrackData> def);
		void setItem(const std::string& key,double val);
		bool hasData1D() const {return (m_defs_1D.empty() ? false : true);}
		bool hasData2D() const {return (m_defs_2D.empty() ? false : true);}
		bool hasDataTrack() const {return (m_tracker_defs.empty() ? false : true);}
    private:

        vec1D m_defs_1D;
        vec2D m_defs_2D;
        vecTrack m_tracker_defs;
};

}