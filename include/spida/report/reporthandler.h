
#ifndef REPORTHANDLER_H_
#define REPORTHANDLER_H_ 

#include <vector> 
#include <filesystem>
#include <memory>
#include <pwutils/report/basedata.hpp>

namespace spida{

class ReportHandler{
    public:
        ReportHandler() {}; 
        void report1D(const std::filesystem::path& dir_path,int repNum) const;
        void report2D(const std::filesystem::path& dir_path,int repNum) const;
        void reportData(const std::filesystem::path& dir_path,int repNum) const;
        void reportTrack(const std::filesystem::path& dir_path,double t) const;
		void addReport(std::unique_ptr<pw::ReportData1D> def);
		void addReport(std::unique_ptr<pw::ReportData2D> def);
		void addReport(std::unique_ptr<pw::TrackData> def);
		void setItem(const std::string& key,double val);
    private:
        std::vector<std::unique_ptr<pw::ReportData1D>> m_defs_1D;
        std::vector<std::unique_ptr<pw::ReportData2D>> m_defs_2D;
        std::vector<std::unique_ptr<pw::TrackData>> m_tracker_defs;
};


}

#endif



