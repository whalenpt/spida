
#ifndef REPORTHANDLER_H_
#define REPORTHANDLER_H_ 

#include <string> 
#include <vector> 
#include <memory>
#include <filesystem>
#include <pwutils/report/basedata.hpp>

namespace spida{

class ReportHandler{
    public:
        ReportHandler() {}; 
        void report(const std::filesystem::path& dir_path,int repNum,double t) const;
        void report1D(const std::filesystem::path& dir_path,int repNum,double t) const;
        void report2D(const std::filesystem::path& dir_path,int repNum,double t) const;
		void addReport(std::unique_ptr<pw::ReportData1D> def);
//		void addReport(std::unique_ptr<pw::VBReportData2D> def);
		void addReport(std::unique_ptr<pw::TrackData> def);

    private:
        std::vector<std::unique_ptr<pw::ReportData1D>> m_defs_1D;
 //       std::vector<std::unique_ptr<pw::VBReportData2D>> m_defs_2D;
        std::vector<std::unique_ptr<pw::TrackData>> m_tracker_defs;
};


}

#endif



