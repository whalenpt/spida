
#include "spida/report/reporthandler.h"
#include <pwutils/report/basereport.h>
#include <memory>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <filesystem>

namespace spida{

void ReportHandler::addReport(std::unique_ptr<pw::ReportData1D> def) 
{
    m_defs_1D.push_back(std::move(def));
}

/*
void ReportHandler::addReport(std::unique_ptr<pw::VBReportData2D> def) 
{
    m_defs_2D.push_back(std::move(def));
}
*/

void ReportHandler::addReport(std::unique_ptr<pw::TrackData> def) 
{
    m_tracker_defs.push_back(std::move(def));
}

#include <iostream>

void ReportHandler::report1D(const std::filesystem::path& dir_path,int repNum,double t) const
{
    std::ofstream os;
    std::vector<std::unique_ptr<pw::ReportData1D>>::const_iterator it; 
    for(it = m_defs_1D.begin(); it != m_defs_1D.end(); it++){
        std::filesystem::path file_path = (*it)->generatePath(dir_path,repNum);
        (*it)->setItem("sim_duration",t);
        try{
            os.open(file_path.string());
            (*it)->report(os);
            os.close();
        } catch(...){
            std::cout << "Failed to report file: " << file_path.string() << std::endl;
        }
    }

    std::vector<std::unique_ptr<pw::TrackData>>::const_iterator tit; 
    for(tit = m_tracker_defs.begin(); tit != m_tracker_defs.end(); tit++){
        std::filesystem::path file_path = (*tit)->generatePath(dir_path);
        (*tit)->updateTracker(t);
        (*tit)->setItem("sim_duration",t);
        try{
            os.open(file_path.string());
            (*tit)->report(os);
            os.close();
        } catch(...){
            std::cout << "Failed to report file: " << file_path.string() << std::endl;
        }
    }
}

void ReportHandler::report2D(const std::filesystem::path& dir_path,int repNum,double t) const
{
    /*
    std::ofstream os;
    std::vector<std::unique_ptr<pw::VBReportData2D>>::const_iterator it; 
    for(it = m_defs_2D.begin(); it != m_defs_2D.end(); it++){
        std::filesystem::path file_path = (*it)->filePath(dir_path,repNum);
        (*it)->setItem("time_elapsed",t);
        try{
            os.open(file_path.string());
            os << *it;
            os.close();
        } catch(...){
            std::cout << "Failed to report file: " << file_path.string() << std::endl;
        }
    }
    */
}

void ReportHandler::report(const std::filesystem::path& dir_path,int repNum,double t) const
{
  report1D(dir_path,repNum,t);
  report2D(dir_path,repNum,t);
}





}

