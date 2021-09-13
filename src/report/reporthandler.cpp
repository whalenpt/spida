
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

void ReportHandler::addReport(std::unique_ptr<pw::ReportData2D> def) 
{
    m_defs_2D.push_back(std::move(def));
}

void ReportHandler::addReport(std::unique_ptr<pw::TrackData> def) 
{
    m_tracker_defs.push_back(std::move(def));
}

void ReportHandler::setItem(const std::string& key,double val)
{
    std::vector<std::unique_ptr<pw::ReportData1D>>::const_iterator it; 
    for(it = m_defs_1D.begin(); it != m_defs_1D.end(); it++)
        (*it)->setItem(key,val);
    std::vector<std::unique_ptr<pw::ReportData2D>>::const_iterator it2; 
    for(it2 = m_defs_2D.begin(); it2 != m_defs_2D.end(); it2++)
        (*it2)->setItem(key,val);
    std::vector<std::unique_ptr<pw::TrackData>>::const_iterator it3; 
    for(it3 = m_tracker_defs.begin(); it3 != m_tracker_defs.end(); it3++)
        (*it3)->setItem(key,val);
}

void ReportHandler::report1D(const std::filesystem::path& dir_path,int repNum) const
{
    std::ofstream os;
    std::vector<std::unique_ptr<pw::ReportData1D>>::const_iterator it; 
    for(it = m_defs_1D.begin(); it != m_defs_1D.end(); it++){
        std::filesystem::path file_path = (*it)->generatePath(dir_path,repNum);
        try{
            os.open(file_path.string());
            (*it)->report(os);
            os.close();
        } catch(...){
            std::cout << "Failed to report file: " << file_path.string() << std::endl;
        }
    }
}

void ReportHandler::report2D(const std::filesystem::path& dir_path,int repNum) const
{
    std::ofstream os;
    std::vector<std::unique_ptr<pw::ReportData2D>>::const_iterator it; 
    for(it = m_defs_2D.begin(); it != m_defs_2D.end(); it++){
        std::filesystem::path file_path = (*it)->generatePath(dir_path,repNum);
        try{
            os.open(file_path.string());
            os << *it;
            os.close();
        } catch(...){
            std::cout << "Failed to report file: " << file_path.string() << std::endl;
        }
    }
}

void ReportHandler::reportTrack(const std::filesystem::path& dir_path,double t) const
{
    std::ofstream os;
    std::vector<std::unique_ptr<pw::TrackData>>::const_iterator it; 
    for(it = m_tracker_defs.begin(); it != m_tracker_defs.end(); it++){
        std::filesystem::path file_path = (*it)->generatePath(dir_path);
        (*it)->updateTracker(t);
        try{
            os.open(file_path.string());
            (*it)->report(os);
            os.close();
        } catch(...){
            std::cout << "Failed to report file: " << file_path.string() << std::endl;
        }
    }
}

void ReportHandler::reportData(const std::filesystem::path& dir_path,int repNum) const
{
  report1D(dir_path,repNum);
  report2D(dir_path,repNum);
}





}

