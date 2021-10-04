
#include "spida/propagator/reporthandler.h"
#include <fstream>
#include <iostream>

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
    for(auto it = m_defs_1D.begin(); it != m_defs_1D.end(); it++)
        (*it)->setItem(key,val);
    for(auto it = m_defs_2D.begin(); it != m_defs_2D.end(); it++)
        (*it)->setItem(key,val);
    for(auto it = m_tracker_defs.begin(); it != m_tracker_defs.end(); it++)
        (*it)->setItem(key,val);
}

void ReportHandler::report1D(const std::filesystem::path& dir_path,int rep_num) const
{
    std::ofstream os;
    for(vec1D::const_iterator it = m_defs_1D.cbegin(); it != m_defs_1D.cend(); it++){
        (*it)->setDirPath(dir_path);
        try{
            (*it)->report(os,rep_num);
            os.close();
        } catch(...){
            std::cout << "Failed to report file: " << (*it)->path(rep_num) << std::endl;
        }
    }
}

void ReportHandler::report2D(const std::filesystem::path& dir_path,int rep_num) const
{
    std::ofstream os;
    for(vec2D::const_iterator it = m_defs_2D.cbegin(); it != m_defs_2D.cend(); it++){
        (*it)->setDirPath(dir_path);
        try{
            (*it)->report(os,rep_num);
            os.close();
        } catch(...){
            std::cout << "Failed to report file: " << (*it)->path(rep_num) << std::endl;
        }
    }
}

void ReportHandler::reportTrack(const std::filesystem::path& dir_path,double t) const
{
    std::ofstream os;
    for(vecTrack::const_iterator it = m_tracker_defs.cbegin(); it != m_tracker_defs.cend(); it++){
        (*it)->setDirPath(dir_path);
        (*it)->updateTracker(t);
        try{
            (*it)->report(os);
            os.close();
        } catch(...){
            std::cout << "Failed to report file: " << (*it)->path() << std::endl;
        }
    }
}

void ReportHandler::reportData(const std::filesystem::path& dir_path,int repNum) const
{
    report1D(dir_path,repNum);
    report2D(dir_path,repNum);
}





}

