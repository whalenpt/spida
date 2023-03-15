#include <fstream>
#include <iostream>
#include "spida/propagator/reporthandler.h"


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
    for(auto const& def : m_defs_1D){
        def->setItem(key,val);
    }
    for(auto const& def : m_defs_2D){
        def->setItem(key,val);
    }
    for(auto const& def : m_tracker_defs){
        def->setItem(key,val);
    }
}

void ReportHandler::report1D(const std::filesystem::path& dir_path,int rep_num) const
{
    std::ofstream os;
    for(auto const& def : m_defs_1D){
        def->setDirPath(dir_path);
        try{
            def->report(os,rep_num);
            os.close();
        } catch(...) {
            std::cout << "Failed to report file: " << def->path(rep_num) << std::endl;
        }
    }
}

void ReportHandler::report2D(const std::filesystem::path& dir_path,int rep_num) const
{
    std::ofstream os;
    for(auto const& def : m_defs_2D){
        def->setDirPath(dir_path);
        try{
            def->report(os,rep_num);
            os.close();
        } catch(...) {
            std::cout << "Failed to report file: " << def->path(rep_num) << std::endl;
        }
    }
}

void ReportHandler::reportTrack(const std::filesystem::path& dir_path,double t) const
{
    std::ofstream os;
    for(auto const& def : m_tracker_defs){
        def->setDirPath(dir_path);
        def->updateTracker(t);
        try{
            def->report(os);
            os.close();
        } catch(...){
            std::cout << "Failed to report file: " << def->path() << std::endl;
        }
    }
}

void ReportHandler::reportData(const std::filesystem::path& dir_path,int repNum) const
{
    report1D(dir_path,repNum);
    report2D(dir_path,repNum);
}

}