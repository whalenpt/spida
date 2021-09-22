// reportcenter.h
#pragma once

#include <vector> 
#include <ctime>
#include <filesystem>
#include <memory>
#include <pwutils/pwstats.h>

namespace spida{

class Propagator;

class ReportCenter{
    public:
        ReportCenter(Propagator* pr,const std::filesystem::path& dir_path,\
                unsigned max_report);
        virtual ~ReportCenter() {}; 
        std::filesystem::path dirPath() const {return m_dir_path;}
        void setDirPath(const std::filesystem::path& dir_path);
        virtual bool stepUpdate(double dist) = 0;
        virtual void report(double dist) = 0;
        virtual void setStepsPerOutput1D(unsigned val) = 0;
        virtual void setStepsPerOutput2D(unsigned val) {};
        void setPropagator(Propagator* pr);
        Propagator* propagator() {return m_pr;}
        bool fileReportOn() const {return (m_pr ? true : false);}
        void setMaxReports(unsigned val) {m_max_reports = val;}
        void setLogProgress(bool val) {m_log_progress = val;}
    protected:
        Propagator* m_pr;
        pw::StatCenter m_stat;
        unsigned m_steps_taken;
        unsigned m_max_reports;
    	std::filesystem::path m_dir_path;
        bool m_log_progress;
};

class ReportCenter1D : public ReportCenter
{
    public:
        ReportCenter1D(Propagator* pr,const std::filesystem::path& dir_path,\
                unsigned step_per1D = 1,
                unsigned max_report = 1000); 
        virtual ~ReportCenter1D() {}; 
        bool stepUpdate(double dist);
        void report(double dist);
        void setStepsPerOutput1D(unsigned val) {m_stepsPerOutput1D = val;}

    protected:
        void report1D(double dist);
        unsigned m_reportCount1D;
        unsigned m_stepsPerOutput1D;
};

class ReportCenter2D : public ReportCenter1D
{
    public:
        ReportCenter2D(Propagator* pr,const std::filesystem::path& dir_path,\
                unsigned step_per1D = 1,\
                unsigned step_per2D = 1,\
                unsigned max_report = 250); 
        bool stepUpdate(double dist);
        void report(double dist);
        virtual ~ReportCenter2D() {}; 
        void setStepsPerOutput2D(unsigned val) {m_stepsPerOutput2D = val;}
    protected:
        void report2D(double dist);
        void report1D2D(double dist);
        unsigned m_reportCount2D;
        unsigned m_stepsPerOutput2D;
};



}



