
#ifndef REPORTCENTER_H_
#define REPORTCENTER_H_ 

#include "spida/propagator/propagator.h"
#include <complex> 
#include <string> 
#include <vector> 
#include <ctime>
#include <filesystem>
#include <pwutils/pwstats.h>

namespace spida{

class ReportCenter{
    public:
        ReportCenter(const std::filesystem::path& dir_path,int max_report = 200);
        virtual ~ReportCenter() {}; 
        virtual bool stepUpdate(double dist) = 0;
        virtual void report(double dist) = 0;
        virtual void setStepsPerOutput1D(int val) = 0;
        virtual void setStepsPerOutput2D(int val) = 0;

        void addPropagator(Propagators* pr) {m_prs = pr;}
        void setMaxReports(int val) {m_max_reports = val;}
        void setLogProgress(bool val) {m_log_progress = val;}
        std::string dirName() {return m_dir_path.string();}
        std::filesystem::path dirPath() {return m_dir_path;}

    protected:
        pw::StatCenter m_stat;
        int m_steps_taken;
        int m_max_reports;
        Propagators* m_prs;
    	const std::filesystem::path m_dir_path;
        bool m_log_progress;
};

class ReportCenter1D : public ReportCenter
{
    public:
        ReportCenter1D(const std::filesystem::path& dir_path,int step_per1D = 1,
                 int max_report = 1000); 
        virtual ~ReportCenter1D() {}; 
        bool stepUpdate(double dist);
        void report(double dist);
        virtual void setStepsPerOutput1D(int val) {m_stepsPerOutput1D = val;}
        virtual void setStepsPerOutput2D(int val) {};

    protected:
        void report1D(double dist);
        int m_reportCount1D;
        int m_stepsPerOutput1D;
};

class ReportCenter2D : public ReportCenter1D
{
    public:
        ReportCenter2D(const std::filesystem::path& dir_path,int step_per1D = 1,int step_per2D=1,\
                int max_report = 250); 
        bool stepUpdate(double dist);
        void report(double dist);
        virtual ~ReportCenter2D() {}; 
        void setStepsPerOutput1D(int val) {m_stepsPerOutput1D = val;}
        void setStepsPerOutput2D(int val) {m_stepsPerOutput2D = val;}
    protected:
        void report2D(double dist);
        void report1D2D(double dist);
        int m_reportCount2D;
        int m_stepsPerOutput2D;
};



}

#endif


