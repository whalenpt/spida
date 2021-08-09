
#include "spida/report/reportcenter.h"
#include "spida/propagator/propagator.h"

namespace spida{

ReportCenter::ReportCenter(const std::filesystem::path& dir_path,int max_report) :
    m_dir_path(dir_path)
{
    m_prs = nullptr;
    m_steps_taken = 0;
    m_max_reports = max_report;
    m_stat.setHeader("REPORT STATS");
}


ReportCenter1D::ReportCenter1D(const std::filesystem::path& dir_path,int step_per1D,int max_report) :
    ReportCenter(dir_path,max_report)
{
    m_reportCount1D = 0;
    m_stepsPerOutput1D = step_per1D;
    m_stat.addCounter("Number Reports1D",1);
    m_stat.addTimer("Time Reports1D");
    m_stat.addTracker("Time Elapsed",0.0);
}

bool ReportCenter1D::stepUpdate(double t)
{
    m_steps_taken++;
    m_stat.updateTracker("Time Elapsed",t);
    if(!(m_steps_taken % m_stepsPerOutput1D))
        report1D(t);
    if(m_reportCount1D >= m_max_reports)
        return false;
    else 
        return true;
}

void ReportCenter1D::report1D(double t) 
{
    m_prs->updateFields(t);
    m_stat.startTimer("Time Reports1D");
    m_prs->reportHandler().report1D(m_dir_path,m_reportCount1D,t);
    m_stat.endTimer("Time Reports1D");
    m_reportCount1D++;
    m_stat.incrementCounter("Number Reports1D");
    if(m_log_progress)
        m_stat.report();
}

void ReportCenter1D::report(double t)
{
    report1D(t);
}


ReportCenter2D::ReportCenter2D(const std::filesystem::path& dir_path,int step_per1D,int step_per2D,int max_report) :
    ReportCenter1D(dir_path,step_per1D,max_report)
{
    m_reportCount2D = 0;
    m_stepsPerOutput2D = step_per2D;
    m_stat.addCounter("Number Reports2D",1);
    m_stat.addTimer("Time Reports2D");
}

bool ReportCenter2D::stepUpdate(double t)
{
    m_steps_taken++;
    m_stat.updateTracker("Time Elapsed",t);
    if(!(m_steps_taken % m_stepsPerOutput1D) && !(m_steps_taken % m_stepsPerOutput2D))
        report1D2D(t);
    else if(!(m_steps_taken % m_stepsPerOutput1D) && (m_steps_taken % m_stepsPerOutput2D) )
        report1D(t);
    else if((m_steps_taken % m_stepsPerOutput1D) && !(m_steps_taken % m_stepsPerOutput2D) )
        report2D(t);
    if(m_reportCount1D >= m_max_reports || m_reportCount2D >= m_max_reports)
        return false;
    else 
        return true;
}

void ReportCenter2D::report(double t) 
{
    report1D2D(t);
}

void ReportCenter2D::report2D(double t) 
{
    m_prs->updateFields(t);
    m_stat.startTimer("Time Reports2D");
    m_prs->reportHandler().report2D(m_dir_path,m_reportCount2D,t);
    m_stat.endTimer("Time Reports2D");
    m_reportCount2D++;
    m_stat.incrementCounter("Number Reports2D");
    if(m_log_progress)
        m_stat.report();
}


void ReportCenter2D::report1D2D(double t) 
{
    m_prs->updateFields(t);

    m_stat.startTimer("Time Reports1D");
    m_prs->reportHandler().report1D(m_dir_path,m_reportCount1D,t);
    m_stat.endTimer("Time Reports1D");

    m_stat.startTimer("Time Reports2D");
    m_prs->reportHandler().report2D(m_dir_path,m_reportCount2D,t);
    m_stat.endTimer("Time Reports2D");

    m_reportCount1D++;
    m_reportCount2D++;
    m_stat.incrementCounter("Number Reports1D");
    m_stat.incrementCounter("Number Reports2D");
    if(m_log_progress)
        m_stat.report();
}

}



