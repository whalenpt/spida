
#include "spida/report/reportcenter.h"
#include "spida/propagator/propagator.h"

namespace spida{

ReportCenter::ReportCenter(Propagator* pr,\
        const std::filesystem::path& dir_path,\
        unsigned max_report) :
    m_pr(pr),
    m_dir_path(dir_path)
{
    m_steps_taken = 0;
    m_max_reports = max_report;
    m_stat.setHeader("REPORT STATS");
}

void ReportCenter::setDirPath(const std::filesystem::path& dir_path)
{
    m_dir_path = dir_path;
}

void ReportCenter::setPropagator(Propagator* pr)
{
    m_pr = pr;
}

ReportCenter1D::ReportCenter1D(Propagator* pr,const std::filesystem::path& dir_path,\
        unsigned step_per1D,\
        unsigned max_report) :
    ReportCenter(pr,dir_path,max_report)
{
    m_reportCount1D = 0;
    m_stepsPerOutput1D = step_per1D;
    m_stat.addCounter("Number Report1D",1);
    m_stat.addTimer("Time Report1D");
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
    if(!m_pr)
        return;

    m_pr->updateFields(t);
    m_stat.startTimer("Time Report1D");
    m_pr->reportHandler().setItem("sim_duration",t);
    m_pr->reportHandler().report1D(m_dir_path,m_reportCount1D);
    m_pr->reportHandler().reportTrack(m_dir_path,t);
    m_stat.endTimer("Time Report1D");
    m_reportCount1D++;
    m_stat.incrementCounter("Number Report1D");
    if(m_log_progress)
        m_stat.report();
}

void ReportCenter1D::report(double t)
{
    report1D(t);
}


ReportCenter2D::ReportCenter2D(Propagator* pr,const std::filesystem::path& dir_path,\
        unsigned step_per1D,\
        unsigned step_per2D,
        unsigned max_report) :
    ReportCenter1D(pr,dir_path,step_per1D,max_report)
{
    m_reportCount2D = 0;
    m_stepsPerOutput2D = step_per2D;
    m_stat.addCounter("Number Report2D",1);
    m_stat.addTimer("Time Report2D");
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
    if(!m_pr)
        return;

    m_pr->updateFields(t);
    m_stat.startTimer("Time Report2D");
    m_pr->reportHandler().setItem("sim_duration",t);
    m_pr->reportHandler().report2D(m_dir_path,m_reportCount2D);
    m_stat.endTimer("Time Report2D");
    m_reportCount2D++;
    m_stat.incrementCounter("Number Report2D");
    if(m_log_progress)
        m_stat.report();
}


void ReportCenter2D::report1D2D(double t) 
{
    if(!m_pr)
        return;

    m_pr->updateFields(t);

    m_pr->reportHandler().setItem("sim_duration",t);
    m_stat.startTimer("Time Report1D");
    m_pr->reportHandler().report1D(m_dir_path,m_reportCount1D);
    m_pr->reportHandler().reportTrack(m_dir_path,t);
    m_stat.endTimer("Time Report1D");

    m_stat.startTimer("Time Report2D");
    m_pr->reportHandler().report2D(m_dir_path,m_reportCount2D);
    m_stat.endTimer("Time Report2D");

    m_reportCount1D++;
    m_reportCount2D++;
    m_stat.incrementCounter("Number Report1D");
    m_stat.incrementCounter("Number Report2D");
    if(m_log_progress)
        m_stat.report();
}

}



