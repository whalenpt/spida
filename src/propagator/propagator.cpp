
#include <stdexcept>
#include <filesystem>
#include "spida/propagator/propagator.h"

namespace spida{

BasePropagator::BasePropagator(const std::filesystem::path& dir_path) :
    m_dir_path(dir_path),
    m_steps_taken(0),
    m_report_count1D(0),
    m_steps_per_out1D(1),
    m_max_reports1D(DEFAULT_MAX_REPORTS_1D),
    m_report_count2D(0),
    m_steps_per_out2D(1),
    m_max_reports2D(DEFAULT_MAX_REPORTS_2D),
    m_steps_per_track(1),
    m_log_progress(false),
    m_log_freq(1)
{
    m_stat.setHeader("REPORT STATS");
    m_stat.addTracker("t",0.0);
}

void BasePropagator::setStepsPerOutput(unsigned val)
{
    setStepsPerOutput1D(val);
    setStepsPerOutput2D(val);
    setStepsPerOutputTrack(val);
}

void BasePropagator::setStepsPerOutput1D(unsigned val)
{
    if(val < 1)
        throw std::invalid_argument("BasePropagator::setStepsPerOutput1D requires a value \
                greater than zero.");
    m_steps_per_out1D = val;
}

void BasePropagator::setStepsPerOutput2D(unsigned val)
{
    if(val < 1)
        throw std::invalid_argument("BasePropagator::setStepsPerOutput2D requires a value \
                greater than zero.");
    m_steps_per_out2D = val;
}

void BasePropagator::setStepsPerOutputTrack(unsigned val)
{
    if(val < 1)
        throw std::invalid_argument("BasePropagator::setStepsPerOutputTrack requires a value \
                greater than zero.");
    m_steps_per_track = val;
}

void BasePropagator::setMaxReports(unsigned val)
{
    setMaxReports1D(val);
    setMaxReports2D(val);
}

void BasePropagator::setMaxReports1D(unsigned val)
{
    if(val < 1)
        throw std::invalid_argument("BasePropagator::setMaxReports1D requires a value \
                greater than zero.");
    m_max_reports1D = val;
}

void BasePropagator::setMaxReports2D(unsigned val)
{
    if(val < 1)
        throw std::invalid_argument("BasePropagator::setMaxReports2D requires a value \
                greater than zero.");
    m_max_reports2D = val;
}


bool BasePropagator::readyForReport() const
{
    if(!(m_steps_taken % m_steps_per_out1D) && m_report_handler.hasData1D()) 
        return true;
    if(!(m_steps_taken % m_steps_per_out2D) && m_report_handler.hasData2D())
        return true;
    if(!(m_steps_taken % m_steps_per_track) && m_report_handler.hasDataTrack())
        return true;
    return false;
}

bool BasePropagator::maxReportReached() const
{
    if(m_report_count1D >= m_max_reports1D)
        return true;
    if(m_report_count2D >= m_max_reports2D)
        return true;
    return false;
}

bool BasePropagator::stepUpdate(double t)
{
    m_steps_taken++;
    m_stat.updateTracker("t",t);
    if(readyForReport())
        updateFields(t);
    if(!(m_steps_taken % m_steps_per_out1D) && m_report_handler.hasData1D())
        report1D(t);
    if(!(m_steps_taken % m_steps_per_out2D) && m_report_handler.hasData2D()) 
        report2D(t);
    if(!(m_steps_taken % m_steps_per_track) && m_report_handler.hasDataTrack())
        reportTrack(t);

    if(m_log_progress && !(m_steps_taken % m_log_freq))
        reportStats();
    if(maxReportReached())
        return false;
    return true;
}

void BasePropagator::reportStats() const
{
    m_stat.report();
}

void BasePropagator::report(double t)
{
    report1D(t);
    report2D(t);
    reportTrack(t);
}

void BasePropagator::report1D(double t) 
{
    m_stat.startTimer("Time Reporting 1D");
    m_report_handler.setItem("t",t);
    m_report_handler.report1D(m_dir_path,m_report_count1D);
    m_stat.endTimer("Time Reporting 1D");
    m_report_count1D++;
    m_stat.incrementCounter("Number Reports 1D");
}

void BasePropagator::report2D(double t) 
{
    m_stat.startTimer("Time Reporting 2D");
    m_report_handler.setItem("t",t);
    m_report_handler.report2D(m_dir_path,m_report_count2D);
    m_stat.endTimer("Time Reporting 2D");
    m_report_count2D++;
    m_stat.incrementCounter("Number Reports 2D");
}

void BasePropagator::reportTrack(double t) 
{
    m_stat.startTimer("Time Reporting Trackers");
    m_report_handler.setItem("t",t);
    m_report_handler.reportTrack(m_dir_path,t);
    m_stat.endTimer("Time Reporting Trackers");
}



}



