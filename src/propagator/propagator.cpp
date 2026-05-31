
#include "spida/propagator/propagator.h"

#include "propagator/reporthandler.h"

#include <filesystem>
#include <format>
#include <stdexcept>

#include <pwutils/pwstats.h>
#include <pwutils/report.hpp>

namespace spida {

namespace {
void validatePositive(std::size_t val, const char* name)
{
    if (val == 0)
        throw std::invalid_argument(
            std::format("BasePropagator::{} requires a value greater than zero.", name));
}
} // namespace

BasePropagator::~BasePropagator() = default;

bool BasePropagator::hasData1D() const
{
    return this->m_report_handler->hasData1D();
}

void BasePropagator::addReport(std::unique_ptr<pw::ReportData1D> def)
{
    this->m_report_handler->addReport(std::move(def));
}

void BasePropagator::addReport(std::unique_ptr<pw::ReportData2D> def)
{
    this->m_report_handler->addReport(std::move(def));
}

void BasePropagator::addReport(std::unique_ptr<pw::TrackData> def)
{
    this->m_report_handler->addReport(std::move(def));
}

BasePropagator::BasePropagator(const std::filesystem::path& dir_path)
    : m_report_handler(std::make_unique<ReportHandler>()),
      m_dir_path(dir_path),
      m_stat(std::make_unique<pw::StatCenter>())
{
    this->m_stat->setHeader("REPORT STATS");
    this->m_stat->addTracker("t", 0.0);
    this->m_stat->addTimer("Time Reporting 1D");
    this->m_stat->addTimer("Time Reporting 2D");
    this->m_stat->addTimer("Time Reporting Trackers");
    this->m_stat->addCounter("Number Reports 1D");
    this->m_stat->addCounter("Number Reports 2D");
    this->m_stat->addCounter("Number Reports Track");
}

void BasePropagator::setStepsPerOutput(std::size_t val)
{
    this->setStepsPerOutput1D(val);
    this->setStepsPerOutput2D(val);
    this->setStepsPerOutputTrack(val);
}

void BasePropagator::setStepsPerOutput1D(std::size_t val)
{
    validatePositive(val, "setStepsPerOutput1D");
    this->m_steps_per_out1D = val;
}

void BasePropagator::setStepsPerOutput2D(std::size_t val)
{
    validatePositive(val, "setStepsPerOutput2D");
    this->m_steps_per_out2D = val;
}

void BasePropagator::setStepsPerOutputTrack(std::size_t val)
{
    validatePositive(val, "setStepsPerOutputTrack");
    this->m_steps_per_track = val;
}

void BasePropagator::setMaxReports(std::size_t val)
{
    this->setMaxReports1D(val);
    this->setMaxReports2D(val);
}

void BasePropagator::setMaxReports1D(std::size_t val)
{
    validatePositive(val, "setMaxReports1D");
    this->m_max_reports1D = val;
}

void BasePropagator::setMaxReports2D(std::size_t val)
{
    validatePositive(val, "setMaxReports2D");
    this->m_max_reports2D = val;
}

void BasePropagator::setLogFrequency(std::size_t val)
{
    validatePositive(val, "setLogFrequency");
    this->m_log_freq = val;
}

bool BasePropagator::ready1D(std::size_t step) const
{
    return !(step % this->m_steps_per_out1D) && this->m_report_handler->hasData1D();
}

bool BasePropagator::ready2D(std::size_t step) const
{
    return !(step % this->m_steps_per_out2D) && this->m_report_handler->hasData2D();
}

bool BasePropagator::readyTrack(std::size_t step) const
{
    return !(step % this->m_steps_per_track) && this->m_report_handler->hasDataTrack();
}

bool BasePropagator::readyForReport(std::size_t step) const
{
    return this->ready1D(step) || this->ready2D(step) || this->readyTrack(step);
}

bool BasePropagator::maxReportReached() const
{
    return this->m_report_count1D >= this->m_max_reports1D ||
           this->m_report_count2D >= this->m_max_reports2D;
}

bool BasePropagator::stepUpdate(double t)
{
    const auto step = this->m_steps_taken + 1;
    this->m_stat->updateTracker("t", t);
    if (this->readyForReport(step))
        this->updateFields(t);
    if (this->ready1D(step))
        this->report1D(t);
    if (this->ready2D(step))
        this->report2D(t);
    if (this->readyTrack(step))
        this->reportTrack(t);
    if (this->m_log_progress && !(step % this->m_log_freq))
        this->reportStats();
    this->m_steps_taken = step;
    return !this->maxReportReached();
}

void BasePropagator::reportStats() const
{
    this->m_stat->report();
}

void BasePropagator::report(double t)
{
    if (this->m_report_handler->hasData1D())
        this->report1D(t);
    if (this->m_report_handler->hasData2D())
        this->report2D(t);
    if (this->m_report_handler->hasDataTrack())
        this->reportTrack(t);
}

void BasePropagator::report1D(double t)
{
    if (!this->m_report_handler->hasData1D())
        return;
    this->m_stat->startTimer("Time Reporting 1D");
    this->m_report_handler->setItem("t", t);
    this->m_report_handler->report1D(this->m_dir_path, this->m_report_count1D);
    this->m_stat->endTimer("Time Reporting 1D");
    this->m_report_count1D++;
    this->m_stat->incrementCounter("Number Reports 1D");
}

void BasePropagator::report2D(double t)
{
    if (!this->m_report_handler->hasData2D())
        return;
    this->m_stat->startTimer("Time Reporting 2D");
    this->m_report_handler->setItem("t", t);
    this->m_report_handler->report2D(this->m_dir_path, this->m_report_count2D);
    this->m_stat->endTimer("Time Reporting 2D");
    this->m_report_count2D++;
    this->m_stat->incrementCounter("Number Reports 2D");
}

void BasePropagator::reportTrack(double t)
{
    if (!this->m_report_handler->hasDataTrack())
        return;
    this->m_stat->startTimer("Time Reporting Trackers");
    this->m_report_handler->setItem("t", t);
    this->m_report_handler->reportTrack(this->m_dir_path, t);
    this->m_stat->endTimer("Time Reporting Trackers");
    this->m_stat->incrementCounter("Number Reports Track");
}

} // namespace spida
