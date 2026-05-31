#include "propagator/reporthandler.h"

#include <format>
#include <fstream>
#include <iostream>

#include <pwutils/pwstats.h>
#include <pwutils/report.hpp>

namespace spida {

ReportHandler::~ReportHandler() = default;

void ReportHandler::addReport(std::unique_ptr<pw::ReportData1D> def)
{
    this->m_defs_1D.push_back(std::move(def));
}

void ReportHandler::addReport(std::unique_ptr<pw::ReportData2D> def)
{
    this->m_defs_2D.push_back(std::move(def));
}

void ReportHandler::addReport(std::unique_ptr<pw::TrackData> def)
{
    this->m_tracker_defs.push_back(std::move(def));
}

void ReportHandler::setItem(std::string_view key, double val)
{
    for (auto const& def : this->m_defs_1D)
        def->setItem(std::string(key), val);
    for (auto const& def : this->m_defs_2D)
        def->setItem(std::string(key), val);
    for (auto const& def : this->m_tracker_defs)
        def->setItem(std::string(key), val);
}

void ReportHandler::report1D(const std::filesystem::path& dir_path, std::size_t rep_num)
{
    for (auto const& def : this->m_defs_1D) {
        def->setDirPath(dir_path);
        std::ofstream os;
        try {
            def->report(os, static_cast<unsigned>(rep_num));
        }
        catch (const std::exception& e) {
            std::cerr << std::format("Failed to report file: {} — {}\n",
                                     def->path(static_cast<int>(rep_num)).string(),
                                     e.what());
        }
        catch (...) {
            std::cerr << std::format("Failed to report file: {} — unknown error\n",
                                     def->path(static_cast<int>(rep_num)).string());
        }
    }
}

void ReportHandler::report2D(const std::filesystem::path& dir_path, std::size_t rep_num)
{
    for (auto const& def : this->m_defs_2D) {
        def->setDirPath(dir_path);
        std::ofstream os;
        try {
            def->report(os, static_cast<unsigned>(rep_num));
        }
        catch (const std::exception& e) {
            std::cerr << std::format("Failed to report file: {} — {}\n",
                                     def->path(static_cast<int>(rep_num)).string(),
                                     e.what());
        }
        catch (...) {
            std::cerr << std::format("Failed to report file: {} — unknown error\n",
                                     def->path(static_cast<int>(rep_num)).string());
        }
    }
}

void ReportHandler::reportTrack(const std::filesystem::path& dir_path, double t)
{
    for (auto const& def : this->m_tracker_defs) {
        def->setDirPath(dir_path);
        def->updateTracker(t);
        std::ofstream os;
        try {
            def->report(os);
        }
        catch (const std::exception& e) {
            std::cerr << std::format(
                "Failed to report file: {} — {}\n", def->path().string(), e.what());
        }
        catch (...) {
            std::cerr << std::format("Failed to report file: {} — unknown error\n",
                                     def->path().string());
        }
    }
}

void ReportHandler::reportData(const std::filesystem::path& dir_path, std::size_t rep_num)
{
    this->report1D(dir_path, rep_num);
    this->report2D(dir_path, rep_num);
}

} // namespace spida
