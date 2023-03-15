/*
 * Helper class for numerical rkstiff methods. Can be used by rkstiff client code
 * to replace manual file reporting of data.
 */
#pragma once

#include <filesystem>
#include <memory>
#include <string>
#include <vector>
#include <pwutils/pwstats.h>
#include <pwutils/report/basedata.hpp>
#include "spida/helper/constants.h"
#include "spida/propagator/reporthandler.h"

namespace spida{


class BasePropagator
{
    static const int DEFAULT_MAX_REPORTS_1D = 500;
    static const int DEFAULT_MAX_REPORTS_2D = 200;

    public:
        explicit BasePropagator(const std::filesystem::path& dir_path);
        virtual ~BasePropagator() = default;
        // updateFields will update all reporting vectors prior to file output 
        // i.e. transform a vector from spectral space to physical space, or perform some
        // other operation
        virtual void updateFields(double t) = 0;

        ReportHandler& reportHandler() {return m_report_handler;}
        void setDirPath(const std::filesystem::path& dir_path) {m_dir_path = dir_path;}
        void setLogProgress(bool val) {m_log_progress = val;}
        void setLogFrequency(unsigned val) {m_log_freq = val;}
        bool logProgress() const {return m_log_progress;}

        void setStepsPerOutput(unsigned val);
        void setStepsPerOutput1D(unsigned val); 
        void setStepsPerOutput2D(unsigned val);
        void setStepsPerOutputTrack(unsigned val);
        void setMaxReports(unsigned val);
        void setMaxReports1D(unsigned val);
        void setMaxReports2D(unsigned val);

        void addReport(std::unique_ptr<pw::ReportData1D> def){m_report_handler.addReport(std::move(def));}
        void addReport(std::unique_ptr<pw::ReportData2D> def){m_report_handler.addReport(std::move(def));}
        void addReport(std::unique_ptr<pw::TrackData> def) {m_report_handler.addReport(std::move(def));}

        std::filesystem::path dirPath() const {return m_dir_path;}
        void report1D(double t);
        void report2D(double t);
        void reportTrack(double t);
        void report(double t);
        void reportStats() const;
        bool stepUpdate(double t);

    private:
        ReportHandler m_report_handler;
        bool readyForReport() const;
        bool maxReportReached() const;
  
        std::filesystem::path m_dir_path;
        unsigned m_steps_taken{0};
        unsigned m_report_count1D{0};
        unsigned m_steps_per_out1D{1};
        unsigned m_max_reports1D{DEFAULT_MAX_REPORTS_1D};
  
        unsigned m_report_count2D{0};
        unsigned m_steps_per_out2D{1};
        unsigned m_max_reports2D{DEFAULT_MAX_REPORTS_2D};
  
        unsigned m_steps_per_track{1};
  
        bool m_log_progress{false};
        unsigned m_log_freq{1};
        pw::StatCenter m_stat;
};

class PropagatorCV : public BasePropagator
{
    public:
        using BasePropagator::BasePropagator;
        virtual std::vector<dcmplx>& propagator() = 0;
};

}