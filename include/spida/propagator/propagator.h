
#ifndef PROPAGATOR_H_
#define PROPAGATOR_H_

#include <string>
#include <vector>
#include <memory>
#include <pwutils/report/basedata.hpp>
#include "spida/helper/constants.h"
#include "spida/report/reporthandler.h"

namespace spida{

class Propagator
{
  public:
      Propagator() {} 
      virtual ~Propagator() {}; 
      virtual void updateFields(double t) = 0;
      ReportHandler& reportHandler() {return m_report_handler;}
      void addReport(std::unique_ptr<pw::ReportData1D> def){m_report_handler.addReport(std::move(def));}
      void addReport(std::unique_ptr<pw::ReportData2D> def){m_report_handler.addReport(std::move(def));}
      void addReport(std::unique_ptr<pw::TrackData> def) {m_report_handler.addReport(std::move(def));}
  private:
      virtual void initReport() = 0;
      ReportHandler m_report_handler;
};


class PropagatorCV : public Propagator
{
  public:
      PropagatorCV() {} 
      virtual ~PropagatorCV() {}; 
      virtual void updateFields(double t) = 0;
      virtual std::vector<dcmplx>& propagator() = 0;
  private:
      virtual void initReport() = 0;
};

}


#endif






