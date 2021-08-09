
#ifndef PROPAGATORS_H_
#define PROPAGATORS_H_

#include "spida/report/reporthandler.h"
#include "spida/constants.h"
#include <string>
#include <vector>


namespace spida{

class Propagators
{
  public:
      Propagators() {} 
      virtual ~Propagators() {}; 
      virtual void updateFields(double t) = 0;
      ReportHandler& reportHandler() {return m_report_handler;}
  private:
      virtual void initReport() = 0;
      ReportHandler m_report_handler;
};

class PropagatorsDC : public Propagators
{
  public:
      PropagatorsDC() {} 
      virtual ~PropagatorsDC() {}; 
      virtual std::vector<double>& real() = 0;
      virtual std::vector<dcmplx>& spec() = 0;
      virtual void updateFields(double t) = 0;
  private:
      virtual void initReport() = 0;
};

}


#endif






