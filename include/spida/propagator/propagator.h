
#ifndef PROPAGATOR_H_
#define PROPAGATOR_H_

#include <string>
#include <vector>
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
      virtual dcmplx* propagator() = 0;
  private:
      virtual void initReport() = 0;
};

}


#endif






