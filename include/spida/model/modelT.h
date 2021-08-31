
#ifndef MODELT_H_
#define MODELT_H_

#include "spida/model/model.h"

namespace spida{

class PeriodicBLT;
class UniformGridT;

class ModelPeriodicBLT : public ModelCV
{
  public:
      ModelPeriodicBLT(const UniformGridT& grid); 
      virtual ~ModelPeriodicBLT() {}
      virtual const std::vector<dcmplx>& linOp() = 0;
      virtual void nonLinResponse(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) = 0;
      //spida::PeriodicBLT& spida() {return *m_spi;}
      spida::PeriodicBLT& spida();
      virtual Dimension dimension() {return Dimension::D1;}
  private:
      std::unique_ptr<spida::PeriodicBLT> m_spi;
};

}


#endif




