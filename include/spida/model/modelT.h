
#ifndef MODELT_H_
#define MODELT_H_

#include "spida/model/model.h"
#include "spida/PeriodicBLT.h"
#include "spida/grid/uniformT.h"

namespace spida{

class ModelPeriodicBLT : public ModelCV
{
  public:
      ModelPeriodicBLT(const UniformGridT& grid); 
      virtual ~ModelPeriodicBLT();
      virtual const std::vector<dcmplx>& linOp() = 0;
      virtual void nonLinResponse(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) = 0;
      PeriodicBLT& spida();
      virtual Dimension dimension() {return Dimension::D1;}
  private:
      std::unique_ptr<spida::PeriodicBLT> m_spi;
};

}


#endif




