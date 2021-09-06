#ifndef MODELRT_H_
#define MODELRT_H_

#include <memory>
#include "spida/model/model.h"
#include "spida/SpidaRBLT.h"
#include "spida/grid/uniformT.h"
#include "spida/grid/besselR.h"

namespace spida{

class ModelRBLT : public ModelCV
{
  public:
      ModelRBLT(const BesselRootGridR& gridR,const UniformGridT& gridT,unsigned int threads=1) :
          ModelCV(threads), m_spi(std::make_unique<SpidaRBLT>(gridR,gridT,threads)) {}
      virtual ~ModelRBLT() {}
      virtual const std::vector<dcmplx>& linOp() = 0;
      virtual void nonLinResponse(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) = 0;
      SpidaRBLT& spida() {return *m_spi;}
      virtual Dimension dimension() {return Dimension::D2;}
  private:
      std::unique_ptr<spida::SpidaRBLT> m_spi;
};

}


#endif




