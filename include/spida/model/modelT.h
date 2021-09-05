
#ifndef MODELT_H_
#define MODELT_H_

#include <memory>
#include "spida/model/model.h"
#include "spida/SpidaBLT.h"
#include "spida/grid/uniformT.h"

namespace spida{

class ModelBLT : public ModelCV
{
  public:
      ModelBLT(const UniformGridT& grid) :
          ModelCV(1), m_spi(std::make_unique<SpidaBLT>(grid)) {}
      virtual ~ModelBLT() {}
      virtual const std::vector<dcmplx>& linOp() = 0;
      virtual void nonLinResponse(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) = 0;
      SpidaBLT& spida() {return *m_spi;}
      virtual Dimension dimension() {return Dimension::D1;}
  private:
      std::unique_ptr<spida::SpidaBLT> m_spi;
};

}


#endif




