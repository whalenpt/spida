
#ifndef MODELT_H_
#define MODELT_H_

#include "spida/model/model.h"
#include "spida/spidaT.h"
#include "spida/grid/uniformT.h"
#include <memory>

namespace spida{

class ModelT : public ModelDC
{
  public:
      ModelT(int nt,double minT,double maxT,double minST,double maxST,int c_nthreads = 1) :
          ModelDC(c_nthreads), m_spi(new PeriodicT(UniformGridT(nt,minT,maxT,minST,maxST))) {}
      ModelT(const UniformGridT& grid,int c_nthreads = 1) :
          ModelDC(c_nthreads), m_spi(new PeriodicT(grid)) {}
      virtual ~ModelT();
      virtual const std::vector<dcmplx>& linOperator() = 0;
      virtual void nonLinResponse(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) = 0;
      int numDimensions() const {return 1;}
      int specDimSize() const {return m_spi->getNst();}
      spida::PeriodicT& spida() {return *m_spi;}
  private:
      std::unique_ptr<spida::PeriodicT> m_spi;
};

}


#endif




