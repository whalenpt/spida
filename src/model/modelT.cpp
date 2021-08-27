
#include "spida/PeriodicBLT.h"
#include "spida/model/modelT.h"
#include "spida/grid/uniformT.h"

namespace spida{

ModelPeriodicBLT::ModelPeriodicBLT(const UniformGridT& grid,int c_nthreads) :
      ModelCV(c_nthreads), 
      m_spi(new PeriodicBLT(grid)) {}


}






