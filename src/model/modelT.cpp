
#include "spida/PeriodicBLT.h"
#include "spida/model/modelT.h"
#include "spida/grid/uniformT.h"

namespace spida{

ModelPeriodicBLT::ModelPeriodicBLT(const UniformGridT& grid) :
      ModelCV(1), 
      m_spi(new PeriodicBLT(grid)) {}


}






