

#include "spida/propagator/propagatorT.h"
#include "spida/model/modelT.h"
#include "spida/shape/shapeT.h"
#include "spida/PeriodicBLT.h"
#include "spida/grid/uniformT.h"

namespace spida{

PropagatorBLT::PropagatorBLT(ModelPeriodicBLT& md) : 
    m_spi(md.spida()), 
    m_spectral_field(m_spi.getGridT().getNst()),
    m_real_field(m_spi.getGridT().getNt()) {}


void PropagatorBLT::initFields(const spida::ShapeT& shape)
{
    shape.shapeReal(m_real_field);
    m_spi.T_To_ST(m_real_field,m_spectral_field);
}

PeriodicBLT& PropagatorBLT::spida() {return m_spi;}

}

