

#include "spida/propagator/propagatorT.h"
#include "spida/model/modelT.h"
#include "spida/shape/shapeT.h"
#include "spida/PeriodicBLT.h"
#include "spida/grid/uniformT.h"

namespace spida{

PropagatorBLT::PropagatorBLT(ModelPeriodicBLT& md) : 
    m_md(md),
    m_spectral_field(md.spida().getGridT().getNst()),
    m_real_field(md.spida().getGridT().getNt()) {}


void PropagatorBLT::initFields(const spida::ShapeT& shape)
{
    shape.shapeReal(m_real_field);
    spida().T_To_ST(m_real_field,m_spectral_field);
}

PeriodicBLT& PropagatorBLT::spida() {return m_md.spida();}

}

