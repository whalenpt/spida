

#include "spida/propagator/propagatorT.h"
#include "spida/model/modelT.h"
#include "spida/shape/shapeT.h"

namespace spida{

void PropagatorsT::initFields(const spida::ShapeT& shape)
{
    shape.shapeReal(m_U);
    m_spi.T_To_ST(m_U,m_V);
}


}

