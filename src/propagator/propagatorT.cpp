

#include "spida/propagator/propagatorT.h"
#include "spida/model/modelT.h"
#include "spida/shape/shapeT.h"

namespace spida{

void PropagatorsT::initFields(const spida::ShapeT& shape)
{
    const std::vector<double>& t = m_spi.getT();
    shape.computeReal(t,m_U);
    m_spi.T_To_ST(m_U,m_V);
}


}

