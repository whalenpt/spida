
#include <cmath>
#include "spida/shape/shapeR.h"

namespace spida{

ShapeR::ShapeR(const GridR& grid,double A,double w0) :
            m_r(grid.getR()),
            m_A(A),
            m_w0(w0)
    {}

void ShapeR::shape(std::vector<double>& v) const
{
    v.clear();
    v.resize(m_r.size());
    for(auto i = 0; i < m_r.size(); i++)
        v[i] = compute(m_r[i]);
}

double GaussR::compute(double r) const
{
    return ShapeR::amplitude()*exp(-pow(r/ShapeR::width(),2));
}

}








