
#include <cmath>
#include "spida/shape/shapeX.h"

namespace spida{

void ShapeX::shape(std::vector<double>& v) const
{
    v.clear();
    v.resize(m_grid.getNx());
    const std::vector<double>& x = m_grid.getX();
    for(auto i = 0; i < m_grid.getNx(); i++)
        v[i] = compute(x[i]);
}

double ExpX::compute(double x) const
{
    return ShapeX::amplitude()*exp(ShapeX::width()*(x-ShapeX::offset())+ShapeX::phase());
}

double CosX::compute(double x) const
{
    return ShapeX::amplitude()*cos(ShapeX::width()*(x-ShapeX::offset())+ShapeX::phase());
}

double SinX::compute(double x) const
{
    return ShapeX::amplitude()*sin(ShapeX::width()*(x-ShapeX::offset())+ShapeX::phase());
}



}








