
#include <cmath>
#include "spida/shape/shapeX.h"


namespace spida{

ShapeX::ShapeX(const GridX& grid,double A,double w0) :
            Shape(grid),
            m_x(grid.getX()),
            m_A(A),
            m_w0(w0) {}


std::vector<dcmplx> ShapeX::shapeCV() const
{
    std::vector<dcmplx> v(m_x.size());
    for(size_t i = 0; i < m_x.size(); i++)
        v[i] = m_A*compute((m_x[i]-m_offset)/m_w0)*exp(ii*m_phi0);
    return v;
}

std::vector<double> ShapeX::shapeRV() const
{
    std::vector<double> v(m_x.size());
    for(size_t i = 0; i < m_x.size(); i++)
        v[i] = m_A*compute((m_x[i]-m_offset)/m_w0)*cos(m_phi0);
    return v;
}

double ExpX::compute(double x) const
{
    return exp(x);
}

double CosX::compute(double x) const
{
    return cos(x);
}

double SinX::compute(double x) const
{
    return sin(x);
}

}