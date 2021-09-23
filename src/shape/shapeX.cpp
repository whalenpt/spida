
#include <cmath>
#include "spida/shape/shapeX.h"


namespace spida{

ShapeX::ShapeX(const GridX& grid,double A,double w0) :
            Shape(grid),
            m_x(grid.getX()),
            m_A(A),
            m_w0(w0),
            m_offset(0.0),
            m_phi0(0.0)
    {}


dcmplx ShapeX::shapeCV(double x) const {
    return m_A*compute(x)*phaseFactor(r);
}

std::vector<dcmplx> ShapeR::shapeCV() const
{
    std::vector<dcmplx> v;
    shapeCV(v);
    return v;
}

std::vector<double> ShapeR::shapeRV() const
{
    std::vector<double> v;
    shapeRV(v);
    return v;
}

void ShapeR::shapeCV(std::vector<dcmplx>& v) const
{
    v.clear();
    v.resize(m_r.size());
    for(auto i = 0; i < m_r.size(); i++)
        v[i] = shapeCV(m_r[i]);
}

void ShapeR::shapeRV(std::vector<double>& v) const
{
    v.clear();
    v.resize(m_r.size());
    for(auto i = 0; i < m_r.size(); i++)
        v[i] = shapeRV(m_r[i]);
}

double GaussR::compute(double r) const
{
    return exp(-pow(r/ShapeR::width(),2));
}

}









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








