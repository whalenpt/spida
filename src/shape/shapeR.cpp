
#include <cmath>
#include "spida/shape/shapeR.h"

namespace spida{

ShapeR::ShapeR(const GridR& grid,double A,double w0) :
            Shape(grid),
            m_r(grid.getR()),
            m_A(A),
            m_w0(w0),
            m_bool_focus(false),
            m_f(0.0)
    {}

void ShapeR::setFocus(double f)
{
    if(f <= 0.0)
        throw std::domain_error("ShapeR::setFocus error: focus must be a positive value");
    m_f = f;
    m_bool_focus = true;
}

dcmplx ShapeR::focusPhaseFactor(double r) const {
    if(m_bool_focus)
        return exp(-ii*pow(r,2)/m_f);
    return 1.0;
}

dcmplx ShapeR::computeShape(double r) const {
    if(m_bool_focus)
        return m_A*compute(r)*focusPhaseFactor(r);
    return m_A*compute(r);
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
        v[i] = computeShape(m_r[i]);
}

void ShapeR::shapeRV(std::vector<double>& v) const
{
    v.clear();
    v.resize(m_r.size());
    for(auto i = 0; i < m_r.size(); i++)
        v[i] = computeShapeReal(m_r[i]);
}

double GaussR::compute(double r) const
{
    return exp(-pow(r/ShapeR::width(),2));
}

}








