#include <cmath>
#include <vector>
#include "spida/shape/shapeR.h"

namespace spida{

ShapeR::ShapeR(const GridR& grid,double A,double w0) :
            Shape(grid),
            m_r(grid.getR()),
            m_A(A),
            m_w0(w0) {}

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

std::vector<dcmplx> ShapeR::shapeCV() const
{
    std::vector<dcmplx> v(m_r.size());
    for(size_t i = 0; i < m_r.size(); i++)
        v[i] = m_A*this->compute(m_r[i]/m_w0);
    if(m_bool_focus){
        for(size_t i = 0; i < m_r.size(); i++){
            v[i] = v[i]*this->focusPhaseFactor(m_r[i]);
        }
    }
    return v;
}

std::vector<double> ShapeR::shapeRV() const
{
    auto cv = this->shapeCV();
    std::vector<double> v(m_r.size());
    for(size_t i = 0; i < m_r.size(); i++)
        v[i] = cv[i].real();
    return v;
}

double GaussR::compute(double r) const
{
    return exp(-pow(r,2));
}

}