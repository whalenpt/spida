

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include "spida/helper/constants.h"
#include "spida/shape/shapeT.h"
#include <boost/math/special_functions/airy.hpp>
#include <boost/math/special_functions/bessel.hpp>

namespace spida{


std::vector<dcmplx> ShapeT::shapeCV() const
{
    std::vector<dcmplx> v(m_t.size());
    auto v_env = this->envelope();
    for(size_t i = 0; i < m_t.size(); i++)
        v[i] = v_env[i]*fastPhaseFactor(m_t[i]-m_offset);
    return v;
}

std::vector<double> ShapeT::shapeRV() const
{
    std::vector<double> v(m_t.size());
    auto v_env = this->envelope();
    for(size_t i = 0; i < m_t.size(); i++)
        v[i] = (v_env[i]*fastPhaseFactor(m_t[i]-m_offset)).real();
    return v;
}

std::vector<dcmplx> ShapeT::envelope() const
{
    std::vector<dcmplx> v(m_t.size());
    for(size_t i = 0; i < m_t.size(); i++)
        v[i] = m_A*compute((m_t[i]-m_offset)/m_tp)*slowPhaseFactor(m_t[i]-m_offset);
    return v;
}

dcmplx ShapeT::slowPhaseFactor(double t) const {
    return exp(-ii*m_chirp*pow(t,2)+ii*m_slow_phase); 
}

dcmplx ShapeT::fastPhaseFactor(double t) const {
    return exp(-ii*m_omega0*t);
}

double GaussT::compute(double t) const
{
    return exp(-pow(t,2));
}

double SechT::compute(double t) const
{
    return 1.0/cosh(t);
}

double SuperGaussT::compute(double t) const {
    return exp(-pow(t,2*m_M));
}

double AiryT::compute(double t) const
{
	double airy = boost::math::airy_ai<double>(t);
    double apodization = exp(-pow(m_apod*t,2));
    return airy*apodization;
}

double  BesselT::compute(double t) const
{
    double bessel = boost::math::cyl_bessel_j<double>(0,m_j1*fabs(t));
    double apodization = exp(-pow(m_apod*t,2));
    return bessel*apodization;
}

}