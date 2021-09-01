

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>
#include "spida/shape/shapeT.h"
#include "spida/helper/constants.h"
#include <boost/math/special_functions/airy.hpp>
#include <boost/math/special_functions/bessel.hpp>

namespace spida{

ShapeT::ShapeT(const GridT& grid,double A,double tp) :
            Shape(grid),
            m_t(grid.getT()),
            m_A(A),m_tp(tp),
            m_offset(0.0),
            m_chirp(0.0),
            m_slow_phase(0.0),
            m_omega0(0.0) {}

std::vector<dcmplx> ShapeT::shapeCV() const
{
    std::vector<dcmplx> v;
    shapeCV(v);
    return v;
}

std::vector<double> ShapeT::shapeRV() const
{
    std::vector<double> v;
    shapeRV(v);
    return v;
}

std::vector<dcmplx> ShapeT::envelope() const
{
    std::vector<dcmplx> v;
    envelope(v);
    return v;
}

void ShapeT::shapeCV(std::vector<dcmplx>& v) const
{
    v.clear();
    v.resize(m_t.size());
    for(auto i = 0; i < m_t.size(); i++)
        v[i] = computeShape(m_t[i]);
}

void ShapeT::shapeRV(std::vector<double>& v) const
{
    v.clear();
    v.resize(m_t.size());
    for(auto i = 0; i < m_t.size(); i++)
        v[i] = computeShapeReal(m_t[i]);
}

void ShapeT::envelope(std::vector<dcmplx>& v) const
{
    v.clear();
    v.resize(m_t.size());
    for(auto i = 0; i < m_t.size(); i++)
        v[i] = computeEnvelope(m_t[i]);
}

dcmplx ShapeT::slowPhaseFactor(double t) const {
    return exp(-ii*m_chirp*pow((t-m_offset),2)+ii*m_slow_phase); 
}

dcmplx ShapeT::fastPhaseFactor(double t) const {
    return exp(-ii*m_omega0*(t-m_offset));
}

double GaussT::compute(double t) const
{
    return ShapeT::amplitude()*exp(-pow((t-ShapeT::offset())/ShapeT::width(),2));
}

double SechT::compute(double t) const
{
    return ShapeT::amplitude()*(1.0/cosh((t-ShapeT::offset()\
                        )/ShapeT::width()));
}

double SuperGaussT::compute(double t) const {
    return ShapeT::amplitude()*exp(-pow((t-ShapeT::offset())/ShapeT::width(),2*m_M));
}

double AiryT::compute(double t) const
{
	double airy = boost::math::airy_ai<double>((t-ShapeT::offset())/ShapeT::width());
    double apodization = exp(-pow(m_apod*(t-ShapeT::offset())/ShapeT::width(),2));
    return ShapeT::amplitude()*airy*apodization;
}

BesselT::BesselT(const GridT& grid,double A,double tp,double apod) : 
    ShapeT(grid,A,tp), 
    m_apod(apod), 
    m_j1(boost::math::cyl_bessel_j_zero<double>(0,1)) {}


double  BesselT::compute(double t) const
{
    double bessel = boost::math::cyl_bessel_j<double>(0,m_j1*fabs(t-ShapeT::offset())/ShapeT::width());
    double apodization = exp(-pow(m_apod*(t-ShapeT::offset())/ShapeT::width(),2));
    return ShapeT::amplitude()*bessel*apodization;
}




}



