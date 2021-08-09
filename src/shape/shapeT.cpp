

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>
#include "spida/shape/shapeT.h"
#include "spida/constants.h"
#include <boost/math/special_functions/airy.hpp>
#include <boost/math/special_functions/bessel.hpp>

namespace spida{

dcmplx ShapeT::computePhaseFactor(double t) const {
    return exp(-ii*m_omega0*(t-m_offset)-ii*m_chirp*pow((t-m_offset),2)+ii*m_slow_phase); 
}


void ShapeT::computeReal(const std::vector<double>& t,std::vector<double>& y) const {
    assert(t.size() == y.size());
    for(auto i = 0; i < t.size(); i++) y[i] = computeReal(t[i]);
}

void ShapeT::envelope(const std::vector<double>& t,std::vector<dcmplx>& y) const
{
    assert(t.size() == y.size());
    for(auto i = 0; i < t.size(); i++) y[i] = envelope(t[i]);
}

dcmplx GaussT::compute(double t) const
{
    return ShapeT::amplitude()*exp(-pow((t-ShapeT::offset())/ShapeT::width(),2))\
        *ShapeT::computePhaseFactor(t);
}

dcmplx SechT::compute(double t) const
{
    return ShapeT::amplitude()*(1.0/cosh((t-ShapeT::offset()\
                        )/ShapeT::width()))*ShapeT::computePhaseFactor(t);
}

dcmplx SuperGaussT::compute(double t) const {
    return ShapeT::amplitude()*exp(-pow((t-ShapeT::offset())/ShapeT::width(),2*m_M))\
        *ShapeT::computePhaseFactor(t);
}

dcmplx AiryT::compute(double t) const
{
	double airy = boost::math::airy_ai<double>((t-ShapeT::offset())/ShapeT::width());
    double apodization = exp(-pow(m_apod*(t-ShapeT::offset())/ShapeT::width(),2));
    return ShapeT::amplitude()*airy*apodization*ShapeT::computePhaseFactor(t);
}

BesselT::BesselT(double A,double tp,double omega0) : ShapeT(A,tp,omega0),
                m_j1(boost::math::cyl_bessel_j_zero<double>(0,1)),
                m_apod(0.25) {}

dcmplx BesselT::compute(double t) const
{
    double bessel = boost::math::cyl_bessel_j<double>(0,m_j1*fabs(t-ShapeT::offset())/ShapeT::width());
    double apodization = exp(-pow(m_apod*(t-ShapeT::offset())/ShapeT::width(),2));
    return ShapeT::amplitude()*bessel*apodization*ShapeT::computePhaseFactor(t);
}




}



