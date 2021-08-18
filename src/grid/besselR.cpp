
#include <cmath>
#include <string>
#include <stdexcept>
#include <boost/math/special_functions/bessel.hpp>
#include "spida/grid/besselR.h"

namespace spida{


BesselRootGridR::BesselRootGridR(int nr,double maxr) : GridR(nr,maxr),
    m_r(nr), m_sr(nr)
{
    //OutputIterator cyl_bessel_j_zero(
    //                 T v,                       // Floating-point value for Jv.
    //                 int start_index,           // 1-based index of first zero.
    //                 unsigned number_of_roots,  // How many roots to generate.
    //                 OutputIterator out_it);

    // Want J0 -> v = 0, starting with first root -> start_index=1
    boost::math::cyl_bessel_j_zero<double>(0,1,nr,std::back_inserter(m_roots));
    // 1-based index of zero (use nr+1 for m_jN rather than nr)
    m_jN = boost::math::cyl_bessel_j_zero<double>(0,nr+1);
    // Set physical grid
    for(auto i = 0; i < nr; i++)
        m_r[i] = m_roots[i]*GridR::getMaxR()/m_jN;
    // Set spectral grid
    for(auto i = 0; i < nr; i++)
        m_sr[i] = m_roots[i]/GridR::getMaxR();
}


}






