
#include "spida/grid/besselR.h"

#include <cmath>
#include <stdexcept>
#include <string>

#include <boost/math/special_functions/bessel.hpp>

namespace spida {

BesselRootGridR::BesselRootGridR(unsigned nr, double maxr) : GridR(nr, maxr), m_r(nr), m_sr(nr)
{
    boost::math::cyl_bessel_j_zero<double>(
        0.0, 1, static_cast<int>(nr), std::back_inserter(m_roots));
    this->m_jN = boost::math::cyl_bessel_j_zero<double>(0.0, static_cast<int>(nr) + 1);
    for (unsigned i = 0; i < nr; i++)
        this->m_r[i] = this->m_roots[i] * GridR::getMaxR() / this->m_jN;
    for (unsigned i = 0; i < nr; i++)
        this->m_sr[i] = this->m_roots[i] * this->getMaxSR() / this->m_jN;
}

} // namespace spida