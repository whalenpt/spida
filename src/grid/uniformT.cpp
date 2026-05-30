#include "spida/grid/uniformT.h"

#include "spida/grid/gridT.h"

#include <stdexcept>
#include <string>

namespace spida {

UniformGridT::UniformGridT(unsigned nt, double minT, double maxT) : GridT(nt, minT, maxT), m_t(nt)
{
    if (nt < 2)
        throw std::invalid_argument("UniformGridT(nt,minT,maxT) error: nt must be >= 2.");
    if (minT >= maxT) {
        std::string msg = "UniformGridT(nt,minT,maxT) error: minT must be less than maxT.";
        throw std::invalid_argument(msg);
    }
    double dt = (maxT - minT) / static_cast<double>(nt);
    for (unsigned i = 0; i < nt; i++)
        m_t[i] = minT + i * dt;
}

} // namespace spida