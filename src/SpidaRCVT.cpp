
#include "spida/grid/besselR.h"
#include "spida/grid/uniformCVT.h"
#include "spida/RCVT.h"
#include "spida/transform/hankelfftRCVT.h"

#include <cmath>
#include <memory>

#include <utils/except.h>

namespace spida {

SpidaRCVT::SpidaRCVT(const BesselRootGridR& gridR, const UniformGridCVT& gridT, unsigned threads)
    : m_gridR(std::make_unique<BesselRootGridR>(gridR)),
      m_gridT(std::make_unique<UniformGridCVT>(gridT)),
      m_tr(std::make_unique<HankelFFTRCVT>(gridR, gridT, threads))
{
}

const BesselRootGridR& SpidaRCVT::getGridR() const
{
    return *m_gridR;
}

const UniformGridCVT& SpidaRCVT::getGridT() const
{
    return *m_gridT;
}

const HankelFFTRCVT& SpidaRCVT::getTransformRT() const
{
    return *m_tr;
}

const std::vector<double>& SpidaRCVT::getR() const
{
    return m_gridR->getR();
}

const std::vector<double>& SpidaRCVT::getSR() const
{
    return m_gridR->getSR();
}

const std::vector<double>& SpidaRCVT::getT() const
{
    return m_gridT->getT();
}

const std::vector<double>& SpidaRCVT::getST() const
{
    return m_gridT->getST();
}

void SpidaRCVT::RT_To_SRST(const std::vector<dcmplx>& in, std::vector<dcmplx>& out)
{
    m_tr->RT_To_SRST(in, out);
}

void SpidaRCVT::SRST_To_RT(const std::vector<dcmplx>& in, std::vector<dcmplx>& out)
{
    m_tr->SRST_To_RT(in, out);
}

unsigned SpidaRCVT::spectralSize() const
{
    return m_gridR->getNr() * m_gridT->getNt();
}

unsigned SpidaRCVT::physicalSize() const
{
    return m_gridR->getNsr() * m_gridT->getNst();
}

void SpidaRCVT::mirrorR(const std::vector<dcmplx>& in, std::vector<dcmplx>& out) const
{
    auto nt = m_gridT->getNt();
    auto nr = m_gridR->getNr();
    for (unsigned j = 0; j < nt; j++) {
        // Reversed half (indices [0, nr)): out[k] = in[nr-1-k] for k = 0..nr-1 --
        // rewritten from a decrementing "for (i = nr-1; i >= 0; i--)" loop, which
        // is an unsigned-underflow bug: i >= 0 is always true for an unsigned
        // type, so i wraps past 0 to UINT_MAX instead of stopping, running
        // in/out far out of bounds (the segfault this fixes). Equivalent by
        // substituting k = (nr-1)-i, verified to produce the identical
        // assignments in the opposite iteration order -- each iteration writes
        // a disjoint output index, so order doesn't affect the result.
        for (unsigned k = 0; k < nr; k++)
            out[k * nt + j] = in[((nr - 1) - k) * nt + j];
        // Original half (indices [nr, 2*nr)): unchanged, already forward.
        for (unsigned i = 0; i < nr; i++)
            out[(i + nr) * nt + j] = in[i * nt + j];
    }
}

} // namespace spida
