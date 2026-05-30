

#include "spida/grid/uniformRVT.h"

#include <cmath>
#include <format>
#include <stdexcept>
#include <string>

#include <pwutils/pwindexing.hpp>

namespace spida {

std::vector<double> UniformGridRVT::constructGridST(unsigned nt, double minT, double maxT) const
{
    std::vector<double> st(nt / 2 + 1);
    double dst = 2.0 * PI / (maxT - minT);
    for (unsigned i = 0; i <= nt / 2; i++)
        st[i] = i * dst;
    return st;
}

UniformGridRVT::UniformGridRVT(unsigned nt, double minT, double maxT, double minST, double maxST)
    : UniformGridT(nt, minT, maxT)
{
    this->verifyFrequencyRange(minST, maxST);
    std::vector<double> full_st = this->constructGridST(nt, minT, maxT);
    this->m_minI = pw::nearestIndex(full_st, minST);
    this->m_maxI = pw::nearestIndex(full_st, maxST);
    this->m_nst = this->m_maxI - this->m_minI + 1;
    this->m_st.resize(this->m_nst, 0.0);
    for (auto i = this->m_minI; i <= this->m_maxI; i++)
        this->m_st[i - this->m_minI] = full_st[i];
    this->m_minST = this->m_st[0];
    this->m_maxST = this->m_st.back();
}

UniformGridRVT::UniformGridRVT(unsigned nt, double minT, double maxT)
    : UniformGridT(nt, minT, maxT),
      m_st(constructGridST(nt, minT, maxT)),
      m_nst(nt / 2 + 1),
      m_minST(indxToFreq(0)),
      m_maxST(indxToFreq(nt / 2)),
      m_maxI(nt / 2)
{
}

double UniformGridRVT::maxPossibleFreq() const
{
    return this->indxToFreq(this->getNt() / 2);
}

unsigned UniformGridRVT::freqToIndx(double omeg) const
{
    if (omeg < 0.0)
        throw std::domain_error("freqToIndx only processes positive frequencies");
    if (omeg > this->maxPossibleFreq())
        throw std::domain_error(
            std::format("freqToIndx(double omeg) error: The provided frequency of {} "
                        "is bigger than the maximum specified frequency of {}",
                        omeg,
                        this->maxPossibleFreq()));
    double dst = 2.0 * PI / this->getLT();
    return static_cast<unsigned>(std::round(omeg / dst));
}

double UniformGridRVT::indxToFreq(unsigned indx) const
{
    if (indx > this->getNt() / 2)
        throw std::domain_error(std::format(
            "indxToFreq(indx) error: Specified index of {} which is bigger than nt/2 of {}",
            indx,
            this->getNt() / 2));
    double dst = 2.0 * PI / this->getLT();
    return dst * indx;
}

void UniformGridRVT::verifyFrequencyRange(double minST, double maxST) const
{
    if (minST >= maxST)
        throw std::domain_error(std::format(
            "UniformGridRVT minST of {} is greater than or equal to the specified maxST of {}",
            minST,
            maxST));
    if (minST < 0.0)
        throw std::domain_error(std::format("UniformGridRVT minST of {} is less than zero, please "
                                            "increase to a non-negative value.",
                                            minST));
    if (maxST > this->maxPossibleFreq()) {
        auto str2 =
            std::to_string(static_cast<int>((this->getMaxT() - this->getMinT()) * maxST / PI) + 1);
        throw std::domain_error(
            std::format("Temporal grid spacing not small enough to accomodate the maximum "
                        "grid frequency of {}. Increase NT to more than {}",
                        this->maxPossibleFreq(),
                        str2));
    }
}

} // namespace spida
