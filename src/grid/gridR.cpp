
#include "spida/grid/gridR.h"

#include <vector>

namespace spida {

namespace {
template <typename T> void mirrorImpl(unsigned nr, const T* in, T* out, bool signReverse)
{
    for (unsigned i = 0; i < nr; i++)
        out[i] = signReverse ? -in[nr - 1 - i] : in[nr - 1 - i];
    for (unsigned i = 0; i < nr; i++)
        out[i + nr] = in[i];
}
} // namespace

std::vector<double> GridR::mirrorGrid(const std::vector<double>& in, bool signReverse) const
{
    std::vector<double> out(2 * in.size());
    this->mirrorGrid(in, out, signReverse);
    return out;
}

std::vector<dcmplx> GridR::mirrorGrid(const std::vector<dcmplx>& in, bool signReverse) const
{
    std::vector<dcmplx> out(2 * in.size());
    this->mirrorGrid(in, out, signReverse);
    return out;
}

void GridR::mirrorGrid(const std::vector<double>& in,
                       std::vector<double>& out,
                       bool signReverse) const
{
    if (out.size() != 2 * in.size())
        out.resize(2 * in.size());
    this->mirrorGrid(in.data(), out.data(), signReverse);
}

void GridR::mirrorGrid(const std::vector<dcmplx>& in,
                       std::vector<dcmplx>& out,
                       bool signReverse) const
{
    if (out.size() != 2 * in.size())
        out.resize(2 * in.size());
    this->mirrorGrid(in.data(), out.data(), signReverse);
}

void GridR::mirrorGrid(const double* in, double* out, bool signReverse) const
{
    mirrorImpl(this->m_nr, in, out, signReverse);
}

void GridR::mirrorGrid(const dcmplx* in, dcmplx* out, bool signReverse) const
{
    mirrorImpl(this->m_nr, in, out, signReverse);
}

} // namespace spida
