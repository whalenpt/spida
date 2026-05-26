
#include <cassert>
#include "spida/grid/uniformCVX.h"
#include "spida/helper/constants.h"

namespace spida{

UniformGridCVX::UniformGridCVX(unsigned nx,double minX,double maxX) :
    UniformGridX(nx,minX,maxX),
    m_sx(nx)
{
    double L = maxX - minX;
    double dsx = 2.0*PI/L;
    for(unsigned i = 0; i <= nx/2; i++)
        this->m_sx[i] = i*dsx;
    for(unsigned i = nx/2+1; i < nx; i++)
        this->m_sx[i] = -static_cast<double>(nx-i)*dsx;
}

std::vector<double> UniformGridCVX::freqshift(const std::vector<double>& in) const
{
    auto nx = this->getNx();
    assert(nx >= 2);
    std::vector<double> out(in.size());
    for(unsigned i = nx/2+1; i < nx; i++)
        out[i-(nx/2+1)] = in[i];
    for(unsigned i = 0; i <= nx/2; i++)
        out[i+(nx/2-1)] = in[i];
    return out;
}

std::vector<dcmplx> UniformGridCVX::freqshift(const std::vector<dcmplx>& in) const
{
    auto sz = static_cast<unsigned>(in.size());
    assert(sz >= 2);
    std::vector<dcmplx> out(in.size());
    for(unsigned i = sz/2+1; i < sz; i++)
        out[i-(sz/2+1)] = in[i];
    for(unsigned i = 0; i <= sz/2; i++)
        out[i+(sz/2-1)] = in[i];
    return out;
}

void UniformGridCVX::freqshift(const std::vector<double>& in,std::vector<double>& out) const
{
    auto nx = this->getNx();
    assert(nx >= 2);
    for(unsigned i = nx/2+1; i < nx; i++)
        out[i-(nx/2+1)] = in[i];
    for(unsigned i = 0; i <= nx/2; i++)
        out[i+(nx/2-1)] = in[i];
}

void UniformGridCVX::freqshift(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) const
{
    auto nsx = this->getNsx();
    assert(nsx >= 2);
    for(unsigned i = nsx/2+1; i < nsx; i++)
        out[i-(nsx/2+1)] = in[i];
    for(unsigned i = 0; i <= nsx/2; i++)
        out[i+(nsx/2-1)] = in[i];
}

}
