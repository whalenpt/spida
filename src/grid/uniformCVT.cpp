#include <stdexcept>
#include <string>
#include "spida/grid/uniformCVT.h"
#include "spida/helper/constants.h"

namespace spida{

UniformGridCVT::UniformGridCVT(unsigned nt,double minT,double maxT) : 
    UniformGridT(nt,minT,maxT),
    m_st(nt)
{
    if(nt % 2 != 0)
        throw std::invalid_argument("UniformGridCVT(nt,minT,maxT) error: nt must be divisible by 2");
    double L = maxT - minT;
    double dst = 2.0*PI/L;
    for(auto i = 0; i <= nt/2; i++)
        m_st[i] = i*dst;
    for(auto i = nt/2+1; i < nt; i++)
        m_st[i] = -static_cast<double>(nt-i)*dst;
}

std::vector<double> UniformGridCVT::freqshift(const std::vector<double>& in) const
{
    auto sz = in.size();
    std::vector<double> out(sz);
    for(auto i = sz/2+1; i < sz; i++)
        out[i-(sz/2+1)] = in[i];
    for(auto i = 0; i <= sz/2; i++)
        out[i+(sz/2-1)] = in[i];
    return out;
}

std::vector<dcmplx> UniformGridCVT::freqshift(const std::vector<dcmplx>& in) const
{
    auto sz = in.size();
    std::vector<dcmplx> out(sz);
    for(auto i = sz/2+1; i < sz; i++)
        out[i-(sz/2+1)] = in[i];
    for(auto i = 0; i <= sz/2; i++)
        out[i+(sz/2-1)] = in[i];
    return out;
}

void UniformGridCVT::freqshift(const std::vector<double>& in,std::vector<double>& out) const
{
    auto sz = in.size();
    for(auto i = sz/2+1; i < sz; i++)
        out[i-(sz/2+1)] = in[i];
    for(auto i = 0; i <= sz/2; i++)
        out[i+(sz/2-1)] = in[i];
}

void UniformGridCVT::freqshift(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) const
{
    auto sz = in.size();
    for(auto i = sz/2+1; i < sz; i++)
        out[i-(sz/2+1)] = in[i];
    for(auto i = 0; i <= sz/2; i++)
        out[i+(sz/2-1)] = in[i];
}

}