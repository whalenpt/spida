
#include <cmath>
#include <string>
#include <stdexcept>
#include "spida/grid/uniformX.h"
#include "spida/helper/constants.h"
#include <iostream>

namespace spida{

std::vector<double> buildUniformX(unsigned nx,double minX,double maxX)
{
    if(minX >= maxX){
        std::string msg = "Error in setUniformX: minX must be less than maxX.";
        throw std::invalid_argument(msg);
    }
    std::vector<double> x(nx);
    double dx = (maxX - minX)/nx;
    for(int i = 0; i < nx; i++) x[i] = minX + i*dx; 
    return x;
}

std::vector<double> buildUniformSX(unsigned nx,double minX,double maxX)
{
    if(minX >= maxX){
        std::string msg = "Error in setUniformSX: minX must be less than maxX.";
        throw std::invalid_argument(msg);
    }
    std::vector<double> sx(nx);
    double L = maxX - minX;
    double dsx = 2.0*PI/L;

    for(auto i = 0; i <= nx/2; i++)
        sx[i] = i*dsx;
    for(auto i = nx/2+1; i < nx; i++)
        sx[i] = -static_cast<double>(nx-i)*dsx;
    return sx;
}


UniformGridX::UniformGridX(unsigned nx,double minX,double maxX) : GridX(nx,minX,maxX),
    m_x(nx), m_sx(nx)
{
    m_x = buildUniformX(nx,minX,maxX);
    m_sx = buildUniformSX(nx,minX,maxX);
}

std::vector<double> UniformGridX::freqshift(const std::vector<double>& in) const
{
    std::vector<double> out(in.size());
    unsigned nx = getNx();
    for(auto i = nx/2+1; i < nx; i++)
        out[i-(nx/2+1)] = in[i];
    for(auto i = 0; i <= nx/2; i++)
        out[i+(nx/2-1)] = in[i];
    return out;
}

std::vector<dcmplx> UniformGridX::freqshift(const std::vector<dcmplx>& in) const
{
    std::vector<dcmplx> out(in.size());
    unsigned sz = in.size();
    for(auto i = sz/2+1; i < sz; i++)
        out[i-(sz/2+1)] = in[i];
    for(auto i = 0; i <= sz/2; i++)
        out[i+(sz/2-1)] = in[i];
    return out;
}

void UniformGridX::freqshift(const std::vector<double>& in,std::vector<double>& out) const
{
    unsigned nx = getNx();
    for(auto i = nx/2+1; i < nx; i++)
        out[i-(nx/2+1)] = in[i];
    for(auto i = 0; i <= nx/2; i++)
        out[i+(nx/2-1)] = in[i];
}

void UniformGridX::freqshift(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) const
{
    unsigned nsx = getNsx();
    for(auto i = nsx/2+1; i < nsx; i++)
        out[i-(nsx/2+1)] = in[i];
    for(auto i = 0; i <= nsx/2; i++)
        out[i+(nsx/2-1)] = in[i];
}


}



