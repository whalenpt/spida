
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
    double dsx = 2.0*spida::PI/(maxX - minX);
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

}

