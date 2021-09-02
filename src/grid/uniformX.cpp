
#include <cmath>
#include <string>
#include <stdexcept>
#include "spida/grid/uniformX.h"
#include "spida/helper/constants.h"

namespace spida{

std::vector<double> buildUniformX(unsigned int nx,double minX,double maxX)
{
    if(minX >= maxX){
        std::string msg = "Error in setUniformX: minX must be less than maxX.";
        throw std::invalid_argument(msg);
    }
    std::vector<double> x(nx);
    double dx = (maxX - minX)/(nx-1);
    for(int i = 0; i < nx; i++) x[i] = minX + i*dx; 
    return x;
}

std::vector<double> buildUniformSX(unsigned int nx,double minX,double maxX)
{
    if(minX >= maxX){
        std::string msg = "Error in setUniformSX: minX must be less than maxX.";
        throw std::invalid_argument(msg);
    }
    std::vector<double> sx(nx);
    double dsx = 2.0*spida::PI/(maxX - minX);
    for(int k = 0; k < nx; k++) sx[k] = (k-nx/2.0)*dsx; 
    return sx;
}

UniformGridX::UniformGridX(unsigned int nx,double minX,double maxX) : GridX(nx,minX,maxX),
    m_x(nx), m_sx(nx)
{
    m_x = buildUniformX(nx,minX,maxX);
    m_sx = buildUniformSX(nx,minX,maxX);
}

}

