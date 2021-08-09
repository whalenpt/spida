
#include <cmath>
#include <string>
#include <stdexcept>
#include "spida/grid/uniformX.h"
#include "spida/constants.h"

namespace spida{

void setUniformX(double minX,double maxX,std::vector<double>& x)
{
    if(minX >= maxX){
        std::string msg = "Error in setUniformX: minX must be less than maxX.";
        throw std::invalid_argument(msg);
    }
    int nx = x.size();
    double dx = (maxX - minX)/(nx-1);
    for(int i = 0; i < nx; i++) x[i] = minX + i*dx; 
}

void setUniformSX(double minX,double maxX,std::vector<double>& sx)
{
    if(minX >= maxX){
        std::string msg = "Error in setUniformSX: minX must be less than maxX.";
        throw std::invalid_argument(msg);
    }
    int nsx = sx.size();
    double dsx = 2.0*spida::PI/(maxX - minX);
    for(int k = 0; k < nsx; k++) sx[k] = (k-nsx/2.0)*dsx; 
}

UniformGridX::UniformGridX(int nx,double min,double max) : GridX(nx,min,max),
    m_x(nx), m_sx(nx)
{
    setUniformX(min,max,m_x);
    setUniformSX(min,max,m_sx);
}

}

