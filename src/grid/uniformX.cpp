#include <stdexcept>
#include <string>
#include "spida/grid/uniformX.h"

namespace spida{

UniformGridX::UniformGridX(unsigned nx,double minX,double maxX) : 
    GridX(nx,minX,maxX),
    m_x(nx)
{
    if(minX >= maxX){
        std::string msg = "Error in setUniformX: minX must be less than maxX.";
        throw std::invalid_argument(msg);
    }
    double dx = (maxX - minX)/static_cast<double>(nx);
    for(auto i = 0; i < nx; i++) m_x[i] = minX + i*dx; 
}

}