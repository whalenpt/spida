#include <stdexcept>
#include <string>
#include "spida/grid/gridT.h"
#include "spida/grid/uniformT.h"

namespace spida{

UniformGridT::UniformGridT(unsigned nt,double minT,double maxT) : 
    GridT(nt,minT,maxT),
    m_t(nt)
{
    if(minT >= maxT){
        std::string msg = "UniformGridT(nt,minT,maxT) error: minT must be less than maxT.";
        throw std::invalid_argument(msg);
    }
    double dt = (maxT - minT)/static_cast<double>(nt);
    for(unsigned i = 0; i < nt; i++) m_t[i] = minT + i*dt; 
}

}