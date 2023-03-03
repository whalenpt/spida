#include "spida/grid/uniformRVX.h"
#include "spida/helper/constants.h"

namespace spida{

UniformGridRVX::UniformGridRVX(unsigned nx,double minX,double maxX) : 
    UniformGridX(nx,minX,maxX),
    m_sx(nx/2+1)
{
    double L = maxX - minX;
    double dsx = 2.0*PI/L;
    for(auto i = 0; i <= nx/2; i++)
        m_sx[i] = i*dsx;
}

}