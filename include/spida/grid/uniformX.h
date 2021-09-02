
#ifndef SPIDA_GRID_UNIFORMX_H_
#define SPIDA_GRID_UNIFORMX_H_

#include <vector>
#include "spida/grid/gridX.h"
#include "spida/helper/constants.h"

namespace spida{

std::vector<double> buildUniformX(unsigned int nx,double minX,double maxX);
std::vector<double> buildUniformSX(unsigned int nx,double minX,double maxX);

class UniformGridX : public GridX
{
    public:
        UniformGridX(unsigned int nx,double min,double max);
        ~UniformGridX() {}; 
        const std::vector<double>& getX() const {return m_x;}
        const std::vector<double>& getSX() const {return m_sx;}
        double getMinSX() const {return -getNx()*spida::PI/getLX();} 
        double getMaxSX() const {return getNx()*spida::PI/getLX();}
    private:
        std::vector<double> m_x;
        std::vector<double> m_sx;
};

}

#endif


