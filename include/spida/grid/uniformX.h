
#ifndef SPIDA_GRID_UNIFORMX_H_
#define SPIDA_GRID_UNIFORMX_H_

#include <vector>
#include "spida/grid/gridX.h"
#include "spida/helper/constants.h"

namespace spida{

void setUniformX(double minX,double maxX,std::vector<double>& x);
void setUniformSX(double minSX,double maxSX,std::vector<double>& sx);

class UniformGridX : public GridX
{
    public:
        UniformGridX(int nx,double min,double max);
        ~UniformGridX() {}; 
        const std::vector<double>& getX() const {return m_x;}
        const std::vector<double>& getSX() const {return m_sx;}
        double getL() const {return GridX::getMaxX()-GridX::getMinX();}
        double getMinSX() const {return -getNx()*spida::PI/getL();} 
        double getMaxSX() const {return getNx()*spida::PI/getL();}
    private:
        std::vector<double> m_x;
        std::vector<double> m_sx;
};

}

#endif


