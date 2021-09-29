// uniformRVX.h
#pragma once

#include <vector>
#include "spida/grid/uniformX.h"
#include "spida/helper/constants.h"

namespace spida{

class UniformGridRVX : public UniformGridX
{
    public:
        UniformGridRVX(unsigned nx,double min,double max);
        ~UniformGridRVX() {}; 
        const std::vector<double>& getSX() const {return m_sx;}
        double getMinSX() const {return 0.0;}
        double getMaxSX() const {return getNx()*spida::PI/getLX();}
        unsigned getNsx() const {return getNx()/2+1;}
    private:
        std::vector<double> m_sx;
};



}



