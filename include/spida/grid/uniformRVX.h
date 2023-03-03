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
        ~UniformGridRVX() final = default;
        const std::vector<double>& getSX() const final {return m_sx;}
        double getMinSX() const final {return 0.0;}
        double getMaxSX() const final {return getNx()*spida::PI/getLX();}
        unsigned getNsx() const final {return getNx()/2+1;}
    private:
        std::vector<double> m_sx;
};

}