// uniformX.h
#pragma once

#include <vector>
#include "spida/grid/gridX.h"

namespace spida{


class UniformGridX : public GridX
{
    public:
        UniformGridX(unsigned nx,double min,double max);
        ~UniformGridX() {}; 
        const std::vector<double>& getX() const {return m_x;}
        virtual const std::vector<double>& getSX() const = 0;
        virtual double getMinSX() const = 0;
        virtual double getMaxSX() const = 0;
    private:
        std::vector<double> m_x;
};

}



