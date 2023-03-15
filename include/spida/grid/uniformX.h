// uniformX.h
#pragma once

#include <vector>
#include "spida/grid/gridX.h"

namespace spida{

class UniformGridX : public GridX
{
    public:
        UniformGridX(unsigned nx,double min,double max);
        ~UniformGridX() override = default; 
        const std::vector<double>& getX() const override {return m_x;}
    private:
        std::vector<double> m_x;
};

}