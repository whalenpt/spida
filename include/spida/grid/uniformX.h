// uniformX.h
#pragma once

#include <vector>
#include "spida/grid/gridX.h"
#include "spida/helper/constants.h"

namespace spida{

std::vector<double> buildUniformX(unsigned int nx,double minX,double maxX);
std::vector<double> buildUniformSX(unsigned int nx,double minX,double maxX);

class UniformGridX : public GridX
{
    public:
        UniformGridX(unsigned nx,double min,double max);
        ~UniformGridX() {}; 
        const std::vector<double>& getX() const {return m_x;}
        const std::vector<double>& getSX() const {return m_sx;}
        double getMinSX() const {return -getNx()*spida::PI/getLX();} 
        double getMaxSX() const {return getNx()*spida::PI/getLX();}
        std::vector<double> freqshift(const std::vector<double>& in) const;
        std::vector<dcmplx> freqshift(const std::vector<dcmplx>& in) const;
        void freqshift(const std::vector<double>& in,std::vector<double>& out) const;
        void freqshift(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) const;
    private:
        std::vector<double> m_x;
        std::vector<double> m_sx;
};

}



