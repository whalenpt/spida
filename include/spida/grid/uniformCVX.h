// uniformCVX.h
#pragma once

#include <vector>
#include "spida/grid/uniformX.h"
#include "spida/helper/constants.h"

namespace spida{

class UniformGridCVX : public UniformGridX
{
    public:
        UniformGridCVX(unsigned nx,double min,double max);
        ~UniformGridCVX() {}; 
        const std::vector<double>& getSX() const {return m_sx;}
        double getMinSX() const {return -getNx()*spida::PI/getLX();} 
        double getMaxSX() const {return getNx()*spida::PI/getLX();}
        std::vector<double> freqshift(const std::vector<double>& in) const;
        std::vector<dcmplx> freqshift(const std::vector<dcmplx>& in) const;
        void freqshift(const std::vector<double>& in,std::vector<double>& out) const;
        void freqshift(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) const;
    private:
        std::vector<double> m_sx;
};


}



