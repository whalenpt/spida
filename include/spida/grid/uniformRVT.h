// uniformRVT.h
#pragma once

#include <vector>
#include "spida/grid/uniformT.h"
#include "spida/helper/constants.h"

namespace spida{

class UniformGridRVT : public UniformGridT
{
    public:
        explicit UniformGridRVT(unsigned nt,double minT,double maxT,double minST,double maxST); 
        explicit UniformGridRVT(unsigned nt,double minT,double maxT); 
        ~UniformGridRVT() final = default; 
        const std::vector<double>& getST() const final {return m_st;}
        unsigned getNst() const final {return m_nst;}
        double getMinST() const final {return m_minST;}
        double getMaxST() const final {return m_maxST;}
        unsigned getMinI() const {return m_minI;}
        unsigned getMaxI() const {return m_maxI;}
        double getDST() const {return getLST()/(getNst()-1);}
        double maxPossibleFreq() const;
        unsigned freqToIndx(double omeg) const;
        double indxToFreq(unsigned indx) const;
    private:
        std::vector<double> m_st;
        unsigned m_nst;
        double m_minST;
        double m_maxST;
        unsigned m_minI{0};
        unsigned m_maxI;
        void verifyFrequencyRange(double minST,double maxST) const;
        std::vector<double> constructGridST(unsigned nt,double minT,double maxT) const;
};

}