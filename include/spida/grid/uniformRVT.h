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
        explicit UniformGridRVT(const UniformGridRVT& grid);
        ~UniformGridRVT() {}; 
        const std::vector<double>& getST() const {return m_st;}
        unsigned getNst() const {return m_nst;}
        double getMinST() const {return m_minST;}
        double getMaxST() const {return m_maxST;}
        unsigned getMinI() const {return m_minI;}
        unsigned getMaxI() const {return m_maxI;}
        double getDST() const {return getLST()/(getNst()-1);}
        double getDT() const {return getLT()/(getNt()-1);}
        double getLST() const {return m_maxST-m_minST;}
        double maxPossibleFreq() const;
        unsigned freqToIndx(double omeg) const;
        double indxToFreq(unsigned indx) const;
    private:
        std::vector<double> m_st;
        unsigned m_nst;
        double m_minST;
        double m_maxST;
        unsigned m_minI;
        unsigned m_maxI;
        unsigned resFreqToIndxT(double omeg) const;
        std::vector<double> constructGridST(unsigned nt,double minST,double maxST);
        void verifyFrequencyRange(double minST,double maxST) const;
};

}



