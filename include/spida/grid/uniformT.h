// uniformT.h
#pragma once

#include <vector>
#include "spida/grid/gridT.h"

namespace spida{

class UniformGridT : public GridT
{
  public:
      explicit UniformGridT(unsigned nt,double minT,double maxT); 
      explicit UniformGridT(unsigned nt,double minT,double maxT,double minST,double maxST); 
      explicit UniformGridT(const UniformGridT& grid);
      ~UniformGridT() {}; 
      const std::vector<double>& getT() const {return m_t;}
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
      std::vector<double> constructGridT(unsigned int nt,double minT,double maxT) const;
      std::vector<double> constructGridST(unsigned int nt,double minT,double maxT) const;
  private:
      unsigned int m_nst;
      double m_minST;
      double m_maxST;
      unsigned m_minI;
      unsigned m_maxI;
      std::vector<double> m_t;
      std::vector<double> m_st;

      void constructGridST(unsigned int nt,double minT,double maxT,double minST,double maxST);
      unsigned resFreqToIndxT(double omeg) const;
      void verifyFrequencyRange(double minST,double maxST) const;
};




}



