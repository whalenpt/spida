
#ifndef SPIDA_GRID_UNIFORMT_H_
#define SPIDA_GRID_UNIFORMT_H_

#include <vector>
#include "spida/grid/gridT.h"

namespace spida{

std::vector<double> buildUniformT(unsigned int nt,double minT,double maxT);
std::vector<double> buildUniformST(unsigned int nt,double minT,double maxT);

class UniformGridT : public GridT
{
  public:
      explicit UniformGridT(unsigned int nt,double minT,double maxT); 
      explicit UniformGridT(unsigned int nt,double minT,double maxT,double minST,double maxST); 
      explicit UniformGridT(const UniformGridT& grid);
      ~UniformGridT() {}; 
      const std::vector<double>& getT() const {return m_t;}
      const std::vector<double>& getST() const {return m_st;}
      unsigned int getNst() const {return m_nst;}
      double getMinST() const {return m_st[0];}
      double getMaxST() const {return m_st.back();}
      unsigned int getMinI() const {return m_minI;}
      unsigned int getMaxI() const {return m_maxI;}
      double getDST() const {return getLST()/(getNst()-1);}
      double getDT() const {return getLT()/(getNt()-1);}
      double getLST() const {return m_maxST-m_minST;}
      double maxPossibleFreq() const;
      unsigned int freqToIndx(double omeg) const;
      double indxToFreq(unsigned int indx) const;
  private:
      unsigned int m_nst;
      double m_minST;
      double m_maxST;
      unsigned int m_minI;
      unsigned int m_maxI;
      std::vector<double> m_t;
      std::vector<double> m_st;
      unsigned int resFreqToIndxT(double omeg) const;
      void verifyFrequencyRange(double minST,double maxST) const;
};




}

#endif


