
#ifndef SPIDA_GRID_UNIFORMT_H_
#define SPIDA_GRID_UNIFORMT_H_

#include <vector>

namespace spida{

void setUniformT(double minT,double maxT,std::vector<double>& t);
void setUniformST(double minST,double maxST,std::vector<double>& st);

class UniformGridT 
{
  public:
      UniformGridT(int nt,double minT,double maxT); 
      UniformGridT(int nt,double minT,double maxT,double minST,double maxST); 
      UniformGridT(const UniformGridT& grid);
      ~UniformGridT() {}; 
      const std::vector<double>& getT() const {return m_t;}
      const std::vector<double>& getST() const {return m_st;}
      int getNt() const {return m_t.size();}
      int getNst() const {return m_st.size();}
      double getMinT() const {return m_t[0];}
      double getMinST() const {return m_st[0];}
      double getMaxT() const {return m_t.back();}
      double getMaxST() const {return m_st.back();}
      double getDT() const {return (m_t.back()-m_t[0])/(m_t.size()-1);}
      double getDST() const {return (m_st.back()-m_st[0])/(m_st.size()-1);}
      double getLT() const {return m_t.back() - m_t[0];}
      double getLST() const {return m_st.back() - m_st[0];}
      int convertFreqToIndx(double omeg) const;
      double convertIndxToFreq(int indx) const;
      int getMinI() const {return m_minI;}
      int getMaxI() const {return m_maxI;}
      void checkGrid();
      void checkMinFreq(double minST);
      void checkMaxFreq(double maxST);
      void checkRangeFreq(double minST,double maxST);
      void checkNumSpectralPoints(int nt,int nST);
  private:
      std::vector<double> m_t;
      std::vector<double> m_st;
      int m_minI; 
      int m_maxI; 
      int resFreqToIndxT(double omeg) const;
};




}

#endif


