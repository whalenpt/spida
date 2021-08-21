#ifndef SPIDA_GRIDT_H_
#define SPIDA_GRIDT_H_

#include <vector>

namespace spida{

class GridT{
  public:
      explicit GridT(int nt,double minT,double maxT)  : 
          m_nt(nt),m_minT(minT),m_maxT(maxT) {}
      virtual ~GridT() {}
      virtual const std::vector<double>& getT() const = 0;
      virtual const std::vector<double>& getST() const = 0;
      int getNt() const {return m_nt;}
      virtual int getNst() const = 0; 
      double getMinT() const {return m_minT;}
      double getMaxT() const {return m_maxT;}
  private:
      int m_nt;
      double m_minT;
      double m_maxT;
};



}

#endif

