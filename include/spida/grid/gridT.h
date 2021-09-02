#ifndef SPIDA_GRIDT_H_
#define SPIDA_GRIDT_H_

#include <vector>
#include "spida/grid/grid.h"

namespace spida{

class GridT : public Grid
{
  public:
      explicit GridT(unsigned int nt,double minT,double maxT)  : 
          m_nt(nt),m_minT(minT),m_maxT(maxT) {}
      virtual ~GridT() {}
      virtual const std::vector<double>& getT() const = 0;
      virtual const std::vector<double>& getST() const = 0;
      unsigned int getNt() const {return m_nt;}
      virtual unsigned int getNst() const = 0; 
      double getMinT() const {return m_minT;}
      double getMaxT() const {return m_maxT;}
      double getLT() const {return m_maxT-m_minT;}
      double getDT() const {return getLT()/(getNt()-1);}

  private:
      unsigned int m_nt;
      double m_minT;
      double m_maxT;
};



}

#endif


