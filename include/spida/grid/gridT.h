#pragma once

#include <vector>
#include "spida/grid/grid.h"

namespace spida{

class GridT : public Grid
{
  public:
      explicit GridT(unsigned nt,double minT,double maxT)  : 
          m_nt(nt),m_minT(minT),m_maxT(maxT) {}
      ~GridT() override = default;
      virtual const std::vector<double>& getT() const = 0;
      virtual const std::vector<double>& getST() const = 0;
      unsigned getNt() const {return m_nt;}
      virtual unsigned getNst() const = 0; 
      double getMinT() const {return m_minT;}
      double getMaxT() const {return m_maxT;}
      double getLT() const {return m_maxT-m_minT;}
  private:
      unsigned m_nt;
      double m_minT;
      double m_maxT;
};

}