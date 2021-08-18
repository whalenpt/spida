
#ifndef SPIDA_GRIDR_H_
#define SPIDA_GRIDR_H_

#include <vector>

namespace spida{

class GridR{
  public:
      GridR(int nr,double rmax)  : 
          m_nr(nr),m_rmax(rmax) {}
      virtual ~GridR() {}
      virtual const std::vector<double>& getR() const = 0;
      virtual const std::vector<double>& getSR() const = 0;
      int getNr() const {return m_nr;}
      int getNsr() const {return m_nr;}
      double getMaxR() const {return m_rmax;}
  private:
      int m_nr;
      double m_rmax;
};



}

#endif


