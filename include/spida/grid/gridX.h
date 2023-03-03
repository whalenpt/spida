#pragma once

#include <vector>
#include "spida/grid/grid.h"

namespace spida{

class GridX : public Grid
{
  public:
      GridX(unsigned nx,double minX,double maxX)  : 
          m_nx(nx),m_minX(minX),m_maxX(maxX) {}
      ~GridX() override = default;
      virtual const std::vector<double>& getX() const = 0;
      virtual const std::vector<double>& getSX() const = 0;
      unsigned getNx() const {return m_nx;}
      double getLX() const {return m_maxX-m_minX;}
      double getMinX() const {return m_minX;}
      double getMaxX() const {return m_maxX;}
      virtual unsigned getNsx() const {return m_nx;}
      virtual double getMinSX() const = 0;
      virtual double getMaxSX() const = 0;
      unsigned int getIndexX(double xval);
  private:
      unsigned m_nx;
      double m_minX;
      double m_maxX;
};

}