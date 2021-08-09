 
#ifndef SPIDA_GRIDX_H_
#define SPIDA_GRIDX_H_

#include <vector>
int indexFromVal(double val,std::vector<double> vec);

namespace spida{

class GridX{
  public:
      GridX(int nx,double minX,double maxX)  : 
          m_nx(nx),m_minX(minX),m_maxX(maxX) {}
      virtual ~GridX() {}
      virtual const std::vector<double>& getX() const = 0;
      virtual const std::vector<double>& getSX() const = 0;
      int getNx() const {return m_nx;}
      int getNsx() const {return m_nx;}
      double getMinX() const {return m_minX;}
      double getMaxX() const {return m_maxX;}
      virtual double getMinSX() const = 0;
      virtual double getMaxSX() const = 0;
      int getIndexX(double xval);
  private:
      int m_nx;
      double m_minX;
      double m_maxX;
};



}

#endif


