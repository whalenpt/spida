// gridR.h
#pragma once

#include <vector>
#include "spida/grid/grid.h"
#include "spida/helper/constants.h"

namespace spida{

class GridR : public Grid
{
  public:
      GridR(unsigned nr,double rmax)  : 
          m_nr(nr),m_rmax(rmax) {}
      virtual ~GridR() {}
      virtual const std::vector<double>& getR() const = 0;
      virtual const std::vector<double>& getSR() const = 0;
      int getNr() const {return m_nr;}
      int getNsr() const {return m_nr;}
      double getMaxR() const {return m_rmax;}

      std::vector<double> mirrorGrid(const std::vector<double>& in,\
              bool sign_reverse=false) const;
      std::vector<dcmplx> mirrorGrid(const std::vector<dcmplx>& in,\
              bool sign_reverse=false) const;

      void mirrorGrid(const std::vector<double>& in,std::vector<double>& out,\
              bool sign_reverse=false) const;
      void mirrorGrid(const std::vector<dcmplx>& in,std::vector<dcmplx>& out,\
              bool sign_reverse=false) const;
      void mirrorGrid(const double* in,double* out,bool sign_reverse=false) const;
      void mirrorGrid(const dcmplx* in,dcmplx* out,bool sign_reverse=false) const;


  private:
      unsigned m_nr;
      double m_rmax;
};



}



