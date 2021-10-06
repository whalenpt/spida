// uniformT.h
#pragma once

#include <vector>
#include "spida/helper/constants.h" // definition of dcmplx
#include "spida/grid/uniformT.h"
#include "spida/grid/uniformCVT.h"

namespace spida{


class UniformGridCVT : public UniformGridT
{
public:
      explicit UniformGridCVT(unsigned nt,double minT,double maxT);
      explicit UniformGridCVT(const UniformGridCVT& grid);
      ~UniformGridCVT() {}; 
      const std::vector<double>& getST() const {return m_st;}
      unsigned getNst() const {return getNt();}
      double getMinST() const {return -getNt()*spida::PI/getLT();} 
      double getMaxST() const {return getNt()*spida::PI/getLT();}
      double getDST() const {return getLST()/(getNst()-1);}
      double getDT() const {return getLT()/(getNt()-1);}
      double getLST() const {return getMaxST()-getMinST();}
      std::vector<double> freqshift(const std::vector<double>& in) const;
      std::vector<dcmplx> freqshift(const std::vector<dcmplx>& in) const;
      void freqshift(const std::vector<double>& in,std::vector<double>& out) const;
      void freqshift(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) const;
private:
      std::vector<double> m_st;
};



}


