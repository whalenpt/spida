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
      ~UniformGridCVT() override = default; 
      const std::vector<double>& getST() const final {return m_st;}
      unsigned getNst() const final {return getNt();}
      double getMinST() const final {return -(getNt()*spida::PI)/getLT();} 
      double getMaxST() const final {return getNt()*spida::PI/getLT();}
      double getDST() const {return getLST()/(getNst()-1);}
      std::vector<double> freqshift(const std::vector<double>& in) const;
      std::vector<dcmplx> freqshift(const std::vector<dcmplx>& in) const;
      void freqshift(const std::vector<double>& in,std::vector<double>& out) const;
      void freqshift(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) const;
private:
      std::vector<double> m_st;
};

}