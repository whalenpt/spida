// uniformT.h
#pragma once

#include <vector>
#include "spida/grid/gridT.h"

namespace spida{

class UniformGridT : public GridT
{
  public:
      explicit UniformGridT(unsigned nt,double minT,double maxT); 
      ~UniformGridT() override = default;
      const std::vector<double>& getT() const final {return m_t;}
      virtual double getMinST() const = 0;
      virtual double getMaxST() const = 0;
      double getLST() const {return getMaxST()-getMinST();}
      double getDT() const {return GridT::getLT()/(GridT::getNt()-1);}
  private:
      std::vector<double> m_t;
};


}



