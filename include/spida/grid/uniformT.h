// uniformT.h
#pragma once

#include <vector>
#include "spida/grid/gridT.h"

namespace spida{

class UniformGridT : public GridT
{
  public:
      explicit UniformGridT(unsigned nt,double minT,double maxT); 
      ~UniformGridT() {}; 
      const std::vector<double>& getT() const {return m_t;}
      virtual const std::vector<double>& getST() const = 0;
      virtual unsigned getNst() const = 0;
      virtual double getMinST() const = 0;
      virtual double getMaxST() const = 0;
      double getDT() const {return GridT::getLT()/(GridT::getNt()-1);}
  private:
      std::vector<double> m_t;
};


}



