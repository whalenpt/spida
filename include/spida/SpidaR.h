
#pragma once

#include <vector>
#include "spida/helper/constants.h"

namespace spida{
  class BesselRootGridR;
  class HankelTransformR;

  // Assumes real data
  class SpidaR 
  {
    public:
      SpidaR(const BesselRootGridR& gridR);
      SpidaR() = delete;
      ~SpidaR();
      const BesselRootGridR& getGridR() const;
      const HankelTransformR& getTransformR() const;
      const std::vector<double>& getR() const;
      const std::vector<double>& getSR() const;
      void R_To_SR(const std::vector<double>& in,std::vector<double>& out); 
      void SR_To_R(const std::vector<double>& in,std::vector<double>& out); 
      void R_To_SR(const std::vector<dcmplx>& in,std::vector<dcmplx>& out); 
      void SR_To_R(const std::vector<dcmplx>& in,std::vector<dcmplx>& out); 
    private:
      std::unique_ptr<BesselRootGridR> m_gridR;
      std::unique_ptr<HankelTransformR> m_tr;
  };
}



