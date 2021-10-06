// SpidaRCVT.h
#pragma once

#include <vector>
#include "spida/helper/constants.h"

namespace spida{
  class UniformGridCVT;
  class BesselRootGridR;
  class HankelFFTRCVT;

  // Assumes real data
  class SpidaRCVT
  {
    public:
      SpidaRCVT(const BesselRootGridR& gridR,const UniformGridCVT& gridT,unsigned threads=1);
      SpidaRCVT() = delete;
      ~SpidaRCVT();
      const BesselRootGridR& getGridR() const;
      const UniformGridCVT& getGridT() const;
      const HankelFFTRCVT& getTransformRT() const;
      const std::vector<double>& getR() const;
      const std::vector<double>& getSR() const;
      const std::vector<double>& getT() const;
      const std::vector<double>& getST() const;
      void RT_To_SRST(const std::vector<dcmplx>& in,std::vector<dcmplx>& out); 
      void SRST_To_RT(const std::vector<dcmplx>& in,std::vector<dcmplx>& out); 
      unsigned spectralSize() const;
      unsigned physicalSize() const;
      void mirrorR(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) const;
//      void mirrorGrid(

    private:
      std::unique_ptr<BesselRootGridR> m_gridR;
      std::unique_ptr<UniformGridCVT> m_gridT;
      std::unique_ptr<HankelFFTRCVT> m_tr;
  };
}


