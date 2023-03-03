// RRVT.h
#pragma once

#include <vector>
#include "spida/grid/besselR.h"
#include "spida/grid/uniformRVT.h"
#include "spida/helper/constants.h"
#include "spida/transform/hankelfftRRVT.h"

namespace spida{
  // Assumes real data
  class SpidaRRVT 
  {
    public:
      SpidaRRVT(const BesselRootGridR& gridR,const UniformGridRVT& gridT,unsigned threads=1);
      SpidaRRVT() = delete;
      ~SpidaRRVT() = default;
      const BesselRootGridR& getGridR() const;
      const UniformGridRVT& getGridT() const;
      const HankelFFTRRVT& getTransformRT() const;
      const std::vector<double>& getR() const;
      const std::vector<double>& getSR() const;
      const std::vector<double>& getT() const;
      const std::vector<double>& getST() const;
      void RT_To_SRST(const std::vector<double>& in,std::vector<dcmplx>& out); 
      void SRST_To_RT(const std::vector<dcmplx>& in,std::vector<double>& out); 
      void CVRT_To_SRST(const std::vector<dcmplx>& in,std::vector<dcmplx>& out); 
      void SRST_To_CVRT(const std::vector<dcmplx>& in,std::vector<dcmplx>& out); 

    private:
      std::unique_ptr<BesselRootGridR> m_gridR;
      std::unique_ptr<UniformGridRVT> m_gridT;
      std::unique_ptr<HankelFFTRRVT> m_tr;
  };

}