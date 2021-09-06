
#ifndef SPIDARBLT_H_
#define SPIDARBLT_H_

#include <vector>
#include "spida/helper/constants.h"

namespace spida{
  class UniformGridT;
  class BesselRootGridR;
  class HankelFFTRBLT;

  // Assumes real data
  class SpidaRBLT 
  {
    public:
      SpidaRBLT(const BesselRootGridR& gridR,const UniformGridT& gridT,unsigned int threads=1);
      SpidaRBLT() = delete;
      ~SpidaRBLT() {};
      const BesselRootGridR& getGridR() const;
      const UniformGridT& getGridT() const;
      const HankelFFTRBLT& getTransformRT() const;
      const std::vector<double>& getR() const;
      const std::vector<double>& getSR() const;
      const std::vector<double>& getT() const;
      const std::vector<double>& getST() const;
      void RT_To_SRST(const std::vector<double>& in,std::vector<dcmplx>& out); 
      void SRST_To_RT(const std::vector<dcmplx>& in,std::vector<double>& out); 
//      void dT(const std::vector<double>& in,std::vector<double>& out,unsigned int n = 1); 

    private:
      std::unique_ptr<BesselRootGridR> m_gridR;
      std::unique_ptr<UniformGridT> m_gridT;
      std::unique_ptr<HankelFFTRBLT> m_tr;
//      std::vector<dcmplx> m_vs;
  };
}

#endif


