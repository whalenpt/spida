
#include <cmath>
#include <memory>
#include "spida/RRVT.h"
#include "spida/grid/uniformRVT.h"
#include "spida/grid/besselR.h"
#include "spida/transform/hankelfftRRVT.h"
#include <pwutils/pwexcept.h>

namespace spida{

  SpidaRRVT::SpidaRRVT(const BesselRootGridR& gridR,\
          const UniformGridRVT& gridT,unsigned threads) :
      m_gridR(std::make_unique<BesselRootGridR>(gridR)),
      m_gridT(std::make_unique<UniformGridRVT>(gridT)),
      m_tr(std::make_unique<HankelFFTRRVT>(gridR,gridT,threads)) { }
  
  const BesselRootGridR& SpidaRRVT::getGridR() const { return *m_gridR; }
  const UniformGridRVT& SpidaRRVT::getGridT() const { return *m_gridT; }
  const HankelFFTRRVT& SpidaRRVT::getTransformRT() const { return *m_tr; }
  const std::vector<double>& SpidaRRVT::getR() const {return m_gridR->getR();}
  const std::vector<double>& SpidaRRVT::getSR() const {return m_gridR->getSR();}
  const std::vector<double>& SpidaRRVT::getT() const {return m_gridT->getT();}
  const std::vector<double>& SpidaRRVT::getST() const {return m_gridT->getST();}
  void SpidaRRVT::RT_To_SRST(const std::vector<double>& in,std::vector<dcmplx>& out) {m_tr->RT_To_SRST(in,out);} 
  void SpidaRRVT::SRST_To_RT(const std::vector<dcmplx>& in,std::vector<double>& out) {m_tr->SRST_To_RT(in,out);} 

}
