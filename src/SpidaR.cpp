
#include <cmath>
#include <memory>
#include "spida/R.h"
#include "spida/helper/constants.h"
#include "spida/grid/besselR.h"
#include "spida/transform/hankelR.h"

namespace spida{


  SpidaR::SpidaR(const BesselRootGridR& gridR) :
      m_gridR(std::make_unique<BesselRootGridR>(gridR)),
      m_tr(std::make_unique<HankelTransformR>(gridR)) { }
  
  const BesselRootGridR& SpidaR::getGridR() const { return *m_gridR; }
  const HankelTransformR& SpidaR::getTransformR() const { return *m_tr; }
  const std::vector<double>& SpidaR::getR() const {return m_gridR->getR();}
  const std::vector<double>& SpidaR::getSR() const {return m_gridR->getSR();}
  void SpidaR::R_To_SR(const std::vector<double>& in,std::vector<double>& out) {m_tr->R_To_SR(in,out);} 
  void SpidaR::SR_To_R(const std::vector<double>& in,std::vector<double>& out) {m_tr->SR_To_R(in,out);} 
  void SpidaR::R_To_SR(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) {m_tr->R_To_SR(in,out);} 
  void SpidaR::SR_To_R(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) {m_tr->SR_To_R(in,out);} 
}








