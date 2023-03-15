
#include <cmath>
#include <memory>
#include "spida/RVX.h"
#include "spida/grid/uniformRVX.h"
#include "spida/transform/fftRVX.h"

namespace spida{

  SpidaRVX::SpidaRVX(const UniformGridRVX& grid) :
      m_gr(std::make_unique<UniformGridRVX>(grid)),
      m_tr(std::make_unique<FFTRVX>(grid)),
      m_vs(grid.getNsx()) {}

  const UniformGridRVX& SpidaRVX::getGridX() const { return *m_gr; }
  const FFTRVX& SpidaRVX::getTransformX() const { return *m_tr; }
  const std::vector<double>& SpidaRVX::getX() const {return m_gr->getX();}
  const std::vector<double>& SpidaRVX::getSX() const {return m_gr->getSX();}
  void SpidaRVX::X_To_SX(const std::vector<double>& in,std::vector<dcmplx>& out) noexcept {m_tr->X_To_SX(in,out);} 
  void SpidaRVX::SX_To_X(const std::vector<dcmplx>& in,std::vector<double>& out) noexcept {m_tr->SX_To_X(in,out);} 

  void SpidaRVX::dX(const std::vector<double>& in,std::vector<double>& out,unsigned n) noexcept
  {
    if(n == 0)
        return;
    m_tr->X_To_SX(in,m_vs);
    dSX(m_vs,n);
    m_tr->SX_To_X(m_vs,out);
  }

  void SpidaRVX::dSX(const std::vector<dcmplx>& in,std::vector<dcmplx>& out,unsigned n) noexcept
  {
    std::copy(std::cbegin(in),std::cend(in),std::begin(out));
    dSX(out,n);
  }

  void SpidaRVX::dSX(std::vector<dcmplx>& in,unsigned n) noexcept
  {
    if(n == 0)
        return;
    const auto& kx = m_gr->getSX();
    for(size_t i = 0; i < kx.size(); i++)
        in[i] = std::pow(ii*kx[i],n)*in[i];
  }



}






