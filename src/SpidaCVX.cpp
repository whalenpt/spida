
#include <cmath>
#include <memory>
#include "spida/SpidaX.h"
#include "spida/grid/uniformX.h"
#include "spida/transform/fftX.h"
#include <pwutils/pwexcept.h>

namespace spida{

  SpidaX::SpidaX(const UniformGridX& grid) :
      m_gr(std::make_unique<UniformGridX>(grid)),
      m_tr(std::make_unique<FFTX>(grid)),
      m_vs(grid.getNsx()) {}

  SpidaX::~SpidaX() {}
  const UniformGridX& SpidaX::getGridX() const { return *m_gr; }
  const FFTX& SpidaX::getTransformX() const { return *m_tr; }
  const std::vector<double>& SpidaX::getX() const {return m_gr->getX();}
  const std::vector<double>& SpidaX::getSX() const {return m_gr->getSX();}
  void SpidaX::X_To_SX(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) noexcept {m_tr->X_To_SX(in,out);} 
  void SpidaX::SX_To_X(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) noexcept {m_tr->SX_To_X(in,out);} 

  void SpidaX::dX(const std::vector<dcmplx>& in,std::vector<dcmplx>& out,int n) noexcept
  {
    if(n < 0)
        return;
    m_tr->X_To_SX(in,m_vs);
    const std::vector<double>& kx = m_gr->getSX();
    unsigned nx = m_gr->getNx();
    for(auto i = 0; i < nx; i++)
        m_vs[i] = std::pow(ii*kx[i],n)*m_vs[i];
    m_tr->SX_To_X(m_vs,out);
  }

  void SpidaX::dSX(const std::vector<dcmplx>& in,std::vector<dcmplx>& out,int n) noexcept
  {
    if(n < 0)
        return;
    const std::vector<double>& kx = m_gr->getSX();
    for(auto i = 0; i < kx.size(); i++)
        out[i] = std::pow(ii*kx[i],n)*in[i];
  }

}






