
#include <cmath>
#include <memory>
#include "spida/CVX.h"
#include "spida/grid/uniformCVX.h"
#include "spida/transform/fftCVX.h"

namespace spida{

  SpidaCVX::SpidaCVX(const UniformGridCVX& grid) :
      m_gr(std::make_unique<UniformGridCVX>(grid)),
      m_tr(std::make_unique<FFTCVX>(grid)),
      m_vs(grid.getNsx()) {}

  const UniformGridCVX& SpidaCVX::getGridX() const { return *m_gr; }
  const FFTCVX& SpidaCVX::getTransformX() const { return *m_tr; }
  const std::vector<double>& SpidaCVX::getX() const {return m_gr->getX();}
  const std::vector<double>& SpidaCVX::getSX() const {return m_gr->getSX();}
  void SpidaCVX::X_To_SX(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) noexcept {m_tr->X_To_SX(in,out);} 
  void SpidaCVX::SX_To_X(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) noexcept {m_tr->SX_To_X(in,out);} 

  void SpidaCVX::dX(const std::vector<dcmplx>& in,std::vector<dcmplx>& out,unsigned n) noexcept
  {
    if(n == 0)
        return;
    m_tr->X_To_SX(in,out);
    dSX(out,m_vs,n);
    m_tr->SX_To_X(m_vs,out);
  }

  void SpidaCVX::dSX(const std::vector<dcmplx>& in,std::vector<dcmplx>& out,unsigned n) const noexcept
  {
    if(n == 0)
        return;
    const auto& kx = m_gr->getSX();
    for(size_t i = 0; i < kx.size(); i++)
        out[i] = std::pow(ii*kx[i],n)*in[i];
  }

}






