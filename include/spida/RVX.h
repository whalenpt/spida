/*  
 *  Spectral integration/differentiation algorithms for 
 *  real-valued data on uniform grids with periodic boundaries. 
 *
 *  Transform uses KissFFT, with phase adjustment since the physical grid
 *  is not assumed to run from 0 to 1.
 *
 *  X_To_SX and SX_To_X transform from physical space to spectral space and vice versa.
 *
 *  dX transforms to spectral space, then takes derivative, then transforms back to real space.
 *  dSX takes a derivate in spectral space
 */

#pragma once

#include <vector>
#include "spida/grid/uniformRVX.h"
#include "spida/helper/constants.h"
#include "spida/transform/fftRVX.h"

namespace spida{
  class UniformGridRVX;
  class FFTRVX;

  // Assumes complex data
  class SpidaRVX
  {
    public:
      explicit SpidaRVX(const UniformGridRVX& grid);
      SpidaRVX() = delete;
      ~SpidaRVX() = default;
      void dX(const std::vector<double>& in,std::vector<double>& out,unsigned n = 1) noexcept; 
      void dSX(const std::vector<dcmplx>& in,std::vector<dcmplx>& out,unsigned n = 1) noexcept; 
      void dSX(std::vector<dcmplx>& in,unsigned n = 1) noexcept; 
      const std::vector<double>& getX() const;
      const std::vector<double>& getSX() const;
      void X_To_SX(const std::vector<double>& in,std::vector<dcmplx>& out) noexcept;
      void SX_To_X(const std::vector<dcmplx>& in,std::vector<double>& out) noexcept; 
      const UniformGridRVX& getGridX() const;
      const FFTRVX& getTransformX() const;
    private:
      std::unique_ptr<UniformGridRVX> m_gr;
      std::unique_ptr<FFTRVX> m_tr;
      std::vector<dcmplx> m_vs;
  };
}