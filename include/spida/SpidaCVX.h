// SpidaCVX.h
/*  
 *  Spectral integration/differentiation algorithms for 
 *  complex data on uniform grids with periodic boundaries. 
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
#include "spida/helper/constants.h"

namespace spida{
  class UniformGridCVX;
  class FFTCVX;

  // Assumes complex data
  class SpidaCVX
  {
    public:
      SpidaCVX(const UniformGridCVX& grid);
      SpidaCVX() = delete;
      ~SpidaCVX();
      void dX(const std::vector<dcmplx>& in,std::vector<dcmplx>& out,int n = 1) noexcept; 
      void dSX(const std::vector<dcmplx>& in,std::vector<dcmplx>& out,int n = 1) noexcept; 
      const std::vector<double>& getX() const;
      const std::vector<double>& getSX() const;
      void X_To_SX(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) noexcept;
      void SX_To_X(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) noexcept; 
      const UniformGridCVX& getGridX() const;
      const FFTCVX& getTransformX() const;
    private:
      std::unique_ptr<UniformGridCVX> m_gr;
      std::unique_ptr<FFTCVX> m_tr;
      std::vector<dcmplx> m_vs;
  };
}



