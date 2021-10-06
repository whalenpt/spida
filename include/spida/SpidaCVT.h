// SpidaCVT.h
/*  
 *  Spectral integration/differentiation algorithms for 
 *  complex data on uniform grids with periodic boundaries. 
 *
 *  Transform uses KissFFT, with phase adjustment since the physical grid
 *  is not assumed to run from 0 to 1.
 *
 *  T_To_ST and ST_To_T transform from physical space to spectral space and vice versa.
 *
 *  dT transforms to spectral space, then takes derivative, then transforms back to real space.
 *  dST takes a derivate in spectral space
 */

#pragma once

#include <vector>
#include "spida/helper/constants.h"

namespace spida{
  class UniformGridCVT;
  class FFTCVT;

  // Assumes complex data
  class SpidaCVT
  {
    public:
      SpidaCVT(const UniformGridCVT& grid);
      SpidaCVT() = delete;
      ~SpidaCVT();
      void dT(const std::vector<dcmplx>& in,std::vector<dcmplx>& out,unsigned n = 1) noexcept; 
      void dST(const std::vector<dcmplx>& in,std::vector<dcmplx>& out,unsigned n = 1) noexcept; 
      const std::vector<double>& getT() const;
      const std::vector<double>& getST() const;
      void T_To_ST(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) noexcept;
      void ST_To_T(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) noexcept; 
      const UniformGridCVT& getGridT() const;
      const FFTCVT& getTransformT() const;
    private:
      std::unique_ptr<UniformGridCVT> m_gr;
      std::unique_ptr<FFTCVT> m_tr;
      std::vector<dcmplx> m_vs;
  };
}



