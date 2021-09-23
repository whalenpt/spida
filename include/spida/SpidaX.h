// SpidaX.h
#pragma once

#include <vector>
#include "spida/helper/constants.h"

namespace spida{
  class UniformGridX;
  class FFTX;

  // Assumes complex data
  class SpidaX
  {
    public:
      SpidaX(const UniformGridX& grid);
      SpidaX() = delete;
      ~SpidaX();
      void dX(const std::vector<dcmplx>& in,std::vector<dcmplx>& out,int n = 1) noexcept; 
      const std::vector<double>& getX() const;
      const std::vector<double>& getSX() const;
      void X_To_SX(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) noexcept;
      void SX_To_X(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) noexcept; 
      const UniformGridX& getGridX() const;
      const FFTX& getTransformX() const;
    private:
      std::unique_ptr<UniformGridX> m_gr;
      std::unique_ptr<FFTX> m_tr;
      std::vector<dcmplx> m_vs;
  };
}



