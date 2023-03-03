// RVT.h 
// Class for (spectrally) band-limited temporal functions
#pragma once

#include <vector>
#include "spida/helper/constants.h"
#include "spida/grid/uniformRVT.h"
#include "spida/transform/fftRVT.h"

namespace spida{
  // Assumes real data
  class SpidaRVT 
  {
    public:
      explicit SpidaRVT(const UniformGridRVT& grid);
      SpidaRVT() = delete;
      ~SpidaRVT() = default;
      void dT(const std::vector<double>& in,std::vector<double>& out,int n = 1); 
      const std::vector<double>& getT() const;
      const std::vector<double>& getST() const;
      void T_To_ST(const std::vector<double>& in,std::vector<dcmplx>& out); 
      void ST_To_T(const std::vector<dcmplx>& in,std::vector<double>& out);
      void CVT_To_ST(const std::vector<dcmplx>& in,std::vector<dcmplx>& out);
      void ST_To_CVT(const std::vector<dcmplx>& in,std::vector<dcmplx>& out); 
      const UniformGridRVT& getGridT() const;
      const FFTRVT& getTransformT() const;
    private:
      std::unique_ptr<UniformGridRVT> m_gr;
      std::unique_ptr<FFTRVT> m_tr;
      std::vector<dcmplx> m_vs;
  };

}