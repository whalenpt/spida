// SpidaBLT.h 
// Class for (spectrally) band-limited temporal functions
#pragma once

#include <vector>
#include "spida/helper/constants.h"

namespace spida{
  class UniformGridRVT;
  class FFTBLT;

  // Assumes real data
  class SpidaBLT 
  {
    public:
      //SpidaBLT(int nt,double minT,double maxT,double minST,double maxST);
      SpidaBLT(const UniformGridRVT& grid);
      SpidaBLT() = delete;
      ~SpidaBLT();
      void dT(const std::vector<double>& in,std::vector<double>& out,int n = 1); 
      const std::vector<double>& getT() const;
      const std::vector<double>& getST() const;
      void T_To_ST(const std::vector<double>& in,std::vector<dcmplx>& out); 
      void ST_To_T(const std::vector<dcmplx>& in,std::vector<double>& out);
      void CVT_To_ST(const std::vector<dcmplx>& in,std::vector<dcmplx>& out);
      void ST_To_CVT(const std::vector<dcmplx>& in,std::vector<dcmplx>& out); 
      const UniformGridRVT& getGridT() const;
      const FFTBLT& getTransformT() const;
    private:
      std::unique_ptr<UniformGridRVT> m_gr;
      std::unique_ptr<FFTBLT> m_tr;
      std::vector<dcmplx> m_vs;
  };
}



