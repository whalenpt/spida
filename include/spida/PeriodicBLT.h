
#ifndef SPIDA_PERIODICBLT_H_
#define SPIDA_PERIODICBLT_H_

#include <vector>
#include "spida/helper/constants.h"

namespace spida{
  class UniformGridT;
  class FFTBLT;

  // Assumes real data
  class PeriodicBLT 
  {
    public:
      //PeriodicBLT(int nt,double minT,double maxT,double minST,double maxST);
      PeriodicBLT(const UniformGridT& grid);
      PeriodicBLT() = delete;
      ~PeriodicBLT();
      void dT(const std::vector<double>& in,std::vector<double>& out,int n = 1); 
      const std::vector<double>& getT() const;
      const std::vector<double>& getST() const;

      void T_To_ST(const std::vector<double>& in,std::vector<dcmplx>& out); 
      void ST_To_T(const std::vector<dcmplx>& in,std::vector<double>& out);
      void T_To_ST_c(const std::vector<dcmplx>& in,std::vector<dcmplx>& out);
      void ST_To_T_c(const std::vector<dcmplx>& in,std::vector<dcmplx>& out); 
      const UniformGridT& getGridT() const;
      const FFTBLT& getTransformT() const;
    private:
      UniformGridT* m_gr;
      FFTBLT* m_tr;
      std::vector<dcmplx> m_vs;
  };
}

#endif


