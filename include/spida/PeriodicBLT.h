
#ifndef SPIDA_PERIODICBLT_H_
#define SPIDA_PERIODICBLT_H_

#include <vector>
#include "spida/helper/constants.h"
#include "spida/grid/uniformT.h"
#include "spida/transform/fftBLT.h"

namespace spida{
  // Assumes real data
  class PeriodicBLT 
  {
    public:
      //PeriodicBLT(int nt,double minT,double maxT,double minST,double maxST);
      PeriodicBLT(const UniformGridT& grid);
      PeriodicBLT() = delete;
      ~PeriodicBLT();
      void dT(const std::vector<double>& in,std::vector<double>& out,int n = 1); 

      const std::vector<double>& getT() const {return m_gr->getT();}
      const std::vector<double>& getST() const {return m_gr->getST();}
      void T_To_ST(const std::vector<double>& in,std::vector<dcmplx>& out) {m_tr->T_To_ST(in,out);} 
      void ST_To_T(const std::vector<dcmplx>& in,std::vector<double>& out) {m_tr->ST_To_T(in,out);} 
      void T_To_ST_c(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) {m_tr->T_To_ST_c(in,out);} 
      void ST_To_T_c(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) {m_tr->ST_To_T_c(in,out);} 
      const UniformGridT& getGridT() const {return *m_gr;}
      const FFTBLT& getTransformT() const {return *m_tr;}
    private:
      UniformGridT* m_gr;
      FFTBLT* m_tr;
      std::vector<dcmplx> m_vs;
  };
}

#endif


