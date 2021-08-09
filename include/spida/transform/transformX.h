
#ifndef SPIDA_TRANSFORMSX_H_
#define SPIDA_TRANSFORMSX_H_ 

#include <complex>
#include <vector>
#include "spida/grid/gridX.h"

namespace spida{


//  class PlanXr;
//  class PlanChebTransformXdct;

  // interface class
  class TransformX
  {
      public:
          TransformX(const GridX& grid) {}
          virtual ~TransformX() {}
          TransformX()=delete;
          TransformX(const TransformX&)=delete;
          TransformX& operator=(const TransformX&)=delete;
          virtual void X_To_SX(const std::vector<double>& in,std::vector<double>& out) = 0;
          virtual void SX_To_X(const std::vector<double>& in,std::vector<double>& out) = 0;
  };

//  class ChebTransformX : public TransformX
//  {
//    public:
//        ChebTransformX(const ChebRootGridX& grid);
//        ~ChebTransformX();
//        void X_To_SX(const std::vector<double>& in,std::vector<double>& out);
//        void SX_To_X(const std::vector<double>& in,std::vector<double>& out);
//    private:
//        int m_nx;
//        std::vector<double> m_uFFT;
//  };
//
//  class PlanChebTransformXdct 
//  {
//    public: 
//      PlanChebTransformXdct(int nxx);
//      PlanChebTransformXdct()=delete;
//      PlanChebTransformXdct(const PlanXr& sp)=delete;
//      PlanChebTransformXdct& operator=(const PlanXr& sp)=delete;
//      ~PlanChebTransformXdct();
//      void SX_To_X(const double* in,double* out); 
//      void X_To_SX(const double* in,double* out); 
//      void SpT_To_Re_2D(const double* in,double* out,int nyy=1); 
//      void Re_To_SpT_2D(const double* in,double* out,int nyy=1); 
//      void ReT_To_SpT_2D(const double* in,double* out,int nd2 = 1); 
//      void SpT_To_ReT_2D(const double* in,double* out,int nd2 = 1);
//    private:
//      int sz;
//      int nlogic;
//      int N;
//      double* uFFT;
//      fftw_plan sp_fftw_plan;
//      fftw_plan re_fftw_plan;
//  };
//
//  class PlanXr 
//  {
//    public: 
//      PlanXr(int nxx);
//      PlanXr()=delete;
//      PlanXr(const PlanXr& sp)=delete;
//      PlanXr& operator=(const PlanXr& sp)=delete;
//      ~PlanXr();
//      void SX_To_X(const double* in,double* out); 
//      void X_To_SX(const double* in,double* out); 
//      void SpT_To_Re_2D(const double* in,double* out,int nyy=1); 
//      void Re_To_SpT_2D(const double* in,double* out,int nyy=1); 
//    private:
//      int sz;
//      double* uFFT;
//      fftw_plan sp_fftw_plan;
//      fftw_plan re_fftw_plan;
//  };

}


#endif


