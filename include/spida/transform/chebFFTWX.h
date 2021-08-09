

#ifndef CHEBFFTWX_H_
#define CHEBFFTWX_H_ 

#include <complex>
#include <vector>
#include "spida/transform/transformX.h"
#include <fftw3.h>

namespace spida{


  class PlanXr;
  class PlanChebTransformFFTWXdct;
  class ChebExtremaGridX;

  class ChebTransformFFTWX : public TransformX
  {
    public:
        ChebTransformFFTWX(const ChebExtremaGridX& grid);
        void X_To_SX(const std::vector<double>& in,std::vector<double>& out);
        void SX_To_X(const std::vector<double>& in,std::vector<double>& out);
        ~ChebTransformFFTWX();
    private:
        PlanChebTransformFFTWXdct* m_plan; 
  };

  class PlanChebTransformFFTWXdct 
  {
    public: 
      PlanChebTransformFFTWXdct(int nxx);
      PlanChebTransformFFTWXdct()=delete;
      PlanChebTransformFFTWXdct(const PlanXr& sp)=delete;
      PlanChebTransformFFTWXdct& operator=(const PlanXr& sp)=delete;
      ~PlanChebTransformFFTWXdct();
      void SX_To_X(const double* in,double* out); 
      void X_To_SX(const double* in,double* out); 
      void SpT_To_Re_2D(const double* in,double* out,int nyy=1); 
      void Re_To_SpT_2D(const double* in,double* out,int nyy=1); 
      void ReT_To_SpT_2D(const double* in,double* out,int nd2 = 1); 
      void SpT_To_ReT_2D(const double* in,double* out,int nd2 = 1);
    private:
      int sz;
      int nlogic;
      int N;
      double* uFFT;
      fftw_plan sp_fftw_plan;
      fftw_plan re_fftw_plan;
  };

  class PlanXr 
  {
    public: 
      PlanXr(int nxx);
      PlanXr()=delete;
      PlanXr(const PlanXr& sp)=delete;
      PlanXr& operator=(const PlanXr& sp)=delete;
      ~PlanXr();
      void SX_To_X(const double* in,double* out); 
      void X_To_SX(const double* in,double* out); 
      void SpT_To_Re_2D(const double* in,double* out,int nyy=1); 
      void Re_To_SpT_2D(const double* in,double* out,int nyy=1); 
    private:
      int sz;
      double* uFFT;
      fftw_plan sp_fftw_plan;
      fftw_plan re_fftw_plan;
  };

}


#endif


