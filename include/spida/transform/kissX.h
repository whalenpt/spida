

#ifndef KISSX_H_
#define KISSX_H_ 

#include <complex>
#include <vector>
#include "kiss_fft.h"
#include "kiss_fftr.h"


namespace spida{


  class PlanXr;
  class PlanChebTransformXdct;

  class ChebTransformX 
  {
    public:
        ChebTransformX(int cnx);
        void X_To_SX(const std::vector<double>& in,std::vector<double>& out);
        void SX_To_X(const std::vector<double>& in,std::vector<double>& out);
        ~ChebTransformX();
    private:
        PlanChebTransformXdct* m_plan; 
  };

  class PlanChebTransformXdct 
  {
    public: 
      PlanChebTransformXdct(int nxx);
      PlanChebTransformXdct()=delete;
      PlanChebTransformXdct(const PlanXr& sp)=delete;
      PlanChebTransformXdct& operator=(const PlanXr& sp)=delete;
      ~PlanChebTransformXdct();
      void SX_To_X(const std::vector<double>& in,std::vector<double>& out); 
      void X_To_SX(const std::vector<double>& in,std::vector<double>& out); 
      void SpT_To_Re_2D(const std::vector<double>& in,std::vector<double>& out,int nyy=1); 
      void Re_To_SpT_2D(const std::vector<double>& in,std::vector<double>& out,int nyy=1); 
      void ReT_To_SpT_2D(const std::vector<double>& in,std::vector<double>& out,int nd2 = 1); 
      void SpT_To_ReT_2D(const std::vector<double>& in,std::vector<double>& out,int nd2 = 1);
    private:
      int m_nx;
      int m_nlogic;
      int m_N;
      std::vector<double> m_uFFT;
      std::vector<dcmplx> m_uFFTs;
      kiss_fftr_cfg m_rcfg_forward;
      kiss_fftr_cfg m_rcfg_reverse;
  };

  class PlanXr 
  {
    public: 
      PlanXr(int nxx);
      PlanXr()=delete;
      PlanXr(const PlanXr& sp)=delete;
      PlanXr& operator=(const PlanXr& sp)=delete;
      ~PlanXr();
      void SX_To_X(const std::vector<double>& in,std::vector<double>& out); 
      void X_To_SX(const std::vector<double>& in,std::vector<double>& out); 
      void SpT_To_Re_2D(const std::vector<double>& in,std::vector<double>& out,int nyy=1); 
      void Re_To_SpT_2D(const std::vector<double>& in,std::vector<double>& out,int nyy=1); 
    private:
      int m_nx;
      std::vector<double> m_uFFT;
      std::vector<dcmplx> m_uFFTs;
      kiss_fftr_cfg m_rcfg_forward;
      kiss_fftr_cfg m_rcfg_reverse;
  };

}


#endif


