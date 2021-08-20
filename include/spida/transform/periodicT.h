
#ifndef SPIDA_TRANSFORMPERIODICT_H_
#define SPIDA_TRANSFORMPERIODICT_H_ 

#include <vector>
#include <thread>
#include "spida/constants.h"
#include "spida/transform/transformT.h"
#include "spida/grid/uniformT.h" 
#include "kiss_fft.h"
#include "kiss_fftr.h"

namespace spida{

class PeriodicTransformT : public TransformsT
{
    public:
        explicit PeriodicTransformT(const UniformGridT& grid);
        ~PeriodicTransformT();
        PeriodicTransformT()=delete;
        PeriodicTransformT(const PeriodicTransformT& sp)=delete;
        PeriodicTransformT& operator=(const PeriodicTransformT& sp)=delete;
        void T_To_ST(const std::vector<double>& in,std::vector<dcmplx>& out); 
        void T_To_ST(const double* in,dcmplx* out); 
        void ST_To_T(const std::vector<dcmplx>& in,std::vector<double>& out); 
        void ST_To_T(const dcmplx* in,double* out); 
        void ST_To_T_c(const std::vector<dcmplx>& in,std::vector<dcmplx>& out); 
        void T_To_ST_c(const std::vector<dcmplx>& in,std::vector<dcmplx>& out); 
    private:
        int m_nt;
        int m_nst;
        int m_minI;
        int m_maxI;
        std::vector<double> m_rFFTr;
        std::vector<dcmplx> m_rFFTs;
        std::vector<dcmplx> m_cFFT;

        void execute_forward();
        void execute_backward();

        kiss_fft_cfg m_cfg_forward;
        kiss_fft_cfg m_cfg_reverse;
        kiss_fftr_cfg m_rcfg_forward;
        kiss_fftr_cfg m_rcfg_reverse;
};

}


#endif


