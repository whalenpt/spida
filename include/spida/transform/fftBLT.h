// fftBLT.h
#pragma once

#include <vector>
#include <thread>
#include "spida/helper/constants.h"
#include "spida/grid/uniformRVT.h" 
#include "kiss_fft.h"
#include "kiss_fftr.h"

namespace spida{

class FFTBLT 
{
    public:
        explicit FFTBLT(const UniformGridRVT& grid);
        ~FFTBLT();
        FFTBLT()=delete;
        FFTBLT(const FFTBLT& sp)=delete;
        FFTBLT& operator=(const FFTBLT& sp)=delete;

        void T_To_ST(const std::vector<double>& in,std::vector<dcmplx>& out); 
        void ST_To_T(const std::vector<dcmplx>& in,std::vector<double>& out); 
        void ST_To_CVT(const std::vector<dcmplx>& in,std::vector<dcmplx>& out); 
        void CVT_To_ST(const std::vector<dcmplx>& in,std::vector<dcmplx>& out); 

        void T_To_ST(const double* in,dcmplx* out); 
        void ST_To_T(const dcmplx* in,double* out); 
        void ST_To_CVT(const dcmplx* in,dcmplx* out); 
        void CVT_To_ST(const dcmplx* in,dcmplx* out); 
    private:
        unsigned m_nt;
        unsigned m_nst;
        unsigned m_minI;
        unsigned m_maxI;
        std::vector<double> m_rFFTr;
        std::vector<dcmplx> m_rFFTs;
        std::vector<dcmplx> m_cFFT;
        std::vector<double> m_omega;
        double m_mint;
        double m_L;

        kiss_fft_cfg m_cfg_forward;
        kiss_fft_cfg m_cfg_reverse;
        kiss_fftr_cfg m_rcfg_forward;
        kiss_fftr_cfg m_rcfg_reverse;
};

}




