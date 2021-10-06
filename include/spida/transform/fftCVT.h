// fftCVT.h
#pragma once

#include <vector>
#include <thread>
#include "spida/helper/constants.h"
#include "spida/grid/uniformCVT.h" 
#include "kiss_fft.h"

namespace spida{

class FFTCVT 
{
    public:
        explicit FFTCVT(const UniformGridCVT& grid);
        ~FFTCVT();
        FFTCVT()=delete;
        FFTCVT(const FFTCVT& sp)=delete;
        FFTCVT& operator=(const FFTCVT& sp)=delete;
        void T_To_ST(const dcmplx* in,dcmplx* out) noexcept;
        void ST_To_T(const dcmplx* in,dcmplx* out) noexcept;
        void T_To_ST(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) noexcept {\
            T_To_ST(in.data(),out.data());}
        void ST_To_T(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) noexcept {\
            ST_To_T(in.data(),out.data());}
    private:
        unsigned m_nt;
        double m_mint;
        std::vector<double> m_omega;
        std::vector<dcmplx> m_temp;

        kiss_fft_cfg m_cfg_forward;
        kiss_fft_cfg m_cfg_reverse;
};

}




