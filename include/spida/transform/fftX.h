// fftX.h
#pragma once

#include <vector>
#include <thread>
#include "spida/helper/constants.h"
#include "spida/grid/uniformX.h" 
#include "kiss_fft.h"

namespace spida{

class FFTX
{
    public:
        explicit FFTX(const UniformGridX& grid);
        ~FFTX();
        FFTX()=delete;
        FFTX(const FFTX& sp)=delete;
        FFTX& operator=(const FFTX& sp)=delete;

        void X_To_SX(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) noexcept; 
        void SX_To_X(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) noexcept; 
        void X_To_SX(const dcmplx* in,dcmplx* out) noexcept; 
        void SX_To_X(const dcmplx* in,dcmplx* out) noexcept; 
    private:
        int m_nx;
        double m_minx;
        double m_L;
        kiss_fft_cfg m_cfg_forward;
        kiss_fft_cfg m_cfg_reverse;
};

}




