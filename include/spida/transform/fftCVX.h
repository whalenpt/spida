// fftCVX.h
/*
 * Wrappers for Fourier transform applied on a UniformGridX grid.
 * Transforms are non-unitary and assume angular frequency (kx = 2pi/L where L is grid length)
 */

#pragma once

#include <vector>
#include <thread>
#include "spida/helper/constants.h"
#include "spida/grid/uniformCVX.h" 
#include "kiss_fft.h"

namespace spida{

class FFTCVX 
{
    public:
        explicit FFTCVX(const UniformGridX& grid);
        ~FFTCVX();
        FFTCVX()=delete;
        FFTCVX(const FFTCVX& sp)=delete;
        FFTCVX& operator=(const FFTCVX& sp)=delete;
        void X_To_SX(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) noexcept; 
        void SX_To_X(const std::vector<dcmplx>& in,std::vector<dcmplx>& out) noexcept; 
        void X_To_SX(const dcmplx* in,dcmplx* out) noexcept; 
        void SX_To_X(const dcmplx* in,dcmplx* out) noexcept; 
    private:
        int m_nx;
        double m_minx;
        double m_L;
        std::vector<double> m_kx;
        std::vector<dcmplx> m_temp;
        kiss_fft_cfg m_cfg_forward;
        kiss_fft_cfg m_cfg_reverse;
};

}




