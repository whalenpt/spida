// fftRVX.h
/*
 * Wrappers for Fourier transform applied on a UniformGridRVX grid.
 * Transforms are non-unitary and assume angular frequency (kx = 2pi/L where L is grid length)
 */

#pragma once

#include <vector>
#include <thread>
#include "spida/helper/constants.h"
#include "spida/grid/uniformRVX.h" 
#include "kiss_fftr.h"

namespace spida{

// interface class
class FFTRVX
{
    public:
        explicit FFTRVX(const UniformGridRVX& grid);
        ~FFTRVX();
        FFTRVX()=delete;
        FFTRVX(const FFTRVX& sp)=delete;
        FFTRVX& operator=(const FFTRVX& sp)=delete;

        void X_To_SX(const std::vector<double>& in,std::vector<dcmplx>& out) noexcept; 
        void SX_To_X(const std::vector<dcmplx>& in,std::vector<double>& out) noexcept; 
        void X_To_SX(const double* in,dcmplx* out) noexcept; 
        void SX_To_X(const dcmplx* in,double* out) noexcept; 
    private:
        unsigned m_nx;
        double m_minx;
        double m_L;
        std::vector<double> m_kx;
        std::vector<dcmplx> m_temp;
        kiss_fftr_cfg m_rcfg_forward;
        kiss_fftr_cfg m_rcfg_reverse;
};

}




