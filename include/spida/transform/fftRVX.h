/*
 * Wrappers for Fourier transform applied on a UniformGridRVX grid.
 * Transforms are non-unitary and assume angular frequency (kx = 2pi/L where L is grid length)
 */

#pragma once

#include "spida/grid/uniformRVX.h"
#include "spida/helper/constants.h"

#include <vector>

// kiss_fftr_cfg is kissfft's own opaque pointer typedef (struct
// kiss_fftr_state*, never defined even by kiss_fftr.h itself) — forward
// declared here so this public header doesn't leak a kissfft header into
// every consumer's translation unit. fftRVX.cpp includes the real
// kiss_fftr.h for the definition it actually needs.
struct kiss_fftr_state;
using kiss_fftr_cfg = kiss_fftr_state*;

namespace spida {

// interface class
class FFTRVX {
public:
    explicit FFTRVX(const UniformGridRVX& grid);
    ~FFTRVX();
    FFTRVX() = delete;
    FFTRVX(const FFTRVX& sp) = delete;
    FFTRVX& operator=(const FFTRVX& sp) = delete;

    void X_To_SX(const std::vector<double>& in, std::vector<dcmplx>& out) noexcept;
    void SX_To_X(const std::vector<dcmplx>& in, std::vector<double>& out) noexcept;
    void X_To_SX(const double* in, dcmplx* out) noexcept;
    void SX_To_X(const dcmplx* in, double* out) noexcept;

private:
    unsigned m_nx;
    double m_minx;
    double m_L;
    std::vector<double> m_kx;
    std::vector<dcmplx> m_temp;
    kiss_fftr_cfg m_rcfg_forward;
    kiss_fftr_cfg m_rcfg_reverse;
};

} // namespace spida