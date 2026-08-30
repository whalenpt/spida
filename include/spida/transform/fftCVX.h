/*
 * Wrappers for Fourier transform applied on a UniformGridX grid.
 * Transforms are non-unitary and assume angular frequency (kx = 2pi/L where L is grid length)
 */

#pragma once

#include "spida/grid/uniformCVX.h"
#include "spida/helper/constants.h"

#include <vector>

// kiss_fft_cfg is kissfft's own opaque pointer typedef (struct
// kiss_fft_state*, never defined even by kiss_fft.h itself) — forward
// declared here so this public header doesn't leak a kissfft header into
// every consumer's translation unit. fftCVX.cpp includes the real
// kiss_fft.h for the definition it actually needs.
struct kiss_fft_state;
using kiss_fft_cfg = kiss_fft_state*;

namespace spida {

class FFTCVX {
public:
    explicit FFTCVX(const UniformGridCVX& grid);
    ~FFTCVX();
    FFTCVX() = delete;
    FFTCVX(const FFTCVX& sp) = delete;
    FFTCVX& operator=(const FFTCVX& sp) = delete;
    void X_To_SX(const std::vector<dcmplx>& in, std::vector<dcmplx>& out) noexcept;
    void SX_To_X(const std::vector<dcmplx>& in, std::vector<dcmplx>& out) noexcept;
    void X_To_SX(const dcmplx* in, dcmplx* out) noexcept;
    void SX_To_X(const dcmplx* in, dcmplx* out) noexcept;

private:
    unsigned m_nx;
    double m_minx;
    double m_L;
    std::vector<double> m_kx;
    std::vector<dcmplx> m_temp;
    kiss_fft_cfg m_cfg_forward;
    kiss_fft_cfg m_cfg_reverse;
};

} // namespace spida