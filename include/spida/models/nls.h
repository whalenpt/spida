#pragma once

// Radial nonlinear Schrödinger equation (cubic Kerr nonlinearity):
//   dz A = -i kr^2 A + i*gamma*|A|^2 A
// promoted from demos/NLSR.cpp — see docs/adr/0003-worker-relocation-and-
// cooperative-cancellation.md for the promotion pattern this follows
// exactly (spida/models/{burgers,kdv,ks}.h). Complex-valued physical field
// on a BesselRootGridR (Hankel transform, not FFT) — the first model wired
// that isn't real-valued-on-a-uniform-periodic-grid, exercising both
// PropagatorCV's complex-field path and GridKind::bessel_root_r for the
// first time (see docs/adr/0001-spida-console-backend-groundwork.md's
// Phase C addendum).
//
// Numerical verification: cubic NLS with a purely dispersive linear
// operator and a pointwise Kerr nonlinearity conserves the spectral-space
// L2 norm sum_k |A_k|^2 exactly for the CONTINUOUS equation — the linear
// step is a per-mode phase rotation (|A_k| invariant), and the nonlinear
// step dA/dz = i*gamma*|A|^2*A preserves |A| pointwise (d/dz|A|^2 =
// 2 Re(A* dA/dz) = 2 gamma |A|^4 Re(i) = 0). The discrete ETD scheme only
// approximates this; checked empirically (not just asserted) against this
// NlsR/NlsRPropagator pair — see test/config_tests.cpp's
// NLS_R_POWER_IS_CONSERVED for the actual measured error and how it was
// confirmed to be a shrinking time-integration effect, not a fixed
// transform-normalization offset.

#include "spida/R.h"
#include "spida/grid/besselR.h"
#include "spida/helper/constants.h"
#include "spida/propagator/propagator.h"

#include <cmath>
#include <complex>
#include <filesystem>
#include <functional>
#include <memory>
#include <vector>

#include <utils/report.hpp>

namespace spida::models {

/// dz A = -i kr^2 A + i*gamma*|A|^2 A. L(k) is purely dispersive
/// (imaginary); NL is a pointwise Kerr nonlinearity evaluated in physical
/// (r) space via the Hankel transform. `gamma` is the one PDE coefficient
/// exposed via modelParams (default 2.0, matching demos/NLSR.cpp).
class NlsR {
public:
    NlsR(const spida::BesselRootGridR& grid, double gamma)
        : m_spi(grid), m_uphys(grid.getNr()), m_L(grid.getNsr())
    {
        const auto& kr = grid.getSR();
        for (std::size_t i = 0; i < kr.size(); i++)
            m_L[i] = spida::dcmplx(0.0, -std::pow(kr[i], 2));
        m_NL = [this, gamma](const std::vector<spida::dcmplx>& in, std::vector<spida::dcmplx>& out) {
            m_spi.SR_To_R(in, m_uphys);
            for (auto& v : m_uphys)
                v = spida::dcmplx(0.0, gamma) * v * std::norm(v);
            m_spi.R_To_SR(m_uphys, out);
        };
    }

    [[nodiscard]] std::vector<spida::dcmplx>& L()
    {
        return m_L;
    }

    [[nodiscard]] std::function<void(const std::vector<spida::dcmplx>&, std::vector<spida::dcmplx>&)>&
    NL()
    {
        return m_NL;
    }

    [[nodiscard]] spida::SpidaR& spida()
    {
        return m_spi;
    }

private:
    spida::SpidaR m_spi;
    std::vector<spida::dcmplx> m_uphys;
    std::vector<spida::dcmplx> m_L;
    std::function<void(const std::vector<spida::dcmplx>&, std::vector<spida::dcmplx>&)> m_NL;
};

/// Gaussian initial condition A(r,0) = amplitude * exp(-r^2) — matches
/// demos/NLSR.cpp's profile (its `a` width parameter fixed at 1.0; only
/// amplitude is exposed via modelParams here, default 2.0 matching the
/// demo). Reports on the grid's own (non-uniform) r/kr coordinates, not a
/// synthesized linear axis — matching the proposal's "gridCoords never
/// assumed uniform" principle without needing any extra bookkeeping.
class NlsRPropagator : public spida::PropagatorCV {
public:
    NlsRPropagator(const std::filesystem::path& path, NlsR& model, double amplitude)
        : PropagatorCV(path),
          m_spi(model.spida()),
          m_usp(model.spida().getGridR().getNsr(), 0.0),
          m_uphys(model.spida().getGridR().getNr(), 0.0)
    {
        constexpr double width = 1.0; // matches demos/NLSR.cpp's `a`
        const auto& r = m_spi.getR();
        for (std::size_t i = 0; i < r.size(); i++)
            m_uphys[i] = amplitude * std::exp(-width * r[i] * r[i]);
        m_spi.R_To_SR(m_uphys, m_usp);
        initReport();
    }

    ~NlsRPropagator() override = default;

    void updateFields(double) override
    {
        m_spi.SR_To_R(m_usp, m_uphys);
    }

    [[nodiscard]] std::vector<spida::dcmplx>& propagator() override
    {
        return m_usp;
    }

    bool checkDiverged(double) override
    {
        for (const auto& v : m_uphys) {
            if (!std::isfinite(v.real()) || !std::isfinite(v.imag()))
                return true;
        }
        return false;
    }

private:
    void initReport()
    {
        const auto& r = m_spi.getGridR().getR();
        auto report = std::make_unique<spida::ReportComplex1D<double, double>>("R", r, m_uphys);
        PropagatorCV::addReport(std::move(report));
        const auto& sr = m_spi.getGridR().getSR();
        auto reportsp = std::make_unique<spida::ReportComplex1D<double, double>>("SR", sr, m_usp);
        PropagatorCV::addReport(std::move(reportsp));
    }

    spida::SpidaR& m_spi;
    std::vector<spida::dcmplx> m_usp;
    std::vector<spida::dcmplx> m_uphys;
};

} // namespace spida::models
