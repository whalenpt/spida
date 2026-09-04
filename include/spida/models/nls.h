#pragma once

// Radial nonlinear Schrödinger equation (cubic Kerr nonlinearity):
//   dz A = -i kr^2 A + i*gamma*|A|^2 A
// promoted from demos/NLSR.cpp — see docs/adr/0002-worker-model-coverage-
// and-config-registry.md for the promotion pattern this follows exactly
// (spida/models/{burgers,kdv,ks}.h). Complex-valued physical field on a
// BesselRootGridR (Hankel transform, not FFT) — the first model wired
// that isn't real-valued-on-a-uniform-periodic-grid, exercising both
// PropagatorCV's complex-field path and GridKind::bessel_root_r for the
// first time (see docs/adr/0002-worker-model-coverage-and-config-registry.md).
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

// Also holds NlsRt/NlsRtPropagator below (ModelKind::nls_rt) — the 2D
// radial + time/frequency cubic NLS, promoted from demos/NLSRT.cpp. Needs
// TWO independent grids (BesselRootGridR for r, UniformGridCVT for t) via
// SpidaRCVT — representing a second grid dimension is a SimulationConfig-
// level wire-contract question (see docs/adr/0002-worker-model-coverage-
// and-config-registry.md for how SimulationConfig.gridT resolved it).

#include "spida/R.h"
#include "spida/RCVT.h"
#include "spida/grid/besselR.h"
#include "spida/grid/uniformCVT.h"
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

/// 2D (r, t) cubic NLS: dz A = (-i kr^2 + i*sigma*omega^2) A + i*gamma*|A|^2 A,
/// sigma fixed at 0.5 (matching demos/NLSRT.cpp, no free coefficient exposed
/// for it — same precedent as Ks). `gamma` is the Kerr coefficient, exposed
/// via modelParams like NlsR's. L(k) is purely dispersive (imaginary) here
/// too, so the same spectral-L2-norm-conservation verification argument
/// NlsR's header comment gives applies unchanged.
class NlsRt {
public:
    NlsRt(const spida::BesselRootGridR& gridR,
          const spida::UniformGridCVT& gridT,
          double gamma,
          unsigned threads = 1)
        : m_spi(gridR, gridT, threads),
          m_uphys(gridR.getNr() * gridT.getNt()),
          m_L(gridR.getNsr() * gridT.getNst())
    {
        constexpr double sigma = 0.5; // matches demos/NLSRT.cpp
        const auto& kr = gridR.getSR();
        const auto& omega = gridT.getST();
        for (std::size_t i = 0; i < kr.size(); i++)
            for (std::size_t j = 0; j < omega.size(); j++)
                m_L[i * omega.size() + j] =
                    spida::dcmplx(0.0, -std::pow(kr[i], 2) + sigma * std::pow(omega[j], 2));
        m_NL = [this, gamma](const std::vector<spida::dcmplx>& in, std::vector<spida::dcmplx>& out) {
            m_spi.SRST_To_RT(in, m_uphys);
            for (auto& v : m_uphys)
                v = spida::dcmplx(0.0, gamma) * v * std::norm(v);
            m_spi.RT_To_SRST(m_uphys, out);
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

    [[nodiscard]] spida::SpidaRCVT& spida()
    {
        return m_spi;
    }

private:
    spida::SpidaRCVT m_spi;
    std::vector<spida::dcmplx> m_uphys;
    std::vector<spida::dcmplx> m_L;
    std::function<void(const std::vector<spida::dcmplx>&, std::vector<spida::dcmplx>&)> m_NL;
};

/// Gaussian pulse A(r,t,0) = amplitude * exp(-(r/w0)^2 - (t/tp)^2), w0/tp
/// fixed at 1.0/0.5 (matching demos/NLSRT.cpp); only amplitude is exposed
/// via modelParams (default 4.0, matching the demo). Reports the model's
/// first 2D data — "RT" (physical field) and "SR" (spectral) as
/// ReportComplex2D — on the grids' own (non-uniform in r) coordinates
/// directly, no mirroring/downsampling: demos/NLSRT.cpp mirrors r and
/// strides both axes purely for nicer plots, which is presentation
/// concern, not something a config-driven backend model should bake in
/// (see NlsRPropagator's header comment for the same reasoning re:
/// mirroring). This is the first model to actually exercise
/// ReportingConfig's stepsPerOutput2D/maxReports2D fields and the worker
/// manifest's "field2d" classification, both already in place before any
/// 2D model existed — see this file's own header comment and
/// docs/adr/0002-worker-model-coverage-and-config-registry.md.
class NlsRtPropagator : public spida::PropagatorCV {
public:
    NlsRtPropagator(const std::filesystem::path& path, NlsRt& model, double amplitude)
        : PropagatorCV(path),
          m_spi(model.spida()),
          m_usp(model.spida().spectralSize(), 0.0),
          m_uphys(model.spida().physicalSize(), 0.0)
    {
        constexpr double w0 = 1.0; // matches demos/NLSRT.cpp
        constexpr double tp = 0.5; // matches demos/NLSRT.cpp
        const auto& r = m_spi.getGridR().getR();
        const auto& t = m_spi.getGridT().getT();
        for (std::size_t i = 0; i < r.size(); i++)
            for (std::size_t j = 0; j < t.size(); j++)
                m_uphys[i * t.size() + j] =
                    amplitude * std::exp(-std::pow(r[i] / w0, 2) - std::pow(t[j] / tp, 2));
        m_spi.RT_To_SRST(m_uphys, m_usp);
        initReport();
    }

    ~NlsRtPropagator() override = default;

    void updateFields(double) override
    {
        m_spi.SRST_To_RT(m_usp, m_uphys);
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
        const auto& t = m_spi.getGridT().getT();
        auto report =
            std::make_unique<spida::ReportComplex2D<double, double, double>>("RT", r, t, m_uphys);
        PropagatorCV::addReport(std::move(report));
        const auto& kr = m_spi.getGridR().getSR();
        const auto& omega = m_spi.getGridT().getST();
        auto reportsp = std::make_unique<spida::ReportComplex2D<double, double, double>>(
            "SR", kr, omega, m_usp);
        PropagatorCV::addReport(std::move(reportsp));
    }

    spida::SpidaRCVT& m_spi;
    std::vector<spida::dcmplx> m_usp;
    std::vector<spida::dcmplx> m_uphys;
};

} // namespace spida::models
