#pragma once

// Kuramoto-Sivashinsky (u_t + u u_x + u_xx + u_xxxx = 0), promoted from the
// spida-worker implementation originally authored in spida-console (see
// docs/adr/0002-worker-model-coverage-and-config-registry.md) — verified
// there against linear stability theory (growth/decay rate at two Fourier
// modes, both exact to 6 decimal places against the analytic L(k)) rather
// than a closed-form nonlinear solution, since KS has none. Mirrors
// spida/models/burgers.h's structure exactly.

#include "spida/RVX.h"
#include "spida/grid/uniformRVX.h"
#include "spida/helper/constants.h"
#include "spida/propagator/propagator.h"

#include <cmath>
#include <filesystem>
#include <functional>
#include <memory>
#include <vector>

#include <utils/report.hpp>

namespace spida::models {

/// u_t + u u_x + u_xx + u_xxxx = 0, periodic on [a, b] — standard
/// normalization, no free coefficient at all (unlike Burgers' mu or KdV's
/// implicit soliton speed). Rearranged to u_t = -u_xx - u_xxxx - u u_x and
/// transformed the same way as Burgers'/Kdv's L() (derivative -> (ik)^n):
///   L(k)  = k^2 - k^4   (-u_xx -> -((ik)^2)u_hat = k^2 u_hat;
///                        -u_xxxx -> -((ik)^4)u_hat = -k^4 u_hat)
///   NL(u) = -u * u_x    (same form as Burgers' NL)
/// L(k) = k^2(1-k^2) is positive (unstable, growing) for 0<k<1 — maximum
/// +0.25 at k=1/sqrt(2) — and negative (decaying) for k>1: the mechanism
/// behind KS's spatiotemporal chaos.
class Ks {
public:
    explicit Ks(const spida::UniformGridRVX& grid)
        : m_grid(grid),
          m_spi(grid),
          m_uphys(grid.getNx()),
          m_uxphys(grid.getNx()),
          m_uxsp(grid.getNsx()),
          m_L(grid.getNsx())
    {
        const auto& sx = grid.getSX();
        for (std::size_t i = 0; i < sx.size(); i++)
            m_L[i] = std::pow(sx[i], 2) - std::pow(sx[i], 4);
        m_NL = [this](const std::vector<spida::dcmplx>& in, std::vector<spida::dcmplx>& out) {
            m_spi.SX_To_X(in, m_uphys);
            m_spi.dSX(in, m_uxsp);
            m_spi.SX_To_X(m_uxsp, m_uxphys);
            for (std::size_t i = 0; i < m_grid.getNx(); i++)
                m_uphys[i] = -m_uphys[i] * m_uxphys[i];
            m_spi.X_To_SX(m_uphys, out);
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

    [[nodiscard]] spida::SpidaRVX& spida()
    {
        return m_spi;
    }

private:
    spida::UniformGridRVX m_grid;
    spida::SpidaRVX m_spi;
    std::vector<double> m_uphys;
    std::vector<double> m_uxphys;
    std::vector<spida::dcmplx> m_uxsp;
    std::vector<spida::dcmplx> m_L;
    std::function<void(const std::vector<spida::dcmplx>&, std::vector<spida::dcmplx>&)> m_NL;
};

/// Kassam-Trefethen reference IC: u(x,0) = cos(x/16)(1 + sin(x/16)) — the
/// standard citable KS demo (Kassam & Trefethen, "Fourth-Order Time-Stepping
/// for Stiff PDEs", 2005), used specifically because KS has no closed-form
/// nonlinear solution to build an exact IC from instead. Not parameterized —
/// this equation's standard normalization has no free coefficient, and
/// modelParams for "ks" is deliberately empty in SimulationConfig. Only
/// reproduces the literature reference pattern exactly on the reference
/// domain [0, 32*pi]; evaluated on whatever [a, b] the caller passes
/// otherwise, same as Burgers'/Kdv's fixed/parameterized ICs.
class KsPropagator : public spida::PropagatorCV {
public:
    KsPropagator(const std::filesystem::path& path, Ks& model)
        : PropagatorCV(path),
          m_spi(model.spida()),
          m_usp(model.spida().getGridX().getNsx(), 0.0),
          m_uphys(model.spida().getGridX().getNx(), 0.0)
    {
        const auto& x = m_spi.getX();
        for (std::size_t i = 0; i < x.size(); i++)
            m_uphys[i] = std::cos(x[i] / 16.0) * (1.0 + std::sin(x[i] / 16.0));
        m_spi.X_To_SX(m_uphys, m_usp);
        initReport();
    }

    ~KsPropagator() override = default;

    void updateFields(double) override
    {
        m_spi.SX_To_X(m_usp, m_uphys);
    }

    [[nodiscard]] std::vector<spida::dcmplx>& propagator() override
    {
        return m_usp;
    }

    bool checkDiverged(double) override
    {
        for (double v : m_uphys) {
            if (!std::isfinite(v))
                return true;
        }
        return false;
    }

private:
    void initReport()
    {
        const auto& x = m_spi.getGridX().getX();
        auto report = std::make_unique<spida::Report1D<double, double>>("X", x, m_uphys);
        PropagatorCV::addReport(std::move(report));
        const auto& sx = m_spi.getGridX().getSX();
        auto reportsp =
            std::make_unique<spida::ReportComplex1D<double, double>>("SX", sx, m_usp);
        PropagatorCV::addReport(std::move(reportsp));
    }

    spida::SpidaRVX& m_spi;
    std::vector<spida::dcmplx> m_usp;
    std::vector<double> m_uphys;
};

} // namespace spida::models
