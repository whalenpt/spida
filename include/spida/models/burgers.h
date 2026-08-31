#pragma once

// Reusable, parameterized version of the Burgers' equation model demos/burgers.cpp
// hand-wires with mu hardcoded in main(). Promoted to a public header so
// spida/config/simulationbuilder.h (and, eventually, other consumers such as a
// spida-console worker linking SPIDA::spida directly) can construct it from data
// instead of copy-pasting the demo's classes. demos/burgers.cpp itself is left
// untouched to avoid any regression risk to a working, tested demo.

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

/// Burgers' equation: du/dt = -u*du/dx + mu*d^2u/dx^2, real-valued physical
/// field over a periodic UniformGridRVX. `mu` is the one physics parameter
/// demos/burgers.cpp hardcodes (0.0005) — here it's a constructor argument so
/// it can come from a SimulationConfig.
class Burgers {
public:
    Burgers(const spida::UniformGridRVX& grid, double mu)
        : m_grid(grid),
          m_spi(grid),
          m_uphys(grid.getNx()),
          m_uxphys(grid.getNx()),
          m_uxsp(grid.getNsx()),
          m_L(grid.getNsx())
    {
        const auto& sx = grid.getSX();
        for (std::size_t i = 0; i < sx.size(); i++)
            m_L[i] = -mu * std::pow(sx[i], 2);
        m_NL = [this](const std::vector<spida::dcmplx>& in, std::vector<spida::dcmplx>& out) {
            m_spi.SX_To_X(in, m_uphys);
            // dSX takes derivative in spectral space, i.e. multiply by i*kx
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

/// Propagator pairing a Burgers model with the same initial condition and
/// "X"/"SX" report definitions demos/burgers.cpp registers by hand, plus a
/// checkDiverged() override (see spida::BasePropagator) that flags a
/// non-finite physical-space value — a concrete example of the divergence
/// hook added for the future job-service backend's error taxonomy.
class BurgersPropagator : public spida::PropagatorCV {
public:
    BurgersPropagator(const std::filesystem::path& path, Burgers& model)
        : PropagatorCV(path),
          m_spi(model.spida()),
          m_usp(model.spida().getGridX().getNsx(), 0.0),
          m_uphys(model.spida().getGridX().getNx(), 0.0)
    {
        const auto& x = m_spi.getX();
        for (std::size_t i = 0; i < x.size(); i++)
            m_uphys[i] = std::exp(-10.0 * std::pow(std::sin(x[i] / 2.0), 2));
        m_spi.X_To_SX(m_uphys, m_usp);
        initReport();
    }

    ~BurgersPropagator() override = default;

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
