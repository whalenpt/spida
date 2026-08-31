#pragma once

// Builds a runnable grid/model/propagator/solver quartet from a
// spida::config::SimulationConfig — the factory a future spida-worker would
// call after parsing a request, replacing the hand-assembly every demo/usage
// example does today (see e.g. demos/burgers.cpp, usage/main.cpp).
//
// Scope: only ModelKind::burgers + GridKind::uniform_rvx are wired end to
// end, matching the SPIDA Console architecture proposal's own recommended
// Phase 1 pilot ("prove the loop... one model only — Burgers or KdV, closest
// to an existing demo"). Extending to kdv/ks/nls means: add that model's
// params struct to simulationconfig.h, add its spida::models:: type (see
// spida/models/burgers.h for the pattern), and add a case in buildRun()
// below — the grid/solver factories are already model-agnostic.

#include <spida/config/simulationconfig.h>
#include <spida/grid/uniformRVX.h>
#include <spida/models/burgers.h>
#include <spida/propagator/propagator.h>
#include <spida/rkstiff/ETDAS.h>
#include <spida/rkstiff/IFAS.h>
#include <spida/rkstiff/solver.h>

#include <memory>
#include <stdexcept>
#include <utility>

namespace spida::config {

[[nodiscard]] inline std::unique_ptr<spida::SolverCV_AS>
buildSolver(const SolverConfig& cfg, const spida::LinOp& L, const spida::NLfunc& NL)
{
    std::unique_ptr<spida::SolverCV_AS> solver;
    switch (cfg.kind) {
    case SolverKind::etd35:
        solver = std::make_unique<spida::ETD35>(L, NL);
        break;
    case SolverKind::etd34:
        solver = std::make_unique<spida::ETD34>(L, NL);
        break;
    case SolverKind::if34:
        solver = std::make_unique<spida::IF34>(L, NL);
        break;
    case SolverKind::if45dp:
        solver = std::make_unique<spida::IF45DP>(L, NL);
        break;
    default:
        throw std::invalid_argument("buildSolver: unhandled SolverKind");
    }
    solver->setEpsRel(cfg.epsRel);
    solver->setLogProgress(cfg.logProgress);
    if (cfg.logProgress)
        solver->setLogFrequency(cfg.logFrequency);
    return solver;
}

/// Owns the grid/model/propagator/solver built from a SimulationConfig and
/// exposes the same evolve() call demos/usage examples make by hand. Only
/// the Burgers pilot (see file header) is wired; the constructor throws
/// std::invalid_argument — the config_validation failure category from the
/// SPIDA Console proposal's error taxonomy — for anything else.
class SimulationRun {
public:
    explicit SimulationRun(const SimulationConfig& cfg)
        : m_cfg(cfg), m_grid(cfg.grid.n, cfg.grid.a, cfg.grid.b), m_model(m_grid, cfg.burgersParams.mu)
    {
        if (cfg.model != ModelKind::burgers)
            throw std::invalid_argument(
                "SimulationRun: only ModelKind::burgers is currently wired "
                "(see simulationbuilder.h header comment)");
        if (cfg.grid.kind != GridKind::uniform_rvx)
            throw std::invalid_argument(
                "SimulationRun: only GridKind::uniform_rvx is currently wired");

        m_propagator = std::make_unique<spida::models::BurgersPropagator>(
            std::filesystem::path("sim_" + cfg.name), m_model);
        m_propagator->setStepsPerOutput1D(cfg.reporting.stepsPerOutput1D);
        m_propagator->setMaxReports1D(cfg.reporting.maxReports1D);
        m_propagator->setLogProgress(cfg.solver.logProgress);
        if (cfg.solver.logProgress)
            m_propagator->setLogFrequency(cfg.solver.logFrequency);

        m_solver = buildSolver(cfg.solver, m_model.L(), m_model.NL());
    }

    /// Same bool contract as spida::SolverCV::evolve(): false means a step
    /// itself failed at the solver level. Check propagator().stopReason()
    /// either way to see whether/why the run stopped early (cancelled,
    /// diverged, or hit its report budget — see spida::StopReason).
    [[nodiscard]] bool run()
    {
        return m_solver->evolve(*m_propagator, m_cfg.solver.t0, m_cfg.solver.tf, m_cfg.solver.hInit);
    }

    [[nodiscard]] spida::BasePropagator& propagator()
    {
        return *m_propagator;
    }

    [[nodiscard]] const spida::BasePropagator& propagator() const
    {
        return *m_propagator;
    }

private:
    SimulationConfig m_cfg;
    spida::UniformGridRVX m_grid;
    spida::models::Burgers m_model;
    std::unique_ptr<spida::models::BurgersPropagator> m_propagator;
    std::unique_ptr<spida::SolverCV_AS> m_solver;
};

} // namespace spida::config
