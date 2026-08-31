#pragma once

// Builds a runnable grid/model/propagator/solver quartet from a
// spida::config::SimulationConfig — the factory spida-worker (worker/) calls
// after parsing a job's config.json, replacing the hand-assembly every
// demo/usage example still does. See docs/adr/0003-worker-relocation-and-
// cooperative-cancellation.md for how this came to cover three real models
// instead of just the original Burgers pilot.
//
// Scope: burgers, kdv_rv, and ks are wired end to end — the three models
// spida-worker already implemented and numerically verified before this
// factory existed (see spida/models/{burgers,kdv,ks}.h's header comments for
// what was checked). kdv_cv/nls_r/nls_rt remain in ModelKind for wire-format
// stability but are not implemented; the constructor throws
// std::invalid_argument for them (config_validation, in the job-service
// error taxonomy). Only GridKind::uniform_rvx is wired — every model here
// happens to use it, not a real per-model choice yet.

#include <spida/config/simulationconfig.h>
#include <spida/grid/uniformRVX.h>
#include <spida/models/burgers.h>
#include <spida/models/kdv.h>
#include <spida/models/ks.h>
#include <spida/propagator/propagator.h>
#include <spida/rkstiff/ETDAS.h>
#include <spida/rkstiff/IFAS.h>
#include <spida/rkstiff/solver.h>

#include <filesystem>
#include <memory>
#include <stdexcept>
#include <variant>

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
    return solver;
}

/// Owns the grid/model/propagator/solver built from a SimulationConfig and
/// exposes the same evolve() call demos/usage examples make by hand. The
/// model object's concrete type varies by ModelKind (spida::models::Burgers/
/// Kdv/Ks each hold their own L()/NL() storage, and BurgersPropagator/etc.
/// hold a live reference to the model's SpidaRVX transform — so the model
/// must outlive the propagator/solver, hence the variant member rather than
/// a temporary local).
class SimulationRun {
public:
    /// outDir: where report files land (BasePropagator's dir_path). Defaults
    /// to "sim_<name>" (relative to CWD) for demo/test convenience when the
    /// caller doesn't care — a real caller managing job directories itself
    /// (e.g. worker/src/main.cpp, given <output-dir> on argv) should always
    /// pass one explicitly. Passing a default-constructed empty path is
    /// equivalent to omitting the argument.
    explicit SimulationRun(const SimulationConfig& cfg, std::filesystem::path outDir = {})
        : m_cfg(cfg), m_grid(cfg.grid.n, cfg.grid.a, cfg.grid.b)
    {
        if (outDir.empty())
            outDir = std::filesystem::path("sim_" + cfg.name);

        if (cfg.grid.kind != GridKind::uniform_rvx)
            throw std::invalid_argument(
                "SimulationRun: only GridKind::uniform_rvx is currently wired");

        switch (cfg.model) {
        case ModelKind::burgers: {
            const double mu = cfg.modelParams.value("mu", 0.0005);
            auto& model = m_model.emplace<spida::models::Burgers>(m_grid, mu);
            m_propagator = std::make_unique<spida::models::BurgersPropagator>(outDir, model);
            m_solver = buildSolver(cfg.solver, model.L(), model.NL());
            break;
        }
        case ModelKind::kdv_rv: {
            const double solitonSpeed = cfg.modelParams.value("solitonSpeed", 1.0);
            auto& model = m_model.emplace<spida::models::Kdv>(m_grid);
            m_propagator =
                std::make_unique<spida::models::KdvPropagator>(outDir, model, solitonSpeed);
            m_solver = buildSolver(cfg.solver, model.L(), model.NL());
            break;
        }
        case ModelKind::ks: {
            // No modelParams read here — ks's standard normalization has no
            // free coefficient (see spida/models/ks.h's header comment).
            auto& model = m_model.emplace<spida::models::Ks>(m_grid);
            m_propagator = std::make_unique<spida::models::KsPropagator>(outDir, model);
            m_solver = buildSolver(cfg.solver, model.L(), model.NL());
            break;
        }
        default:
            throw std::invalid_argument(
                "SimulationRun: ModelKind not yet wired (only burgers/kdv_rv/ks are) — "
                "see simulationbuilder.h's header comment");
        }

        m_propagator->setStepsPerOutput1D(cfg.reporting.stepsPerOutput1D);
        m_propagator->setMaxReports1D(cfg.reporting.maxReports1D);
        m_propagator->setLogProgress(cfg.solver.logProgress);
        if (cfg.solver.logProgress)
            m_propagator->setLogFrequency(cfg.reporting.logFrequency);
        // cfg.solver.tf is known up front (it's part of the submitted config,
        // not something evolve() discovers) — set it so every
        // ProgressSnapshot a caller observes via setProgressObserver() has
        // tf populated, not just t. Previously left unset: nothing here ever
        // called setFinalTime(), so every progress event's tf was null even
        // though the value was available the whole time.
        m_propagator->setFinalTime(cfg.solver.tf);
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
    std::variant<std::monostate, spida::models::Burgers, spida::models::Kdv, spida::models::Ks> m_model;
    std::unique_ptr<spida::PropagatorCV> m_propagator;
    std::unique_ptr<spida::SolverCV_AS> m_solver;
};

} // namespace spida::config
