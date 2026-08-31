#pragma once

// Builds a runnable grid/model/propagator/solver quartet from a
// spida::config::SimulationConfig — the factory spida-worker (worker/) calls
// after parsing a job's config.json, replacing the hand-assembly every
// demo/usage example still does. See docs/adr/0003-worker-relocation-and-
// cooperative-cancellation.md for how this came to cover three real models
// instead of just the original Burgers pilot, and docs/adr/0001-spida-
// console-backend-groundwork.md's Phase C addendum for nls_r/bessel_root_r.
//
// Scope: burgers/kdv_rv/ks (GridKind::uniform_rvx), kdv_cv
// (GridKind::uniform_cvx), and nls_r (GridKind::bessel_root_r) are wired
// end to end — see spida/models/{burgers,kdv,ks,nls}.h's header comments
// for what was checked. nls_rt remains in ModelKind for wire-format
// stability but is not implemented — it needs a second, independent grid
// dimension SimulationConfig can't represent yet (SpidaRCVT combines a
// BesselRootGridR and a UniformGridCVT), a bigger wire-contract change
// than a per-model case here; validate() (validation.h) rejects a nls_rt
// request before this class is ever constructed (config_validation).
//
// The grid itself is a per-case LOCAL variable in the constructor below,
// not a member: Burgers/Kdv/KdvCv/Ks/NlsR all copy the grid into their own
// internal storage (see e.g. spida/models/kdv.h's `m_grid(grid)`), so
// nothing needs it to outlive construction — which is what lets each
// ModelKind's case use whichever concrete grid type it actually needs
// (UniformGridRVX / UniformGridCVX / BesselRootGridR) without a grid-side
// variant.

#include <spida/config/simulationconfig.h>
#include <spida/config/validation.h>
#include <spida/grid/besselR.h>
#include <spida/grid/uniformCVX.h>
#include <spida/grid/uniformRVX.h>
#include <spida/models/burgers.h>
#include <spida/models/kdv.h>
#include <spida/models/ks.h>
#include <spida/models/nls.h>
#include <spida/propagator/propagator.h>
#include <spida/rkstiff/ETDAS.h>
#include <spida/rkstiff/IFAS.h>
#include <spida/rkstiff/solver.h>

#include <filesystem>
#include <memory>
#include <stdexcept>
#include <string>
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
    ///
    /// Validates cfg via spida::config::validate() before touching the
    /// grid/model/solver at all (see the delegating constructor below) —
    /// this class stays a caller of validate(), not a second, possibly
    /// drifting copy of what makes a config valid. Throws
    /// std::invalid_argument (joining every field/message pair) if cfg is
    /// invalid, same contract as before this split.
    explicit SimulationRun(const SimulationConfig& cfg, std::filesystem::path outDir = {})
        : SimulationRun(cfg, std::move(outDir), requireValid(cfg))
    {
    }

    /// Same bool contract as spida::SolverCV::evolve(): false means a step
    /// itself failed at the solver level. Check propagator().stopReason()
    /// either way to see whether/why the run stopped early (cancelled,
    /// diverged, or hit its report budget — see spida::StopReason).
    [[nodiscard]] bool run()
    {
        return m_solver->evolve(*m_propagator, m_cfg.solver.t0, m_cfg.solver.tf, m_cfg.solver.hInit);
    }

    /// Returns PropagatorCV&, not just BasePropagator& — m_propagator is
    /// always a PropagatorCV in practice (every wired model's Propagator
    /// class derives from it), so this is a strictly more informative type
    /// than what every caller so far has actually needed (stopReason(),
    /// requestCancel(), etc., all declared on BasePropagator). Widened so a
    /// caller that DOES need the raw field (PropagatorCV::propagator(),
    /// e.g. a numerical-verification check reading the propagated
    /// spectral-space array directly) doesn't have to downcast.
    [[nodiscard]] spida::PropagatorCV& propagator()
    {
        return *m_propagator;
    }

    [[nodiscard]] const spida::PropagatorCV& propagator() const
    {
        return *m_propagator;
    }

private:
    // Private tag type — only requireValid() can produce one, so the
    // three-argument constructor below can only run after validate(cfg) has
    // already passed. Argument expressions (requireValid(cfg) here) are
    // evaluated before the delegated-to constructor's member-initializer
    // list runs, so this guarantees validation happens before the
    // constructor body (and the per-case grid it constructs) even starts.
    struct Validated {};

    [[nodiscard]] static Validated requireValid(const SimulationConfig& cfg)
    {
        auto errors = validate(cfg);
        if (!errors.empty()) {
            std::string msg =
                "SimulationRun: invalid config (" + std::to_string(errors.size()) + " error(s)): ";
            for (std::size_t i = 0; i < errors.size(); ++i) {
                if (i != 0)
                    msg += "; ";
                msg += errors[i].field + ": " + errors[i].message;
            }
            throw std::invalid_argument(msg);
        }
        return {};
    }

    SimulationRun(const SimulationConfig& cfg, std::filesystem::path outDir, Validated)
        : m_cfg(cfg)
    {
        if (outDir.empty())
            outDir = std::filesystem::path("sim_" + cfg.name);

        // model/grid.kind (and their pairing) are already guaranteed valid
        // by requireValid() above — the switch's default case below is
        // unreachable in practice, kept only as a defensive fallback.
        switch (cfg.model) {
        case ModelKind::burgers: {
            spida::UniformGridRVX grid(cfg.grid.n, cfg.grid.a, cfg.grid.b);
            const double mu = cfg.modelParams.value("mu", 0.0005);
            auto& model = m_model.emplace<spida::models::Burgers>(grid, mu);
            m_propagator = std::make_unique<spida::models::BurgersPropagator>(outDir, model);
            m_solver = buildSolver(cfg.solver, model.L(), model.NL());
            break;
        }
        case ModelKind::kdv_rv: {
            spida::UniformGridRVX grid(cfg.grid.n, cfg.grid.a, cfg.grid.b);
            const double solitonSpeed = cfg.modelParams.value("solitonSpeed", 1.0);
            auto& model = m_model.emplace<spida::models::Kdv>(grid);
            m_propagator =
                std::make_unique<spida::models::KdvPropagator>(outDir, model, solitonSpeed);
            m_solver = buildSolver(cfg.solver, model.L(), model.NL());
            break;
        }
        case ModelKind::ks: {
            // No modelParams read here — ks's standard normalization has no
            // free coefficient (see spida/models/ks.h's header comment).
            spida::UniformGridRVX grid(cfg.grid.n, cfg.grid.a, cfg.grid.b);
            auto& model = m_model.emplace<spida::models::Ks>(grid);
            m_propagator = std::make_unique<spida::models::KsPropagator>(outDir, model);
            m_solver = buildSolver(cfg.solver, model.L(), model.NL());
            break;
        }
        case ModelKind::kdv_cv: {
            // No modelParams read here — the 5-soliton initial condition
            // is fixed, matching demos/kdvCV.cpp (see
            // spida/models/kdv.h's KdvCvPropagator header comment).
            spida::UniformGridCVX grid(cfg.grid.n, cfg.grid.a, cfg.grid.b);
            auto& model = m_model.emplace<spida::models::KdvCv>(grid);
            m_propagator = std::make_unique<spida::models::KdvCvPropagator>(outDir, model);
            m_solver = buildSolver(cfg.solver, model.L(), model.NL());
            break;
        }
        case ModelKind::nls_r: {
            spida::BesselRootGridR grid(cfg.grid.n, cfg.grid.rMax);
            const double gamma = cfg.modelParams.value("gamma", 2.0);
            const double amplitude = cfg.modelParams.value("amplitude", 2.0);
            auto& model = m_model.emplace<spida::models::NlsR>(grid, gamma);
            m_propagator =
                std::make_unique<spida::models::NlsRPropagator>(outDir, model, amplitude);
            m_solver = buildSolver(cfg.solver, model.L(), model.NL());
            break;
        }
        default:
            throw std::invalid_argument(
                "SimulationRun: ModelKind not yet wired (only burgers/kdv_rv/ks/kdv_cv/nls_r "
                "are) — see simulationbuilder.h's header comment");
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

    SimulationConfig m_cfg;
    std::variant<std::monostate,
                spida::models::Burgers,
                spida::models::Kdv,
                spida::models::Ks,
                spida::models::KdvCv,
                spida::models::NlsR>
        m_model;
    std::unique_ptr<spida::PropagatorCV> m_propagator;
    std::unique_ptr<spida::SolverCV_AS> m_solver;
};

} // namespace spida::config
