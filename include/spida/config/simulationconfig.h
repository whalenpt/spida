#pragma once

// JSON-serializable configuration mirroring the SPIDA Console architecture
// proposal's domain.ts SimulationConfig (discriminated by ModelKind/GridKind/
// SolverKind) — the C++-side source of truth a future spida-worker would parse
// a request into, and demos/usage examples build grids/models/solvers from
// today by hand. burgers/kdv_rv/ks (GridKind::uniform_rvx) and nls_r
// (GridKind::bessel_root_r) are wired end to end by simulationbuilder.h so
// far (see its header comment); the enums below are intentionally larger
// than that to keep the on-disk/wire shape stable as more models and grids
// come online.

#include <nlohmann/json.hpp>

#include <string>

namespace spida::config {

enum class ModelKind { burgers, kdv_cv, kdv_rv, ks, nls_r, nls_rt };
NLOHMANN_JSON_SERIALIZE_ENUM(ModelKind,
                             {
                                 {ModelKind::burgers, "burgers"},
                                 {ModelKind::kdv_cv, "kdv_cv"},
                                 {ModelKind::kdv_rv, "kdv_rv"},
                                 {ModelKind::ks, "ks"},
                                 {ModelKind::nls_r, "nls_r"},
                                 {ModelKind::nls_rt, "nls_rt"},
                             })

// uniform_cvx/uniform_cvt are NOT in the architecture proposal's own
// domain.ts GridKind union — added because kdv_cv and nls_rt (both IN the
// proposal's ModelKind list) genuinely need grid types the proposal never
// enumerated: uniform_cvx is a full-complex-FFT grid distinct from
// kdv_rv's real-optimized uniform_rvx transform; uniform_cvt is the
// time/frequency dimension nls_rt needs ALONGSIDE bessel_root_r (see
// SimulationConfig.gridT below). See docs/adr/0001-spida-console-backend-
// groundwork.md's Phase D/E addenda.
enum class GridKind { uniform_rvx, uniform_rvt, uniform_cvx, uniform_cvt, bessel_root_r, cheb_x };
NLOHMANN_JSON_SERIALIZE_ENUM(GridKind,
                             {
                                 {GridKind::uniform_rvx, "uniform_rvx"},
                                 {GridKind::uniform_rvt, "uniform_rvt"},
                                 {GridKind::uniform_cvx, "uniform_cvx"},
                                 {GridKind::uniform_cvt, "uniform_cvt"},
                                 {GridKind::bessel_root_r, "bessel_root_r"},
                                 {GridKind::cheb_x, "cheb_x"},
                             })

enum class SolverKind { etd35, etd34, if34, if45dp };
NLOHMANN_JSON_SERIALIZE_ENUM(SolverKind,
                             {
                                 {SolverKind::etd35, "etd35"},
                                 {SolverKind::etd34, "etd34"},
                                 {SolverKind::if34, "if34"},
                                 {SolverKind::if45dp, "if45dp"},
                             })

/// n/a/b are UniformGridRVX/UniformGridCVX/UniformGridCVT-shaped (all
/// three take the same (n, min, max) constructor signature); rMax is
/// BesselRootGridR-shaped (only meaningful when kind == bessel_root_r —
/// see simulationbuilder.h's nls_r case). cheb_x's own params remain left
/// for whoever wires it next.
struct GridConfig {
    GridKind kind{GridKind::uniform_rvx};
    unsigned n{256};
    double a{0.0};
    double b{1.0};
    double rMax{5.0};
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE_WITH_DEFAULT(GridConfig, kind, n, a, b, rMax)

struct SolverConfig {
    SolverKind kind{SolverKind::etd35};
    double epsRel{1e-8};
    double t0{0.0};
    double tf{1.0};
    double hInit{0.01};
    bool logProgress{false};
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE_WITH_DEFAULT(SolverConfig, kind, epsRel, t0, tf, hInit, logProgress)

/// Mirrors BasePropagator's setters (include/spida/propagator/propagator.h)
/// and the proposal's own domain.ts ReportingConfig field-for-field.
/// logFrequency lives here, not on SolverConfig, matching the real worker
/// wire contract (spida-console's apps/web sends it nested under
/// "reporting", not "solver" — see ADR-0003).
///
/// stepsPerOutput2D/maxReports2D/stepsPerOutputTrack are defaulted and
/// currently unused by simulationbuilder.h — no model wired today
/// (burgers/kdv_rv/ks) reports 2D or track data. Kept present rather than
/// added later specifically to freeze the wire shape now: ADR-0003 already
/// had to realign this struct's shape once (modelParams, logFrequency's
/// nesting) when the real worker/frontend contract turned out to differ
/// from what was speculatively designed — adding these fields after a 2D
/// model exists would risk the same kind of break a second time.
struct ReportingConfig {
    unsigned stepsPerOutput1D{5};
    unsigned stepsPerOutput2D{1};
    unsigned stepsPerOutputTrack{1};
    unsigned maxReports1D{500};
    unsigned maxReports2D{200};
    unsigned logFrequency{200};
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE_WITH_DEFAULT(ReportingConfig,
                                                stepsPerOutput1D,
                                                stepsPerOutput2D,
                                                stepsPerOutputTrack,
                                                maxReports1D,
                                                maxReports2D,
                                                logFrequency)

struct SimulationConfig {
    std::string name{"run"};
    ModelKind model{ModelKind::burgers};
    GridConfig grid;
    /// A second, independent grid dimension — needed only by nls_rt
    /// (SpidaRCVT combines a BesselRootGridR from `grid` with a
    /// UniformGridCVT from `gridT`; see simulationbuilder.h's nls_rt
    /// case). Absent/default for every other model, which is why this is
    /// a sibling field rather than turning `grid` itself into an array —
    /// keeps `grid`'s existing single-dimension meaning untouched for
    /// every model that doesn't need a second one. Not in the
    /// architecture proposal's own domain.ts at all; see docs/adr/0001-
    /// spida-console-backend-groundwork.md's Phase E addendum for why.
    GridConfig gridT;
    SolverConfig solver;
    /// Model-specific parameters, generic rather than a per-model typed
    /// struct — matches the real wire contract exactly (a single nested
    /// "modelParams" object whose shape depends on "model"): {"mu": ...}
    /// for burgers, {"solitonSpeed": ...} for kdv_rv, {} for ks. See
    /// simulationbuilder.h's per-model construction for what each model
    /// reads out of it, mirroring spida-worker's own
    /// modelParams.value(key, default) pattern exactly.
    // Direct construction, not brace-init: nlohmann::json's initializer_list
    // constructor would treat {nlohmann::json::object()} as a one-element
    // JSON *array* containing an empty object ([{}]), not the empty object
    // itself ({}) — a real bug caught by SimulationRunTest.
    // BURGERS_PILOT_RUNS_TO_COMPLETION exercising the never-overridden
    // default (SimulationRun's cfg.modelParams.value("mu", ...) throws
    // "cannot use value() with array" against the array form).
    nlohmann::json modelParams = nlohmann::json::object();
    ReportingConfig reporting;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE_WITH_DEFAULT(
    SimulationConfig, name, model, grid, gridT, solver, modelParams, reporting)

} // namespace spida::config
