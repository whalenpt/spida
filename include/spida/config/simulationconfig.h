#pragma once

// JSON-serializable configuration mirroring the SPIDA Console architecture
// proposal's domain.ts SimulationConfig (discriminated by ModelKind/GridKind/
// SolverKind) — the C++-side source of truth a future spida-worker would parse
// a request into, and demos/usage examples build grids/models/solvers from
// today by hand. Only ModelKind::burgers and GridKind::uniform_rvx are wired
// end to end by simulationbuilder.h so far (see its header comment); the enums
// below are intentionally larger than that to keep the on-disk/wire shape
// stable as more models and grids come online.

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

enum class GridKind { uniform_rvx, uniform_rvt, bessel_root_r, cheb_x };
NLOHMANN_JSON_SERIALIZE_ENUM(GridKind,
                             {
                                 {GridKind::uniform_rvx, "uniform_rvx"},
                                 {GridKind::uniform_rvt, "uniform_rvt"},
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

/// n/a/b are UniformGridRVX-shaped (the only kind simulationbuilder.h wires
/// today); rMax/etc. for the other GridKinds are left for whoever wires them.
struct GridConfig {
    GridKind kind{GridKind::uniform_rvx};
    unsigned n{256};
    double a{0.0};
    double b{1.0};
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE_WITH_DEFAULT(GridConfig, kind, n, a, b)

struct SolverConfig {
    SolverKind kind{SolverKind::etd35};
    double epsRel{1e-8};
    double t0{0.0};
    double tf{1.0};
    double hInit{0.01};
    bool logProgress{false};
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE_WITH_DEFAULT(SolverConfig, kind, epsRel, t0, tf, hInit, logProgress)

/// Mirrors the subset of BasePropagator's setters simulationbuilder.h wires;
/// see include/spida/propagator/propagator.h for the full set (2D/track
/// cadence are omitted here since every model wired today only reports 1D).
/// logFrequency lives here, not on SolverConfig, matching the real worker
/// wire contract (spida-console's apps/web sends it nested under
/// "reporting", not "solver" — see ADR-0003).
struct ReportingConfig {
    unsigned stepsPerOutput1D{5};
    unsigned maxReports1D{500};
    unsigned logFrequency{200};
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE_WITH_DEFAULT(
    ReportingConfig, stepsPerOutput1D, maxReports1D, logFrequency)

struct SimulationConfig {
    std::string name{"run"};
    ModelKind model{ModelKind::burgers};
    GridConfig grid;
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
    SimulationConfig, name, model, grid, solver, modelParams, reporting)

} // namespace spida::config
