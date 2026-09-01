#pragma once

// Single source of truth for the static, per-ModelKind facts that
// validation.h, simulationbuilder.h, and capabilities.h each used to encode
// by hand, independently: which GridKind(s) a model needs, its modelParams
// schema (name/type/default/description), a numerically-sensible default
// SimulationConfig for it, and the report series (name/kind/valueType) it
// registers. All three consumers now read from modelRegistry()/describe()
// instead of maintaining their own copy.
//
// Why this exists: the previous three-hand-copies design had already
// drifted for real once by the time this was written — capabilities.h and
// simulationbuilder.h each independently hardcoded burgers' "mu" default as
// the literal 0.0005, and three different files each hardcoded their own
// "only burgers/kdv_rv/ks/kdv_cv/nls_r/nls_rt are wired" list. describe()
// returning nullptr for an unrecognized ModelKind is now the ONE definition
// of "wired" — nothing else should hardcode that list again.
//
// Default grid/solver/reporting values below are sourced directly from
// worker/README.md's own worked, numerically-verified examples (see that
// file's "Numerical verification" section for what was checked) — not
// invented here.

#include <spida/config/simulationconfig.h>

#include <array>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace spida::config {

struct ParamSpec {
    std::string name;
    std::string type; // "number" for every modelParams entry today
    double defaultValue{};
    std::string description;
};

struct SeriesSpec {
    std::string name;      // report label, e.g. "X", "SX", "RT"
    std::string kind;      // "field1d" | "field2d" -- matches manifest.json/ResultSeriesDescriptor
    std::string valueType; // "real" | "complex"
    std::string description;
};

struct ModelDescriptor {
    ModelKind model;
    std::string description;
    GridKind gridKind;
    std::optional<GridKind> gridTKind; // set only for nls_rt today
    std::vector<ParamSpec> modelParams;
    GridConfig defaultGrid;
    std::optional<GridConfig> defaultGridT; // set only for nls_rt today
    SolverConfig defaultSolver;
    ReportingConfig defaultReporting;
    std::vector<SeriesSpec> series;
};

/// Order matches CAPABILITIES_TEST.DESCRIBES_ALL_SIX_WIRED_MODELS's
/// expected ordering (test/config_tests.cpp) -- kept stable deliberately,
/// not just an accident of declaration order, since describeCapabilities()
/// iterates this array directly.
[[nodiscard]] inline const std::array<ModelDescriptor, 6>& modelRegistry()
{
    static const std::array<ModelDescriptor, 6> registry{{
        {
            .model = ModelKind::burgers,
            .description = "u_t + u u_x = mu u_xx, real-valued, uniform periodic grid.",
            .gridKind = GridKind::uniform_rvx,
            .gridTKind = std::nullopt,
            .modelParams = {{.name = "mu",
                             .type = "number",
                             .defaultValue = 0.0005,
                             .description = "diffusion coefficient"}},
            .defaultGrid = {.kind = GridKind::uniform_rvx,
                            .n = 8192,
                            .a = -3.14159265358979,
                            .b = 3.14159265358979},
            .defaultGridT = std::nullopt,
            .defaultSolver = {.kind = SolverKind::etd35,
                              .epsRel = 1e-8,
                              .t0 = 0.0,
                              .tf = 1.0,
                              .hInit = 0.5},
            .defaultReporting = {.stepsPerOutput1D = 5, .maxReports1D = 500, .logFrequency = 200},
            .series = {{.name = "X",
                       .kind = "field1d",
                       .valueType = "real",
                       .description = "physical-space field u(x,t)"},
                      {.name = "SX",
                       .kind = "field1d",
                       .valueType = "complex",
                       .description = "spectral-space field (real-optimized transform)"}},
        },
        {
            .model = ModelKind::kdv_rv,
            .description = "u_t + 6 u u_x + u_xxx = 0 (standard normalization), real-valued, "
                           "uniform periodic grid; initial condition is a single soliton.",
            .gridKind = GridKind::uniform_rvx,
            .gridTKind = std::nullopt,
            .modelParams = {{.name = "solitonSpeed",
                             .type = "number",
                             .defaultValue = 1.0,
                             .description =
                                 "translation speed used to build the fixed initial condition"}},
            .defaultGrid = {.kind = GridKind::uniform_rvx, .n = 1024, .a = -20.0, .b = 20.0},
            .defaultGridT = std::nullopt,
            .defaultSolver = {.kind = SolverKind::etd35,
                              .epsRel = 1e-10,
                              .t0 = 0.0,
                              .tf = 5.0,
                              .hInit = 0.01},
            .defaultReporting = {.stepsPerOutput1D = 20, .maxReports1D = 500, .logFrequency = 500},
            .series = {{.name = "X",
                       .kind = "field1d",
                       .valueType = "real",
                       .description = "physical-space field u(x,t)"},
                      {.name = "SX",
                       .kind = "field1d",
                       .valueType = "complex",
                       .description = "spectral-space field"}},
        },
        {
            .model = ModelKind::ks,
            .description = "Kuramoto-Sivashinsky: u_t + u u_x + u_xx + u_xxxx = 0, real-valued, "
                           "uniform periodic grid; chaotic attractor on a wide-enough domain.",
            .gridKind = GridKind::uniform_rvx,
            .gridTKind = std::nullopt,
            .modelParams = {},
            .defaultGrid =
                {.kind = GridKind::uniform_rvx, .n = 256, .a = 0.0, .b = 100.530964914873},
            .defaultGridT = std::nullopt,
            .defaultSolver = {.kind = SolverKind::etd35,
                              .epsRel = 1e-8,
                              .t0 = 0.0,
                              .tf = 150.0,
                              .hInit = 0.01},
            .defaultReporting = {.stepsPerOutput1D = 5, .maxReports1D = 2000, .logFrequency = 500},
            .series = {{.name = "X",
                       .kind = "field1d",
                       .valueType = "real",
                       .description = "physical-space field u(x,t)"},
                      {.name = "SX",
                       .kind = "field1d",
                       .valueType = "complex",
                       .description = "spectral-space field"}},
        },
        {
            .model = ModelKind::kdv_cv,
            .description = "Same PDE as kdv_rv, solved on a full-complex-FFT grid (uniform_cvx) "
                           "instead of the real-optimized half-spectrum transform; fixed "
                           "5-soliton initial condition, no free modelParams.",
            .gridKind = GridKind::uniform_cvx,
            .gridTKind = std::nullopt,
            .modelParams = {},
            .defaultGrid = {.kind = GridKind::uniform_cvx, .n = 512, .a = -150.0, .b = 150.0},
            .defaultGridT = std::nullopt,
            .defaultSolver = {.kind = SolverKind::etd35,
                              .epsRel = 1e-4,
                              .t0 = 0.0,
                              .tf = 600.0,
                              .hInit = 0.1},
            .defaultReporting = {.stepsPerOutput1D = 16, .maxReports1D = 500, .logFrequency = 16},
            .series = {{.name = "X",
                       .kind = "field1d",
                       .valueType = "complex",
                       .description = "physical-space field (complex-typed; stays real-valued "
                                      "to machine precision -- see KdvCvPropagator)"},
                      {.name = "SX",
                       .kind = "field1d",
                       .valueType = "complex",
                       .description = "spectral-space field"}},
        },
        {
            .model = ModelKind::nls_r,
            .description = "Radial cubic NLS: dz A = -i kr^2 A + i gamma |A|^2 A, "
                           "complex-valued, on a non-uniform Hankel-transform grid.",
            .gridKind = GridKind::bessel_root_r,
            .gridTKind = std::nullopt,
            .modelParams = {{.name = "gamma",
                             .type = "number",
                             .defaultValue = 2.0,
                             .description = "Kerr nonlinearity coefficient"},
                            {.name = "amplitude",
                             .type = "number",
                             .defaultValue = 2.0,
                             .description = "initial Gaussian pulse's peak amplitude"}},
            .defaultGrid = {.kind = GridKind::bessel_root_r, .n = 100, .rMax = 5.0},
            .defaultGridT = std::nullopt,
            .defaultSolver = {.kind = SolverKind::etd35,
                              .epsRel = 1e-8,
                              .t0 = 0.0,
                              .tf = 0.8,
                              .hInit = 0.01},
            .defaultReporting = {.stepsPerOutput1D = 10, .maxReports1D = 500, .logFrequency = 500},
            .series = {{.name = "R",
                       .kind = "field1d",
                       .valueType = "complex",
                       .description = "physical-space field A(r)"},
                      {.name = "SR",
                       .kind = "field1d",
                       .valueType = "complex",
                       .description = "spectral-space field"}},
        },
        {
            .model = ModelKind::nls_rt,
            .description = "2D radial + time/frequency cubic NLS: dz A = (-i kr^2 + "
                           "i*0.5*omega^2) A + i gamma |A|^2 A, complex-valued, 2D -- needs BOTH "
                           "grid (bessel_root_r) and gridT (uniform_cvt).",
            .gridKind = GridKind::bessel_root_r,
            .gridTKind = GridKind::uniform_cvt,
            .modelParams = {{.name = "gamma",
                             .type = "number",
                             .defaultValue = 2.0,
                             .description = "Kerr nonlinearity coefficient"},
                            {.name = "amplitude",
                             .type = "number",
                             .defaultValue = 4.0,
                             .description = "initial pulse's peak amplitude"}},
            .defaultGrid = {.kind = GridKind::bessel_root_r, .n = 80, .rMax = 4.0},
            .defaultGridT = GridConfig{.kind = GridKind::uniform_cvt, .n = 512, .a = -6.0, .b = 6.0},
            .defaultSolver = {.kind = SolverKind::etd35,
                              .epsRel = 1e-8,
                              .t0 = 0.0,
                              .tf = 0.3,
                              .hInit = 0.01},
            .defaultReporting =
                {.stepsPerOutput2D = 4, .maxReports2D = 100, .logFrequency = 12},
            .series = {{.name = "RT",
                       .kind = "field2d",
                       .valueType = "complex",
                       .description = "physical-space field A(r,t)"},
                      {.name = "SR",
                       .kind = "field2d",
                       .valueType = "complex",
                       .description = "spectral-space field (kr, omega)"}},
        },
    }};
    return registry;
}

/// nullptr means "not wired" -- the one definition of that concept, used by
/// validate() (validation.h) and SimulationRun (simulationbuilder.h) alike.
[[nodiscard]] inline const ModelDescriptor* describe(ModelKind model)
{
    for (auto const& d : modelRegistry()) {
        if (d.model == model)
            return &d;
    }
    return nullptr;
}

/// Looks up a modelParams default by name. Throws std::logic_error (a
/// programmer error, not a user-input error) if `name` isn't in `d`'s
/// modelParams -- every call site passes a name it already knows is there,
/// so this should never actually throw; it exists to fail loudly rather
/// than silently return 0.0 if a future edit to the registry and its
/// caller ever fall out of sync.
[[nodiscard]] inline double paramDefault(const ModelDescriptor& d, std::string_view name)
{
    for (auto const& p : d.modelParams) {
        if (p.name == name)
            return p.defaultValue;
    }
    throw std::logic_error("paramDefault: \"" + std::string(name) +
                           "\" is not in this ModelDescriptor's modelParams");
}

} // namespace spida::config
