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

// ---- Axis/evolution semantics (docs/adr/0004-self-describing-result-
// metadata.md) --------------------------------------------------------
// A report series' bare numeric arrays ("x", "yr"/"yi", ...) don't say what
// an axis *means* -- nls_r's "R" series' "x" is a non-uniform radial
// (Bessel) coordinate, byte-for-byte indistinguishable on the wire from a
// Cartesian "x". These types let modelRegistry() say so explicitly, once,
// so capabilities.h and spida-worker's ManifestBuilder can both serialize
// it instead of the frontend guessing from a series' name (see the ADR's
// Context for the real bug that guessing caused).

/// Geometry/domain of one axis.
enum class AxisCoordinate { cartesian, radial, spectral };
NLOHMANN_JSON_SERIALIZE_ENUM(AxisCoordinate,
                             {
                                 {AxisCoordinate::cartesian, "cartesian"},
                                 {AxisCoordinate::radial, "radial"},
                                 {AxisCoordinate::spectral, "spectral"},
                             })

/// Whether an axis's grid points are evenly spaced. Bessel-root/Chebyshev
/// grids are "nonuniform"; every uniform_* GridKind is "uniform".
enum class AxisSpacing { uniform, nonuniform };
NLOHMANN_JSON_SERIALIZE_ENUM(AxisSpacing,
                             {
                                 {AxisSpacing::uniform, "uniform"},
                                 {AxisSpacing::nonuniform, "nonuniform"},
                             })

/// Only meaningful when coordinate == spectral -- which integral transform
/// produced the axis. Orthogonal to AxisOrdering below: a real-optimized
/// half-spectrum FFT is "fourier" AND "ascending" (burgers'/kdv_rv's/ks's
/// "SX"); a full-complex FFT is "fourier" AND "fftNatural" (kdv_cv's "SX");
/// a Hankel transform is "hankel" AND always "ascending" (nls_r/nls_rt's
/// "SR").
enum class AxisTransform { fourier, hankel };
NLOHMANN_JSON_SERIALIZE_ENUM(AxisTransform,
                             {
                                 {AxisTransform::fourier, "fourier"},
                                 {AxisTransform::hankel, "hankel"},
                             })

/// Whether a spectral axis's wire values are already sorted ascending, or
/// arrive in an FFT's *natural* order -- 0, dk, ..., kNyquist, -kNyquist,
/// ..., -dk -- which wraps from the most positive value straight to the
/// most negative partway through the array. Only a full-complex FFT
/// (uniform_cvx/uniform_cvt grids) does this; a real-optimized
/// half-spectrum FFT and a Hankel transform are always non-negative and
/// already ascending. A frontend line renderer that requires ascending x
/// (e.g. uPlot) needs "fftNatural" to know it must reorder the axis (and
/// its per-point data, in lockstep) before plotting -- exactly a
/// numpy.fft.fftshift.
enum class AxisOrdering { ascending, fftNatural };
NLOHMANN_JSON_SERIALIZE_ENUM(AxisOrdering,
                             {
                                 {AxisOrdering::ascending, "ascending"},
                                 {AxisOrdering::fftNatural, "fftNatural"},
                             })

/// The scrubbed coordinate a model marches along. spida propagates some
/// PDE models in time (t) and others in propagation distance (z); the
/// on-wire report meta key is always literally "t" either way (see
/// BasePropagator::report1D/report2D/reportTrack) -- EvolutionMeta::label
/// below is what tells a caller to display it as "z" instead.
enum class EvolutionKind { time, space };
NLOHMANN_JSON_SERIALIZE_ENUM(EvolutionKind,
                             {
                                 {EvolutionKind::time, "time"},
                                 {EvolutionKind::space, "space"},
                             })

/// Describes one axis of one series. Every field is optional/best-effort:
/// an empty string or unset std::optional means "not specified" and is
/// omitted from the serialized JSON (see capabilities.h and
/// worker/src/main.cpp's ManifestBuilder) rather than emitted as null.
struct AxisMeta {
    std::string label;    // display label, e.g. "r", "x", "kx"
    std::string units;    // e.g. "um", "1/um"; empty when dimensionless/unspecified
    std::string quantity; // e.g. "radius", "position", "wavenumber"
    std::optional<AxisCoordinate> coordinate;
    std::optional<AxisSpacing> spacing;
    std::optional<AxisTransform> transform; // spectral only
    std::optional<AxisOrdering> ordering;   // spectral only
};

/// Describes the evolution (marching) coordinate for a whole model -- the
/// same for every series/frame within one run, so this lives on
/// ModelDescriptor rather than per-series.
struct EvolutionMeta {
    std::string label; // "t" or "z"
    std::string units;
    EvolutionKind quantity;
};

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
    /// [x] for kind == "field1d"/"track"; [x, y] for kind == "field2d" --
    /// same order as manifest.json's gridCoords/gridCoordsY. Empty means
    /// "no axis metadata specified" (every field optional; see AxisMeta).
    std::vector<AxisMeta> axes;
    std::string valueLabel; // label for the plotted quantity, e.g. "u", "A"
    std::string valueUnits;
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
    std::optional<EvolutionMeta> evolution;
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
            .defaultSolver =
                {.kind = SolverKind::etd35, .epsRel = 1e-8, .t0 = 0.0, .tf = 1.0, .hInit = 0.5},
            .defaultReporting = {.stepsPerOutput1D = 5, .maxReports1D = 500, .logFrequency = 200},
            .series = {{.name = "X",
                        .kind = "field1d",
                        .valueType = "real",
                        .description = "physical-space field u(x,t)",
                        .axes = {{.label = "x",
                                  .quantity = "position",
                                  .coordinate = AxisCoordinate::cartesian,
                                  .spacing = AxisSpacing::uniform}},
                        .valueLabel = "u"},
                       {.name = "SX",
                        .kind = "field1d",
                        .valueType = "complex",
                        .description = "spectral-space field (real-optimized transform)",
                        .axes = {{.label = "kx",
                                  .quantity = "wavenumber",
                                  .coordinate = AxisCoordinate::spectral,
                                  .transform = AxisTransform::fourier,
                                  .ordering = AxisOrdering::ascending}},
                        .valueLabel = "SX"}},
            .evolution = EvolutionMeta{.label = "t", .quantity = EvolutionKind::time},
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
            .defaultSolver =
                {.kind = SolverKind::etd35, .epsRel = 1e-10, .t0 = 0.0, .tf = 5.0, .hInit = 0.01},
            .defaultReporting = {.stepsPerOutput1D = 20, .maxReports1D = 500, .logFrequency = 500},
            .series = {{.name = "X",
                        .kind = "field1d",
                        .valueType = "real",
                        .description = "physical-space field u(x,t)",
                        .axes = {{.label = "x",
                                  .quantity = "position",
                                  .coordinate = AxisCoordinate::cartesian,
                                  .spacing = AxisSpacing::uniform}},
                        .valueLabel = "u"},
                       {.name = "SX",
                        .kind = "field1d",
                        .valueType = "complex",
                        .description = "spectral-space field",
                        .axes = {{.label = "kx",
                                  .quantity = "wavenumber",
                                  .coordinate = AxisCoordinate::spectral,
                                  .transform = AxisTransform::fourier,
                                  .ordering = AxisOrdering::ascending}},
                        .valueLabel = "SX"}},
            .evolution = EvolutionMeta{.label = "t", .quantity = EvolutionKind::time},
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
            .defaultSolver =
                {.kind = SolverKind::etd35, .epsRel = 1e-8, .t0 = 0.0, .tf = 150.0, .hInit = 0.01},
            .defaultReporting = {.stepsPerOutput1D = 5, .maxReports1D = 2000, .logFrequency = 500},
            .series = {{.name = "X",
                        .kind = "field1d",
                        .valueType = "real",
                        .description = "physical-space field u(x,t)",
                        .axes = {{.label = "x",
                                  .quantity = "position",
                                  .coordinate = AxisCoordinate::cartesian,
                                  .spacing = AxisSpacing::uniform}},
                        .valueLabel = "u"},
                       {.name = "SX",
                        .kind = "field1d",
                        .valueType = "complex",
                        .description = "spectral-space field",
                        .axes = {{.label = "kx",
                                  .quantity = "wavenumber",
                                  .coordinate = AxisCoordinate::spectral,
                                  .transform = AxisTransform::fourier,
                                  .ordering = AxisOrdering::ascending}},
                        .valueLabel = "SX"}},
            .evolution = EvolutionMeta{.label = "t", .quantity = EvolutionKind::time},
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
            .defaultSolver =
                {.kind = SolverKind::etd35, .epsRel = 1e-4, .t0 = 0.0, .tf = 600.0, .hInit = 0.1},
            .defaultReporting = {.stepsPerOutput1D = 16, .maxReports1D = 500, .logFrequency = 16},
            .series = {{.name = "X",
                        .kind = "field1d",
                        .valueType = "complex",
                        .description = "physical-space field (complex-typed; stays real-valued "
                                       "to machine precision -- see KdvCvPropagator)",
                        .axes = {{.label = "x",
                                  .quantity = "position",
                                  .coordinate = AxisCoordinate::cartesian,
                                  .spacing = AxisSpacing::uniform}},
                        .valueLabel = "u"},
                       {.name = "SX",
                        .kind = "field1d",
                        .valueType = "complex",
                        .description = "spectral-space field",
                        // Full-complex FFT (uniform_cvx), unlike burgers'/
                        // kdv_rv's/ks's real-optimized half-spectrum "SX" --
                        // arrives in the FFT's natural order (0, dk, ...,
                        // kNyquist, -kNyquist, ..., -dk), which wraps rather
                        // than ascending. See docs/adr/0004's Context for
                        // the real rendering bug this distinction fixes.
                        .axes = {{.label = "kx",
                                  .quantity = "wavenumber",
                                  .coordinate = AxisCoordinate::spectral,
                                  .transform = AxisTransform::fourier,
                                  .ordering = AxisOrdering::fftNatural}},
                        .valueLabel = "SX"}},
            .evolution = EvolutionMeta{.label = "t", .quantity = EvolutionKind::time},
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
            .defaultSolver =
                {.kind = SolverKind::etd35, .epsRel = 1e-8, .t0 = 0.0, .tf = 0.8, .hInit = 0.01},
            .defaultReporting = {.stepsPerOutput1D = 10, .maxReports1D = 500, .logFrequency = 500},
            .series = {{.name = "R",
                        .kind = "field1d",
                        .valueType = "complex",
                        .description = "physical-space field A(r)",
                        // Bessel-root grid: non-uniform radial coordinate,
                        // indistinguishable from a Cartesian "x" on the wire
                        // without this -- see docs/adr/0004's Context.
                        .axes = {{.label = "r",
                                  .quantity = "radius",
                                  .coordinate = AxisCoordinate::radial,
                                  .spacing = AxisSpacing::nonuniform}},
                        .valueLabel = "A"},
                       {.name = "SR",
                        .kind = "field1d",
                        .valueType = "complex",
                        .description = "spectral-space field",
                        .axes = {{.label = "kr",
                                  .quantity = "wavenumber",
                                  .coordinate = AxisCoordinate::spectral,
                                  .transform = AxisTransform::hankel,
                                  .ordering = AxisOrdering::ascending}},
                        .valueLabel = "SR"}},
            // dz A = ... -- this model marches in propagation distance,
            // not time, even though the on-wire report meta key is still
            // literally "t" (BasePropagator hardcodes that key name for
            // every model). This is what tells a caller to relabel it "z".
            .evolution = EvolutionMeta{.label = "z", .quantity = EvolutionKind::space},
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
            .defaultGridT =
                GridConfig{.kind = GridKind::uniform_cvt, .n = 512, .a = -6.0, .b = 6.0},
            .defaultSolver =
                {.kind = SolverKind::etd35, .epsRel = 1e-8, .t0 = 0.0, .tf = 0.3, .hInit = 0.01},
            .defaultReporting = {.stepsPerOutput2D = 4, .maxReports2D = 100, .logFrequency = 12},
            .series = {{.name = "RT",
                        .kind = "field2d",
                        .valueType = "complex",
                        .description = "physical-space field A(r,t)",
                        // Two independent axes: x = non-uniform radial
                        // (bessel_root_r), y = the LOCAL envelope-time
                        // coordinate (uniform_cvt) -- distinct from
                        // `evolution` below, which is the propagation
                        // distance z this whole run marches along. Do not
                        // conflate the two just because both use "t"-ish
                        // labels; see docs/adr/0004's Decision section.
                        .axes = {{.label = "r",
                                  .quantity = "radius",
                                  .coordinate = AxisCoordinate::radial,
                                  .spacing = AxisSpacing::nonuniform},
                                 {.label = "t",
                                  .quantity = "time",
                                  .coordinate = AxisCoordinate::cartesian,
                                  .spacing = AxisSpacing::uniform}},
                        .valueLabel = "A"},
                       {.name = "SR",
                        .kind = "field2d",
                        .valueType = "complex",
                        .description = "spectral-space field (kr, omega)",
                        .axes = {{.label = "kr",
                                  .quantity = "wavenumber",
                                  .coordinate = AxisCoordinate::spectral,
                                  .transform = AxisTransform::hankel,
                                  .ordering = AxisOrdering::ascending},
                                 // omega rides on gridT's uniform_cvt --
                                 // a full-complex FFT, same fftNatural
                                 // wraparound as kdv_cv's "SX" above.
                                 {.label = "omega",
                                  .quantity = "frequency",
                                  .coordinate = AxisCoordinate::spectral,
                                  .transform = AxisTransform::fourier,
                                  .ordering = AxisOrdering::fftNatural}},
                        .valueLabel = "SR"}},
            .evolution = EvolutionMeta{.label = "z", .quantity = EvolutionKind::space},
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
