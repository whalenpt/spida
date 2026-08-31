#include "spida/config/capabilities.h"
#include "spida/config/simulationbuilder.h"
#include "spida/config/simulationconfig.h"
#include "spida/config/validation.h"

#include <algorithm>
#include <filesystem>
#include <vector>

#include <gtest/gtest.h>

namespace fs = std::filesystem;
using namespace spida::config;

// ============================================================
//  SimulationConfig <-> JSON round trip
// ============================================================

TEST(SIMULATION_CONFIG_TEST, DEFAULT_ROUND_TRIPS_THROUGH_JSON)
{
    SimulationConfig cfg;
    nlohmann::json j = cfg;
    SimulationConfig parsed = j.get<SimulationConfig>();

    EXPECT_EQ(parsed.name, cfg.name);
    EXPECT_EQ(parsed.model, cfg.model);
    EXPECT_EQ(parsed.grid.kind, cfg.grid.kind);
    EXPECT_EQ(parsed.grid.n, cfg.grid.n);
    EXPECT_DOUBLE_EQ(parsed.solver.epsRel, cfg.solver.epsRel);
    EXPECT_EQ(parsed.solver.kind, cfg.solver.kind);
    EXPECT_EQ(parsed.modelParams, cfg.modelParams);
    EXPECT_EQ(parsed.reporting.maxReports1D, cfg.reporting.maxReports1D);
    EXPECT_EQ(parsed.reporting.logFrequency, cfg.reporting.logFrequency);
}

TEST(SIMULATION_CONFIG_TEST, PARSES_PARTIAL_JSON_WITH_DEFAULTS)
{
    // Only override a couple of fields — everything else should fall back
    // to SimulationConfig's own defaults (NLOHMANN_DEFINE_TYPE_..._WITH_DEFAULT).
    nlohmann::json j = {{"name", "sweep-01"}, {"modelParams", {{"mu", 0.001}}}};
    SimulationConfig cfg = j.get<SimulationConfig>();

    EXPECT_EQ(cfg.name, "sweep-01");
    EXPECT_DOUBLE_EQ(cfg.modelParams.value("mu", -1.0), 0.001);
    EXPECT_EQ(cfg.model, ModelKind::burgers); // default
    EXPECT_EQ(cfg.grid.n, 256u);              // default
}

TEST(SIMULATION_CONFIG_TEST, MODEL_KIND_SERIALIZES_TO_EXPECTED_STRING)
{
    nlohmann::json j = ModelKind::burgers;
    EXPECT_EQ(j.get<std::string>(), "burgers");
}

// ============================================================
//  validate() (Phase B: structured, field-keyed config validation)
// ============================================================

TEST(VALIDATE_TEST, DEFAULT_CONFIG_IS_VALID)
{
    SimulationConfig cfg;
    EXPECT_TRUE(validate(cfg).empty());
}

TEST(VALIDATE_TEST, REJECTS_UNWIRED_MODEL_KIND_WITH_FIELD)
{
    // All six ModelKind enumerators are wired as of nls_rt (Phase E) --
    // nothing legitimate is left to construct here. Exercises validate()'s
    // defensive default: branch instead, via a value outside the
    // enumerator set (well-defined: ModelKind's underlying type can
    // represent it, only the switch's named cases can't match it).
    SimulationConfig cfg;
    cfg.model = static_cast<ModelKind>(-1);
    auto errors = validate(cfg);
    ASSERT_EQ(errors.size(), 1u);
    EXPECT_EQ(errors[0].field, "model");
}

TEST(VALIDATE_TEST, REJECTS_UNWIRED_GRID_KIND_WITH_FIELD)
{
    SimulationConfig cfg;
    cfg.grid.kind = GridKind::cheb_x;
    auto errors = validate(cfg);
    ASSERT_EQ(errors.size(), 1u);
    EXPECT_EQ(errors[0].field, "grid.kind");
}

TEST(VALIDATE_TEST, REJECTS_ZERO_GRID_N)
{
    SimulationConfig cfg;
    cfg.grid.n = 0;
    auto errors = validate(cfg);
    ASSERT_EQ(errors.size(), 1u);
    EXPECT_EQ(errors[0].field, "grid.n");
}

TEST(VALIDATE_TEST, REJECTS_INVERTED_GRID_BOUNDS)
{
    SimulationConfig cfg;
    cfg.grid.a = 1.0;
    cfg.grid.b = -1.0;
    auto errors = validate(cfg);
    ASSERT_EQ(errors.size(), 1u);
    EXPECT_EQ(errors[0].field, "grid.b");
}

TEST(VALIDATE_TEST, REJECTS_NEGATIVE_EPS_REL)
{
    // Regression check: this is exactly the case that used to fall through
    // to spida-worker's generic "runtime_exception" catch instead of
    // "config_validation" (Control::setEpsRel throws detail::Exception, not
    // std::invalid_argument) — see validation.h's header comment.
    SimulationConfig cfg;
    cfg.solver.epsRel = -1.0;
    auto errors = validate(cfg);
    ASSERT_EQ(errors.size(), 1u);
    EXPECT_EQ(errors[0].field, "solver.epsRel");
}

TEST(VALIDATE_TEST, REJECTS_NON_POSITIVE_H_INIT)
{
    SimulationConfig cfg;
    cfg.solver.hInit = 0.0;
    auto errors = validate(cfg);
    ASSERT_EQ(errors.size(), 1u);
    EXPECT_EQ(errors[0].field, "solver.hInit");
}

TEST(VALIDATE_TEST, REJECTS_TF_NOT_AFTER_T0)
{
    SimulationConfig cfg;
    cfg.solver.t0 = 1.0;
    cfg.solver.tf = 1.0;
    auto errors = validate(cfg);
    ASSERT_EQ(errors.size(), 1u);
    EXPECT_EQ(errors[0].field, "solver.tf");
}

TEST(VALIDATE_TEST, REJECTS_ZERO_STEPS_PER_OUTPUT_1D)
{
    SimulationConfig cfg;
    cfg.reporting.stepsPerOutput1D = 0;
    auto errors = validate(cfg);
    ASSERT_EQ(errors.size(), 1u);
    EXPECT_EQ(errors[0].field, "reporting.stepsPerOutput1D");
}

TEST(VALIDATE_TEST, REJECTS_ZERO_MAX_REPORTS_1D)
{
    SimulationConfig cfg;
    cfg.reporting.maxReports1D = 0;
    auto errors = validate(cfg);
    ASSERT_EQ(errors.size(), 1u);
    EXPECT_EQ(errors[0].field, "reporting.maxReports1D");
}

TEST(VALIDATE_TEST, REJECTS_ZERO_LOG_FREQUENCY_ONLY_WHEN_LOG_PROGRESS_ENABLED)
{
    SimulationConfig cfg;
    cfg.reporting.logFrequency = 0;
    EXPECT_TRUE(validate(cfg).empty()); // logProgress defaults to false

    cfg.solver.logProgress = true;
    auto errors = validate(cfg);
    ASSERT_EQ(errors.size(), 1u);
    EXPECT_EQ(errors[0].field, "reporting.logFrequency");
}

TEST(VALIDATE_TEST, REJECTS_MISMATCHED_MODEL_GRID_PAIRING)
{
    // The proposal's own error-taxonomy example: a model paired with a
    // grid.kind it doesn't use. nls_r requires bessel_root_r, not
    // uniform_rvx (the default) -- distinct from grid.kind being entirely
    // unwired (REJECTS_UNWIRED_GRID_KIND_WITH_FIELD, above).
    SimulationConfig cfg;
    cfg.model = ModelKind::nls_r;
    cfg.grid.kind = GridKind::uniform_rvx;
    auto errors = validate(cfg);
    ASSERT_EQ(errors.size(), 1u);
    EXPECT_EQ(errors[0].field, "grid.kind");
    EXPECT_NE(errors[0].message.find("bessel_root_r"), std::string::npos);
}

TEST(VALIDATE_TEST, REJECTS_ZERO_R_MAX_FOR_BESSEL_ROOT_R)
{
    SimulationConfig cfg;
    cfg.model = ModelKind::nls_r;
    cfg.grid.kind = GridKind::bessel_root_r;
    cfg.grid.rMax = 0.0;
    auto errors = validate(cfg);
    ASSERT_EQ(errors.size(), 1u);
    EXPECT_EQ(errors[0].field, "grid.rMax");
}

TEST(VALIDATE_TEST, NLS_R_WITH_BESSEL_ROOT_R_IS_VALID)
{
    SimulationConfig cfg;
    cfg.model = ModelKind::nls_r;
    cfg.grid.kind = GridKind::bessel_root_r;
    cfg.grid.rMax = 5.0;
    EXPECT_TRUE(validate(cfg).empty());
}

TEST(VALIDATE_TEST, REJECTS_KDV_CV_WITH_UNIFORM_RVX)
{
    // kdv_cv requires uniform_cvx, not the default uniform_rvx -- distinct
    // from grid.kind being entirely unwired.
    SimulationConfig cfg;
    cfg.model = ModelKind::kdv_cv;
    cfg.grid.kind = GridKind::uniform_rvx;
    auto errors = validate(cfg);
    ASSERT_EQ(errors.size(), 1u);
    EXPECT_EQ(errors[0].field, "grid.kind");
    EXPECT_NE(errors[0].message.find("uniform_cvx"), std::string::npos);
}

TEST(VALIDATE_TEST, KDV_CV_WITH_UNIFORM_CVX_IS_VALID)
{
    SimulationConfig cfg;
    cfg.model = ModelKind::kdv_cv;
    cfg.grid.kind = GridKind::uniform_cvx;
    cfg.grid.a = -150.0;
    cfg.grid.b = 150.0;
    EXPECT_TRUE(validate(cfg).empty());
}

TEST(VALIDATE_TEST, REJECTS_NLS_RT_MISSING_GRID_T)
{
    // grid.kind correct (bessel_root_r) but gridT.kind left at its default
    // (uniform_rvx, meaningless for nls_rt) -- the two-grid pairing check.
    SimulationConfig cfg;
    cfg.model = ModelKind::nls_rt;
    cfg.grid.kind = GridKind::bessel_root_r;
    auto errors = validate(cfg);
    ASSERT_EQ(errors.size(), 1u);
    EXPECT_EQ(errors[0].field, "gridT.kind");
}

TEST(VALIDATE_TEST, REJECTS_NLS_RT_WITH_WRONG_GRID_KIND)
{
    SimulationConfig cfg;
    cfg.model = ModelKind::nls_rt;
    cfg.grid.kind = GridKind::uniform_rvx; // wrong -- needs bessel_root_r
    cfg.gridT.kind = GridKind::uniform_cvt;
    cfg.gridT.a = -6.0;
    cfg.gridT.b = 6.0;
    auto errors = validate(cfg);
    ASSERT_EQ(errors.size(), 1u);
    EXPECT_EQ(errors[0].field, "grid.kind");
}

TEST(VALIDATE_TEST, REJECTS_INVERTED_GRID_T_BOUNDS)
{
    SimulationConfig cfg;
    cfg.model = ModelKind::nls_rt;
    cfg.grid.kind = GridKind::bessel_root_r;
    cfg.gridT.kind = GridKind::uniform_cvt;
    cfg.gridT.a = 6.0;
    cfg.gridT.b = -6.0;
    auto errors = validate(cfg);
    ASSERT_EQ(errors.size(), 1u);
    EXPECT_EQ(errors[0].field, "gridT.b");
}

TEST(VALIDATE_TEST, NLS_RT_WITH_BOTH_GRIDS_IS_VALID)
{
    SimulationConfig cfg;
    cfg.model = ModelKind::nls_rt;
    cfg.grid.kind = GridKind::bessel_root_r;
    cfg.grid.rMax = 4.0;
    cfg.gridT.kind = GridKind::uniform_cvt;
    cfg.gridT.a = -6.0;
    cfg.gridT.b = 6.0;
    EXPECT_TRUE(validate(cfg).empty());
}

TEST(VALIDATE_TEST, REJECTS_ZERO_STEPS_PER_OUTPUT_2D)
{
    // Checked unconditionally now that simulationbuilder.h calls
    // setStepsPerOutput2D()/setMaxReports2D() for every model, not just
    // nls_rt -- see validation.h's own comment.
    SimulationConfig cfg;
    cfg.reporting.stepsPerOutput2D = 0;
    auto errors = validate(cfg);
    ASSERT_EQ(errors.size(), 1u);
    EXPECT_EQ(errors[0].field, "reporting.stepsPerOutput2D");
}

TEST(VALIDATE_TEST, REJECTS_ZERO_MAX_REPORTS_2D)
{
    SimulationConfig cfg;
    cfg.reporting.maxReports2D = 0;
    auto errors = validate(cfg);
    ASSERT_EQ(errors.size(), 1u);
    EXPECT_EQ(errors[0].field, "reporting.maxReports2D");
}

TEST(VALIDATE_TEST, REPORTS_MULTIPLE_ERRORS_AT_ONCE)
{
    SimulationConfig cfg;
    cfg.grid.n = 0;
    cfg.solver.hInit = -1.0;
    EXPECT_EQ(validate(cfg).size(), 2u);
}

TEST(VALIDATE_TEST, VALIDATION_ERROR_ROUND_TRIPS_THROUGH_JSON)
{
    ValidationError err{"grid.n", "grid point count must be greater than zero"};
    nlohmann::json j = err;
    EXPECT_EQ(j.at("field"), "grid.n");
    ValidationError parsed = j.get<ValidationError>();
    EXPECT_EQ(parsed.field, err.field);
    EXPECT_EQ(parsed.message, err.message);
}

TEST(SIMULATION_RUN_CONSTRUCTOR_TEST, INVALID_CONFIG_STILL_THROWS_INVALID_ARGUMENT)
{
    // SimulationRun's constructor calls validate() itself (see
    // simulationbuilder.h) -- same std::invalid_argument contract as
    // before the split, now driven by validate() instead of a duplicate
    // set of checks.
    SimulationConfig cfg;
    cfg.solver.epsRel = -1.0;
    EXPECT_THROW(SimulationRun run(cfg), std::invalid_argument);
}

// ============================================================
//  Capability discovery (Phase B)
// ============================================================

TEST(CAPABILITIES_TEST, DESCRIBES_ALL_SIX_WIRED_MODELS)
{
    auto caps = spida::config::describeCapabilities();
    ASSERT_TRUE(caps.contains("models"));
    std::vector<std::string> names;
    for (auto const& m : caps.at("models"))
        names.push_back(m.at("model").get<std::string>());
    EXPECT_EQ(names,
              (std::vector<std::string>{"burgers", "kdv_rv", "ks", "kdv_cv", "nls_r", "nls_rt"}));
}

TEST(CAPABILITIES_TEST, LISTS_ALL_FOUR_SOLVER_KINDS)
{
    auto caps = spida::config::describeCapabilities();
    ASSERT_TRUE(caps.contains("solvers"));
    EXPECT_EQ(caps.at("solvers").size(), 4u);
}

// ============================================================
//  SimulationRun (config -> grid/model/propagator/solver -> evolve)
// ============================================================

class SimulationRunTest : public ::testing::Test {
protected:
    void TearDown() override
    {
        // SimulationRun writes report files under "sim_<name>" in the CWD;
        // clean up so repeated test runs don't accumulate directories.
        std::error_code ec;
        fs::remove_all(fs::path("sim_" + m_cfg.name), ec);
    }

    SimulationConfig m_cfg;
};

TEST_F(SimulationRunTest, REJECTS_UNWIRED_MODEL_KIND)
{
    // All six ModelKind enumerators are wired as of nls_rt (Phase E) --
    // see VALIDATE_TEST.REJECTS_UNWIRED_MODEL_KIND_WITH_FIELD for why this
    // now uses an out-of-range value instead of a real enumerator.
    m_cfg.model = static_cast<ModelKind>(-1);
    EXPECT_THROW(SimulationRun run(m_cfg), std::invalid_argument);
}

TEST_F(SimulationRunTest, KDV_RV_PILOT_RUNS_TO_COMPLETION)
{
    m_cfg.name = "run_kdv_rv";
    m_cfg.model = ModelKind::kdv_rv;
    m_cfg.modelParams = {{"solitonSpeed", 1.0}};
    m_cfg.grid.n = 128;
    m_cfg.grid.a = -20.0;
    m_cfg.grid.b = 20.0;
    m_cfg.solver.epsRel = 1e-10;
    m_cfg.solver.t0 = 0.0;
    m_cfg.solver.tf = 0.05;
    m_cfg.solver.hInit = 0.01;
    m_cfg.reporting.stepsPerOutput1D = 1;

    SimulationRun run(m_cfg);
    EXPECT_TRUE(run.run());
    EXPECT_EQ(run.propagator().stopReason(), spida::StopReason::None);
}

TEST_F(SimulationRunTest, KS_PILOT_RUNS_TO_COMPLETION)
{
    m_cfg.name = "run_ks";
    m_cfg.model = ModelKind::ks;
    m_cfg.grid.n = 64;
    m_cfg.grid.a = 0.0;
    m_cfg.grid.b = 100.530964914873; // 32*pi, the reference domain
    m_cfg.solver.t0 = 0.0;
    m_cfg.solver.tf = 0.5;
    m_cfg.solver.hInit = 0.01;
    m_cfg.reporting.stepsPerOutput1D = 1;

    SimulationRun run(m_cfg);
    EXPECT_TRUE(run.run());
    EXPECT_EQ(run.propagator().stopReason(), spida::StopReason::None);
}

TEST_F(SimulationRunTest, KDV_CV_PILOT_RUNS_TO_COMPLETION)
{
    m_cfg.name = "run_kdv_cv";
    m_cfg.model = ModelKind::kdv_cv;
    m_cfg.grid.kind = GridKind::uniform_cvx;
    m_cfg.grid.n = 256;
    m_cfg.grid.a = -150.0;
    m_cfg.grid.b = 150.0;
    m_cfg.solver.t0 = 0.0;
    m_cfg.solver.tf = 1.0;
    m_cfg.solver.hInit = 0.1;
    m_cfg.reporting.stepsPerOutput1D = 1;

    SimulationRun run(m_cfg);
    EXPECT_TRUE(run.run());
    EXPECT_EQ(run.propagator().stopReason(), spida::StopReason::None);
}

TEST_F(SimulationRunTest, KDV_CV_STAYS_REAL_VALUED)
{
    // Both the PDE (real coefficients) and the 5-soliton initial condition
    // are real-valued in PHYSICAL space, so u(x,t) should stay real for
    // all t under exact evolution -- a standard fact about real-
    // coefficient PDEs. PropagatorCV::propagator() exposes the SPECTRAL
    // array, not the physical field directly (an earlier version of this
    // test checked Im(usp) against Re(usp) directly and failed loudly --
    // comparable magnitude, not noise -- because a real signal's DFT is
    // generally complex at nonzero frequencies; only its Hermitian
    // symmetry, usp[N-k] == conj(usp[k]), reflects physical-space
    // realness). Checked as that spectral-domain invariant instead, at two
    // points, both measured empirically before picking tolerances:
    //  - t0 (right after construction, before any evolution): relative
    //    asymmetry ~5e-16 -- exact machine precision, confirming the
    //    initial condition and transform themselves introduce no
    //    asymmetry at all.
    //  - tf=0.01: relative asymmetry grows to ~2.6e-6 -- a real,
    //    legitimate accumulation from adaptive time-stepping (same
    //    category of effect as NLS_R_POWER_IS_CONSERVED's conservation
    //    drift), not a structural bug -- confirmed by the t0 check above
    //    ruling out an IC/transform-level cause.
    m_cfg.name = "run_kdv_cv_real";
    m_cfg.model = ModelKind::kdv_cv;
    m_cfg.grid.kind = GridKind::uniform_cvx;
    m_cfg.grid.n = 256;
    m_cfg.grid.a = -150.0;
    m_cfg.grid.b = 150.0;
    m_cfg.solver.t0 = 0.0;
    m_cfg.solver.tf = 0.01;
    m_cfg.solver.hInit = 0.005;
    m_cfg.reporting.stepsPerOutput1D = 1;

    auto maxRelativeAsymmetry = [](const std::vector<spida::dcmplx>& usp) {
        const std::size_t N = usp.size();
        double maxMag = 0.0;
        double maxAsymmetry = std::abs(usp[0].imag()); // usp[0] (DC) must itself be real
        for (std::size_t k = 1; k < N; ++k) {
            maxMag = std::max(maxMag, std::abs(usp[k]));
            maxAsymmetry = std::max(maxAsymmetry, std::abs(usp[k] - std::conj(usp[N - k])));
        }
        return maxMag > 0.0 ? maxAsymmetry / maxMag : maxAsymmetry;
    };

    SimulationRun run(m_cfg);
    EXPECT_LT(maxRelativeAsymmetry(run.propagator().propagator()), 1e-10); // t0: near machine precision

    EXPECT_TRUE(run.run());
    EXPECT_LT(maxRelativeAsymmetry(run.propagator().propagator()), 1e-4); // tf: ~38x margin over measured ~2.6e-6
}

TEST_F(SimulationRunTest, NLS_R_PILOT_RUNS_TO_COMPLETION)
{
    m_cfg.name = "run_nls_r";
    m_cfg.model = ModelKind::nls_r;
    m_cfg.grid.kind = GridKind::bessel_root_r;
    m_cfg.grid.n = 64;
    m_cfg.grid.rMax = 5.0;
    m_cfg.modelParams = {{"gamma", 2.0}, {"amplitude", 2.0}};
    m_cfg.solver.epsRel = 1e-8;
    m_cfg.solver.t0 = 0.0;
    m_cfg.solver.tf = 0.05;
    m_cfg.solver.hInit = 0.001;
    m_cfg.reporting.stepsPerOutput1D = 1;

    SimulationRun run(m_cfg);
    EXPECT_TRUE(run.run());
    EXPECT_EQ(run.propagator().stopReason(), spida::StopReason::None);
}

TEST_F(SimulationRunTest, NLS_R_POWER_IS_CONSERVED)
{
    // Cubic NLS with a purely dispersive linear operator conserves the
    // spectral-space L2 norm sum_k |A_k|^2 exactly for the CONTINUOUS
    // equation (see spida/models/nls.h's header comment for the argument)
    // -- the discrete ETD scheme only approximates that, with error
    // shrinking as tf/step size shrink (confirmed empirically before
    // picking this tolerance: ~6e-6 relative at tf=0.01/epsRel=1e-10,
    // ~1.4e-6 at tf=0.005, ~3e-4 at tf=0.05 -- roughly quadratic in tf,
    // consistent with a genuine, shrinking time-integration truncation
    // effect rather than a fixed transform-normalization mismatch, which
    // was checked too: the error is unchanged between grid.n=64 and
    // grid.n=256, ruling out a spatial-resolution cause). 1e-4 relative
    // leaves a ~16x margin over the measured value at these settings.
    m_cfg.name = "run_nls_r_power";
    m_cfg.model = ModelKind::nls_r;
    m_cfg.grid.kind = GridKind::bessel_root_r;
    m_cfg.grid.n = 64;
    m_cfg.grid.rMax = 5.0;
    m_cfg.modelParams = {{"gamma", 2.0}, {"amplitude", 2.0}};
    m_cfg.solver.epsRel = 1e-10;
    m_cfg.solver.t0 = 0.0;
    m_cfg.solver.tf = 0.01;
    m_cfg.solver.hInit = 0.0005;
    m_cfg.reporting.stepsPerOutput1D = 1;

    SimulationRun run(m_cfg);
    auto power = [](const std::vector<spida::dcmplx>& v) {
        double p = 0.0;
        for (auto const& c : v)
            p += std::norm(c);
        return p;
    };
    const double initialPower = power(run.propagator().propagator());
    ASSERT_GT(initialPower, 0.0); // sanity: the initial condition isn't all-zero

    EXPECT_TRUE(run.run());
    const double finalPower = power(run.propagator().propagator());

    EXPECT_NEAR(finalPower, initialPower, 1e-4 * initialPower);
}

TEST_F(SimulationRunTest, NLS_RT_PILOT_RUNS_TO_COMPLETION)
{
    m_cfg.name = "run_nls_rt";
    m_cfg.model = ModelKind::nls_rt;
    m_cfg.grid.kind = GridKind::bessel_root_r;
    m_cfg.grid.n = 32;
    m_cfg.grid.rMax = 4.0;
    m_cfg.gridT.kind = GridKind::uniform_cvt;
    m_cfg.gridT.n = 64;
    m_cfg.gridT.a = -6.0;
    m_cfg.gridT.b = 6.0;
    m_cfg.modelParams = {{"gamma", 2.0}, {"amplitude", 4.0}};
    m_cfg.solver.epsRel = 1e-8;
    m_cfg.solver.t0 = 0.0;
    m_cfg.solver.tf = 0.01;
    m_cfg.solver.hInit = 0.001;
    m_cfg.reporting.stepsPerOutput1D = 1;
    m_cfg.reporting.stepsPerOutput2D = 1;

    SimulationRun run(m_cfg);
    EXPECT_TRUE(run.run());
    EXPECT_EQ(run.propagator().stopReason(), spida::StopReason::None);
}

TEST_F(SimulationRunTest, NLS_RT_POWER_IS_CONSERVED)
{
    // Same conservation argument as NLS_R_POWER_IS_CONSERVED above: L(k) is
    // purely dispersive here too (see spida/models/nls.h's NlsRt comment),
    // so the spectral-space L2 norm should be approximately conserved,
    // with error shrinking as tf/step size shrink. Measured empirically at
    // these settings before picking the tolerance below: ~6.35e-4 relative
    // deviation (tf=0.005, epsRel=1e-10) — 2e-3 leaves ~3x margin.
    m_cfg.name = "run_nls_rt_power";
    m_cfg.model = ModelKind::nls_rt;
    m_cfg.grid.kind = GridKind::bessel_root_r;
    m_cfg.grid.n = 32;
    m_cfg.grid.rMax = 4.0;
    m_cfg.gridT.kind = GridKind::uniform_cvt;
    m_cfg.gridT.n = 64;
    m_cfg.gridT.a = -6.0;
    m_cfg.gridT.b = 6.0;
    m_cfg.modelParams = {{"gamma", 2.0}, {"amplitude", 4.0}};
    m_cfg.solver.epsRel = 1e-10;
    m_cfg.solver.t0 = 0.0;
    m_cfg.solver.tf = 0.005;
    m_cfg.solver.hInit = 0.0005;
    m_cfg.reporting.stepsPerOutput1D = 1;
    m_cfg.reporting.stepsPerOutput2D = 1;

    SimulationRun run(m_cfg);
    auto power = [](const std::vector<spida::dcmplx>& v) {
        double p = 0.0;
        for (auto const& c : v)
            p += std::norm(c);
        return p;
    };
    const double initialPower = power(run.propagator().propagator());
    ASSERT_GT(initialPower, 0.0);

    EXPECT_TRUE(run.run());
    const double finalPower = power(run.propagator().propagator());

    EXPECT_NEAR(finalPower, initialPower, 2e-3 * initialPower);
}

TEST_F(SimulationRunTest, REJECTS_UNWIRED_GRID_KIND)
{
    m_cfg.grid.kind = GridKind::cheb_x;
    EXPECT_THROW(SimulationRun run(m_cfg), std::invalid_argument);
}

TEST_F(SimulationRunTest, BURGERS_PILOT_RUNS_TO_COMPLETION)
{
    m_cfg.name = "run_to_completion";
    m_cfg.grid.n = 64;
    m_cfg.grid.a = -spida::PI;
    m_cfg.grid.b = spida::PI;
    m_cfg.solver.t0 = 0.0;
    m_cfg.solver.tf = 0.05;
    m_cfg.solver.hInit = 0.01;
    m_cfg.reporting.stepsPerOutput1D = 1;
    m_cfg.reporting.maxReports1D = 500;

    SimulationRun run(m_cfg);
    EXPECT_TRUE(run.run());
    EXPECT_EQ(run.propagator().stopReason(), spida::StopReason::None);
}

TEST_F(SimulationRunTest, PROGRESS_SNAPSHOT_TF_IS_POPULATED)
{
    // SimulationRun's constructor should call setFinalTime(cfg.solver.tf)
    // so every ProgressSnapshot a caller observes has tf set, not just t —
    // previously left unset (see simulationbuilder.h's header comment).
    m_cfg.name = "run_progress_tf";
    m_cfg.grid.n = 64;
    m_cfg.grid.a = -spida::PI;
    m_cfg.grid.b = spida::PI;
    m_cfg.solver.t0 = 0.0;
    m_cfg.solver.tf = 0.05;
    m_cfg.solver.hInit = 0.01;
    m_cfg.reporting.stepsPerOutput1D = 1;

    SimulationRun run(m_cfg);
    std::vector<spida::ProgressSnapshot> seen;
    run.propagator().setProgressObserver(
        [&seen](const spida::ProgressSnapshot& s) { seen.push_back(s); });

    EXPECT_TRUE(run.run());
    ASSERT_FALSE(seen.empty());
    ASSERT_TRUE(seen.front().tf.has_value());
    EXPECT_DOUBLE_EQ(*seen.front().tf, m_cfg.solver.tf);
}

TEST_F(SimulationRunTest, CANCEL_REQUESTED_BEFORE_RUN_STOPS_EARLY)
{
    m_cfg.name = "run_cancelled";
    m_cfg.grid.n = 64;
    m_cfg.grid.a = -spida::PI;
    m_cfg.grid.b = spida::PI;
    m_cfg.solver.t0 = 0.0;
    m_cfg.solver.tf = 1.0; // long enough that cancellation, not tf, ends it
    m_cfg.solver.hInit = 0.001;
    m_cfg.reporting.stepsPerOutput1D = 1;
    m_cfg.reporting.maxReports1D = 500;

    SimulationRun run(m_cfg);
    run.propagator().requestCancel();
    EXPECT_TRUE(run.run()); // solver-level evolve() still returns true (no step failed)
    EXPECT_EQ(run.propagator().stopReason(), spida::StopReason::CancelRequested);
}
