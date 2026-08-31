#pragma once

// Structured, field-keyed config validation — split out of SimulationRun's
// constructor (simulationbuilder.h) so a caller (spida-worker today; a
// future api-server's POST /configs dry-run eventually) can find out which
// field of a SimulationConfig was wrong without parsing a free-text
// exception message, and without constructing a grid/model/solver at all.
// This directly serves the proposal's own UX principle ("inline field
// errors on the config form; no run is created" — see docs/adr/0001-spida-
// console-backend-groundwork.md's Consequences) which plain exceptions
// can't: one exception carries one free-text message, not a list of
// (field, message) pairs a form could highlight.
//
// Each check below mirrors an exception that already fires somewhere
// deeper in construction (SimulationRun's own throws, UniformGridX's
// minX>=maxX check, Control::setEpsRel, BasePropagator's validatePositive())
// — see each check's comment for exactly where. validate() doesn't replace
// those checks; SimulationRun's constructor calls validate() itself first
// (see simulationbuilder.h) and throws std::invalid_argument if it's
// non-empty, so this stays the single source of truth for what makes a
// config valid rather than a second opinion layered on top.
//
// One real bug this closes: Control::setEpsRel throws a detail::Exception
// (src/utils/except.h), not std::invalid_argument — spida-worker's catch
// chain (worker/src/main.cpp) previously only mapped std::invalid_argument/
// nlohmann::json::exception to failureReason "config_validation", so a
// negative epsRel fell through to the generic "runtime_exception" catch
// instead. validate() rejects it before SimulationRun (and therefore
// Control::setEpsRel) is ever reached.

#include <spida/config/simulationconfig.h>

#include <nlohmann/json.hpp>

#include <string>
#include <vector>

namespace spida::config {

struct ValidationError {
    std::string field;
    std::string message;
};
NLOHMANN_DEFINE_TYPE_NON_INTRUSIVE(ValidationError, field, message)

/// Pure and side-effect-free — never constructs a grid/model/solver, so
/// it's safe to call before deciding whether a run should even be queued.
/// Empty return means cfg is valid (as far as this function checks; see
/// the header comment above for what it deliberately does NOT check, e.g.
/// per-model modelParams ranges).
[[nodiscard]] inline std::vector<ValidationError> validate(const SimulationConfig& cfg)
{
    std::vector<ValidationError> errors;

    // Mirrors SimulationRun's own "not yet wired" throws.
    if (cfg.model != ModelKind::burgers && cfg.model != ModelKind::kdv_rv &&
        cfg.model != ModelKind::ks) {
        errors.push_back(
            {"model", "model kind is not yet implemented (only burgers, kdv_rv, ks are wired)"});
    }
    if (cfg.grid.kind != GridKind::uniform_rvx) {
        errors.push_back(
            {"grid.kind", "grid kind is not yet implemented (only uniform_rvx is wired)"});
    }

    // Mirrors UniformGridX's own throw (src/grid/uniformX.cpp).
    if (cfg.grid.n == 0) {
        errors.push_back({"grid.n", "grid point count must be greater than zero"});
    }
    if (!(cfg.grid.a < cfg.grid.b)) {
        errors.push_back({"grid.b", "grid.a must be less than grid.b"});
    }

    // Mirrors Control::setEpsRel's own throw (src/rkstiff/solver.cpp).
    if (cfg.solver.epsRel < 0.0) {
        errors.push_back({"solver.epsRel", "epsRel must not be negative"});
    }
    if (!(cfg.solver.hInit > 0.0)) {
        errors.push_back({"solver.hInit", "hInit must be greater than zero"});
    }
    if (!(cfg.solver.t0 < cfg.solver.tf)) {
        errors.push_back({"solver.tf", "tf must be greater than t0"});
    }

    // Mirrors BasePropagator::validatePositive()'s throws
    // (src/propagator/propagator.cpp), reached via SimulationRun's calls to
    // setStepsPerOutput1D()/setMaxReports1D()/setLogFrequency().
    if (cfg.reporting.stepsPerOutput1D == 0) {
        errors.push_back(
            {"reporting.stepsPerOutput1D", "stepsPerOutput1D must be greater than zero"});
    }
    if (cfg.reporting.maxReports1D == 0) {
        errors.push_back({"reporting.maxReports1D", "maxReports1D must be greater than zero"});
    }
    if (cfg.solver.logProgress && cfg.reporting.logFrequency == 0) {
        errors.push_back({"reporting.logFrequency",
                          "logFrequency must be greater than zero when logProgress is enabled"});
    }

    return errors;
}

} // namespace spida::config
