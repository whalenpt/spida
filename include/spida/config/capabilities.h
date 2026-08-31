#pragma once

// Capability discovery — what a given spida binary build actually
// supports (ModelKind x GridKind combinations, each model's modelParams
// schema, and the full SolverKind list), exposed as JSON. Written because
// ADR-0003 already hit real drift once between spida-console's
// hand-maintained TS domain.ts and what spida-worker actually implements
// ("the wire shape had to be realigned... rather than the shape ADR-0001
// had speculatively designed") — this lets a frontend/api-server
// introspect a worker binary directly instead of hand-syncing enums across
// the repo boundary. Wired into spida-worker as `spida-worker --describe`
// (worker/src/main.cpp); see docs/api/openapi.yaml's (not-yet-implemented)
// GET /capabilities for the HTTP surface this is meant to eventually back.

#include <nlohmann/json.hpp>

namespace spida::config {

/// Kept by hand in sync with validate() (validation.h) and SimulationRun's
/// own switch statements (simulationbuilder.h) — there's no single
/// generator for all three yet, so extend this in the same commit that
/// wires a new ModelKind/GridKind.
[[nodiscard]] inline nlohmann::json describeCapabilities()
{
    return {
        {"schemaVersion", 1},
        {"models",
         nlohmann::json::array({
             {{"model", "burgers"},
              {"grids", nlohmann::json::array({"uniform_rvx"})},
              {"modelParams",
               nlohmann::json::array(
                   {{{"name", "mu"}, {"type", "number"}, {"default", 0.0005}}})}},
             {{"model", "kdv_rv"},
              {"grids", nlohmann::json::array({"uniform_rvx"})},
              {"modelParams",
               nlohmann::json::array(
                   {{{"name", "solitonSpeed"}, {"type", "number"}, {"default", 1.0}}})}},
             {{"model", "ks"},
              {"grids", nlohmann::json::array({"uniform_rvx"})},
              {"modelParams", nlohmann::json::array()}},
         })},
        {"solvers", nlohmann::json::array({"etd35", "etd34", "if34", "if45dp"})},
    };
}

} // namespace spida::config
