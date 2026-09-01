#pragma once

// Capability discovery — what a given spida binary build actually
// supports (ModelKind x GridKind combinations, each model's modelParams
// schema, sensible default configs, and the report series it produces),
// exposed as JSON. Serialized directly from modelregistry.h's
// ModelDescriptor table — the single source of truth this file, validate()
// (validation.h), and SimulationRun (simulationbuilder.h) all read, rather
// than each independently hardcoding the same per-model facts (which had
// already drifted for real once — see modelregistry.h's own header
// comment). Wired into spida-worker as `spida-worker --describe`
// (worker/src/main.cpp); see docs/api/openapi.yaml's Capabilities schema
// for the HTTP surface this is meant to eventually back.

#include <spida/config/modelregistry.h>

#include <nlohmann/json.hpp>

#include <utility>

namespace spida::config {

[[nodiscard]] inline nlohmann::json describeCapabilities()
{
    nlohmann::json models = nlohmann::json::array();
    for (auto const& d : modelRegistry()) {
        nlohmann::json params = nlohmann::json::array();
        for (auto const& p : d.modelParams) {
            params.push_back({
                {"name", p.name},
                {"type", p.type},
                {"default", p.defaultValue},
                {"description", p.description},
            });
        }
        nlohmann::json series = nlohmann::json::array();
        for (auto const& s : d.series) {
            series.push_back({
                {"name", s.name},
                {"kind", s.kind},
                {"valueType", s.valueType},
                {"description", s.description},
            });
        }
        nlohmann::json entry = {
            {"model", d.model},
            {"description", d.description},
            {"grids", nlohmann::json::array({d.gridKind})},
            {"modelParams", params},
            {"defaultGrid", d.defaultGrid},
            {"defaultSolver", d.defaultSolver},
            {"defaultReporting", d.defaultReporting},
            {"series", series},
        };
        // gridT/defaultGridT present only for models needing a second grid
        // dimension (nls_rt today) — matches SimulationConfig.gridT's own
        // "absent/default for every other model" contract.
        if (d.gridTKind.has_value()) {
            entry["gridT"] = nlohmann::json::array({*d.gridTKind});
            entry["defaultGridT"] = *d.defaultGridT;
        }
        models.push_back(std::move(entry));
    }
    return {
        {"schemaVersion", 2},
        {"models", models},
        {"solvers", nlohmann::json::array({"etd35", "etd34", "if34", "if45dp"})},
    };
}

} // namespace spida::config
