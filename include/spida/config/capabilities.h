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
            nlohmann::json sEntry = {
                {"name", s.name},
                {"kind", s.kind},
                {"valueType", s.valueType},
                {"description", s.description},
            };
            // axes/valueLabel/valueUnits are all optional (docs/adr/0004) --
            // omitted entirely rather than emitted empty/null when a series
            // doesn't set them, matching gridT/defaultGridT's own
            // conditional-emit style just below.
            if (!s.axes.empty()) {
                nlohmann::json axes = nlohmann::json::array();
                for (auto const& a : s.axes) {
                    nlohmann::json aEntry = nlohmann::json::object();
                    if (!a.label.empty())
                        aEntry["label"] = a.label;
                    if (!a.units.empty())
                        aEntry["units"] = a.units;
                    if (!a.quantity.empty())
                        aEntry["quantity"] = a.quantity;
                    if (a.coordinate.has_value())
                        aEntry["coordinate"] = *a.coordinate;
                    if (a.spacing.has_value())
                        aEntry["spacing"] = *a.spacing;
                    if (a.transform.has_value())
                        aEntry["transform"] = *a.transform;
                    if (a.ordering.has_value())
                        aEntry["ordering"] = *a.ordering;
                    axes.push_back(std::move(aEntry));
                }
                sEntry["axes"] = std::move(axes);
            }
            if (!s.valueLabel.empty())
                sEntry["valueLabel"] = s.valueLabel;
            if (!s.valueUnits.empty())
                sEntry["valueUnits"] = s.valueUnits;
            series.push_back(std::move(sEntry));
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
        // evolution (docs/adr/0004) -- the marching coordinate's label/
        // units/quantity ("t"/time vs "z"/space), same for every series in
        // this model. Optional; omitted when not set.
        if (d.evolution.has_value()) {
            nlohmann::json evo = {{"quantity", d.evolution->quantity}};
            if (!d.evolution->label.empty())
                evo["label"] = d.evolution->label;
            if (!d.evolution->units.empty())
                evo["units"] = d.evolution->units;
            entry["evolution"] = std::move(evo);
        }
        models.push_back(std::move(entry));
    }
    return {
        // Bumped 2 -> 3: adds optional per-series axes/valueLabel/
        // valueUnits and per-model evolution (docs/adr/0004). Additive and
        // backward-compatible -- a caller ignoring the new keys sees
        // exactly the same shape as schemaVersion 2 -- bumped anyway so a
        // caller that DOES want the new fields can detect their presence
        // without probing for individual keys.
        {"schemaVersion", 3},
        {"models", models},
        {"solvers", nlohmann::json::array({"etd35", "etd34", "if34", "if45dp"})},
    };
}

} // namespace spida::config
