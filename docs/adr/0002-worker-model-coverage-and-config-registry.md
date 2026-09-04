# ADR-0002: Worker model/grid coverage and the config registry

## Status

Accepted — implemented.

## Context

ADR-0001's seams enable a worker; building one for real means covering the six
`ModelKind` values the frontend's own design already names, each pairing with a
specific `GridKind` (or, for `nls_rt`, two — a second, independent grid dimension).
As models were added one at a time, three separate places ended up each hand-encoding
the same per-model facts: `validate()`'s grid/model pairing checks, `SimulationRun`'s
construction switch (including each model's `modelParams` defaults), and
`describeCapabilities()`'s JSON literal. This had already drifted for real by the
time it was noticed — `capabilities.h` and `simulationbuilder.h` each independently
hardcoded burgers' `mu` default as the literal `0.0005`.

## Decision

1. **`spida-worker` lives in-tree** (`worker/`), a plain CMake executable target
   (`add_executable(spida-worker src/main.cpp)`, `SPIDA_WORKER` CMake option) linking
   `SPIDA::spida` directly — no separate Conan consumer, no pinned installed package.
   Runs exactly one job from a JSON config and exits; see `worker/README.md`.

2. **All six `ModelKind` values are wired end to end**, each promoted into the
   library itself (`include/spida/models/{burgers,kdv,ks,nls}.h`):
   `burgers`/`kdv_rv`/`ks` on `GridKind::uniform_rvx`; `kdv_cv` (same PDE as `kdv_rv`,
   full-complex-FFT pipeline) on `GridKind::uniform_cvx`; `nls_r` (first
   complex-valued field, first non-uniform grid) on `GridKind::bessel_root_r`;
   `nls_rt` (2D, first model reporting `field2d` data) needing **both**
   `GridKind::bessel_root_r` (`grid`) and `GridKind::uniform_cvt` (`gridT`, a sibling
   field on `SimulationConfig` added specifically for this — `grid` itself stays a
   single-dimension field for every other model). Full per-model scope, worked
   config examples, and numerical verification methodology (closed-form soliton
   checks, conservation laws, linear stability theory, Hermitian-symmetry checks)
   live in `worker/README.md` — not duplicated here.

3. **`include/spida/config/modelregistry.h`** — a `ModelDescriptor` table that is now
   the single source of truth for the per-model facts `validate()`, `SimulationRun`,
   and `describeCapabilities()` all need: required `grid`/`gridT` kind, `modelParams`
   schema (name/type/default/description), a numerically-sensible default
   grid/gridT/solver/reporting config (sourced from `worker/README.md`'s own verified
   worked examples, not the generic `GridConfig{}` default, which is wrong for most
   models — e.g. a soliton needs room to translate without periodic wraparound), and
   the report series (name/kind/valueType/description) the model registers.
   `describe(ModelKind) -> const ModelDescriptor*` returning `nullptr` is now the one
   definition of "not wired" — the three files above call it instead of each keeping
   their own copy of the wired-model list. Closes the real `mu`-default drift
   mentioned above, not just documents the risk of it.

4. **Structured, field-keyed validation** — `spida::config::validate()`
   (`include/spida/config/validation.h`) returns `std::vector<ValidationError>`
   (`{field, message}`), not a free-text exception, so a config form can highlight
   the specific field that was wrong. `spida-worker` calls it before constructing
   anything or writing `status: "running"`.

5. **Capability introspection** — `spida::config::describeCapabilities()`
   (`spida-worker --describe`) serializes the registry directly: every wired
   `ModelKind`×`GridKind` pairing, each model's `modelParams`/default config/report
   series, and the full `SolverKind` list — as JSON, so a frontend/api-server can
   introspect a worker binary instead of hand-syncing enums against this repo.

## Consequences

- **Positive**: the variant-based construction pattern in `SimulationRun` genuinely
  generalizes to a complex-valued field on a non-uniform grid, and to a second
  independent grid dimension, without structural rework. The registry means a config
  form can pre-fill sensible per-model defaults and discover report series before
  any run exists, and a new model/grid pairing can't drift between validation,
  construction, and capability discovery the way it already did once.
- **Negative / accepted limitation**: `modelregistry.h` is still hand-maintained (no
  generator) — extend it in the same commit that wires a new `ModelKind`/`GridKind`.
- **Notable pitfalls closed along the way** (condensed; full detail in
  `worker/README.md` and this repo's own commit history): `SimulationConfig`'s
  default `modelParams` field was originally brace-initialized, which
  `nlohmann::json`'s initializer-list constructor treats as a one-element JSON
  *array*, not an empty object — crashed any code path that never overrode it; fixed
  with direct construction. `ProgressSnapshot.tf` was silently always `null` because
  nothing called `setFinalTime()`, even though `cfg.solver.tf` is known at parse
  time. `SimulationRun`'s constructor originally hardcoded report output to a path
  relative to CWD, ignoring the caller's real output directory. `kdv_cv`'s "stays
  real-valued" check originally compared `Im`/`Re` of the *spectral* array directly
  and failed loudly — a real signal's DFT is generally complex at nonzero
  frequencies; only the spectral array's Hermitian symmetry reflects physical-space
  realness, not `Im()` itself.
