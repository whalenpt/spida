# spida-worker

One executable, links `SPIDA::spida` directly (in-tree — see this
directory's `CMakeLists.txt`), runs exactly one job, exits. See
`src/main.cpp`'s header comment for what it reads/writes, and
`docs/adr/0003-worker-relocation-and-cooperative-cancellation.md` for how
this came to live here instead of as a separate Conan consumer in
spida-console, and exactly what changed along the way.

## Scope

Three models are implemented, all real-valued on a uniform periodic grid
with `ETD35` by default (any of `spida::config::SolverKind` is honored —
see below), promoted into the library itself at
`include/spida/models/{burgers,kdv,ks}.h`:

- `"model": "burgers"` — `u_t + u u_x = mu u_xx`. `modelParams.mu` is the
  diffusion coefficient (default `0.0005`).
- `"model": "kdv_rv"` — `u_t + 6 u u_x + u_xxx = 0`, the standard
  normalization with no free PDE coefficient. `modelParams.solitonSpeed`
  (default `1.0`) sets the initial condition instead: a single soliton
  whose exact solution translates at speed `c` with no change of shape —
  see `spida::models::KdvPropagator`'s header comment.
- `"model": "ks"` — Kuramoto-Sivashinsky, `u_t + u u_x + u_xx + u_xxxx = 0`.
  No `modelParams` — this equation's standard normalization has no free
  coefficient.

Anything else in `config.json`'s `"model"` field is rejected with
`status: "failed", failureReason: "config_validation"`. NLS/`kdv_cv` come
later — `spida::config::ModelKind` already has room for them so the wire
shape won't need to change again when they're wired.

Construction goes through `spida::config::SimulationRun`
(`include/spida/config/simulationbuilder.h`), which also means any of the
four `SolverKind`s (`etd35`/`etd34`/`if34`/`if45dp`) works via
`config.json`'s `"solver": {"kind": ...}` — the original worker always
hardcoded `ETD35` regardless of what a request asked for.

## Building

In-tree — no separate Conan consumer, no pinned `spida/0.1.1` sitting in a
local Conan cache, no `conan-recipes/spida` GitHub-tarball recipe. Just the
same sequence as every other target in this repo (see the top-level
`CLAUDE.md`):

```bash
conan install . -of build/Release --build=missing -s build_type=Release
cmake --preset conan-release -DSPIDA_TEST=ON -DSPIDA_DEMOS=ON -DSPIDA_WORKER=ON
cmake --build --preset conan-release --parallel
```

Binary lands at `build/Release/worker/spida-worker`.

## Usage

```bash
./build/Release/worker/spida-worker config.json /path/to/output-dir
```

`config.json` shape (`spida::config::SimulationConfig`, see
`include/spida/config/simulationconfig.h` — same shape
`@spida-console/domain`'s `SimulationConfig` sends):

```json
{
  "name": "run-1",
  "model": "burgers",
  "modelParams": { "mu": 0.0005 },
  "grid": { "n": 8192, "a": -3.14159265358979, "b": 3.14159265358979 },
  "solver": { "kind": "etd35", "epsRel": 1e-8, "t0": 0.0, "tf": 1.0, "hInit": 0.5 },
  "reporting": { "stepsPerOutput1D": 5, "maxReports1D": 500, "logFrequency": 200 }
}
```

```json
{
  "name": "run-2",
  "model": "kdv_rv",
  "modelParams": { "solitonSpeed": 1.0 },
  "grid": { "n": 1024, "a": -20.0, "b": 20.0 },
  "solver": { "epsRel": 1e-10, "t0": 0.0, "tf": 5.0, "hInit": 0.01 },
  "reporting": { "stepsPerOutput1D": 20, "maxReports1D": 500, "logFrequency": 500 }
}
```

KdV needs a much wider domain than Burgers' `[-π, π]` — the soliton needs
room to translate without wrapping around and re-interacting with its own
periodic image before `tf`.

```json
{
  "name": "run-3",
  "model": "ks",
  "grid": { "n": 256, "a": 0.0, "b": 100.530964914873 },
  "solver": { "epsRel": 1e-8, "t0": 0.0, "tf": 150.0, "hInit": 0.01 },
  "reporting": { "stepsPerOutput1D": 5, "maxReports1D": 2000, "logFrequency": 500 }
}
```

`ks` needs both a much wider domain (`b = 32π`, the standard chaotic-
attractor reference domain) and a much longer `tf` than either other model
— `maxReports1D` is raised well above the shared `500` default accordingly,
or the run truncates with `stopReason: "max_reports_reached"` before `tf`.

## Cancellation

`SIGTERM` now triggers real cooperative cancellation
(`BasePropagator::requestCancel()`), checked at the run's next report
checkpoint — not an immediate hard kill. The process exits normally (code
`0`) and writes its own `status.json` (`status: "completed"`,
`stopReason: "cancel_requested"`) once it reaches that checkpoint, the same
way a run that finishes any other way does. See
`docs/adr/0003-worker-relocation-and-cooperative-cancellation.md` for the
api-server-side consequence this has (not addressed here — api-server lives
in spida-console).

## Numerical verification

Carried over from the original implementation (see the ADR for the git
history this moved from) rather than re-run here:

- **Burgers**: built and run for real (`n=256`, `tf=0.2`), producing the
  exact files this README's shapes above describe.
- **`kdv_rv`**: checked against the exact closed-form soliton solution.
  With `solitonSpeed: 1.0`, `n=1024`, `a=-20, b=20`, `tf=5.0`: first-report
  peak amplitude `0.5000` (`= c/2`, exact); last-report peak translated by
  exactly `5.0000` (`= c·tf`) with the same amplitude — a soliton's
  defining property. Max absolute error against the shifted closed form:
  `0.0024`, RMS `0.0007` — consistent with `n=1024`'s grid spacing, not
  numerical drift.
- **`ks`**: no closed-form nonlinear solution exists, so verified against
  linear stability theory instead — growth rate at the max-growth mode
  (`k=1/√2`) and decay rate at `k=√2`, both matching the analytic `L(k)` to
  6 decimal places, plus a conservation-law check (`∫u dx` stayed at
  `~1e-13`–`1e-12`, i.e. floating-point noise, over a full `tf=150` chaotic
  run). **A caveat, not a regression**: KS on this domain has a positive
  leading Lyapunov exponent — two runs of the identical config will diverge
  exponentially in their exact pointwise values after enough Lyapunov
  times, even though both are equally valid trajectories on the same
  chaotic attractor. The growth/decay-rate check above, not bitwise
  reproducibility, is what actually verifies this model.
