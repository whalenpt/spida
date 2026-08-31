# spida-worker

One executable, links `SPIDA::spida` directly (in-tree — see this
directory's `CMakeLists.txt`), runs exactly one job, exits. See
`src/main.cpp`'s header comment for what it reads/writes, and
`docs/adr/0003-worker-relocation-and-cooperative-cancellation.md` for how
this came to live here instead of as a separate Conan consumer in
spida-console, and exactly what changed along the way.

## Scope

All six `ModelKind` values are implemented, promoted into the library
itself at `include/spida/models/{burgers,kdv,ks,nls}.h`, `ETD35` by
default (any of `spida::config::SolverKind` is honored — see below):

- `"model": "burgers"` — `u_t + u u_x = mu u_xx`, real-valued, uniform
  periodic grid (`"grid": {"kind": "uniform_rvx", ...}`). `modelParams.mu`
  is the diffusion coefficient (default `0.0005`).
- `"model": "kdv_rv"` — `u_t + 6 u u_x + u_xxx = 0`, the standard
  normalization with no free PDE coefficient, same grid kind.
  `modelParams.solitonSpeed` (default `1.0`) sets the initial condition
  instead: a single soliton whose exact solution translates at speed `c`
  with no change of shape — see `spida::models::KdvPropagator`'s header
  comment.
- `"model": "ks"` — Kuramoto-Sivashinsky, `u_t + u u_x + u_xx + u_xxxx = 0`,
  same grid kind. No `modelParams` — this equation's standard normalization
  has no free coefficient.
- `"model": "kdv_cv"` — the *same* PDE as `kdv_rv`, on a
  `"grid": {"kind": "uniform_cvx", "n": ..., "a": ..., "b": ...}` (full
  complex FFT instead of the real-optimized half-spectrum transform — see
  ADR-0001's Phase D addendum for why this exists as a separate ModelKind).
  Fixed 5-soliton initial condition, no `modelParams`; needs a domain at
  least `[-150, 150]` wide for the soliton centers to fit (see
  `spida::models::KdvCvPropagator`'s header comment).
- `"model": "nls_r"` — radial cubic NLS, `dz A = -i kr^2 A + i gamma |A|^2 A`,
  **complex-valued**, on a `"grid": {"kind": "bessel_root_r", "n": ..., "rMax": ...}`
  (Hankel transform, not FFT — the first non-uniform grid wired; see
  ADR-0001's Phase C addendum). `modelParams.gamma` (default `2.0`) is the
  Kerr nonlinearity coefficient; `modelParams.amplitude` (default `2.0`)
  sets the initial Gaussian's peak amplitude.
- `"model": "nls_rt"` — 2D radial + time/frequency cubic NLS,
  `dz A = (-i kr^2 + i*0.5*omega^2) A + i gamma |A|^2 A`, **complex-valued,
  2D**. Needs TWO grids simultaneously: `"grid": {"kind": "bessel_root_r", ...}`
  (radial) AND `"gridT": {"kind": "uniform_cvt", "n": ..., "a": ..., "b": ...}`
  (time/frequency) — see ADR-0001's Phase E addendum for why this needed a
  new `SimulationConfig.gridT` field. `modelParams.gamma` (default `2.0`)
  and `modelParams.amplitude` (default `4.0`, the initial pulse peak). The
  first model reporting 2D data (`ReportComplex2D`, `"RT"`/`"SR"`).

Requesting a wired model with the wrong `grid.kind`/missing `gridT` (e.g.
`"nls_r"` with `"uniform_rvx"`, or `"nls_rt"` without a `gridT`) is
rejected the same way as an entirely unwired model:
`status: "failed", failureReason: "config_validation"`, along with a
structured `"validationErrors"` array (`{field, message}` per problem —
see `include/spida/config/validation.h`) so a caller doesn't have to parse
`"failureDetail"`'s free text to find out which field was wrong. The same
check (`spida::config::validate()`) also rejects bad numeric input —
negative `epsRel`, non-positive `hInit`, `tf <= t0`, zero `grid.n`/`gridT.n`,
`grid.a >= grid.b`/`gridT.a >= gridT.b` (uniform grids), non-positive
`grid.rMax` (bessel_root_r), zero `stepsPerOutput1D`/`maxReports1D`/
`stepsPerOutput2D`/`maxReports2D` — before anything is constructed or
`status: "running"` is even written.

Run `spida-worker --describe` (no config/output-dir needed) to print which
`ModelKind`×`GridKind` combinations and `SolverKind`s this build actually
supports, plus each model's `modelParams` schema, as JSON — see
`include/spida/config/capabilities.h`. Meant to let a frontend/api-server
introspect a worker binary directly instead of hand-syncing enums against
this README.

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

```json
{
  "name": "run-4",
  "model": "kdv_cv",
  "grid": { "kind": "uniform_cvx", "n": 512, "a": -150.0, "b": 150.0 },
  "solver": { "epsRel": 1e-4, "t0": 0.0, "tf": 600.0, "hInit": 0.1 },
  "reporting": { "stepsPerOutput1D": 16, "maxReports1D": 500, "logFrequency": 16 }
}
```

`kdv_cv` needs `"grid": {"kind": "uniform_cvx", ...}` — same PDE as
`kdv_rv`, but solved on a full-complex-FFT grid; requesting it with
`uniform_rvx` fails `config_validation`. The fixed 5-soliton initial
condition needs the domain shown (or wider) for the soliton centers
(spanning `x` in `[-120, 0]`) to fit without immediate periodic wraparound.

```json
{
  "name": "run-5",
  "model": "nls_r",
  "modelParams": { "gamma": 2.0, "amplitude": 2.0 },
  "grid": { "kind": "bessel_root_r", "n": 100, "rMax": 5.0 },
  "solver": { "epsRel": 1e-8, "t0": 0.0, "tf": 0.8, "hInit": 0.01 },
  "reporting": { "stepsPerOutput1D": 10, "maxReports1D": 500, "logFrequency": 500 }
}
```

`nls_r` needs `"grid": {"kind": "bessel_root_r", ...}`, not `uniform_rvx`
— requesting it with the default grid kind fails `config_validation` (see
`ValidationError`'s `"grid.kind"` field). Its `"R"`/`"SR"` reports are
complex (`ReportComplex1D`), unlike every other wired model's real-valued
reports.

```json
{
  "name": "run-6",
  "model": "nls_rt",
  "modelParams": { "gamma": 2.0, "amplitude": 4.0 },
  "grid": { "kind": "bessel_root_r", "n": 80, "rMax": 4.0 },
  "gridT": { "kind": "uniform_cvt", "n": 512, "a": -6.0, "b": 6.0 },
  "solver": { "epsRel": 1e-8, "t0": 0.0, "tf": 0.3, "hInit": 0.01 },
  "reporting": { "stepsPerOutput2D": 4, "maxReports2D": 100, "logFrequency": 12 }
}
```

`nls_rt` needs BOTH `"grid": {"kind": "bessel_root_r", ...}` (radial) AND
`"gridT": {"kind": "uniform_cvt", ...}` (time/frequency) — omitting
`gridT`, or setting its `kind` to anything else, fails `config_validation`
with a `"gridT.kind"` field error. Its `"RT"`/`"SR"` reports are
`ReportComplex2D` — the first (and, today, only) model reporting 2D data;
`reporting.stepsPerOutput2D`/`maxReports2D` (not the `...1D` fields) govern
its cadence.

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
- **`nls_r`**: no closed-form solution for a general Gaussian initial
  condition, so verified against a conservation law instead — cubic NLS
  with a purely dispersive linear operator conserves the spectral-space L2
  norm `sum_k |A_k|^2` exactly for the continuous equation. Measured
  empirically (not just asserted) against the real ETD35 output: ~6e-6
  relative deviation at `tf=0.01`/`epsRel=1e-10`, ~1.4e-6 at `tf=0.005`,
  scaling down roughly quadratically with `tf` — confirmed to be a
  shrinking time-integration effect (as expected for a finite-order
  adaptive scheme), not a fixed transform-normalization bug, by checking
  the error was unchanged between `grid.n=64` and `grid.n=256` (ruling out
  a spatial-resolution cause). See `spida/models/nls.h`'s header comment
  and `test/config_tests.cpp`'s `NLS_R_POWER_IS_CONSERVED`.
- **`kdv_cv`**: same PDE as `kdv_rv`, whose closed-form check already
  covers the physics — this model's own numerical question is whether the
  full-complex-FFT pipeline stays faithful to a real-coefficient PDE with a
  real initial condition, which should keep `u(x,t)` real for all `t`.
  `PropagatorCV::propagator()` exposes the *spectral* array though, so the
  checkable invariant is that array's Hermitian symmetry
  (`usp[N-k] == conj(usp[k])`), not `Im(u)` directly — an earlier version
  of this check compared `Im`/`Re` of the spectral array directly and
  failed loudly (comparable magnitude, not noise), which is how this
  distinction was caught. Measured for real: ~5e-16 relative asymmetry
  right after construction (confirms the initial condition/transform
  introduce none), growing to ~2.6e-6 at `tf=0.01` (a legitimate
  time-stepping accumulation, same category as `nls_r`'s conservation
  drift above — confirmed, not assumed, since the t0 measurement rules out
  an IC/transform-level cause). See `spida/models/kdv.h`'s
  `KdvCvPropagator` header comment and `test/config_tests.cpp`'s
  `KDV_CV_STAYS_REAL_VALUED`.
- **`nls_rt`**: same conservation argument as `nls_r` (L(k) is purely
  dispersive here too), checked the same way — spectral-space L2 norm
  measured at `t0` and `tf`. ~6.35e-4 relative deviation measured at
  `tf=0.005`/`epsRel=1e-10`/`grid.n=32`/`gridT.n=64`; test tolerance set to
  `2e-3` (~3x margin). See `spida/models/nls.h`'s `NlsRt` header comment
  and `test/config_tests.cpp`'s `NLS_RT_POWER_IS_CONSERVED`.
