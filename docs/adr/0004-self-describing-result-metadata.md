# ADR-0004: Self-describing result metadata (axis semantics, labels, units)

## Status

Accepted, implemented on the `spida`/`spida-worker` side (this repo). The
`packages/domain` TS types and `apps/web` consumption (mirror-about-r=0,
log-amplitude default, the `fftshift`-equivalent reorder, deleting the
`radialSeries()` name-guess and `sortAscendingByAxis` runtime sortedness
scan) are `spida-console`'s to build — a separate repo, does not exist in
this workspace — flagged explicitly below wherever this ADR's decisions
imply work there.

## Context

Result frames (`docs/api/binary-frame-spec.md`) carry only `meta:{t}`, a
`type`, and bare numeric arrays (`x`, `yr`, `yi`, …). `manifest.json`'s
`ResultSeriesDescriptor` (`docs/api/openapi.yaml`) adds only
`gridCoords`/`gridCoordsY` — again bare numbers. **Nothing in the contract
says what an axis *means***: `nls_r`'s `"R"` series' `x` is a non-uniform
radial (Bessel) coordinate, byte-for-byte indistinguishable from a
Cartesian `x`. A frontend plot-spec explorer needs this to (a) mirror
radial data about `r=0`, (b) label axes/colorbars, (c) distinguish spatial
from spectral (wavenumber) axes, and (d) label the evolution
coordinate — sometimes time `t` (`burgers`/`kdv_rv`/`ks`/`kdv_cv`) and
sometimes propagation distance `z` (`nls_r`/`nls_rt`), depending on the
model, even though the on-wire report `meta` key is always literally `"t"`
either way (`BasePropagator::report1D`/`report2D`/`reportTrack` hardcode
that key name for every model). Without this, a frontend is left guessing
radial-ness from a series' *name* — a fragile heuristic with no real
signal behind it.

A second gap surfaced independently, from a real bug reproduced against
this repo's own worker: `spida-worker`'s spectral series don't all share
one wire convention. A real-optimized half-spectrum FFT
(`burgers`'/`kdv_rv`'s/`ks`'s `"SX"`) and a Hankel-transform
radial-spectral axis (`nls_r`/`nls_rt`'s `"SR"`) both report non-negative,
already-ascending wavenumbers — but a full-complex FFT (`kdv_cv`'s
`"SX"`, solved on a `uniform_cvx` grid) reports its wavenumber axis in the
FFT's *natural* order (`0, dk, …, kNyquist, -kNyquist, …, -dk`), which
wraps from the most positive value straight to the most negative partway
through the array. Confirmed for real against this repo's own worker (see
"Numerical verification" below): a `kdv_cv` run's `"SX"` `gridCoords` ran
`0, 0.0209, …, 0.670`, then wrapped to `-0.649, …, -0.0209`. A line
renderer requiring ascending x (e.g. uPlot) would silently mis-render this
instead of erroring. A single `coordinate: "spectral"` value can't express
the distinction that causes this: it can't tell a Hankel or half-spectrum
axis (already ascending) from a full-complex-FFT one (isn't, and wraps).

## Decision

Extend the **manifest `ResultSeriesDescriptor`** (static, per-series,
fetched once) with optional, additive, structured metadata describing each
axis, the plotted value, and the evolution coordinate. Per-frame `meta`
stays reserved for dynamic values (`t`, and now `generatedAt` — see below).
`spida-worker`'s `ManifestBuilder` (`worker/src/main.cpp`) derives these
from `include/spida/config/modelregistry.h`'s `ModelDescriptor`/`SeriesSpec`
— the same single source of truth `validate()`, `SimulationRun`, and
`capabilities.h` already read (see that file's own header comment) — and
`spida-worker --describe`'s `Capabilities.models[].series[]` emits the
identical shape. Every field is optional, so old files still parse and
nothing existing needed to change.

Shape added to `modelregistry.h` (mirrors the wire shape documented in
`docs/api/openapi.yaml`'s `ResultSeriesDescriptor`/`Capabilities` schemas):

```cpp
enum class AxisCoordinate { cartesian, radial, spectral };
enum class AxisSpacing    { uniform, nonuniform };

// Only meaningful when coordinate == spectral -- which integral transform
// produced the axis. Orthogonal to `ordering` below: a real-optimized
// half-spectrum FFT is "fourier" AND "ascending" (burgers'/kdv_rv's/ks's
// "SX"); a full-complex FFT is "fourier" AND "fftNatural" (kdv_cv's "SX");
// a Hankel transform is "hankel" AND always "ascending" (nls_r/nls_rt's
// "SR") -- "spectral" alone can't distinguish any of these, which is what
// made the full-complex case a real, reproduced bug (see Context).
enum class AxisTransform { fourier, hankel };

// Whether a spectral axis's wire values are already sorted ascending, or
// arrive in an FFT's *natural* order -- wraps partway through the array.
// Only a full-complex FFT (uniform_cvx/uniform_cvt grids) does this.
enum class AxisOrdering { ascending, fftNatural };

struct AxisMeta {
    std::string label, units, quantity;   // e.g. "r"/"um"/"radius"
    std::optional<AxisCoordinate> coordinate;
    std::optional<AxisSpacing> spacing;
    std::optional<AxisTransform> transform; // spectral only
    std::optional<AxisOrdering> ordering;   // spectral only
};

enum class EvolutionKind { time, space };
struct EvolutionMeta { std::string label, units; EvolutionKind quantity; };

struct SeriesSpec {
    // ...existing name/kind/valueType/description...
    std::vector<AxisMeta> axes;   // [x] for field1d, [x, y] for field2d
    std::string valueLabel, valueUnits;
};

struct ModelDescriptor {
    // ...existing...
    std::optional<EvolutionMeta> evolution;
};
```

Interpretation rules a consumer (e.g. `spida-console`'s `apps/web`) is
meant to follow from this metadata:

- **Radial** (`axes[i].coordinate == "radial"`) → mirror that axis about
  the origin, reflecting **even** (radial field data is always even about
  `r=0`, so the reflection never flips sign — no `symmetry` field needed).
- **Spectral** (`coordinate == "spectral"`) → default to a log-amplitude
  view; the axis label reads as a wavenumber.
- **Fourier-natural ordering** (`ordering == "fftNatural"`) → reorder the
  axis and its per-point data ascending (a `numpy.fft.fftshift`) before
  handing either to a line/waterfall renderer or an export; every other
  `ordering` value (or its absence) is a no-op. `transform` itself drives
  no behavior directly — it's a label/description signal (`"hankel"` vs
  `"fourier"` informs `quantity`/`label` defaults like "kr" vs "kx") that
  `ordering` alone doesn't carry.
- **Evolution** → the frame slider labels itself from `evolution`
  (e.g. `z = 0.183`) instead of assuming `"t"`.
- **Labels/units** → axis titles, the value/colorbar label, and export
  headers come from `axes[].label`/`units`, `valueLabel`/`valueUnits`.

### Also decided here: per-frame `generatedAt` timestamp

Separately motivated (not part of the axis-semantics gap above, but
decided in the same change): every report frame's own `meta` object now
also carries `generatedAt`, an ISO-8601 UTC wall-clock timestamp of when
that frame's JSON was built — distinct from `meta.t`/`meta.z`
(simulation time/propagation distance). Added once, in
`src/utils/report.hpp`'s `detail::buildMeta()`, which every `buildJson()`
override (`Report1D`/`ReportComplex1D`/`Report2D`/`ReportComplex2D`/
`Track`/`TrackComplex`) already funnels through — a single-point change
covering every report kind uniformly. The formatter itself
(`detail::nowIso8601Utc()`, `src/utils/isotime.hpp`) is shared with
`worker/src/main.cpp`'s `SimulationEvent`/`status.json` timestamps rather
than re-derived a second time.

Note: `ReportHandler::report1D()`/`report2D()`/`reportTrack()` call
`toJson()` twice per checkpoint (once for the file write, once for the
sink) — so the on-disk file and the sink-observed copy of the same frame
carry microsecond-apart `generatedAt` values. Cosmetic, not a correctness
issue; not worth restructuring `ReportHandler` to compute it once given
the size of the actual skew.

## Alternatives Considered

### Alternative 1: Put everything in per-frame `meta`
- **Pros**: `meta` is already an open map, so it's zero-contract-change on
  the frame side.
- **Cons**: Repeats static data on every one of N frames; forces a frame
  fetch to learn a series' geometry when the manifest is already fetched
  first.
- **Why not**: Axis semantics are static per series — the manifest is
  their natural, already-fetched home. `meta` remains for genuinely
  per-frame values (`t`/`z`, and now `generatedAt`).

### Alternative 2: Separate per-series sidecar metadata file
- **Pros**: Keeps the manifest small; fully decoupled.
- **Cons**: A third contract surface and an extra round-trip per series,
  duplicating what the manifest already does.
- **Why not**: The manifest already carries per-series static data
  (`gridCoords`); extending it is simpler than a parallel channel.

### Alternative 3: Keep a frontend series-name heuristic
- **Pros**: No backend work.
- **Cons**: No real signal; breaks on any series not matching the guessed
  naming; cannot express units/labels/spectral/evolution at all.
- **Why not**: It is the problem this ADR exists to remove.

### Alternative 4: A single `coordinate: "spectral"` value, no `transform`/`ordering`
- **Pros**: Smaller schema.
- **Cons**: Can't distinguish a Hankel-transform radial-spectral axis
  (`nls_r`/`nls_rt`'s `"SR"` — non-negative, always ascending) from a
  full-complex Fourier transform (`kdv_cv`'s `"SX"` — symmetric ±k,
  arrives in FFT-natural order and wraps). A consumer would still need its
  own runtime sortedness scan to catch the wrapping case.
- **Why not**: This is the concrete gap that surfaced the reproduced bug
  motivating this ADR (see Context) — `"spectral"` alone isn't enough
  signal.

## Consequences

### Positive
- Radial mirroring, axis/colorbar labels, spectral-vs-spatial handling,
  and the evolution-coordinate readout become data-driven, not guessed —
  for whoever builds `spida-console`'s frontend against this contract.
- Additive & optional → backward compatible; old manifests/report files
  and any existing consumer keep working unchanged.
- `spida-worker --describe`'s `Capabilities` schema version bumped
  `2 → 3` (additive; a consumer ignoring the new keys sees the same shape
  as before) so a caller that *does* want the new fields can detect their
  presence without probing for individual keys.

### Negative
- Two metadata locations to reason about (manifest = static, frame
  `meta` = dynamic) — mitigated by the clear rule stated above.
- The exact per-model axis labels/quantities (`modelregistry.h`) are a
  domain judgment call, not derivable from code alone — reviewed by
  Patrick Whalen before landing.

### Risks
- Producer/consumer drift if a future model is wired without populating
  its `axes`/`evolution` — mitigated by every field being optional (a
  silent gap degrades to "no metadata", not a wrong answer) and by
  `modelregistry.h` remaining the one place per-model facts live.

## Numerical verification

Reproduced for real against this repo's own `spida-worker` binary (not
just asserted):

- `kdv_cv` (`n=64`, `a=-150, b=150`, `tf=1.0`): `"SX"`'s `gridCoords` ran
  `0, 0.0209, …, 0.670`, then wrapped to `-0.649, …, -0.0209` — the exact
  full-complex-FFT natural-order wraparound this ADR's `AxisOrdering`
  field exists to flag. `manifest.json` correctly tagged it
  `{"transform": "fourier", "ordering": "fftNatural"}`.
- `nls_r` (`n=32`, `rMax=5.0`): `"R"`'s `gridCoords` started
  `0.1169, 0.2683, 0.4205, …` — non-uniformly spaced, as expected for a
  Bessel-root grid — correctly tagged
  `{"coordinate": "radial", "spacing": "nonuniform"}`, and
  `manifest.json`'s `evolution` correctly read `{"label": "z", "quantity": "space"}`.
- Every report frame from both runs carried a `meta.generatedAt` ISO-8601
  UTC string alongside the existing `meta.t`.

## Implementation order

1. `include/spida/config/modelregistry.h` — add the types above; populate
   `axes`/`valueLabel`/`evolution` per model/series (done here).
2. `include/spida/config/capabilities.h` — serialize the new fields
   (done here).
3. `worker/src/main.cpp`'s `ManifestBuilder` — take the run's
   `ModelDescriptor` and merge its static `axes`/`valueLabel`/
   `valueUnits`/`evolution` into each `manifest.json` series entry, matched
   by series name (done here).
4. `src/utils/report.hpp`/`src/utils/isotime.hpp` — the per-frame
   `generatedAt` timestamp (done here).
5. `docs/api/openapi.yaml`, `docs/api/binary-frame-spec.md`,
   `worker/README.md` — contract docs (done here).
6. `packages/domain` (TS types) and `apps/web` (mirror default,
   log-amplitude default, the `fftshift`-equivalent reorder driven by
   `ordering == "fftNatural"`, evolution-labeled frame slider, deleting
   `radialSeries()`/`sortAscendingByAxis`'s runtime detection) —
   **`spida-console`'s to pick up next; out of scope for this repo.**
