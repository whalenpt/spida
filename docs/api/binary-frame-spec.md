# Binary frame wire format

Phase 0 contract freeze for `GET /simulations/:id/results/:series/frames/:i`
(proposal §1: *"a small JSON header followed by a raw Float64Array buffer...
never a JSON number array"*). Not implemented — no api-server exists yet.

## Why binary at all

A single field snapshot is up to `N = 8192` complex points; a run keeps up
to `DEFAULT_MAX_REPORTS_1D = 500` of those (`DEFAULT_MAX_REPORTS_2D = 200`)
— both are `BasePropagator`'s own compiled-in defaults
(`include/spida/propagator/propagator.h`), not proposal guesses. JSON
number arrays cost roughly 3–4× the bytes and a slow parse for no benefit
at that size.

## What already exists server-side, and what this format adds

`ReportHandler` (via `ReportData1D`/`ReportComplex1D`/etc. in
`src/utils/report.hpp`) already builds a `nlohmann::json` object per report
— `{"type": "xy", "x": [...], "y": [...]}` for real 1D data,
`{"type": "xy_complex", "x": [...], "yr": [...], "yi": [...]}` for complex
— and, since ADR-0001, can push that same object to a `ReportHandler::Sink`
callback instead of (or in addition to) writing it to disk
(`setReportSink()`/`setWriteReportFiles()` on `BasePropagator`).

This binary format is a **transport re-encoding** of that same data, done
by the (not-yet-built) API server after receiving a report event — it does
not require changing anything in `report.hpp` itself:

- `x` (the grid coordinates) is **not** repeated per frame — it's static
  across all frames of a series and is served once via
  `ResultSeriesDescriptor.gridCoords` (`GET /simulations/:id/results`).
- `y` (real) or `yr`/`yi` (complex) become the Float64Array payload below.

## Wire format

```
byte 0..3     uint32, little-endian   H = header length in bytes
byte 4..4+H-1 UTF-8 JSON              header (see below), exactly H bytes
byte 4+H..end raw IEEE-754 float64,   payload (see below)
              little-endian
```

`Content-Type: application/octet-stream`. No compression, no trailing
padding — the payload ends at the end of the HTTP body.

### Header JSON

```json
{
  "valueType": "real",
  "count": 8192,
  "t": 0.0741721
}
```

| Field | Type | Meaning |
|---|---|---|
| `valueType` | `"real"` \| `"complex"` | Matches `ResultSeriesDescriptor.valueType` for this series. |
| `count` | integer | Number of grid points `N` (matches `ResultSeriesDescriptor.gridCoords.length`). |
| `t` | number | The simulation time this frame was reported at (from the report's `"t"` metadata item — see `ReportHandler::setItem`). |

### Payload

| `valueType` | Payload length | Layout |
|---|---|---|
| `real` | `count` × 8 bytes | `y[0], y[1], ..., y[count-1]` |
| `complex` | `count` × 2 × 8 bytes | interleaved: `re[0], im[0], re[1], im[1], ..., re[count-1], im[count-1]` |

Interleaved re/im (not two separate blocks) per the proposal's own
wording, so a client can read one `Float64Array` and stride-2 index it
directly without a second allocation/copy.

## Client-side decode sketch

```ts
function decodeFrame(buf: ArrayBuffer): ResultFrame {
  const view = new DataView(buf);
  const headerLen = view.getUint32(0, /* littleEndian */ true);
  const headerBytes = new Uint8Array(buf, 4, headerLen);
  const header = JSON.parse(new TextDecoder().decode(headerBytes));
  const payloadOffset = 4 + headerLen;

  if (header.valueType === "real") {
    const y = new Float64Array(buf, payloadOffset, header.count);
    return { valueType: "real", y };
  }
  const interleaved = new Float64Array(buf, payloadOffset, header.count * 2);
  const re = new Float64Array(header.count);
  const im = new Float64Array(header.count);
  for (let i = 0; i < header.count; i++) {
    re[i] = interleaved[2 * i];
    im[i] = interleaved[2 * i + 1];
  }
  return { valueType: "complex", re, im };
}
```

Matches the `ResultFrame` discriminated union from the proposal's §6
`domain.ts` (`RealFrame { valueType: "real"; y: Float64Array }` /
`ComplexFrame { valueType: "complex"; re: Float64Array; im: Float64Array }`).

`payloadOffset` (`4 + headerLen`) is not guaranteed 8-byte aligned, so the
naive `new Float64Array(buf, payloadOffset, ...)` above can throw
(`RangeError: start offset ... is not a multiple of 8`) in a strict
runtime. A real implementation should either pad the header to a multiple
of 8 bytes server-side, or `buf.slice(payloadOffset)` before constructing
the typed array. Left as an open implementation detail here rather than
silently baked into the spec.

## 2D frames

Not yet specified — `Report2D`/`ReportComplex2D` produce a `z` matrix
(`x.length × y.length`), and the proposal's `SeriesKind` includes
`field2d`. Extend this doc (row-major vs column-major, and whether `y`
axis coords also move to `ResultSeriesDescriptor` alongside `x`) when a
2D-reporting model is wired — none is today (only the Burgers 1D pilot).
