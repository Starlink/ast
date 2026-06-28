# Transform Output Oracle — Design

Date: 2026-06-27
Status: Approved design, ready for implementation planning

## Motivation

PR #51 experiments with adding SIMD transform kernels to SphMap, PolyMap, and
MatrixMap.
That work revealed that AST has no recorded source-of-truth for the numerical
output of its transforms.
Without one, we cannot tell whether a new SIMD kernel, a refactor, or a vendored
library change has altered transform results, and if so by how much.

The SIMD kernels are deliberately not bit-for-bit identical to the scalar paths.
`libmvec` is accurate to roughly four ULP, so a SIMD result can differ from the
scalar `libm` result by a few ULP.
Different CPU architectures (ARM vs x86) and different `libm` implementations also
differ at the ULP level for transcendental functions.
A few ULP is far more than adequate for astronomical coordinate transforms.

The oracle must therefore be a tolerance-based regression guard, not a bit-exact
lock.
It records the present output of the transforms over a fixed corpus and flags
when future output drifts beyond a numerical tolerance, while passing ULP-scale
noise from SIMD, FMA contraction, and cross-architecture `libm` differences.

## Goals

- Record the current transform output for a representative corpus of mappings and
  input points.
- Store the reference in a human-readable, git-tracked text form so that any
  change is visible and reviewable in a pull request diff.
- Provide an automated check, wired into `ctest`, that flags output drift beyond
  tolerance.
- Give identical pass/fail verdicts on ARM and x86; no CPU-architecture lock-in.
- Cross-check that the unsimplified (`.map`) and simplified (`.simp`) forms of a
  fixture produce the same transform output, validating that simplification
  preserves the math.

## Non-goals

- Bit-for-bit reproducibility of transform output.
- Performance measurement; that is the job of `perf/bench_simd.c` in PR #51.
- Testing the simplification structure itself; the existing `simplify` test
  already compares simplified dumps against `.simp` fixtures.

## Corpus

Two sources, both already present in the tree, kept in two separate oracle files
because they exercise two distinct loader code paths.

### Native-dump corpus

`ast_tester/simplify_fixtures/*.map` and `ast_tester/simplify_fixtures/*.simp`,
roughly 206 stem pairs.
Each file is an AST object dump loaded via an `astChannel` reading from the file.
Both members of every stem pair are recorded as separate corpus entries.
They are mathematically equivalent but run different transform code: the `.map`
form is typically a `CmpMap` chain, and the `.simp` form is the collapsed result
(for example a single `ZoomMap`).
The class mix is dominated by `CmpMap`, with `SlaMap`, `TimeMap`, `SpecMap`,
`MatrixMap`, `PolyMap`, `ChebyMap`, `WcsMap`, and others, plus a few Region
classes (`Box`, `Interval`, `NullRegion`, `PointList`) which also have a
well-defined transform.

### FITS-header corpus

`ast_tester/*.head`, 132 files.
Each file is a set of FITS header cards loaded via an `astFitsChan` into a
`FrameSet`.
The mapping tested is `astGetMapping( fs, AST__BASE, AST__CURRENT )`.
These exercise WCS pipelines, including the sky projections that drive SphMap.

## Input-point sampling

Sampling lives only in the generator.
The checker never samples; it reads input points from the oracle file.

For a fixture with `Nin` input axes the generator produces tens of input points
(target around 20 to 30, tunable) from:

- A deterministic low-discrepancy sequence over a base range.
  The exact sequence (for example Halton or Sobol) and the point count are
  implementation-tunable parameters with defaults fixed during implementation.
- Explicit edge points: axis extremes, zero, and a few points just outside any
  natural domain.

Range selection:

- FITS-header FrameSets use the pixel grid: fractional positions spread across
  `1.0` to `NAXISn`, plus a few off-image points.
- Native-dump mappings have no declared domain, so they use a fixed base range
  (starting point `[-1000, 1000]`, tunable).

Adaptive broadening is a soft heuristic.
If the recorded outputs for a fixture are degenerately clustered, the generator
may widen the range and resample.
This is bounded to a small fixed number of attempts and is never a hard failure,
because some mappings collapse their output by design.
SwitchMap routes points to different sub-mappings or regions and so has expected
dead zones and collapse; constant and degenerate mappings and Region masking
also collapse legitimately.

Input points are sampled once per stem and shared between that stem's `.map` and
`.simp` entries, so their input columns are identical.
This identical-input pairing is what makes the cross-check (check C below)
possible.

All input coordinates are floating-point values, not integers.

## Oracle file format

Plain text, sectioned, one section per corpus entry.
All numbers, inputs and outputs alike, are written at full `double` precision with
`%.17g` (C11 build, where `AST_DBL_DIG` is 17).
Because the comparison is tolerance-based, full-precision decimal is more than
precise enough; the stored value is never the dominant source of error.

`AST__BAD` is written as the token `BAD`.
The checker treats `BAD` as matching `BAD` and as mismatching any finite value.

Two files, matching the two loader paths:

- `ast_tester/simplify_fixtures.oracle` for the native-dump corpus.
- `ast_tester/headers.oracle` for the FITS-header corpus.

A section header carries the fixture's identity and shape, so the checker knows
which file to load, how to load it (by extension: `.map` and `.simp` via a
`Channel`, `.head` via a `FitsChan`), and the input and output arity.

Example layout:

```
# transform oracle — regenerate with: gen_transform_oracle
# tol: rtol=1e-12 atol=1e-12  equiv_rtol=1e-9 equiv_atol=1e-9

[simplify_fixtures/box_self_simplify.map  nin=2 nout=2 dir=forward]
  1.2340000000000000e+01 -5.6780000000000000e+01   1.234...e+00  5.678...e+00
  ...

[simplify_fixtures/box_self_simplify.simp  nin=2 nout=2 dir=forward]
  1.2340000000000000e+01 -5.6780000000000000e+01   1.234...e+00  5.678...e+00
  ...
```

Each data row is the input coordinates followed by the recorded output
coordinates for the recorded direction.

## Directions and inverse handling

- The forward output is recorded as the golden reference wherever `TranForward`
  is defined.
- Where both directions are defined, the checker additionally asserts the
  round-trip property `inverse(forward(P)) ~= P` within tolerance.
  Round-trip results are not stored; the original inputs are the reference.
- For mappings where only the inverse is defined, the inverse output is recorded
  directly as the golden reference instead.

Storing forward output plus a round-trip assertion avoids redundant data: storing
`inverse(forward(P))` would merely re-store the inputs.
The combination is sound.
A forward change is caught directly by the forward golden check.
An inverse change larger than tolerance is caught by the round-trip assertion.
The only blind spot of round-trip alone, a compensating change to both forward
and inverse, is closed because the forward output is independently pinned by its
golden reference.

Round-trip accuracy must not be assumed uniform.
Some mappings have an inherently inexact inverse.
PolyMap, for example, computes its inverse iteratively, so its round-trip error
is bounded by convergence rather than by ULP-scale noise.
The round-trip tolerance is therefore worked out per fixture during
implementation, not set as a single tight blanket value, and not papered over
with a single very loose one.
A grossly failing round-trip is treated as a signal to investigate, either a real
bug or too few inverse iterations, rather than something to hide behind a wide
tolerance.
Cases that turn out to be legitimately inexact are revisited individually.

## Checks

For each `.map`/`.simp` stem pair, three independent checks run.

| Check | Asserts | Catches |
| --- | --- | --- |
| A | `.map` output matches its stored reference | regression in the unsimplified path |
| B | `.simp` output matches its stored reference | regression in the simplified path |
| C | `.map` live output matches `.simp` live output, row for row | a simplification bug |

Checks A and B are golden regression checks against the recorded file.
Check C is an equivalence invariant computed between the two live transforms.
The checker pairs sections by basename stem and, because their inputs are
identical, compares outputs row for row.

Check C survives regeneration.
If a simplification bug is introduced and the oracle file is then regenerated,
checks A and B would bless the wrong `.simp` numbers, but check C compares the two
live transforms to each other and still flags the disagreement.
Check C is therefore a genuine correctness invariant, not just a snapshot.

FITS-header entries have no `.map`/`.simp` pairing, so only the golden check and
the round-trip assertion apply to them.

## Tolerances

Comparison uses a combined relative and absolute bound:
`|got - ref| <= atol + rtol * |ref|`.

Three knobs:

- Golden checks (A and B) and round-trip: starting defaults `rtol = 1e-12`,
  `atol = 1e-12`.
- Equivalence check (C): looser starting defaults `equiv_rtol = 1e-9`,
  `equiv_atol = 1e-9`.
  Check C compares two genuinely different arithmetic sequences, a compound chain
  versus one collapsed operation, which can legitimately differ by more than a
  few ULP from accumulated rounding, so it needs a looser bound than the
  same-path golden checks.

Defaults are recorded in the oracle file header and are overridable via command
line or environment variable.
Final values are an implementation deliverable, chosen from the actual spread of
outputs observed across the corpus, ideally cross-checked on both ARM and x86 so
that the gate passes with comfortable margin everywhere while still failing on
real algorithmic changes.

## Programs

Two standalone C programs using the public `ast.h` API only.
They share a small helper that loads a Mapping given a fixture path; everything
else is distinct.

### `gen_transform_oracle.c`

Human-driven, built but not registered as a `ctest` test.
It scans the corpus tree, samples input points, transforms them, and writes the
two oracle files.
It is the only component that scans the filesystem.

### `check_transform_oracle.c`

Registered as a `ctest` test.
The oracle files are its sole manifest; it does no filesystem scanning beyond
loading the fixtures the oracle names.
For each section it loads the named fixture, feeds the listed inputs, transforms
them, and runs the applicable checks.
It exits non-zero on any out-of-tolerance result and prints, for each failure,
the fixture, row, axis, expected value, actual value, and relative difference.

Both programs are built in C11 mode, consistent with the project rule that tests
involving native-format serialization precision use C11 so `AST_DBL_DIG` is 17.

## Build and regeneration workflow

- Add `ast_add_test(check_transform_oracle)` to `ast_tester/CMakeLists.txt`.
- Build `gen_transform_oracle` as a normal executable target, not a test.
- Apply the project's whole-archive satellite-library linking pattern, as the
  other tests do.
- Update `PLAN.md` per the project's test-tracking convention.

Intentional-change workflow:

1. Make the transform change.
2. Rerun `gen_transform_oracle` to regenerate the two oracle files.
3. Review the resulting diff in the pull request; the text format makes every
   changed number visible.
4. Commit the regenerated oracle files alongside the code change.

Behaviour on corpus changes:

- A fixture deleted from the tree makes the checker fail to load it, a real
  signal that the oracle is stale.
- A fixture added to the tree is invisible to the checker until someone
  regenerates, making regeneration the deliberate step that blesses new fixtures.

## Open implementation parameters

These are deliberately left for the implementation phase, with the defaults above
as starting points:

- The specific low-discrepancy sequence and the number of points per fixture.
- The base sampling range for native-dump mappings and the clustering criterion
  and attempt cap for adaptive broadening.
- The final tolerance values, tuned from the observed corpus spread.
- Per-fixture round-trip tolerances for mappings with inherently inexact
  inverses, and investigation of any inverse that round-trips grossly out of
  range.

## Implementation amendments (2026-06-28)

Implementation surfaced that round-trip and equivalence checks over arbitrary
sampled inputs fail pervasively for legitimate domain reasons (out-of-range
inputs, coordinate normalization wrapping, region masking).
The design below was refined during execution so that every check is robust
to unknown domains.
The original sections above are retained for context; where they differ, this
section governs.

### Corpus: three files, not two

A third corpus was added: `ast_tester/*.ast` FrameSet dumps, written to
`framesets.oracle`.
A FrameSet is itself a Mapping (base GRID frame to current frame), so it is
loaded with `astChannel`/`astRead` and transformed directly; non-FrameSet
`.ast` dumps (KeyMap, STC) are not Mappings and are skipped by the loader.
The three files are `simplify_fixtures.oracle`, `headers.oracle` and
`framesets.oracle`.

### Checks

The golden forward check is unchanged and is the domain-agnostic core: it
records the output for fixed inputs and detects any change, regardless of
whether the inputs are "valid".

- **Golden inverse** replaces the original round-trip-of-everything idea.
  For every fixture whose inverse is defined, an inverse section is recorded
  whose inputs are the forward outputs (genuine output-space coordinates) and
  whose outputs are `inverse(forward(P))`.
  This pins the inverse kernel domain-agnostically; it does not assert pixel
  recovery.
- **Equivalence (check C)** is gated on stored agreement: the live `.map` and
  `.simp` forward outputs are compared only at points where the recorded
  golden values already agree (within the equivalence tolerance).
  Generation-time agreement marks the shared valid domain; out-of-domain
  divergence is excluded while a future simplification regression at a
  previously-agreeing point is still caught.
- **Round-trip accuracy** is applied only to the GRID-domain corpora
  (`.head`, `.ast`), where the input domain (the pixel grid) is known and the
  inverse is expected to recover the pixel.
  It asserts `inverse(forward(P)) ~= P` using an *absolute pixel* tolerance
  (default `0.05` px): a relative bound is meaningless near pixel 1, where a
  0.01 px error reads as a 1% relative error.
  Points whose forward output is `AST__BAD` (e.g. all-sky projection corners)
  are skipped.
  Bare `.map` fixtures have no known domain, so round-trip accuracy is not
  asserted for them.

### Sampling

- FITS headers sample each input axis over `[1, NAXISj]`, reading `NAXISj`
  from the header.
  The mapping's actual `Nin` is the axis count (not the `NAXIS` keyword), so
  `WCSAXES > NAXIS` and degenerate axes are handled; axes with no matching
  `NAXISj` fall back to a small `[1, 2]` range.
- `.ast` FrameSets carry no `NAXIS`, so their GRID axes use a default positive
  pixel range `[1, 1000]`.
- Native dumps use a fixed symmetric range `[-1000, 1000]`.
- Adaptive broadening was not implemented; ranges are fixed (a deterministic,
  machine-independent choice).
  Off-domain points simply record `AST__BAD`, which golden handles.

### Round-trip overrides

`transform_oracle_rtrip_overrides.txt` lists per-fixture round-trip
exceptions: `<relpath> off` disables the check, `<relpath> <rtol> <atol>`
loosens it.
Currently `tnx-cheb.head` and `tsc.head` are `off`, for distinct reasons:
`tnx-cheb` has no inverse at all (its Chebyshev distortion builds a ChebyMap,
which disables PolyMap's iterative inverse because it cannot represent its own
Jacobian), so `sky->pixel` is `BAD`; `tsc` has an accurate inverse except at
the quad-cube cross "gap" corners of a rectangular grid, where one sampled
pixel round-trips ~982 px away (relative error 0.70) because the reverse
projection does not flag the gap as invalid.
Their forward and inverse outputs remain pinned by the golden checks; only the
round-trip assertion is disabled, pending investigation.

### Skips

Fixtures whose dump contains an `IntraMap` (the unregistered
`simplifyidentity` test extension) are skipped, as are sections whose
transform raises an AST error on the sampled domain (e.g. two SplineMap
fixtures that extrapolate out of range).

### Tolerances (final, this machine)

- Golden forward/inverse: `rtol = atol = 1e-9`.
  Golden is exact on recompute on the same build, but comparing across
  toolchains is not: an oracle generated with one compiler and checked with
  another differed by up to ~2e-12 (relative) for transcendental-heavy WCS
  inverses on a single machine.
  `1e-9` absorbs that and the larger differences expected across ARM/x86
  `libm`, while still catching real algorithmic change (orders of magnitude
  larger, `>= ~1e-6`).
- Equivalence: `rtol = atol = 1e-9`.
- Round-trip accuracy: absolute `0.05` px (`rtol = 1e-6`).

These pass the full corpus on arm64 under both the conda (RelWithDebInfo) and
homebrew (Debug + ASan/UBSan) toolchains; the cross-toolchain agreement is the
proxy that sized the golden tolerance.  Cross-architecture (x86-64)
confirmation remains future work, but the margin is set for it.
