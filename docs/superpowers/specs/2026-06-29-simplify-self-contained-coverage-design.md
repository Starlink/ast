# Self-Contained `astSimplify` Coverage — Design

Date: 2026-06-29
Branch: `u/timj/simp-coverage`

## Motivation

The `ast_tester/simplify_fixtures/` set already pins individual `astSimplify`
rules with focused `.map`/`.simp` pairs (see the 2026-05-05 branch-coverage
campaign).
Today, however, a meaningful fraction of the simplify pathways is exercised
only as a side effect of running the broader suite — FitsChan, WCS, and Region
tests drive Mapping pipelines through `astSimplify` and happen to cover merge
branches that no dedicated fixture covers.

The goal of this effort is to make the simplify fixture set **self-contained**:
running only the `simplify_*` tests should cover as much of the simplify
call graph as possible, without relying on any other test or on FitsChan-driven
simplification.
Every branch that is currently covered only because some unrelated test happens
to drive it is a fixture we are missing.

This is an **open-ended, resumable** effort.
Each gap closed is independently valuable; work commits incrementally and can
pause at any point with the repository in a coherent state, then resume later.

The effort is regression-only.
No simplification logic in `src/` is modified.

## Goal

1. Import 13 new fixtures contributed externally (see below) into the
   permanent fixture set.
2. Establish a **gap ledger** — a regenerable, checked-in record of which
   simplify-pathway branches are not covered by the `simplify_*` tests alone.
3. Drive the ledger toward zero **differential** gaps (the primary milestone),
   then opportunistically close or annotate the remaining **absolute** gaps.

## Externally contributed fixtures (import first)

Thirteen new `.map` fixtures, produced by the C AST library and byte-compatible
with the existing fixtures, are staged for import:

```
dssmap_inv_winmap_absorb     dssmap_winmap_absorb      dssmap_zoom_no_merge
neg_box_asymmetric_2d        reconstruct_right_assoc   series_parallel_absorb
sla_run_identity             sla_run_partial           spec_cel_invert
spec_run_inverted_mid        spec_run4_identity        time_run_identity
wcsconv_gappt_iwc_residue
```

The 248 fixtures the contributor's directory shares with
`ast_tester/simplify_fixtures/` are byte-identical to ours and are not
re-imported.

Twelve of the new fixtures are positive (`.map` + `.simp`).
`neg_box_asymmetric_2d` is a negative fixture: by definition its simplified
output equals its input, so it has no `.simp` and is wired with the `.map` in
both the input and reference columns, `skip_string_compare=yes` — the existing
negative-fixture convention.

Import keeps the contributor's filenames (the `wcsconv_*`, `spec_*`, `sla_*`,
`dssmap_*` names are already descriptive and ease cross-referencing).
Each imported fixture gains an `# AST_FIXTURE` header, an inventory row, and a
`simplify_pathways.md` row, exactly like any other rule fixture.

The baseline coverage measurement (below) is taken **after** these 13 are
imported.
There is no pre-import report.

## Coverage tooling and the differential mechanism

### Collection: gcovr; diff/ledger: stdlib Python

The work runs on both macOS and Linux.
Hand-parsing raw `.gcov` text is the one genuinely toolchain-specific problem:
macOS (`llvm-cov gcov`) and Linux (`gcov`) emit subtly different output.
`gcovr` normalizes that, so it is used for **collection**, emitting branch-level
coverage as JSON of identical shape on both platforms.

The **diff and ledger generation** is a stdlib-only Python script
(`ast_tester/coverage/simplify_coverage_gaps.py`) that consumes the gcovr JSON
and imports nothing third-party, so it runs anywhere Python does.
If `gcovr` is not importable/on `PATH`, the script fails with a clear message
telling the operator to activate their gcovr environment rather than silently
producing a wrong report.

### One coverage build, two measured runs

A thin shell wrapper (`ast_tester/coverage/run_simplify_coverage.sh`)
orchestrates:

1. Configure a coverage build:
   `cmake -B build-cov -DCMAKE_BUILD_TYPE=Debug -DCMAKE_C_FLAGS="--coverage -O0"`
   with `--coverage` also on the linker flags. Build once.
2. `ctest --test-dir build-cov -R '^simplify_'` — run the simplify fixtures
   only.
3. `gcovr --json simplify.json --delete ...` — capture `cov_simplify` and
   delete the `.gcda` files so counters reset to zero.
4. `ctest --test-dir build-cov` — run the full suite.
5. `gcovr --json full.json --delete ...` — capture `cov_full`.

Both gcovr invocations `--filter` to the simplify-pathway source files and emit
branch detail.

### Gap classification

For each branch in a target function (defined below), the diff script
classifies:

| Class | Condition | Meaning |
| --- | --- | --- |
| (covered) | hit in `cov_simplify` | already self-covered — no action |
| `differential` | hit in `cov_full`, missed in `cov_simplify` | provably reachable; covered only by other tests; **top priority** |
| `absolute-only` | hit in neither run | candidate dead code, or a gap nothing exercises; investigate |

Closing all `differential` gaps is the literal definition of "self-contained".

## Target-function scope

The "simplify pathway" is an explicit **allowlist** of functions held at the
top of `simplify_coverage_gaps.py` (one place, easy to extend), matched against
the function table gcovr already emits in its JSON.
Branch line numbers are mapped into that function table per run, so function
boundaries are derived fresh every time — no hardcoded line ranges that go
stale.

The starting allowlist:

- `mapping.c` — the `astSimplify` orchestration: the `Simplify` method and the
  fixed-point merge loop plus its local helpers.
- The `MapMerge` static function in each in-scope source file:
  `box.c`, `cmpmap.c`, `dssmap.c`, `grismmap.c`, `interval.c`, `intramap.c`,
  `lutmap.c`, `mathmap.c`, `matrixmap.c`, `normmap.c`, `nullregion.c`,
  `pcdmap.c`, `permmap.c`, `pointlist.c`, `polymap.c`, `prism.c`, `ratemap.c`,
  `selectormap.c`, `shiftmap.c`, `slamap.c`, `specmap.c`, `sphmap.c`,
  `splinemap.c`, `switchmap.c`, `timemap.c`, `tranmap.c`, `unitmap.c`,
  `unitnormmap.c`, `wcsmap.c`, `winmap.c`, `xphmap.c`, `zoommap.c`.
- The merge helpers and their siblings: `WinMat` (winmap.c),
  `MatWin`/`MatWin2` (matrixmap.c), and cmpmap.c's merge-driver/decompose
  helpers (e.g. `CombineMaps`).

The ledger reports per function, so a missing helper is visible (it appears
with no gap rows / fully covered from the start) and is added to the allowlist.
The allowlist is seeded from the merge call graph and grown as the absolute
pass surfaces neighboring helpers.
The implementation plan enumerates the concrete helper set by grepping each
in-scope file for its `MapMerge` and the static functions it calls.

## The gap ledger

### File and format

`ast_tester/coverage/simplify_coverage_gaps.md`, regenerated by the script.
One row per uncovered target-function branch:

| Column | Meaning |
| --- | --- |
| `Location` | `file:line` |
| `Branch` | branch number at that line |
| `Function` | enclosing target function |
| `Class` | `differential` \| `absolute-only` |
| `Status` | human-owned: `open` \| `fixture=<name>` \| `unreachable:<reason>` \| `wontfix:<reason>` |

A "Closed" section at the foot records branches that were once gaps and are now
covered, with the fixture that closed them, for historical traceability.

### Annotation-preserving regeneration (resumability)

Resumability depends on regeneration never discarding human notes.
On each run the script:

1. Reads the existing ledger and keys prior `Status` annotations by
   `(file, function, line, branch#)`.
2. Regenerates the gap list from fresh `simplify.json` / `full.json`.
3. Re-attaches any prior annotation whose key is still a gap.
4. Moves branches that are now covered into the "Closed" section.
5. Lists genuinely new gaps as `open`.

The result is idempotent: re-running after any change updates the picture
without losing prior reasoning.
This is what lets the operator close a few gaps, commit, walk away, and resume
later — on either platform — from an accurate ledger.

## Workflow phases

Each phase ends in a committable state; the effort can stop after any commit.

### P0 — Import the 13 fixtures (once)

Copy the 13 `.map` files (and the 12 `.simp` files) into
`ast_tester/simplify_fixtures/`, add `# AST_FIXTURE` headers, wire rows into
`simplify_tests.txt` (positives as `.map`/`.simp`; `neg_box_asymmetric_2d` as a
`.map`-in-both-columns negative), and add inventory + `simplify_pathways.md`
rows.
Verify each imported fixture passes both the string-diff and astequal tests;
any with ulp-level float variation get `skip_string_compare=yes` (the existing
`rigby` precedent).
Commit.

### P1 — Baseline ledger

Run `run_simplify_coverage.sh` to produce the first ledger from the
post-import fixture set.
Commit the ledger.

### P2 — Differential pass (primary milestone)

For each `differential` gap, author the smallest fixture that triggers the
branch — focused, cascade (capped at three composed mappings), or negative as
appropriate.
Regenerate; the gap should move to "Closed".
Batch commits per class.
The milestone is reached when no `differential` gaps remain: the simplify
fixture set then covers everything the broader suite covers in the target
functions, with no reliance on other tests.

### P3 — Absolute pass (open-ended)

For each `absolute-only` gap, investigate reachability by reading the source.
Reachable ⇒ author a fixture.
Genuinely dead/unreachable ⇒ annotate `unreachable:<reason>` (or `wontfix:`
with justification).
Pursued opportunistically; the repository is coherent and committable after
every regeneration.

## Integration with existing artifacts

- Every new and imported fixture follows the established convention:
  `# AST_FIXTURE id=… polarity=…` header, an inventory row in
  `docs/superpowers/specs/2026-05-05-simplify-branch-coverage-inventory.md`,
  and a row in `ast_tester/simplify_pathways.md`.
- New fixtures are registered in `ast_tester/simplify_tests.txt`; CMake reads
  that file and registers the `simplify_<name>` and `simplify_<name>_astequal`
  tests automatically.
- The `simplify.c` harness, `gen_simplify_fixtures.c`, and `normalize_fixture.c`
  helpers are reused as-is for authoring and verification.

## Verification gate (per batch)

- `ctest -R '^simplify_'` passes under both `Debug` and `Debug` + sanitizers.
- Full `ctest` shows no new regressions relative to `master`.
- `git diff --stat master..HEAD -- src/` is empty — this effort modifies no
  library source.

## Definition of done

"Done" is a milestone, not a finish line.

- **Primary milestone:** zero `differential` gaps — the self-contained goal.
- **Open-ended continuation:** each `absolute-only` gap closed or annotated.

The effort is in a coherent, committable state after every ledger
regeneration, so it can be declared complete-for-now at any milestone and
resumed later.

## Risks and mitigations

- **gcov format differs across platforms.** Mitigated by collecting through
  gcovr (normalized JSON) rather than parsing raw `.gcov`.
- **Coverage build masks optimizer-dependent branches.** Build at `-O0` so
  branch structure matches the source; note this in the wrapper.
- **Ledger regeneration loses human notes.** Mitigated by the
  annotation-preserving merge keyed on `(file, function, line, branch#)`.
- **A target helper is missing from the allowlist.** The per-function ledger
  makes an absent function visible; the allowlist is extended as helpers are
  discovered.
- **Imported fixture serialization varies at ulp level.** Detected at import by
  running both the string-diff and astequal checks; falls back to
  `skip_string_compare=yes`.
- **An `absolute-only` gap is genuinely unreachable.** Annotated in the ledger
  with justification rather than left as a silent uncovered branch.
