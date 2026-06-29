# Simplify coverage gap tracking

`simplify_coverage_gaps.md` records simplify-pathway branches that the
`simplify_*` tests do not cover on their own.

## Regenerate

    source ~/pyenv/bin/activate      # gcovr 8.x
    ast_tester/coverage/run_simplify_coverage.sh

This builds `build-cov` with `--coverage -O0`, runs the simplify tests then
the full suite (capturing branch coverage after each via gcovr), and rewrites
the ledger.
The gcov tool is auto-selected to match the compiler that built `build-cov`
(the `llvm-cov gcov` beside a clang, or `xcrun llvm-cov gcov` for Apple's
toolchain); override with the `GCOV` environment variable if needed.

## Reading the ledger

- **differential** — branch covered by the full suite but not by the simplify
  fixtures alone. Highest priority: author a fixture so the simplify set is
  self-contained.
- **absolute-only** — branch covered by nothing. Investigate: author a fixture
  if reachable, else set Status to `unreachable:<reason>`.

Branches on a **pure `astOK` status guard** (`if ( !astOK ) ...`,
`if ( astOK ) {`, `while ( astOK )`) are filtered out: their uncovered
direction is the error path, which never fires in a passing run, so it is
non-actionable noise. The header reports how many were filtered. Compound
conditions like `while ( astOK && cond )` are kept, since they carry real
logic.

The **Status** column is the only field you edit by hand:
`open` | `fixture=<name>` | `unreachable:<reason>` | `wontfix:<reason>`.
Regeneration preserves your Status notes by `(file, function, line, branch)`.

The **Deferred** section at the foot lists Region geometric `Simplify`
branches (Box/Interval/Prism/etc. self-simplification). These are out of
scope for the merge-engine self-containment goal — they need a different,
geometry-based fixture methodology — and are recorded only as a known
backlog for a possible future effort.

## Capture method (how the `cap_*` fixtures were produced)

The `cap_*` fixtures in `../simplify_fixtures/` (catalogued in
`captured_fixtures.md`) were not hand-authored. They are real top-level
Mappings that the wider AST test suite hands to `astSimplify`, captured and
replayed as self-contained fixtures. This closes differential gaps without
manually reverse-engineering each branch's trigger.

The method, reproducible end to end:

1. **Instrument** `astSimplify_` in `src/mapping.c` to clone each top-level
   input Mapping and, at process exit (via `atexit`), write each clone as a
   native dump. The write must be deferred to exit and use an explicit
   `status` argument — creating a Channel mid-simplify corrupts AST's global
   channel/warning context, and the variadic `astChannel` macro needs the
   status passed explicitly inside the library. The instrumentation is a
   **temporary patch behind `#ifdef AST_SIMPLIFY_CAPTURE`; it is never
   committed** (the campaign is regression-only — no `src` changes ship).
2. Build with `-DAST_SIMPLIFY_CAPTURE` and run the full suite, then dedup the
   dumps by content.
3. `build_capture_matrix.py` measures, per distinct candidate, which target
   branches it covers (resumable; writes `/tmp/matrix.jsonl`).
4. `greedy_select.py` runs set-cover to pick a near-minimal covering subset.
5. The selected `.map` files are installed and their `.simp` references
   regenerated with the clean (un-instrumented) harness.

`build_capture_matrix.py` and `greedy_select.py` are kept here as a record of
the method; rerunning them requires re-applying the temporary instrumentation
described above.
