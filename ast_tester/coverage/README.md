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

The **Status** column is the only field you edit by hand:
`open` | `fixture=<name>` | `unreachable:<reason>` | `wontfix:<reason>`.
Regeneration preserves your Status notes by `(file, function, line, branch)`.
