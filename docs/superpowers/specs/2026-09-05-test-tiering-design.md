# Test Tiering for Distribution — Design

Date: 2026-09-05
Status: Approved design, ready for implementation planning

## Motivation

`ast_tester/fixtures/` is 20 MB, slightly larger than the entire `src/` tree.
All of it currently ships in `make dist`, and the CMake build must work from that
tarball because ctest is the full suite going forward.

The needs of the two audiences for that tarball are distinct.
An end user unpacking a release wants to answer one question: did this build
work on my machine?
A developer preparing a release wants to answer a different one: have I broken
anything that months of fixture work was added to detect?

Serving only the second audience makes every end user pay 20 MB for tests that
tell them nothing they need to know.
The serialization round-trip corpus is the clearest example: it exists to catch
an attribute dropped by a loader, which is a library defect, not a property of
the user's compiler.

## Goals

- A default `make dist` whose fixture tree is a few MB rather than 20 MB.
- Both `make check` and `ctest` run correctly in that tarball, testing what the
  shipped fixtures support and skipping the rest.
- A `make fulldist` carrying everything, for release verification.
- No loss of the existing protection against a manifest typo, or against a
  source file silently falling out of the tarball.
- No curated "representative" list that a human must maintain and that will
  therefore fall behind.

## Non-goals

- Changing which tests exist, or what any test asserts.
- Reducing coverage in a git checkout. A developer working in the repository
  continues to get all 2475 ctest tests.
- Slimming `src/`, the documentation, or anything outside `ast_tester/fixtures/`.

## Fixture inventory

The 20 MB is unevenly distributed, and the blocks map cleanly onto test
families.

| Block | Size | Files | Consumer |
| --- | --- | --- | --- |
| `oracle/simplify_fixtures.oracle` | 3.7M | 1 | `transform_oracle_simplify` (1601 sections) |
| `simplify/` | 4.2M | 828 | 982 `simplify_*`, most of 996 `roundtrip_*` |
| `wcsconv/expected/` | 1.9M | — | 386 `wcsconv_*` |
| `wcsconv/framesets/` | 1.6M | — | `transform_oracle_framesets` (277 of 351 sections) |
| `oracle/framesets.oracle` | 920K | 1 | `transform_oracle_framesets` |
| `wcsconv/inputs/` | 704K | 142 | `transform_oracle_headers`, grid heads, `testfitschan` |
| `oracle/headers.oracle` | 704K | 1 | `transform_oracle_headers` (277 sections) |
| `plot/expected/*.svg` | 2.1M | 20 | 20 `grid_*` |
| `plot/expected/*.ps` | 1.7M | 20 | **nothing** |
| `programs/` | 1.9M | 41 | the C and Fortran test programs |
| `serialisation/` | 160K | 27 | `roundtrip_*`, 54 `framesets.oracle` sections |
| `*/cases.txt` | 84K | 3 | manifests for the data-driven families |

Two findings from the survey shape the design.

The 1.7 MB of `plot/expected/*.ps` is referenced by no `cases.txt` row, no CMake
driver, and no test source.
The PLplot `plotter_*` tests are pure smoke tests that write `.pdf` and compare
against nothing.

They were not dead by accident: they served as a visual reference a human could
open and compare a plot against by eye.
The `.svg` references now do that job, and unlike the `.ps` files they are also
compared automatically by the 20 `grid_*` tests, so the role is filled by files
that a test keeps honest.
The `.ps` files can therefore be deleted from git, not merely excluded from the
tarball, so that they stop costing 1.7 MB in `fulldist` as well.

The `grid_*` tests are not self-contained.
They read their head files from `wcsconv/inputs/` and compare against
`plot/expected/*.svg`, so keeping them is what requires a FITS header subset to
ship.
The two requirements are one requirement, and the 20 grid heads are a subset
that is already justified by a test rather than by judgment.

## The two tiers

### Default tier

Fixture tree of approximately 5.5 MB, carrying:

- `programs/**` — inputs for the C and Fortran test programs.
- `plot/cases.txt` and `plot/expected/*.svg` — the 20 grid tests.
- `wcsconv/inputs/**` — the header oracle corpus, the grid heads, and the SIP
  headers `testfitschan.c` and `testfitschan.f` read directly.
- `oracle/headers.oracle`, `oracle/transform_oracle_overrides.txt`,
  `oracle/README.md`.
- All three `cases.txt` manifests, including those whose fixtures are absent.

Approximately 113 tests run: 51 C programs, 33 Fortran, 20 `grid_*`,
`transform_oracle_headers` covering 277 sections, 5 `roundtrip_*` over the
`.ast` files in `programs/`, and the `transform_oracle_selftest`,
`test_oracle_util` and `compare_dumps_selftest` self-checks.
The 33 Fortran tests assume the user builds the Fortran interface; without it
the figure is 80.
The 20 PLplot `plotter_*` tests also run where PLplot is found, since their head
files come from `wcsconv/inputs/`.

This is a numerical check rather than a smoke test.
The 277-section header oracle verifies that the build turns pixel coordinates
into the right sky coordinates across every projection AST supports, and the
grid tests exercise the whole Plot and `astGrid` path against byte-comparable
output.

Approximately 2362 tests are skipped: 982 `simplify_*`, 991 `roundtrip_*`,
386 `wcsconv_*`, `transform_oracle_simplify`, `transform_oracle_framesets` and
`keymap_mapsz`.

`keymap_mapsz` is deliberately in the skipped set.
It checks that a committed capture still describes this library's KeyMap growth
policy, which is a developer concern and cannot be affected by the user's build.
Omitting `oracle/keymap_mapsz.txt` skips it through the ordinary mechanism, with
no special case.

### Full tier

Everything, as today, via `make fulldist`.
It produces a differently named tarball — `ast-<version>-full.tar.gz` — so that
a full tarball is never mistaken for a release artifact.

## Detection

Tests are gated on whether their fixtures are present, not on a configure
switch.
One mechanism then serves both cases: the cut-down tarball lacks the files, so
those tests skip; a git checkout has everything, so they all run.
There is no flag to set inconsistently, and no way for a tarball to claim a test
it cannot run.

The difficulty is that a missing fixture means two different things.
It may be a typo in a manifest, which must fail loudly — this is exactly what
`ast_require_fixture`'s `FATAL_ERROR` was added for, so that a typo surfaces at
configure time rather than part-way through a ctest run.
Or it may be a fixture deliberately not shipped, which must skip quietly.

Two files distinguish them.

`ast_tester/fixtures/DIST_MANIFEST` is tracked in git and declares the default
tier's keep-set as a list of path patterns.
It is the single source of truth, read by the dist hook and by the fixture
check, rather than a pattern buried in a recipe.

`ast_tester/fixtures/DIST_TIER` is generated into the tarball only, containing
`cutdown` or `full`.
Its absence means a git checkout.

`ast_require_fixture` then resolves a missing fixture as follows.

| `DIST_TIER` | Fixture matches `DIST_MANIFEST` | Verdict |
| --- | --- | --- |
| absent (git checkout) | — | `FATAL_ERROR` — a typo, or a deleted fixture |
| `full` | — | `FATAL_ERROR` — a packaging defect |
| `cutdown` | yes | `FATAL_ERROR` — should have shipped; a packaging defect |
| `cutdown` | no | deliberate omission; report the test as skipped |

Every strictness the build has today is preserved in the case where typos are
actually introduced, which is the git checkout.
The relaxation applies only inside a tarball, and only to files that the
manifest agrees should not be there.

`ast_require_fixture` changes from a void guard to one returning a status, so
that callers can register a skipped test instead of aborting.

## CMake behavior

Skipped tests are registered and reported rather than omitted.
Omitting them would make `ctest` print `113 tests` with nothing to indicate that
2362 more exist, leaving "did I run the full suite?" unanswerable from the
output.
This is why the `cases.txt` manifests ship even when their fixtures do not: the
test rows come from the manifest, so without it CMake cannot know there are 982
simplify tests to report.
At 84 KB for all three, that is a negligible price for an honest test count.

A test whose fixture is absent is registered and marked `DISABLED TRUE`, so
`ctest` names it under "the following tests did not run" and accounts for it
rather than omitting it from the run entirely.
`DISABLED` is chosen over `SKIP_RETURN_CODE` because the latter requires the
test to run and return the code, which a test with no fixtures cannot do without
a stub command; if the word "Skipped" in the output matters more than avoiding
that stub, the choice can be revisited during implementation.

A configure-time summary states the situation once:

```
-- AST test fixtures: cut-down distribution (fixtures/DIST_TIER)
--   running 113 tests; 2362 skipped for fixtures not distributed
--   omitted: simplify, wcsconv expected, framesets oracle, keymap capture
--   use `make fulldist` from a git checkout for the full corpus
```

`roundtrip`'s `file(GLOB_RECURSE)` needs no guard, because a glob finds only
what is present and cannot contain a typo.
Generators that read a manifest — `gen_frameset_fixtures` reads
`wcsconv/cases.txt` — are guarded the same way as the tests.

## Autotools behavior

The `TESTS` list stays static and the scripts self-skip, so `Makefile.am` needs
no new conditionals.

Each `.tap` producer emits `1..0 # SKIP fixtures not distributed` and exits 0
when its fixtures are absent, which is valid TAP that `tap-driver.sh` reports as
a skipped file.
Where only some rows are absent, the producer emits a per-row
`ok N - name # SKIP fixture not distributed`.

`test_keymap_mapsz.sh` exits 77 when its fixture is absent, the idiom
`testhuge.sh` already uses.

CMake therefore gates at configure time and autotools at run time.
The asymmetry is deliberate: each is the idiomatic mechanism for its build
system, and both derive from the same block-level fact.

## `make dist` and `make fulldist`

`fulldist` is a thin wrapper, so that there is one dist path rather than two
that can diverge:

```make
fulldist:
	$(MAKE) AST_DIST_FULL=1 distdir=$(PACKAGE)-$(VERSION)-full dist
```

Overriding `distdir` on the command line is the intended way to get the distinct
`ast-<version>-full.tar.gz` name, but it is not verified: automake computes both
`distdir` and the archive name from `$(PACKAGE)-$(VERSION)`, and whether the
override propagates cleanly through the `dist` recipe needs checking against a
real automake run.
If it does not, the fallback is a `fulldist` rule that runs `dist` and renames
the resulting archive.

The `ast_tester` dist hook branches on `AST_DIST_FULL`.
In full mode it reverts to today's `find fixtures -type f`.
In default mode it copies only what `DIST_MANIFEST` matches, then writes
`DIST_TIER`.

## The oracle filter

The oracle files are section-keyed by fixture path:

```
[programs/testcmpmap/splittest1.ast  nin=4 nout=4 dir=forward]
```

Filtering is therefore a short `awk` pass that copies a section if and only if
that fixture landed in `$(distdir)`.
The rule is derived from what shipped rather than curated separately, so it
cannot drift.

Applied to the default tier this keeps all 277 sections of `headers.oracle`;
`framesets.oracle` and `simplify_fixtures.oracle` are not shipped at all.

Because the shipped oracle names only shipped fixtures, it is internally
consistent and `check_transform_oracle` needs no change.

## The `dist-check-tracked` guard

The guard currently requires every git-tracked file under `cmake/`,
`ast_tester/` and `src/` to be present in the tarball.
That is what caught `testimmutable.c` and `cmake/run_simplify_noop_test.cmake`
before release.
A cut-down dist deliberately omits tracked fixtures, so the guard as written
would fail, and weakening it wholesale would discard exactly the protection it
was built for.

It is therefore split:

- Non-fixture tracked files keep the existing check unchanged, in both tiers.
  A new `.c`, `.tap` or `.cmake` still cannot silently fall out.
- Fixtures are checked against `DIST_MANIFEST`: every shipped fixture must match
  a manifest entry, and every manifest entry must match at least one shipped
  file, so a stale pattern is caught too.
- Under `AST_DIST_FULL` the original all-tracked-files check applies to
  fixtures as well.

This keeps "a file nobody thought of" caught for sources, and makes the fixture
set a stated decision rather than an oversight — the same principle as the
existing `DIST_CHECK_EXCLUDE` comment.

## Testing

- `make distcheck` on the default tarball: approximately 113 tests pass and
  nothing errors.
- `make fulldist` followed by a full `make check` and `ctest` in the unpacked
  tree: all tests run.
- Unpack the default tarball and configure CMake inside it, then run `ctest`.
  This is the requirement driving the whole design, and the autotools CI job
  already builds CMake from the tarball, so the check belongs there.
- A negative test for the guard: delete a fixture named in a shipped manifest
  from a git checkout and confirm configure still fails with `FATAL_ERROR`.

## Deferred

`make distcheck` cannot be run in the current development environment.
`autoreconf` loops in m4 under the system autoconf 2.69 with m4 1.4.19,
generating multi-gigabyte trace files; it reproduces from a bare `configure.ac`
with no Starlink macro stubs, so it is unrelated to this work.
CI uses autoconf 2.71 or later, so the autotools half of this design will be
verified there rather than locally.
