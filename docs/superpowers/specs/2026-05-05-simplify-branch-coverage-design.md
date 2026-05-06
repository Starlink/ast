# `astSimplify` Branch Coverage — Design

Date: 2026-05-05
Branch: `u/timj/simplify-testing`

## Motivation

`astSimplify` is the most critical API in the AST library. It collapses
chains of Mappings into simpler equivalents, and is invoked all over the
codebase and by downstream users. Today the simplification pathways are
exercised mostly by integration-style tests (e.g. FITS, WCS, region tests)
and by three scenario fixtures (`brad`, `lsst1`, `rigby`).

The work on this branch already adds ~50 focused `.map`/`.simp` fixtures
that pin individual rules. This design extends that effort into a
systematic, branch-level coverage of every `MapMerge` pathway, plus a
permanent reference document describing the catalogue of simplifications.

## Goal

Establish a permanent test fixture set covering every code path inside
`MapMerge` methods that `astSimplify` can exercise on a Mapping pipeline,
along with a permanent reference document mapping each pathway to the
fixture that pins it.

The work is intended to be regression-only. No simplification logic is
modified.

## Scope

### In scope

1. Every `MapMerge` branch in pure-Mapping classes:
   `cmpmap.c`, `unitmap.c`, `zoommap.c`, `shiftmap.c`, `winmap.c`,
   `matrixmap.c`, `permmap.c`, `lutmap.c`, `polymap.c`, `mathmap.c`,
   `ratemap.c`, `intramap.c`, `splinemap.c`, `normmap.c`, `unitnormmap.c`,
   `selectormap.c`, `switchmap.c`, `tranmap.c`, `timemap.c`, `slamap.c`,
   `specmap.c`, `wcsmap.c`, `sphmap.c`, `pcdmap.c`, `dssmap.c`,
   `grismmap.c`, `xphmap.c`.

2. Every `MapMerge` branch in Region-as-Mapping classes:
   `box.c`, `interval.c`, `nullregion.c`, `pointlist.c`, `prism.c`.

3. Both polarities — positive (simplification fires) and negative (guard
   rejects).

4. The `astSimplify` orchestration code in `mapping.c` (the wrapper that
   calls `MapMerge` repeatedly until fixed point).

### Out of scope (deferred)

- Region `Simplify` overrides (`region.c`, `circle.c`, `ellipse.c`,
  `polygon.c`, `prism.c`, `cmpregion.c`, `stc.c`). These simplify a
  Region's internal geometric representation rather than its embedding
  in a Mapping pipeline. Different methodology (geometric assertions),
  separate effort.

### Non-goals

- No refactoring of `MapMerge` functions.
- No new simplifications.
- No fixes to existing simplifications. Bugs found during this work get
  filed as follow-ups.

## Granularity and fixture composition

### Branch-level granularity

One fixture per distinct branch inside each `MapMerge`. A branch is
distinct if changing the input shape from "matches branch A" to
"matches branch B" produces a different simplified output structure or
refuses simplification for a different reason. Pure conditional logging
or error paths are not branches. A loop iteration that handles each map
in a list uniformly is one branch, not N.

This level of granularity is required because the stated goal is to
detect regressions in any specific pathway. A coarser rule-level
granularity ("ZoomMap merges with neighbouring maps") cannot pin which
specific neighbour-merge rule broke when a regression appears.

### Fixture types

Each fixture is one of:

- **focused** — exercises a single branch in isolation. The default.
- **cascade** — exercises a branch that fires only after another
  simplification reshapes the input. Required for swap-then-merge style
  rules whose effect is only visible across multiple components.
- **scenario** — pre-existing integration regressions (`brad`, `lsst1`,
  `rigby`). Kept as integration safety nets, not branch coverage.

A branch's input requires a multi-component composition to trigger ⇒ it
is a cascade. Otherwise it is focused. Cascade fixtures are capped at
three composed mappings; if a cascade needs more, document why in the
inventory row.

### Polarity

For each branch with a non-trivial guard condition, an additional
**negative** fixture is added that violates the guard. Negative fixtures
catch the highest-risk class of regression: a code change that
accidentally relaxes a guard, silently collapsing pipelines that should
be left alone. The `.simp` reference for a negative fixture is the input
mapping unchanged in structure (with whatever cosmetic
`IsSimp = 1`/normalization a simplify pass adds, but no structural
collapse).

## Inventory artifact

### File location

Working inventory: `docs/superpowers/specs/2026-05-05-simplify-branch-coverage-inventory.md`,
next to this design doc. Survives the campaign as a historical artifact
because it carries the source line ranges that the permanent reference
document deliberately drops.

### Structure

One `##` section per source class, in the order listed under "In scope"
above. Each section opens with a one-paragraph summary of the class's
`MapMerge` strategy, then a single Markdown table of branches.

### Row schema

| Column | Purpose |
| --- | --- |
| `ID` | Stable identifier `<class>-NN`, e.g. `winmap-04`. Referenced from commits, fixture header comments, and the permanent reference doc. |
| `Fixture` | `.map` filename, or `via <id>` for branches reached only through a cascade. |
| `Type` | `focused` \| `cascade` \| `scenario`. |
| `Polarity` | `positive` \| `negative`. |
| `Lines` | Source line range in the `MapMerge` function (e.g. `winmap.c:2104-2147`). Used by gcov verification. |
| `Description` | One sentence describing what the branch does or refuses to do. |
| `Trigger` | Minimum input shape, e.g. *"CmpMap(ZoomMap, ZoomMap), Series=1, matching Nin"*. |

## Fixture layout and harness

### Subdirectory

All rule-targeted fixtures live in `ast_tester/simplify_fixtures/`.
Reasons:

- `ast_tester/` already mixes ~200 files of various test types; adding
  another ~200 `.map`/`.simp` pairs makes it unnavigable.
- Grouping fixtures makes gcov scoping trivial.
- A `README.md` in the subdirectory points at the inventory and the
  harness, giving anyone landing there immediate orientation.

The three scenario fixtures (`brad`, `lsst1`, `rigby`) stay at
`ast_tester/`. They predate this work and serve as integration
regressions, not rule fixtures.

### `simplify_tests.txt` and harness

The driver file stays at `ast_tester/simplify_tests.txt`. Filename
columns gain a `simplify_fixtures/` prefix where relevant:

```
unit_series_elision | simplify_fixtures/unit_series_elision.map | simplify_fixtures/unit_series_elision.simp |
brad                | brad.map                                  | brad.simp                                 |
```

The `simplify.c` harness needs no changes — it already takes filenames
as arguments and the test driver computes paths.

### Fixture header convention

Each `.map` (and `.simp` for clarity) gets a leading `# AST_FIXTURE`
comment block (the AST native serializer ignores `#` lines):

```
# AST_FIXTURE id=winmap-04 polarity=positive
# AST_FIXTURE desc=Two adjacent ZoomMaps with matching Nin merge to product zoom
Begin CmpMap
   ...
```

This embeds traceability so a future developer reading a fixture can
find its row in the inventory.

## Campaign workflow

The work runs in five sequential phases. Each phase has an explicit
exit criterion before the next can start. All five phases land as one
PR; phases are commits (or small groups of commits) within that PR.

### Phase 1 — Inventory bootstrap

Walk `simplify_tests.txt` and the existing `.map`/`.simp` pairs. For
each existing fixture, identify the branch it targets by reading the
input + output and the relevant `MapMerge` source. Pre-populate the
inventory with one row per existing fixture, marked with the existing
filename, inferred `Type` and `Polarity`, provisional `Lines`, and
`Description`/`Trigger` derived from the source.

**Exit criterion.** Every fixture in `simplify_tests.txt` has a row.

### Phase 2 — Inventory expansion

Dispatch parallel subagents (one per source class) using the prompt:

> Read `src/<file>.c`. Locate every distinct branch inside `MapMerge`
> that produces a different simplified output (or that explicitly
> refuses to simplify due to a guard). For each branch, fill one row of
> the inventory schema. Include negative cases for every non-trivial
> guard. Do not write fixtures; produce only the inventory section.
> Cite line ranges from current `master`.

Review each subagent's output and merge into the inventory.

**Exit criterion.** Every source class listed in scope has a section,
and a spot-check of three classes confirms the line ranges and trigger
shapes are accurate.

### Phase 3 — Fixture authoring

1. Move existing rule-targeted fixtures into
   `ast_tester/simplify_fixtures/`. Rename where the existing name
   doesn't match the `<class>_<short_rule>` convention. (Existing
   fixtures were created on this branch and have no external
   references, so renaming is safe.)
2. Update `simplify_tests.txt` paths.
3. Add `# AST_FIXTURE` headers to the moved fixtures.
4. Author new fixtures for every inventory row missing one. Each
   commit references the inventory ID(s) it covers.

**Exit criterion.** Every inventory row maps to a `.map`/`.simp` pair
(or a cascade reference), and `ctest -R simplify` passes for all of
them.

### Phase 4 — gcov verification

Build with `--coverage`. Run **only** the simplify-driver tests — every
test row in `simplify_tests.txt` runs through the same `simplify.c`
harness, and the test names are the row names. The concrete invocation
is a `ctest -R` regex matching those test names; the implementation
plan picks the regex.

Generate a coverage report restricted to `MapMerge` line ranges in the
in-scope source files. Cross-reference against the inventory.

**Exit criterion.** Every line range in the inventory is hit, or each
unhit range is explicitly annotated as dead/unreachable in the
inventory with justification.

A secondary check: per-fixture coverage flags any inventory row whose
target lines are hit by more than one fixture. The three scenario
fixtures (`brad`, `lsst1`, `rigby`) are excluded from this check — they
are integration regressions that will hit many branches, and their
overlap with focused fixtures is expected, not a defect. For focused
and cascade fixtures, unexpected overlap should be reviewed but is not
a hard gate.

### Phase 5 — Permanent reference document

Convert the inventory into the permanent reference document at
`ast_tester/simplify_pathways.md` (see next section). Delete
`ast_tester/simplify_coverage.md`.

**Exit criterion.** The new reference doc lives at its permanent path,
and the design and inventory under `docs/superpowers/specs/` are
flagged as historical-only.

## Permanent reference document

### Path

`ast_tester/simplify_pathways.md`. Replaces and extends
`ast_tester/simplify_coverage.md`, which is deleted in Phase 5.

The new filename follows local naming convention (lowercase with
underscores, parallels `simplify_tests.txt`) and signals reference
material describing what `astSimplify` does, not scaffolding tracking
gaps.

### Audience

1. *AST developers* modifying a `MapMerge`. They need to know which
   branches their change touches and which fixtures will catch a
   regression.
2. *Library users* trying to predict what `astSimplify` will do to a
   pipeline they've built. They need a readable description of the
   simplification rules.

### Structure

1. **Top matter** — what `astSimplify` is, how the pipeline of
   `MapMerge` calls works, links to SUN/210, SUN/211,
   `ast_tester/simplify_fixtures/`, and this design document. About
   one page.

2. **One section per source class**, same order as the inventory.
   Each class section contains:
   - 1–2 paragraphs describing the class's simplification strategy.
   - A table with the same rows as the inventory, with two
     presentation tweaks:
     - The `Lines` column is dropped (volatile, useful only during the
       campaign).
     - Rows are grouped by polarity (positives first, negatives
       second) within each class.
   - For cascade rules, a short prose paragraph explaining why
     composition is needed.

3. **Cross-class index at the end** — alphabetical list of every
   inventory ID with a one-line description, for quick lookup from a
   `git blame` or commit message reference.

### Maintenance contract

A short section near the top: when adding or modifying a `MapMerge`
branch, the developer is expected to add a row to the relevant section
and a fixture to `ast_tester/simplify_fixtures/`. CI does not enforce
this in the initial campaign, but `simplify_tests.txt` should contain
one fixture per row.

A follow-up CI check that compares row count to fixture count is a
plausible later improvement, but is not in scope here.

## Risks and mitigations

- **Subagent under-counts branches.** Spot-check three classes by
  hand-reading after Phase 2. If the gap exceeds 10% in any class,
  re-dispatch with a stricter prompt.

- **Fixture isolation drift** — multiple fixtures accidentally hit the
  same target branch, weakening bisection. Phase 4 produces per-fixture
  coverage; flag any inventory row whose target lines are hit by more
  than one fixture for review.

- **Cascade fixtures grow large.** Cap cascade fixtures at three
  composed mappings; if a cascade needs more, document why in the
  inventory row.

- **Branches that are genuinely unreachable.** Annotate in the
  inventory with justification; do not silently leave them uncovered.

## Deliverables summary

- `ast_tester/simplify_fixtures/` directory containing all
  rule-targeted `.map`/`.simp` pairs (existing ones moved/renamed,
  new ones authored) with `# AST_FIXTURE` headers.
- `ast_tester/simplify_tests.txt` updated to reference the new
  locations.
- `ast_tester/simplify_pathways.md` as the permanent reference document.
- `ast_tester/simplify_coverage.md` deleted.
- `docs/superpowers/specs/2026-05-05-simplify-branch-coverage-inventory.md`
  as the historical inventory artifact.
- This design document, retained as the historical record of the
  campaign.
