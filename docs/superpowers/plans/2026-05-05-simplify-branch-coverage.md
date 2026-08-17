# `astSimplify` Branch Coverage Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a complete branch-level fixture set for `astSimplify`, plus a permanent reference document mapping every `MapMerge` pathway to the fixture that pins it.

**Architecture:** Five sequential phases — inventory bootstrap, inventory expansion via parallel subagents, fixture authoring, gcov verification, permanent reference document. All work lands in a single PR on branch `u/timj/simplify-testing`. Fixtures live under `ast_tester/simplify_fixtures/`; each is a `.map`/`.simp` pair driven by the existing `simplify.c` harness through `add_simplify_test()` in `ast_tester/CMakeLists.txt`.

**Tech Stack:** AST native serialization format, `simplify.c` harness, CMake/ctest, gcov for coverage verification, Markdown for inventory and reference docs.

**Companion design doc:** `docs/superpowers/specs/2026-05-05-simplify-branch-coverage-design.md`. Read it first.

---

## Important context for the executor

This plan is unusual. It doesn't add new product code; it adds **regression fixtures for existing behavior**. So conventional TDD ("write failing test → write code to make it pass") does not map cleanly. The adapted cycle is:

1. **Pick an inventory row.** Read its trigger shape and description.
2. **Author an input mapping (`.map`).** Construct the smallest pipeline that triggers the target branch.
3. **Run the harness against the input.** This produces the simplified output.
4. **Inspect the output.** Verify it matches the inventory description. If yes → save as `.simp`, commit. If no → either the input shape was wrong (fix and retry), the inventory description was wrong (fix the row), or there's a real simplification bug (file a follow-up issue, document, move on).
5. **Wire into `simplify_tests.txt`.** Run `ctest` to confirm.

There is **no implementation code to write**. The "code" is the AST native dump format. Treat fixture authoring as data engineering, not programming.

**Hand-authoring vs generator code.** Most fixtures are 5–20 line text files and can be hand-authored. For complex inputs (e.g. PolyMap with many coefficients, LutMap with long lookup tables), write a throwaway C program in a scratch directory, run it to emit the `.map` via `astShow`, commit the resulting fixture, and discard the generator. The design doc allows this; the existing `brad`/`lsst1`/`rigby` are likely generated this way.

**Native serialization gotcha.** Tests run against native AST encoding. C99 mode uses `AST_DBL_DIG = 18`; C11 mode uses 17. The harness is built at the project's default C standard. If a fixture's `.simp` has float values that vary at ulp level between platforms, mark the test row in `simplify_tests.txt` with `skip_string_compare = yes` so only the `astequal` check runs (the existing `rigby` row is the precedent).

**`# AST_FIXTURE` headers.** Lines starting with `#` are ignored by the AST channel parser. Add header comments at the top of every `.map` and `.simp` file in `simplify_fixtures/` to record the inventory ID and description. Format:

```
# AST_FIXTURE id=winmap-04 polarity=positive
# AST_FIXTURE desc=Two adjacent ZoomMaps with matching Nin merge to product zoom
```

---

## File structure

**Files created:**
- `docs/superpowers/specs/2026-05-05-simplify-branch-coverage-inventory.md` — working inventory, single source of truth during the campaign.
- `ast_tester/simplify_fixtures/` — directory holding all rule-targeted `.map`/`.simp` pairs.
- `ast_tester/simplify_fixtures/README.md` — short orientation doc pointing at the harness, inventory, and naming convention.
- `ast_tester/simplify_pathways.md` — permanent reference document (replaces `simplify_coverage.md`).
- ~150–250 new `.map`/`.simp` pairs under `ast_tester/simplify_fixtures/`.

**Files modified:**
- `ast_tester/simplify_tests.txt` — update existing rows with new paths; append rows for new fixtures.

**Files moved/renamed:**
- All ~50 existing rule-targeted fixtures (those created on this branch) → `ast_tester/simplify_fixtures/`. Rename when needed to match `<class>_<short_rule>` convention.

**Files deleted at end:**
- `ast_tester/simplify_coverage.md` — replaced by `simplify_pathways.md`.

**Files unchanged:**
- `ast_tester/simplify.c` — harness needs no changes.
- `ast_tester/CMakeLists.txt` — `add_simplify_test()` uses paths from `simplify_tests.txt`, so nothing to modify (already supports relative paths).
- `ast_tester/brad.{map,simp}`, `lsst1.{map,simp}`, `rigby.{map,simp}` — pre-branch scenario fixtures, stay where they are.
- All `src/*.c` — no source code is modified by this campaign.

---

## Phase 1 — Inventory bootstrap

**Goal.** Produce the inventory file with one row per existing fixture in `simplify_tests.txt`, so subagents in Phase 2 see what's already covered and can concentrate on gaps.

### Task 1: Create the inventory skeleton

**Files:**
- Create: `docs/superpowers/specs/2026-05-05-simplify-branch-coverage-inventory.md`

- [ ] **Step 1: Create the inventory file with section headers and schema.**

Write the file with one `##` section per source class listed below. Each section has the placeholder paragraph and an empty schema table. The exact contents:

```markdown
# `astSimplify` Branch Coverage Inventory

Working inventory for the campaign described in
`docs/superpowers/specs/2026-05-05-simplify-branch-coverage-design.md`.

Each row records one branch in a `MapMerge` method (or in the
`astSimplify` orchestration loop). The fixture column links the branch
to its `.map`/`.simp` pair under `ast_tester/simplify_fixtures/` (or
`via <id>` for branches reached only through a cascade).

## Schema

| Column | Meaning |
| --- | --- |
| ID | `<class>-NN` stable identifier |
| Fixture | `.map` filename, or `via <id>` for cascade-only branches |
| Type | `focused` \| `cascade` \| `scenario` |
| Polarity | `positive` \| `negative` |
| Lines | Source line range, e.g. `winmap.c:2104-2147` |
| Description | One sentence describing the branch behavior |
| Trigger | Minimum input shape needed to fire the branch |

---

## mapping.c — astSimplify orchestration

(One paragraph describing how `astSimplify` loops over `MapMerge`.)

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
| --- | --- | --- | --- | --- | --- | --- |

## cmpmap.c

(One paragraph describing CmpMap's MapMerge strategy.)

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
| --- | --- | --- | --- | --- | --- | --- |

## unitmap.c

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
| --- | --- | --- | --- | --- | --- | --- |

## zoommap.c

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
| --- | --- | --- | --- | --- | --- | --- |

## shiftmap.c

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
| --- | --- | --- | --- | --- | --- | --- |

## winmap.c

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
| --- | --- | --- | --- | --- | --- | --- |

## matrixmap.c

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
| --- | --- | --- | --- | --- | --- | --- |

## permmap.c

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
| --- | --- | --- | --- | --- | --- | --- |

## lutmap.c

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
| --- | --- | --- | --- | --- | --- | --- |

## polymap.c

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
| --- | --- | --- | --- | --- | --- | --- |

## mathmap.c

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
| --- | --- | --- | --- | --- | --- | --- |

## ratemap.c

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
| --- | --- | --- | --- | --- | --- | --- |

## intramap.c

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
| --- | --- | --- | --- | --- | --- | --- |

## splinemap.c

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
| --- | --- | --- | --- | --- | --- | --- |

## normmap.c

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
| --- | --- | --- | --- | --- | --- | --- |

## unitnormmap.c

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
| --- | --- | --- | --- | --- | --- | --- |

## selectormap.c

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
| --- | --- | --- | --- | --- | --- | --- |

## switchmap.c

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
| --- | --- | --- | --- | --- | --- | --- |

## tranmap.c

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
| --- | --- | --- | --- | --- | --- | --- |

## timemap.c

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
| --- | --- | --- | --- | --- | --- | --- |

## slamap.c

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
| --- | --- | --- | --- | --- | --- | --- |

## specmap.c

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
| --- | --- | --- | --- | --- | --- | --- |

## wcsmap.c

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
| --- | --- | --- | --- | --- | --- | --- |

## sphmap.c

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
| --- | --- | --- | --- | --- | --- | --- |

## pcdmap.c

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
| --- | --- | --- | --- | --- | --- | --- |

## dssmap.c

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
| --- | --- | --- | --- | --- | --- | --- |

## grismmap.c

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
| --- | --- | --- | --- | --- | --- | --- |

## xphmap.c

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
| --- | --- | --- | --- | --- | --- | --- |

## box.c (Region as Mapping)

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
| --- | --- | --- | --- | --- | --- | --- |

## interval.c (Region as Mapping)

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
| --- | --- | --- | --- | --- | --- | --- |

## nullregion.c (Region as Mapping)

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
| --- | --- | --- | --- | --- | --- | --- |

## pointlist.c (Region as Mapping)

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
| --- | --- | --- | --- | --- | --- | --- |

## prism.c (Region as Mapping)

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
| --- | --- | --- | --- | --- | --- | --- |
```

- [ ] **Step 2: Commit.**

```bash
git add docs/superpowers/specs/2026-05-05-simplify-branch-coverage-inventory.md
git commit -m "Add empty inventory skeleton for simplify branch coverage"
```

### Task 2: Bootstrap rows for existing fixtures

**Files:**
- Modify: `docs/superpowers/specs/2026-05-05-simplify-branch-coverage-inventory.md`

For each row in `ast_tester/simplify_tests.txt`, identify which class section it belongs in, read the `.map` and `.simp` to understand the targeted branch, locate the corresponding code in `src/<class>.c`, and add a row to the inventory.

- [ ] **Step 1: List the existing fixtures and group by class.**

Run:
```bash
cat ast_tester/simplify_tests.txt | grep -v '^#' | grep -v '^$'
```

Group each non-comment row by the inferred class. Examples (from existing names):
- `unit_*` → `unitmap.c`
- `zoom_*` → `zoommap.c`
- `matrix_*` → `matrixmap.c`
- `win_*` → `winmap.c` (and some win_matrix_* might also be a matrixmap.c branch — pick the class that owns the branch the fixture targets)
- `perm_*` → `permmap.c`
- `lut_*` → `lutmap.c`
- `cmpmap_*` → `cmpmap.c`
- `shift_*` → `shiftmap.c`
- `sph_*` → `sphmap.c`
- `wcs_*` → `wcsmap.c`
- `sla_*` → `slamap.c`
- `spec_*` → `specmap.c`
- `time_*` → `timemap.c`
- `pcd_*` → `pcdmap.c`
- `grism_*` → `grismmap.c`
- `poly_*` → `polymap.c`
- `spline_*` → `splinemap.c`
- `math_*` → `mathmap.c`
- `ratemap_*` → `ratemap.c`
- `tranmap_*` → `tranmap.c`
- `intramap_*` → `intramap.c`
- `brad`, `lsst1`, `rigby` → scenario; place under whichever class they exercise primarily (mapping.c is fine — they exercise the orchestration loop). Mark `Type=scenario`.

- [ ] **Step 2: For each fixture, read the `.map`, the `.simp`, and the relevant `MapMerge` source to identify the branch.**

For example, for `zoom_series_merge`:
```bash
cat ast_tester/zoom_series_merge.map
cat ast_tester/zoom_series_merge.simp
grep -n "MapMerge\|series\|Series" src/zoommap.c | head -20
```

Read the surrounding `MapMerge` code to find the branch that handles two adjacent ZoomMaps in series. Note the line range.

- [ ] **Step 3: Add a row to the inventory under the appropriate class.**

Use ID `<class>-01`, `<class>-02`, ... in the order rows are added. Example:

```
| zoommap-01 | zoom_series_merge.map | focused | positive | zoommap.c:1234-1267 | Two adjacent ZoomMaps with matching Nin merge to single ZoomMap with product zoom | CmpMap(ZoomMap, ZoomMap), Series=1, matching Nin |
```

For now, use the existing fixture filename (no path prefix) — it will be updated to `simplify_fixtures/<name>.map` in Phase 3.

- [ ] **Step 4: Repeat for every existing fixture.** Work through the list systematically. If a fixture is unclear (you cannot identify the branch with confidence), add it with a `???` placeholder in `Lines` and a note in `Description`; flag those for review at the end of the task.

- [ ] **Step 5: Verify every existing fixture has exactly one inventory row.**

```bash
# Count fixture rows in simplify_tests.txt
grep -v '^#' ast_tester/simplify_tests.txt | grep -v '^$' | wc -l

# Count inventory rows referencing existing fixtures (those with `.map` filenames not yet prefixed)
grep -E '\| [a-z_0-9]+\.map \|' docs/superpowers/specs/2026-05-05-simplify-branch-coverage-inventory.md | wc -l
```

Both numbers should match.

- [ ] **Step 6: Resolve any `???` rows by hand-reading or adjust the row to read more honestly (e.g. "scenario fixture, hits multiple branches, primary branch unknown — to be refined post-Phase 2").**

- [ ] **Step 7: Commit.**

```bash
git add docs/superpowers/specs/2026-05-05-simplify-branch-coverage-inventory.md
git commit -m "Bootstrap inventory rows for existing simplify fixtures"
```

---

## Phase 2 — Inventory expansion via parallel subagents

**Goal.** For every source class listed in scope, ensure the inventory has a row for every distinct branch in that class's `MapMerge` (and for `mapping.c`'s orchestration loop), including negatives.

### Task 3: Define and document the subagent dispatch template

**Files:**
- Modify: `docs/superpowers/specs/2026-05-05-simplify-branch-coverage-inventory.md` (add a comment block at the top recording the dispatch template — for traceability if a row's provenance is later questioned)

- [ ] **Step 1: Add a dispatch-template comment block under the schema section.**

Insert after the schema table:

```markdown
## Subagent dispatch template

Each `MapMerge` (and `mapping.c`'s `astSimplify` loop) was extended by
dispatching a subagent with the following prompt. Existing rows from
the bootstrap step are pre-populated in each class section so the
subagent does not duplicate them.

> Read `src/<file>.c`. Locate every distinct branch inside `MapMerge`
> (or, for `mapping.c`, inside the `astSimplify` orchestration loop)
> that produces a different simplified output, or that explicitly
> refuses to simplify due to a guard. For each branch, append one row
> to the inventory section for this class, using the schema in this
> file. Do not modify rows already present (those are the bootstrap
> rows for existing fixtures). Include one negative row for every
> non-trivial guard. Cite line ranges from current `master`. Use ID
> `<class>-NN` continuing from the highest existing N. Leave Fixture
> blank for newly-discovered branches; Phase 3 will fill it in.
```

- [ ] **Step 2: Commit.**

```bash
git add docs/superpowers/specs/2026-05-05-simplify-branch-coverage-inventory.md
git commit -m "Document subagent dispatch template in inventory"
```

### Task 4: Dispatch subagents — Family 1 (core algebraic mappings)

**Goal.** Expand inventory sections for: `mapping.c`, `unitmap.c`, `zoommap.c`, `shiftmap.c`, `winmap.c`, `matrixmap.c`, `permmap.c`, `lutmap.c`, `cmpmap.c`.

These are dispatched in parallel because they're independent (each subagent reads one file).

- [ ] **Step 1: Dispatch nine `Explore` subagents in parallel, one per class.**

Use the prompt template above, specialized per class. Example prompt for `winmap.c`:

> Read `/Users/timj/work/starlink-ast/src/winmap.c`. Focus on the `MapMerge` static function. Locate every distinct branch inside `MapMerge` that produces a different simplified output, or that explicitly refuses to simplify due to a guard. For each branch, return one row in the schema:
>
> `| ID | Fixture | Type | Polarity | Lines | Description | Trigger |`
>
> The current inventory already contains these rows for existing fixtures: `<list winmap rows from current inventory state>`. Do not duplicate those — start your IDs from `winmap-N+1`. Include negative rows for non-trivial guards. Cite line ranges from current `master`. Leave Fixture blank for new rows. Return only the rows; do not write fixtures or modify any file.

Run all nine in a single response with parallel `Agent` tool calls. (`Explore` is the right subagent type — read-only, fast.)

- [ ] **Step 2: For each subagent's response, append rows to the appropriate inventory section.**

After all subagents return, edit the inventory file. For each class section, append the new rows under the existing bootstrap rows. Maintain ID continuity per class.

- [ ] **Step 3: Sanity-check counts.**

Compare the row count per class against an order-of-magnitude expectation. For example, `winmap.c` has dozens of `MapMerge` branches based on the `grep -c` from the design phase (42 internal functions/blocks); expect 30–60 rows. If a subagent returns fewer than 10 rows for a complex class, re-dispatch with a stricter prompt asking it to enumerate every conditional in the function.

- [ ] **Step 4: Commit.**

```bash
git add docs/superpowers/specs/2026-05-05-simplify-branch-coverage-inventory.md
git commit -m "Expand inventory: core algebraic mapping classes"
```

### Task 5: Dispatch subagents — Family 2 (compound and specialty)

**Goal.** Expand inventory sections for: `polymap.c`, `mathmap.c`, `ratemap.c`, `intramap.c`, `splinemap.c`, `normmap.c`, `unitnormmap.c`, `selectormap.c`, `switchmap.c`, `tranmap.c`.

- [ ] **Step 1: Dispatch ten `Explore` subagents in parallel using the same prompt template as Task 4, specialized per class.**

- [ ] **Step 2: Append returned rows to the inventory.**

- [ ] **Step 3: Sanity-check counts.**

- [ ] **Step 4: Commit.**

```bash
git add docs/superpowers/specs/2026-05-05-simplify-branch-coverage-inventory.md
git commit -m "Expand inventory: compound and specialty mapping classes"
```

### Task 6: Dispatch subagents — Family 3 (sky, spectral, time, projections)

**Goal.** Expand inventory sections for: `timemap.c`, `slamap.c`, `specmap.c`, `wcsmap.c`, `sphmap.c`, `pcdmap.c`, `dssmap.c`, `grismmap.c`, `xphmap.c`.

- [ ] **Step 1: Dispatch nine `Explore` subagents in parallel.**

- [ ] **Step 2: Append returned rows to the inventory.**

- [ ] **Step 3: Sanity-check counts.**

- [ ] **Step 4: Commit.**

```bash
git add docs/superpowers/specs/2026-05-05-simplify-branch-coverage-inventory.md
git commit -m "Expand inventory: sky, spectral, time, and projection classes"
```

### Task 7: Dispatch subagents — Family 4 (Region-as-Mapping)

**Goal.** Expand inventory sections for: `box.c`, `interval.c`, `nullregion.c`, `pointlist.c`, `prism.c`. The subagent prompt for these is slightly different — `MapMerge` here is the Region's behavior when treated as a Mapping in a pipeline, not the Region's `Simplify` override (which is out of scope).

- [ ] **Step 1: Dispatch five `Explore` subagents in parallel using the modified prompt:**

> Read `/Users/timj/work/starlink-ast/src/<file>.c`. Focus on the `MapMerge` static function (the Region-as-Mapping merge logic — NOT the Region's own `Simplify` override). Locate every distinct branch in `MapMerge` that produces a different simplified output or explicitly refuses to simplify due to a guard. For each branch, return one row in the schema. Include negative rows for non-trivial guards. Cite line ranges from current `master`. Leave Fixture blank. Return only the rows.

- [ ] **Step 2: Append returned rows.**

- [ ] **Step 3: Commit.**

```bash
git add docs/superpowers/specs/2026-05-05-simplify-branch-coverage-inventory.md
git commit -m "Expand inventory: Region-as-Mapping classes"
```

### Task 8: Spot-check inventory accuracy

**Goal.** Verify subagent output is trustworthy before committing to authoring fixtures against it.

- [ ] **Step 1: Pick three classes by hand for spot-checking.** Choose one from each family — e.g. `winmap.c`, `polymap.c`, `wcsmap.c`.

- [ ] **Step 2: For each picked class, hand-read the `MapMerge` and confirm:**

(a) Every distinct branch the subagent identified actually exists at the cited lines.
(b) The trigger shape would actually fire that branch.
(c) No obvious branch is missing (read top to bottom; flag anything not in the inventory).

- [ ] **Step 3: If any class has more than 10% of branches missing or mis-identified, re-dispatch the subagent with a stricter prompt that includes the specific gap as a counter-example.**

- [ ] **Step 4: Iterate on Step 3 until the spot-checked classes are accurate.**

- [ ] **Step 5: Commit any inventory corrections.**

```bash
git add docs/superpowers/specs/2026-05-05-simplify-branch-coverage-inventory.md
git commit -m "Correct inventory rows after spot-check review"
```

---

## Phase 3 — Fixture authoring

**Goal.** Every inventory row maps to a `.map`/`.simp` pair in `ast_tester/simplify_fixtures/` (or to a cascade fixture via the `via <id>` syntax), with `# AST_FIXTURE` headers, registered in `simplify_tests.txt`, passing under `ctest`.

### Task 9: Create the fixtures subdirectory and README

**Files:**
- Create: `ast_tester/simplify_fixtures/`
- Create: `ast_tester/simplify_fixtures/README.md`

- [ ] **Step 1: Create the directory and README.**

```bash
mkdir -p ast_tester/simplify_fixtures
```

Write `ast_tester/simplify_fixtures/README.md`:

```markdown
# `astSimplify` rule-targeted fixtures

Each `.map` is an AST native dump of a Mapping pipeline. The matching
`.simp` is the expected output of `astSimplify` applied to that input.

The harness is `ast_tester/simplify.c`. Tests are registered in
`ast_tester/simplify_tests.txt` and built by `ast_add_simplify_test()`
in `ast_tester/CMakeLists.txt`.

For the catalogue of which fixture pins which simplification rule, see
`ast_tester/simplify_pathways.md`.

Each fixture begins with `# AST_FIXTURE` comment lines recording its
inventory ID, polarity, and a short description. The AST channel
parser ignores `#` lines.

The three pre-existing scenario fixtures `brad`, `lsst1`, and `rigby`
live at `ast_tester/` (one level up). They are integration regressions,
not rule fixtures.
```

- [ ] **Step 2: Commit.**

```bash
git add ast_tester/simplify_fixtures/README.md
git commit -m "Create simplify_fixtures subdirectory with README"
```

### Task 10: Move and rename existing rule-targeted fixtures

**Goal.** Move all existing rule-targeted fixtures created on this branch into `simplify_fixtures/`, renaming any that don't match the `<class>_<short_rule>` convention.

- [ ] **Step 1: Build the move/rename map from the inventory.**

For each existing fixture (rows with a non-empty `Fixture` column from Phase 1), look at the row's class and rule and decide if the fixture name needs to change. Most are already correct (e.g. `zoom_series_merge` → keep). Some may benefit from class-prefix consistency (e.g. `shift_invert_normalize` → `shiftmap_invert_normalize` if you prefer `shiftmap` everywhere). Document each rename in the commit message.

- [ ] **Step 2: Move every existing rule-targeted fixture using `git mv`.**

```bash
# Run from repo root
for name in unit_series_elision unit_parallel_merge zoom_series_merge \
            zoom_series_cancel zoom_parallel_to_matrix \
            matrix_unit_to_unit matrix_diagonal_to_zoom \
            matrix_full_to_diagonal win_to_shift win_to_matrix \
            lut_linear_to_win lut_inverse_cancel shift_invert_normalize \
            perm_series_merge perm_parallel_merge \
            matrix_matrix_series_merge matrix_zoom_series_merge \
            matrix_perm_series_merge win_shift_series_merge \
            win_zoom_series_merge win_matrix_series_merge \
            win_perm_swap_merge cmpmap_parallel_series_components \
            sph_inverse_cancel wcs_inverse_cancel sla_inverse_cancel \
            spec_inverse_cancel time_inverse_cancel \
            cmpmap_nested_parallel_flatten cmpmap_perm_parallel_swap \
            win_parallel_merge sph_matrix_sandwich pcd_zero_to_unit \
            pcd_inverse_cancel grism_inverse_cancel grism_zoom_merge \
            poly_duplicate_terms poly_inverse_cancel \
            spline_inverse_cancel math_inverse_cancel \
            ratemap_inverse_cancel tranmap_equal_components \
            intramap_inverse_cancel; do
    git mv "ast_tester/${name}.map" "ast_tester/simplify_fixtures/${name}.map"
    git mv "ast_tester/${name}.simp" "ast_tester/simplify_fixtures/${name}.simp"
done
```

(Apply renames inline if the inventory specified any. The list above is taken from `simplify_tests.txt` minus the three scenario rows — verify against the current state of that file before running.)

- [ ] **Step 3: Update `simplify_tests.txt` paths.**

Edit `ast_tester/simplify_tests.txt`. For every row except `brad`, `lsst1`, `rigby`, prepend `simplify_fixtures/` to the `.map` and `.simp` filenames:

```
unit_series_elision     | simplify_fixtures/unit_series_elision.map     | simplify_fixtures/unit_series_elision.simp     |
```

- [ ] **Step 4: Build and run the simplify tests to verify the move didn't break anything.**

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Debug 2>&1 | tail -5
cmake --build build 2>&1 | tail -10
ctest --test-dir build -R '^simplify_' --output-on-failure 2>&1 | tail -20
```

Expected: all `simplify_*` tests pass.

If any test fails because of the move, the most likely cause is `simplify_tests.txt` having a typo in a path or a fixture file having been missed in the loop. Fix and rerun until green.

- [ ] **Step 5: Commit.**

```bash
git add ast_tester/simplify_tests.txt ast_tester/simplify_fixtures/
git commit -m "Move existing rule-targeted simplify fixtures to simplify_fixtures/"
```

### Task 11: Add `# AST_FIXTURE` headers to moved fixtures

**Goal.** Each moved fixture's `.map` and `.simp` gains the two-line `# AST_FIXTURE` header recording its inventory ID, polarity, and description.

- [ ] **Step 1: For each moved fixture, look up its inventory row.**

The inventory row gives you `ID`, `Polarity`, `Description`.

- [ ] **Step 2: Prepend the header lines to both files.**

Example for `simplify_fixtures/zoom_series_merge.map`:

```
# AST_FIXTURE id=zoommap-01 polarity=positive
# AST_FIXTURE desc=Two adjacent ZoomMaps with matching Nin merge to single ZoomMap with product zoom
Begin CmpMap
   ...
```

The corresponding `.simp` gets the same headers (the simplify harness preserves leading `#` comments through the Channel; if a build-time check fails because the `.simp` produced by the harness doesn't include the headers verbatim, the headers go in the `.map` only — the comparison is against the harness output which would not include them).

**Verification step before applying everywhere:** apply the headers to *one* fixture (`zoom_series_merge`), rebuild, and run that one test:

```bash
cmake --build build 2>&1 | tail -5
ctest --test-dir build -R '^simplify_zoom_series_merge' --output-on-failure
```

If it passes with headers in both `.map` and `.simp`, apply to all. If it fails because the harness output lacks the headers (likely), apply headers to `.map` only and re-verify.

- [ ] **Step 3: Apply the verified header convention to all moved fixtures.**

Create a one-off script `/tmp/apply_fixture_headers.py` (do not commit):

```python
#!/usr/bin/env python3
"""Read the inventory and prepend AST_FIXTURE headers to each fixture."""
import re, pathlib

INV = pathlib.Path("docs/superpowers/specs/2026-05-05-simplify-branch-coverage-inventory.md")
FIXDIR = pathlib.Path("ast_tester/simplify_fixtures")

# Parse rows of the form: | id | fixture.map | type | polarity | lines | desc | trigger |
row_re = re.compile(r"^\|\s*([a-z0-9_-]+)\s*\|\s*([a-z0-9_]+\.map)\s*\|\s*\w+\s*\|\s*(\w+)\s*\|[^|]*\|\s*([^|]+?)\s*\|", re.M)

for m in row_re.finditer(INV.read_text()):
    inv_id, fname, polarity, desc = m.group(1), m.group(2), m.group(3), m.group(4)
    map_path = FIXDIR / fname
    if not map_path.exists():
        continue
    body = map_path.read_text()
    if body.startswith("# AST_FIXTURE"):
        continue  # already done
    header = f"# AST_FIXTURE id={inv_id} polarity={polarity}\n# AST_FIXTURE desc={desc}\n"
    map_path.write_text(header + body)
```

Run:
```bash
python3 /tmp/apply_fixture_headers.py
```

Apply the same header to the matching `.simp` only if Step 2's verification confirmed the harness preserves them through round-trip. If not, leave `.simp` files unmodified.

- [ ] **Step 4: Rebuild and run all simplify tests.**

```bash
cmake --build build 2>&1 | tail -5
ctest --test-dir build -R '^simplify_' --output-on-failure 2>&1 | tail -20
```

Expected: all pass.

- [ ] **Step 5: Commit.**

```bash
git add ast_tester/simplify_fixtures/
git commit -m "Add AST_FIXTURE traceability headers to moved fixtures"
```

### Task 12: Author new fixtures (the bulk of the work)

**Goal.** Every inventory row that has a blank `Fixture` column gets a `.map`/`.simp` pair authored, registered, and passing.

This task is iterative. Process inventory rows one class at a time, in the same order as the inventory sections. Within a class, do positives before negatives. Within positives, do focused before cascade.

**Cap on cascade fixtures.** A cascade fixture composes at most three component mappings. If an inventory row needs more than three components to trigger, document why in the row's `Trigger` cell and treat it as scenario-class instead of cascade.

**Per-fixture procedure** (apply for every blank-Fixture row in the inventory):

- [ ] **Step 1: Read the inventory row.** Note ID, polarity, description, trigger shape.

- [ ] **Step 2: Construct the input `.map` file.**

Use the trigger shape to write the smallest possible AST native dump that matches. For simple maps, hand-author. For complex maps (PolyMap with many coefficients, large LutMap), write a throwaway C program in `/tmp/` that constructs the mapping, calls `astShow(map)`, and saves the output. Do not commit the generator. Example skeleton:

```c
/* /tmp/gen_polymap.c — discard after use */
#include "ast.h"
#include <stdio.h>
int main(void) {
    int status = 0;
    astWatch(&status);
    AstPolyMap *p = astPolyMap(...);
    FILE *f = fopen("ast_tester/simplify_fixtures/polymap_NN.map", "w");
    /* Use astChannel with SinkFile to write native dump */
    AstChannel *ch = astChannel(NULL, NULL, "SinkFile=ast_tester/simplify_fixtures/polymap_NN.map");
    astWrite(ch, p);
    return 0;
}
```

Add the `# AST_FIXTURE` headers to the generated `.map`.

- [ ] **Step 3: Run the harness to produce candidate `.simp`.**

```bash
cd build  # or wherever the simplify binary lives after cmake --build
./ast_tester/simplify ../ast_tester/simplify_fixtures/<name>.map /tmp/<name>.simp.candidate
cat /tmp/<name>.simp.candidate
```

- [ ] **Step 4: Verify the candidate matches the inventory description.**

Read the candidate output. Confirm:
- For positive fixtures: the structure is the simpler form described in the inventory row.
- For negative fixtures: the structure is essentially unchanged from the input (only `IsSimp = 1` flag may be added).

If the candidate doesn't match: either the input shape was wrong (revise `.map`), the inventory row's description/trigger was wrong (revise the inventory), or there is a real simplification bug. In the third case, document the bug as a follow-up issue, mark the inventory row with a `BUG-NNN` reference, and skip the fixture for now.

- [ ] **Step 5: Save the verified candidate as the reference `.simp`.**

```bash
cp /tmp/<name>.simp.candidate ast_tester/simplify_fixtures/<name>.simp
```

Add the matching `# AST_FIXTURE` headers (or skip them in `.simp` if the verification in Task 11 found that the harness doesn't preserve them).

- [ ] **Step 6: Update the inventory row with the fixture filename.**

Replace the blank `Fixture` cell with `<name>.map`. For cascade fixtures, also update the rows for branches that are reached only through this cascade with `via <cascade-id>`.

- [ ] **Step 7: Add the row to `simplify_tests.txt`.**

```
<name> | simplify_fixtures/<name>.map | simplify_fixtures/<name>.simp |
```

(Add `yes` in the fourth column if floating-point output varies at ulp level across platforms.)

- [ ] **Step 8: Build and run the new test.**

```bash
cmake --build build 2>&1 | tail -5
ctest --test-dir build -R "^simplify_<name>$|^simplify_<name>_astequal$" --output-on-failure
```

Expected: both tests pass.

- [ ] **Step 9: Commit.**

```bash
git add ast_tester/simplify_fixtures/<name>.map \
        ast_tester/simplify_fixtures/<name>.simp \
        ast_tester/simplify_tests.txt \
        docs/superpowers/specs/2026-05-05-simplify-branch-coverage-inventory.md
git commit -m "Add fixture <name> for inventory <ID>"
```

**Batching commits.** A single commit per fixture is heavy. After verifying each fixture works, batch commits by class section (e.g. one commit per class with a message like "Add fixtures for zoommap MapMerge branches"). Per-fixture verification still happens; only the commit boundary changes.

**Tracking progress.** After each class section is complete, run:

```bash
grep -c '^|' docs/superpowers/specs/2026-05-05-simplify-branch-coverage-inventory.md
grep -c '| simplify_fixtures/' docs/superpowers/specs/2026-05-05-simplify-branch-coverage-inventory.md
```

The two numbers should converge as the campaign progresses.

- [ ] **Step 10: When every inventory row has a fixture, run the full simplify suite.**

```bash
ctest --test-dir build -R '^simplify_' --output-on-failure 2>&1 | tail -20
```

All tests must pass. If any fail, fix the failing fixture before proceeding to Phase 4.

---

## Phase 4 — gcov verification

**Goal.** Confirm that every line range in the inventory is hit by at least one fixture, and that no `MapMerge` line in the in-scope source files is unexpectedly uncovered.

### Task 13: Build with coverage instrumentation

**Files:** none modified — uses build flags only.

- [ ] **Step 1: Configure a coverage build.**

```bash
cmake -B build-coverage \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_C_FLAGS="--coverage -O0" \
  -DCMAKE_EXE_LINKER_FLAGS="--coverage" \
  -DCMAKE_SHARED_LINKER_FLAGS="--coverage"
cmake --build build-coverage 2>&1 | tail -10
```

Expected: build completes without errors. `.gcno` files appear next to object files.

- [ ] **Step 2: Run only the simplify-driver tests with coverage data accumulating.**

```bash
ctest --test-dir build-coverage -R '^simplify_' --output-on-failure 2>&1 | tail -20
```

Expected: all `simplify_*` tests pass (same as Phase 3 final run, but now `.gcda` files are produced).

### Task 14: Generate coverage reports per source file

**Files:**
- Create: `build-coverage/coverage_reports/` (working directory; not committed)

- [ ] **Step 1: For each in-scope source file, generate a per-file gcov report.**

First, locate the directory CMake puts AST object files in:

```bash
GCNO_DIR=$(find build-coverage -name "cmpmap.c.gcno" -exec dirname {} \; | head -1)
echo "$GCNO_DIR"
```

Then generate per-file reports:

```bash
mkdir -p build-coverage/coverage_reports
cd build-coverage
for src in cmpmap unitmap zoommap shiftmap winmap matrixmap permmap \
           lutmap polymap mathmap ratemap intramap splinemap normmap \
           unitnormmap selectormap switchmap tranmap timemap slamap \
           specmap wcsmap sphmap pcdmap dssmap grismmap xphmap \
           mapping box interval nullregion pointlist prism; do
    gcov -o "${GCNO_DIR#build-coverage/}" "../src/${src}.c" \
         > "coverage_reports/${src}.gcov.summary" 2>&1
    mv "${src}.c.gcov" "coverage_reports/" 2>/dev/null \
      || echo "WARNING: no gcov output for ${src}.c"
done
cd ..
```

If any file emits "WARNING: no gcov output", inspect with `find build-coverage -name "${src}.c.gcda"` to see if the source was actually exercised.

- [ ] **Step 2: For each `MapMerge` line range in the inventory, check if every line is covered.**

A `.gcov` line looks like:
```
        5:  1234:    if ( astOK ) {
    #####:  1235:        do_something();
        -:  1236:    /* comment */
```

Field 1 is the hit count (`#####` = uncovered, `-` = non-executable like comments/blanks), field 2 is the source line number, the rest is source code. We need to flag lines where field 1 is `#####` and the line number falls within an inventory range.

Use this Python script `/tmp/check_coverage.py` (do not commit):

```python
#!/usr/bin/env python3
"""Cross-reference inventory line ranges against gcov output."""
import re, pathlib, sys

INV = pathlib.Path("docs/superpowers/specs/2026-05-05-simplify-branch-coverage-inventory.md")
GCOV_DIR = pathlib.Path("build-coverage/coverage_reports")

# Parse "lines" field of inventory rows, e.g. "winmap.c:2104-2147"
range_re = re.compile(r"([a-z]+)\.c:(\d+)-(\d+)")
ranges = {}
for m in range_re.finditer(INV.read_text()):
    src, start, end = m.group(1), int(m.group(2)), int(m.group(3))
    ranges.setdefault(src, []).append((start, end))

uncovered = []
for src, line_ranges in ranges.items():
    gcov_file = GCOV_DIR / f"{src}.c.gcov"
    if not gcov_file.exists():
        print(f"MISSING: {gcov_file}", file=sys.stderr)
        continue
    for line in gcov_file.read_text().splitlines():
        parts = line.split(":", 2)
        if len(parts) < 2:
            continue
        count, lineno_str = parts[0].strip(), parts[1].strip()
        try:
            lineno = int(lineno_str)
        except ValueError:
            continue
        if count != "#####":
            continue
        for start, end in line_ranges:
            if start <= lineno <= end:
                uncovered.append((src, lineno, line.split(":", 2)[2]))
                break

for src, lineno, code in uncovered:
    print(f"{src}.c:{lineno}: {code}")
```

Run:
```bash
python3 /tmp/check_coverage.py > /tmp/uncovered_lines.txt
cat /tmp/uncovered_lines.txt
```

- [ ] **Step 3: For every uncovered line in `/tmp/uncovered_lines.txt`:**

(a) Determine whether it's part of an inventory row or not.
(b) If not in any inventory row → the line is in a `MapMerge` but no inventory row claims it. Either it's dead code (annotate with `Lines = ...; status=unreachable; reason=...` and add a comment in the inventory), or it's a missing branch — go back to Phase 2 for that class, add the row, then back to Phase 3 to author the fixture.
(c) If in an inventory row but uncovered → the fixture for that row doesn't actually exercise the claimed lines. Re-author the fixture.

- [ ] **Step 4: Re-run Steps 1–3 after every fix until `/tmp/uncovered_lines.txt` is empty (or only contains lines explicitly annotated as unreachable in the inventory).**

- [ ] **Step 5: Commit any inventory annotations or new fixtures added during this task.**

```bash
git add docs/superpowers/specs/2026-05-05-simplify-branch-coverage-inventory.md \
        ast_tester/simplify_fixtures/ \
        ast_tester/simplify_tests.txt
git commit -m "gcov verification: close coverage gaps and annotate dead branches"
```

### Task 15: Per-fixture overlap check

**Goal.** Flag any focused/cascade inventory row whose target lines are hit by more than one fixture (excluding scenario fixtures).

- [ ] **Step 1: For each `simplify_*` test, run it in isolation with coverage cleared.**

```bash
cd build-coverage
mkdir -p coverage_per_fixture
GCNO_DIR=$(find . -name "cmpmap.c.gcno" -exec dirname {} \; | head -1)
GCNO_REL="${GCNO_DIR#./}"

for fixture in $(ctest -N -R '^simplify_[a-z_0-9]+$' \
                 | grep -oE 'simplify_[a-z_0-9]+' | sort -u); do
    case "$fixture" in
        simplify_brad|simplify_lsst1|simplify_rigby) continue ;;
        *_astequal) continue ;;  # the astequal sibling exercises the same code
    esac
    find . -name "*.gcda" -delete
    ctest -R "^${fixture}$" --output-on-failure > /dev/null

    out_dir="coverage_per_fixture/${fixture}"
    mkdir -p "$out_dir"
    for src in cmpmap unitmap zoommap shiftmap winmap matrixmap permmap \
               lutmap polymap mathmap ratemap intramap splinemap normmap \
               unitnormmap selectormap switchmap tranmap timemap slamap \
               specmap wcsmap sphmap pcdmap dssmap grismmap xphmap \
               mapping box interval nullregion pointlist prism; do
        gcov -o "$GCNO_REL" "../src/${src}.c" > /dev/null 2>&1
        mv "${src}.c.gcov" "$out_dir/" 2>/dev/null || true
    done
done
cd ..
```

This produces `build-coverage/coverage_per_fixture/<fixture>/<src>.c.gcov` — one set of gcov reports per fixture. Subsequent steps cross-reference each inventory row's target line range against these per-fixture reports to detect overlap.

- [ ] **Step 2: For each inventory row that targets a focused or cascade fixture, count how many fixtures' coverage reports show hits in the row's target line range.**

Use this Python script `/tmp/check_overlap.py` (do not commit):

```python
#!/usr/bin/env python3
"""For each inventory row, count how many per-fixture gcov reports hit it."""
import re, pathlib, collections

INV = pathlib.Path("docs/superpowers/specs/2026-05-05-simplify-branch-coverage-inventory.md")
PERFIX = pathlib.Path("build-coverage/coverage_per_fixture")

# Parse rows: id | fixture | type | polarity | lines | desc | trigger
row_re = re.compile(
    r"^\|\s*([a-z0-9_-]+)\s*\|\s*([a-z0-9_]+\.map|via [a-z0-9_-]+|\s*)\s*\|"
    r"\s*(focused|cascade|scenario)\s*\|\s*\w+\s*\|"
    r"\s*([a-z]+)\.c:(\d+)-(\d+)\s*\|", re.M)

# Build a map: (src, start, end) -> list of fixtures whose gcov shows hits
hits = collections.defaultdict(list)
for m in row_re.finditer(INV.read_text()):
    inv_id, fname, ftype, src, start, end = m.groups()
    if ftype == "scenario":
        continue
    start, end = int(start), int(end)
    for fix_dir in PERFIX.iterdir():
        gcov_file = fix_dir / f"{src}.c.gcov"
        if not gcov_file.exists():
            continue
        for line in gcov_file.read_text().splitlines():
            parts = line.split(":", 2)
            if len(parts) < 2:
                continue
            count, lineno_str = parts[0].strip(), parts[1].strip()
            try:
                lineno = int(lineno_str)
            except ValueError:
                continue
            if start <= lineno <= end and count not in ("#####", "-"):
                hits[(inv_id, src, start, end)].append(fix_dir.name)
                break

for key, fix_list in sorted(hits.items()):
    if len(fix_list) > 1:
        inv_id, src, start, end = key
        print(f"OVERLAP: {inv_id} ({src}.c:{start}-{end}) hit by: {', '.join(fix_list)}")
```

Run:
```bash
python3 /tmp/check_overlap.py > /tmp/overlap_warnings.txt
cat /tmp/overlap_warnings.txt
```

Expected: most inventory rows have one fixture; cascades may show two or three (the cascade itself plus the focused fixtures it subsumes).

- [ ] **Step 3: Review `/tmp/overlap_warnings.txt` and decide for each entry:**

(a) Legitimate (cascade naturally subsumes a focused branch) → no action.
(b) Concerning (two focused fixtures unexpectedly overlap) → revise one to be more targeted.

- [ ] **Step 4: Commit any fixture revisions made as a result of overlap review.**

```bash
git add ast_tester/simplify_fixtures/
git commit -m "Refine fixtures to reduce per-fixture coverage overlap"
```

---

## Phase 5 — Permanent reference document

**Goal.** Produce `ast_tester/simplify_pathways.md` from the inventory and delete `ast_tester/simplify_coverage.md`.

### Task 16: Generate the permanent reference document

**Files:**
- Create: `ast_tester/simplify_pathways.md`

- [ ] **Step 1: Write the top-matter section.**

Open `ast_tester/simplify_pathways.md` and write the introduction. About one page covering:
- What `astSimplify` is (one paragraph).
- How the pipeline of `MapMerge` calls works — `astSimplify` repeatedly invokes each component's `MapMerge`, allowing it to merge with neighbors, until no further changes occur (one paragraph).
- Links: SUN/210, SUN/211, `ast_tester/simplify_fixtures/`, the design doc `docs/superpowers/specs/2026-05-05-simplify-branch-coverage-design.md`, the working inventory `docs/superpowers/specs/2026-05-05-simplify-branch-coverage-inventory.md` (one short list).
- The maintenance contract (about a paragraph): when adding/modifying a `MapMerge` branch, add a row to the relevant section here and a fixture to `ast_tester/simplify_fixtures/`.

- [ ] **Step 2: For each class section in the inventory, generate the corresponding section in the reference doc.**

For each `## <class>` section in the inventory:
1. Write 1–2 paragraphs describing the class's simplification strategy. Use the inventory's per-class summary paragraph as the seed; expand to be reader-facing (no jargon like "MapMerge branch"; use phrasings like "this rule fires when…").
2. Write a Markdown table with columns: `ID | Fixture | Type | Polarity | Description | Trigger`. Drop the `Lines` column. Group rows by `Polarity`: positives first, negatives second.
3. For cascade rules, add a short prose paragraph immediately after the table explaining why a multi-component composition is needed.

- [ ] **Step 3: Add a "Cross-class index" section at the end.**

Alphabetical list of every inventory ID with the one-line description, formatted as:

```markdown
## Cross-class index

- `box-01` — Box-Mapping merge with adjacent unit shift collapses to inset Box.
- `cmpmap-01` — Nested parallel CmpMaps flatten to a single parallel list.
...
```

- [ ] **Step 4: Verify every inventory row appears in the reference doc.**

```bash
grep -oE '^\| [a-z]+map-[0-9]+' \
     docs/superpowers/specs/2026-05-05-simplify-branch-coverage-inventory.md \
  | sort -u > /tmp/inventory_ids.txt

grep -oE '`[a-z]+map-[0-9]+`' ast_tester/simplify_pathways.md \
  | tr -d '`' \
  | sort -u > /tmp/refdoc_ids.txt

diff /tmp/inventory_ids.txt /tmp/refdoc_ids.txt
```

Expected: empty diff. If anything is missing from the reference doc, add it.

- [ ] **Step 5: Commit.**

```bash
git add ast_tester/simplify_pathways.md
git commit -m "Add simplify_pathways.md as permanent reference for astSimplify rules"
```

### Task 17: Delete `simplify_coverage.md` and update cross-references

**Files:**
- Delete: `ast_tester/simplify_coverage.md`
- Modify: any file that references `simplify_coverage.md`

- [ ] **Step 1: Find references to `simplify_coverage.md`.**

```bash
grep -rn "simplify_coverage" --exclude-dir=build --exclude-dir=build-coverage \
                              --exclude-dir=.git .
```

- [ ] **Step 2: Update each reference to point at `simplify_pathways.md` (for reference material) or at the design/inventory under `docs/superpowers/specs/` (for historical context).**

- [ ] **Step 3: Delete `simplify_coverage.md`.**

```bash
git rm ast_tester/simplify_coverage.md
```

- [ ] **Step 4: Verify no stale references remain.**

```bash
grep -rn "simplify_coverage" --exclude-dir=build --exclude-dir=build-coverage \
                              --exclude-dir=.git .
```

Expected: no results.

- [ ] **Step 5: Commit.**

```bash
git add -A
git commit -m "Remove simplify_coverage.md, replaced by simplify_pathways.md"
```

### Task 18: Final verification and PR readiness

- [ ] **Step 1: Run the full simplify suite one more time, in a clean Debug build.**

```bash
rm -rf build-final
cmake -B build-final \
  -DCMAKE_BUILD_TYPE=Debug \
  -DAST_ENABLE_WARNINGS=ON \
  -DAST_ENABLE_SANITIZERS=ON
cmake --build build-final 2>&1 | tail -10
ctest --test-dir build-final -R '^simplify_' --output-on-failure 2>&1 | tail -20
```

Expected: all simplify tests pass under sanitizers and warnings.

- [ ] **Step 2: Run the full ctest suite to confirm nothing else broke.**

```bash
ctest --test-dir build-final --output-on-failure 2>&1 | tail -30
```

Expected: same pass/fail status as `master` (no new regressions). If existing tests fail and they were passing on `master` before this branch's changes, investigate and fix.

- [ ] **Step 3: Check no warnings were introduced by this branch's source changes.**

This branch should not modify any `src/*.c`, but verify:

```bash
git diff --stat master..HEAD -- src/ | tail -5
```

Expected: no source files in the diff.

- [ ] **Step 4: Confirm the branch is ready for PR.**

```bash
git log --oneline master..HEAD | wc -l
git status
```

Expected: a clean working tree with all phase commits in `git log`.

- [ ] **Step 5: Push and open the PR.**

(User to perform — do not push without explicit user authorization.)

---

## Self-review (executor checklist before declaring complete)

After Task 18, before declaring the campaign complete:

1. Every row in the inventory has a non-empty `Fixture` column (or an explicit `via <id>` cascade reference).
2. Every fixture file in `ast_tester/simplify_fixtures/` corresponds to an inventory row.
3. Every `simplify_tests.txt` row points to a file under `simplify_fixtures/` (except the three scenario rows).
4. `ctest -R '^simplify_'` is green under both regular Debug and Debug+sanitizers builds.
5. gcov shows every inventory `Lines` range is hit, with explicit annotations for any unreachable ranges.
6. `ast_tester/simplify_pathways.md` includes every inventory ID.
7. `ast_tester/simplify_coverage.md` is gone, and no file references it.
8. `git diff --stat master..HEAD -- src/` shows no source modifications.
