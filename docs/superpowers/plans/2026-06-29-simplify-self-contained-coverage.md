# Self-Contained `astSimplify` Coverage Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the `ast_tester/simplify_fixtures/` set self-contained — running only the `simplify_*` tests covers the simplify call graph without relying on FitsChan/WCS/Region tests — by importing 13 contributed fixtures and then driving a regenerable coverage gap ledger to zero.

**Architecture:** A coverage build is run twice (simplify-only, then full suite); `gcovr` captures normalized branch-level JSON for each; a stdlib-only Python script diffs them into a checked-in Markdown gap ledger that classifies each uncovered branch as `differential` (covered by other tests only — top priority) or `absolute-only` (covered by nothing). Fixtures are authored to close gaps; the ledger regenerates idempotently, preserving human annotations, so the open-ended effort can pause and resume on either macOS or Linux.

**Tech Stack:** C99 AST native serialization fixtures, the existing `simplify.c` harness, CMake/ctest, `gcovr` 8.6 (from `~/pyenv`), `gcov`/`llvm-cov gcov`, Python 3 stdlib (`json`, `unittest`, `bisect`).

**Companion design doc:** `docs/superpowers/specs/2026-06-29-simplify-self-contained-coverage-design.md`. Read it first.

## Global Constraints

- **Regression-only.** No `src/*.c` is modified. `git diff --stat master..HEAD -- src/` must stay empty.
- **Branch:** all work lands on `u/timj/simp-coverage`. Never push (user pushes manually).
- **Fixtures** live in `ast_tester/simplify_fixtures/`, are registered as rows in `ast_tester/simplify_tests.txt`, and carry **no in-band header** (none of the existing 248 do — traceability is by filename through the inventory).
- **Negative fixtures** (`neg_*`, output == input) use the `.map` in **both** the input and reference columns of `simplify_tests.txt` with the skip flag `yes` (astequal-only). They have no `.simp`.
- **`skip_string_compare=yes`** (4th column) is set for any fixture whose serialized output varies at ulp level across platforms; the existing `rigby` row is the precedent.
- **Coverage build is `-O0`** so branch structure matches source; it is a separate build dir from the sanitizer build.
- **`gcovr` is invoked from the user's `~/pyenv`**; the diff/ledger Python script imports nothing third-party.
- **Cascade fixtures** compose at most three mappings.
- **Inventory** = `docs/superpowers/specs/2026-05-05-simplify-branch-coverage-inventory.md` (schema: `ID | Fixture | Type | Polarity | Lines | Description | Trigger`). **Reference doc** = `ast_tester/simplify_pathways.md`. Both gain a row per new/imported fixture.
- **Target functions** (the "simplify pathway") are the set `{MapMerge, MapList, Simplify, RemoveRegions, CombineMaps, Decompose, WinMat, MatWin, MatWin2}` defined within the in-scope source files (the 32 `MapMerge` files plus `mapping.c`).

---

## File structure

**Files created:**
- `ast_tester/coverage/simplify_coverage_gaps.py` — stdlib-only diff + ledger generator. One responsibility: turn two gcovr JSON files (+ the prior ledger) into the new ledger.
- `ast_tester/coverage/test_simplify_coverage_gaps.py` — `unittest` tests for the generator, using synthetic in-memory gcovr JSON.
- `ast_tester/coverage/run_simplify_coverage.sh` — orchestration wrapper: build, run two ctest passes, capture two gcovr JSONs, invoke the generator.
- `ast_tester/coverage/README.md` — orientation: how to regenerate the ledger and what the columns mean.
- `ast_tester/coverage/simplify_coverage_gaps.md` — the committed, regenerable gap ledger.
- 13 fixtures copied into `ast_tester/simplify_fixtures/` (12 `.map`+`.simp`, 1 `.map`-only).
- New `.map`/`.simp` fixtures authored during Tasks 4–5 (count unknown until baseline).

**Files modified:**
- `ast_tester/simplify_tests.txt` — append rows for imported and newly-authored fixtures.
- `docs/superpowers/specs/2026-05-05-simplify-branch-coverage-inventory.md` — append rows.
- `ast_tester/simplify_pathways.md` — append rows.
- `.gitignore` — ignore transient coverage artifacts.

**Files unchanged:** `ast_tester/simplify.c`, `ast_tester/CMakeLists.txt` (reads `simplify_tests.txt` automatically), all `src/*.c`.

**Fixture → inventory-section assignment for the 13 imports** (used in Task 1):

| Fixture | Section | Polarity |
| --- | --- | --- |
| `dssmap_inv_winmap_absorb` | dssmap.c | positive |
| `dssmap_winmap_absorb` | dssmap.c | positive |
| `dssmap_zoom_no_merge` | dssmap.c | negative (name not `neg_`; classify by reading output) |
| `neg_box_asymmetric_2d` | box.c | negative |
| `reconstruct_right_assoc` | cmpmap.c | positive |
| `series_parallel_absorb` | cmpmap.c | positive |
| `sla_run_identity` | slamap.c | positive |
| `sla_run_partial` | slamap.c | positive |
| `spec_cel_invert` | specmap.c | positive |
| `spec_run_inverted_mid` | specmap.c | positive |
| `spec_run4_identity` | specmap.c | positive |
| `time_run_identity` | timemap.c | positive |
| `wcsconv_gappt_iwc_residue` | wcsmap.c | scenario/cascade |

Confirm each polarity/type by reading the `.map` and the harness output during Step 5 below; the table is the starting assignment, not gospel.

---

## Task 1: Import the 13 contributed fixtures (P0)

**Files:**
- Create: `ast_tester/simplify_fixtures/<name>.map` (×13) and `.simp` (×12)
- Modify: `ast_tester/simplify_tests.txt`
- Modify: `docs/superpowers/specs/2026-05-05-simplify-branch-coverage-inventory.md`
- Modify: `ast_tester/simplify_pathways.md`

**Interfaces:**
- Produces: 13 registered `simplify_*` ctest cases that pass; the post-import fixture set that Task 3 baselines against.

- [ ] **Step 1: Copy the 13 `.map` files and 12 `.simp` files into the fixtures dir.**

```bash
cd /Users/timj/work/starlink-ast
SRC=~/work/lsst/starlink-astrs/crates/ast-core/tests/fixtures/native/reference/simplify_fixtures
DST=ast_tester/simplify_fixtures
for n in dssmap_inv_winmap_absorb dssmap_winmap_absorb dssmap_zoom_no_merge \
         neg_box_asymmetric_2d reconstruct_right_assoc series_parallel_absorb \
         sla_run_identity sla_run_partial spec_cel_invert spec_run_inverted_mid \
         spec_run4_identity time_run_identity wcsconv_gappt_iwc_residue; do
    cp "$SRC/$n.map" "$DST/$n.map"
    [ -f "$SRC/$n.simp" ] && cp "$SRC/$n.simp" "$DST/$n.simp"
done
ls "$DST" | grep -E 'dssmap_inv_winmap_absorb|neg_box_asymmetric_2d|wcsconv_gappt' 
```

Expected: the copied files are listed. `neg_box_asymmetric_2d` has only a `.map`.

- [ ] **Step 2: Append rows to `simplify_tests.txt`.**

Add a section divider comment and 13 rows. Positives use `.map`/`.simp`; the negative uses the `.map` in both columns with `yes`:

```
# Imported contributed fixtures (2026-06-29).
dssmap_inv_winmap_absorb  | simplify_fixtures/dssmap_inv_winmap_absorb.map  | simplify_fixtures/dssmap_inv_winmap_absorb.simp  |
dssmap_winmap_absorb      | simplify_fixtures/dssmap_winmap_absorb.map      | simplify_fixtures/dssmap_winmap_absorb.simp      |
dssmap_zoom_no_merge      | simplify_fixtures/dssmap_zoom_no_merge.map      | simplify_fixtures/dssmap_zoom_no_merge.simp      |
neg_box_asymmetric_2d     | simplify_fixtures/neg_box_asymmetric_2d.map     | simplify_fixtures/neg_box_asymmetric_2d.map      | yes
reconstruct_right_assoc   | simplify_fixtures/reconstruct_right_assoc.map   | simplify_fixtures/reconstruct_right_assoc.simp   |
series_parallel_absorb    | simplify_fixtures/series_parallel_absorb.map    | simplify_fixtures/series_parallel_absorb.simp    |
sla_run_identity          | simplify_fixtures/sla_run_identity.map          | simplify_fixtures/sla_run_identity.simp          |
sla_run_partial           | simplify_fixtures/sla_run_partial.map           | simplify_fixtures/sla_run_partial.simp           |
spec_cel_invert           | simplify_fixtures/spec_cel_invert.map           | simplify_fixtures/spec_cel_invert.simp           |
spec_run_inverted_mid     | simplify_fixtures/spec_run_inverted_mid.map     | simplify_fixtures/spec_run_inverted_mid.simp     |
spec_run4_identity        | simplify_fixtures/spec_run4_identity.map        | simplify_fixtures/spec_run4_identity.simp        |
time_run_identity         | simplify_fixtures/time_run_identity.map         | simplify_fixtures/time_run_identity.simp         |
wcsconv_gappt_iwc_residue | simplify_fixtures/wcsconv_gappt_iwc_residue.map | simplify_fixtures/wcsconv_gappt_iwc_residue.simp |
```

- [ ] **Step 3: Configure and build (plain Debug is enough for this task).**

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Debug 2>&1 | tail -3
cmake --build build 2>&1 | tail -5
```

Expected: build succeeds; CMake re-reads `simplify_tests.txt` and registers `simplify_<name>` + `simplify_<name>_astequal` for each new row.

- [ ] **Step 4: Run the 13 new tests.**

```bash
ctest --test-dir build -R '^simplify_(dssmap_inv_winmap_absorb|dssmap_winmap_absorb|dssmap_zoom_no_merge|neg_box_asymmetric_2d|reconstruct_right_assoc|series_parallel_absorb|sla_run_identity|sla_run_partial|spec_cel_invert|spec_run_inverted_mid|spec_run4_identity|time_run_identity|wcsconv_gappt_iwc_residue)(_astequal)?$' --output-on-failure 2>&1 | tail -40
```

Expected: ideally all pass.

- [ ] **Step 5: Resolve any failures by classification.**

For each failing `simplify_<name>` (string-diff) test:
- If the matching `simplify_<name>_astequal` test **passes**, the output differs only at ulp/formatting level. Set the 4th column of that row to `yes` (astequal-only). Re-run.
- If the `_astequal` test **also fails**, the harness output is structurally different from the contributed `.simp`. Do **not** paper over it: read both with `diff -u`, determine whether (a) the contributed `.simp` was generated against different behavior (regenerate `.simp` from our harness: `./build/ast_tester/simplify ast_tester/simplify_fixtures/<name>.map ast_tester/simplify_fixtures/<name>.simp` and note the discrepancy in the commit message), or (b) it reveals a real simplification difference worth a follow-up issue. Stop and report to the user before continuing if it looks like a real difference.

Re-run Step 4 until green.

- [ ] **Step 6: Add inventory and `simplify_pathways.md` rows for the 13.**

For each fixture, open the inventory section named in the assignment table, find the highest existing `<class>-NN` ID in that section, and append a row continuing the numbering. Read the `.map` and harness output to fill `Description`/`Trigger`, and identify the `Lines` range by grepping the class's `MapMerge`/helper for the branch.

Worked example — `sla_run_identity` (slamap.c section; existing rows run to `slamap-10`+, so this is the next free ID, e.g. `slamap-21`):

```
| slamap-21 | sla_run_identity.map | focused | positive | slamap.c:NNNN-MMMM | SlaMap run of conversions that fully cancels to a UnitMap | SlaMap whose forward+inverse step sequence is identity |
```

Add the parallel row to `simplify_pathways.md` (same columns minus `Lines`, grouped by polarity within the class section). Repeat for all 13. For `wcsconv_gappt_iwc_residue`, if it composes >1 component and exercises orchestration, mark `Type=cascade` (or `scenario` if >3 components) and place it under wcsmap.c.

- [ ] **Step 7: Verify all simplify tests still pass and commit.**

```bash
ctest --test-dir build -R '^simplify_' --output-on-failure 2>&1 | tail -10
git add ast_tester/simplify_fixtures/ ast_tester/simplify_tests.txt \
        docs/superpowers/specs/2026-05-05-simplify-branch-coverage-inventory.md \
        ast_tester/simplify_pathways.md
git commit -m "Import 13 contributed simplify fixtures"
```

Expected: full `simplify_*` suite green; commit created.

---

## Task 2: Build the coverage gap generator (stdlib Python, TDD)

**Files:**
- Create: `ast_tester/coverage/simplify_coverage_gaps.py`
- Test: `ast_tester/coverage/test_simplify_coverage_gaps.py`

**Interfaces:**
- Produces (consumed by Task 3's wrapper):
  - `load_coverage(path: str) -> CovData` where `CovData` is a dict
    `{"branches": {(file, line, idx): count}, "functions": {file: [(lineno, name), ...sorted]}}`.
  - `function_for(file: str, line: int, functions: dict) -> str | None` — name of the enclosing function (greatest `lineno <= line`), or `None`.
  - `classify_gaps(simp: CovData, full: CovData, allowlist: set[str], inscope: set[str]) -> list[Gap]` where `Gap = (file, line, idx, function, cls)` and `cls in {"differential", "absolute-only"}`.
  - `parse_prior(md_text: str) -> dict[(file, function, line, idx), status]`.
  - `render_ledger(gaps, prior, full: CovData, simp: CovData) -> str`.
  - CLI: `python3 simplify_coverage_gaps.py --simplify simplify.json --full full.json --ledger simplify_coverage_gaps.md`.

- [ ] **Step 1: Write the failing test.**

Create `ast_tester/coverage/test_simplify_coverage_gaps.py`:

```python
import unittest
import simplify_coverage_gaps as g

# Minimal gcovr-8.x JSON: one file, one function, three branch lines.
def cov(branch_counts):
    # branch_counts: {line: [counts...]}
    lines = []
    for line, counts in branch_counts.items():
        lines.append({"line_number": line, "count": sum(counts) or 0,
                      "branches": [{"count": c} for c in counts]})
    return {"files": [{"file": "src/winmap.c",
                       "lines": lines,
                       "functions": [{"name": "MapMerge", "lineno": 900},
                                     {"name": "WinMat", "lineno": 2793}]}]}

ALLOW = {"MapMerge", "WinMat"}
INSCOPE = {"winmap.c"}

class ClassifyTest(unittest.TestCase):
    def test_differential_and_absolute(self):
        simp = g.load_coverage_obj(cov({905: [5, 0], 910: [0, 0]}))
        full = g.load_coverage_obj(cov({905: [5, 3], 910: [0, 0]}))
        gaps = g.classify_gaps(simp, full, ALLOW, INSCOPE)
        keys = {(gp[1], gp[2], gp[4]) for gp in gaps}  # (line, idx, cls)
        # line 905 branch idx 1: full-hit, simp-miss -> differential
        self.assertIn((905, 1, "differential"), keys)
        # line 910 both branches hit by nobody -> absolute-only
        self.assertIn((910, 0, "absolute-only"), keys)
        self.assertIn((910, 1, "absolute-only"), keys)
        # line 905 branch idx 0: simp-hit -> not a gap
        self.assertNotIn((905, 0, "differential"), keys)
        self.assertNotIn((905, 0, "absolute-only"), keys)

    def test_function_for(self):
        functions = {"src/winmap.c": [(900, "MapMerge"), (2793, "WinMat")]}
        self.assertEqual(g.function_for("src/winmap.c", 905, functions), "MapMerge")
        self.assertEqual(g.function_for("src/winmap.c", 2800, functions), "WinMat")
        self.assertIsNone(g.function_for("src/winmap.c", 10, functions))

    def test_prior_annotation_preserved(self):
        prior = g.parse_prior(
            "| `src/winmap.c:905` | 1 | MapMerge | differential | fixture=win_x |\n")
        self.assertEqual(prior[("src/winmap.c", "MapMerge", 905, 1)], "fixture=win_x")

if __name__ == "__main__":
    unittest.main()
```

- [ ] **Step 2: Run the test to confirm it fails.**

```bash
cd ast_tester/coverage && python3 -m unittest test_simplify_coverage_gaps -v
```

Expected: FAIL / ERROR — `module 'simplify_coverage_gaps' has no attribute ...`.

- [ ] **Step 3: Implement `simplify_coverage_gaps.py`.**

```python
#!/usr/bin/env python3
"""Diff two gcovr JSON coverage files into a simplify-pathway gap ledger.

Stdlib only. Reads gcovr 8.x JSON (the `--json` format): each file has
`lines` (with `line_number`, `count`, `branches:[{count}]`) and `functions`
(with `name`, `lineno`). Classifies each uncovered branch in a target
function as `differential` (hit by the full suite but not by simplify-only)
or `absolute-only` (hit by neither).
"""
import argparse
import bisect
import json
import os
import re

ALLOWLIST = {"MapMerge", "MapList", "Simplify", "RemoveRegions",
             "CombineMaps", "Decompose", "WinMat", "MatWin", "MatWin2"}

INSCOPE = {
    "mapping.c", "box.c", "cmpmap.c", "dssmap.c", "grismmap.c", "interval.c",
    "intramap.c", "lutmap.c", "mathmap.c", "matrixmap.c", "normmap.c",
    "nullregion.c", "pcdmap.c", "permmap.c", "pointlist.c", "polymap.c",
    "prism.c", "ratemap.c", "selectormap.c", "shiftmap.c", "slamap.c",
    "specmap.c", "sphmap.c", "splinemap.c", "switchmap.c", "timemap.c",
    "tranmap.c", "unitmap.c", "unitnormmap.c", "wcsmap.c", "winmap.c",
    "xphmap.c", "zoommap.c",
}


def load_coverage_obj(obj):
    """Turn a parsed gcovr JSON object into a CovData dict."""
    branches = {}
    functions = {}
    for f in obj.get("files", []):
        path = f["file"]
        functions[path] = sorted(
            (fn["lineno"], fn["name"]) for fn in f.get("functions", []))
        for ln in f.get("lines", []):
            line = ln["line_number"]
            for idx, br in enumerate(ln.get("branches", [])):
                branches[(path, line, idx)] = br.get("count", 0)
    return {"branches": branches, "functions": functions}


def load_coverage(path):
    with open(path) as fh:
        return load_coverage_obj(json.load(fh))


def function_for(file, line, functions):
    table = functions.get(file, [])
    linenos = [t[0] for t in table]
    pos = bisect.bisect_right(linenos, line) - 1
    if pos < 0:
        return None
    return table[pos][1]


def classify_gaps(simp, full, allowlist=ALLOWLIST, inscope=INSCOPE):
    gaps = []
    for key, full_count in full["branches"].items():
        file, line, idx = key
        if os.path.basename(file) not in inscope:
            continue
        fn = function_for(file, line, full["functions"])
        if fn not in allowlist:
            continue
        simp_count = simp["branches"].get(key, 0)
        if simp_count > 0:
            continue  # already self-covered
        cls = "differential" if full_count > 0 else "absolute-only"
        gaps.append((file, line, idx, fn, cls))
    gaps.sort()
    return gaps


_ROW = re.compile(
    r"\|\s*`([^:]+):(\d+)`\s*\|\s*(\d+)\s*\|\s*(\w+)\s*\|\s*[\w-]+\s*\|\s*([^|]*?)\s*\|")


def parse_prior(md_text):
    prior = {}
    for m in _ROW.finditer(md_text):
        file, line, idx, fn, status = m.groups()
        status = status.strip()
        if status and status != "open":
            prior[(file, fn, int(line), int(idx))] = status
    return prior


def render_ledger(gaps, prior, full, simp):
    lines = ["# `astSimplify` Self-Contained Coverage Gap Ledger", "",
             "Regenerated by `run_simplify_coverage.sh`. Do not hand-edit",
             "rows except the Status column. See `README.md`.", ""]
    open_rows, closed_rows = [], []
    live_keys = set()
    for file, line, idx, fn, cls in gaps:
        key = (file, fn, line, idx)
        live_keys.add(key)
        status = prior.get(key, "open")
        open_rows.append(
            f"| `{file}:{line}` | {idx} | {fn} | {cls} | {status} |")
    # Closed: anything previously annotated that is no longer a gap.
    for (file, fn, line, idx), status in sorted(prior.items()):
        if (file, fn, line, idx) not in live_keys:
            closed_rows.append(
                f"| `{file}:{line}` | {idx} | {fn} | {status} |")
    diff_n = sum(1 for r in open_rows if " differential " in r)
    abs_n = sum(1 for r in open_rows if " absolute-only " in r)
    lines += [f"**Open differential gaps:** {diff_n}  ",
              f"**Open absolute-only gaps:** {abs_n}", "",
              "## Open gaps", "",
              "| Location | Branch | Function | Class | Status |",
              "| --- | --- | --- | --- | --- |"]
    lines += open_rows or ["| _none_ |  |  |  |  |"]
    lines += ["", "## Closed", "",
              "| Location | Branch | Function | Status |",
              "| --- | --- | --- | --- |"]
    lines += closed_rows or ["| _none_ |  |  |  |"]
    return "\n".join(lines) + "\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--simplify", required=True)
    ap.add_argument("--full", required=True)
    ap.add_argument("--ledger", required=True)
    args = ap.parse_args()
    simp = load_coverage(args.simplify)
    full = load_coverage(args.full)
    gaps = classify_gaps(simp, full)
    prior = {}
    if os.path.exists(args.ledger):
        with open(args.ledger) as fh:
            prior = parse_prior(fh.read())
    with open(args.ledger, "w") as fh:
        fh.write(render_ledger(gaps, prior, full, simp))
    diff_n = sum(1 for gp in gaps if gp[4] == "differential")
    print(f"Wrote {args.ledger}: {len(gaps)} open gaps "
          f"({diff_n} differential).")


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run the test to confirm it passes.**

```bash
cd ast_tester/coverage && python3 -m unittest test_simplify_coverage_gaps -v
```

Expected: PASS (3 tests).

- [ ] **Step 5: Commit.**

```bash
git add ast_tester/coverage/simplify_coverage_gaps.py \
        ast_tester/coverage/test_simplify_coverage_gaps.py
git commit -m "Add stdlib simplify coverage gap generator with tests"
```

---

## Task 3: Coverage wrapper and baseline ledger (P1)

**Files:**
- Create: `ast_tester/coverage/run_simplify_coverage.sh`
- Create: `ast_tester/coverage/README.md`
- Create: `ast_tester/coverage/simplify_coverage_gaps.md` (generated)
- Modify: `.gitignore`

**Interfaces:**
- Consumes: `simplify_coverage_gaps.py` CLI from Task 2.
- Produces: the committed baseline ledger that Tasks 4–5 work from.

- [ ] **Step 1: Write the wrapper script.**

Create `ast_tester/coverage/run_simplify_coverage.sh`:

```bash
#!/usr/bin/env bash
# Regenerate the simplify coverage gap ledger.
# Requires gcovr on PATH (e.g. `source ~/pyenv/bin/activate`).
set -euo pipefail
cd "$(dirname "$0")/../.."          # repo root
ROOT="$PWD"
BUILD="build-cov"
COVDIR="ast_tester/coverage"

if ! command -v gcovr >/dev/null 2>&1; then
    echo "ERROR: gcovr not found. Activate your gcovr env, e.g.:" >&2
    echo "  source ~/pyenv/bin/activate" >&2
    exit 1
fi

# Pick a gcov compatible with the compiler (clang on macOS needs llvm-cov).
if [ -z "${GCOV:-}" ]; then
    if [ "$(uname)" = "Darwin" ]; then GCOV="xcrun llvm-cov gcov"; else GCOV="gcov"; fi
fi

# In-scope file filter (32 MapMerge files + mapping.c).
FILTER='src/(mapping|box|cmpmap|dssmap|grismmap|interval|intramap|lutmap|mathmap|matrixmap|normmap|nullregion|pcdmap|permmap|pointlist|polymap|prism|ratemap|selectormap|shiftmap|slamap|specmap|sphmap|splinemap|switchmap|timemap|tranmap|unitmap|unitnormmap|wcsmap|winmap|xphmap|zoommap)\.c$'

cmake -B "$BUILD" -DCMAKE_BUILD_TYPE=Debug \
      -DCMAKE_C_FLAGS="--coverage -O0" \
      -DCMAKE_EXE_LINKER_FLAGS="--coverage" \
      -DCMAKE_SHARED_LINKER_FLAGS="--coverage" >/dev/null
cmake --build "$BUILD" >/dev/null

run_gcovr() {  # $1 = output json
    gcovr --root "$ROOT" --search-path "$BUILD" \
          --gcov-executable "$GCOV" \
          --filter "$FILTER" --json --output "$1" --delete >/dev/null
}

echo ">> simplify-only pass"
ctest --test-dir "$BUILD" -R '^simplify_' --output-on-failure >/dev/null
run_gcovr "$COVDIR/simplify.json"

echo ">> full-suite pass"
ctest --test-dir "$BUILD" >/dev/null || true   # coverage even if some non-simplify test fails
run_gcovr "$COVDIR/full.json"

python3 "$COVDIR/simplify_coverage_gaps.py" \
    --simplify "$COVDIR/simplify.json" \
    --full "$COVDIR/full.json" \
    --ledger "$COVDIR/simplify_coverage_gaps.md"
```

```bash
chmod +x ast_tester/coverage/run_simplify_coverage.sh
```

- [ ] **Step 2: Confirm the gcovr JSON field names match the parser.**

```bash
source ~/pyenv/bin/activate 2>/dev/null || true
cmake -B build-cov -DCMAKE_BUILD_TYPE=Debug -DCMAKE_C_FLAGS="--coverage -O0" \
      -DCMAKE_EXE_LINKER_FLAGS="--coverage" -DCMAKE_SHARED_LINKER_FLAGS="--coverage" >/dev/null
cmake --build build-cov >/dev/null
ctest --test-dir build-cov -R '^simplify_zoom_series_merge$' >/dev/null
GCOV=$([ "$(uname)" = Darwin ] && echo "xcrun llvm-cov gcov" || echo gcov)
gcovr --root "$PWD" --search-path build-cov --gcov-executable "$GCOV" \
      --filter 'src/zoommap\.c$' --json --output /tmp/probe.json --delete >/dev/null
python3 -c "import json;d=json.load(open('/tmp/probe.json'));f=d['files'][0];print('FILE KEYS',list(f));print('LINE KEYS',list(f['lines'][0]));print('FUNC KEYS',list(f['functions'][0]))"
```

Expected: `LINE KEYS` includes `line_number`, `count`, `branches`; `FUNC KEYS` includes `name`, `lineno`. If a key name differs in this gcovr build, adjust `load_coverage_obj`/`load_coverage` in `simplify_coverage_gaps.py` and re-run Task 2 Step 4 before proceeding.

- [ ] **Step 3: Add `.gitignore` entries for transient artifacts.**

Append to `.gitignore`:

```
/build-cov/
ast_tester/coverage/simplify.json
ast_tester/coverage/full.json
```

- [ ] **Step 4: Write `ast_tester/coverage/README.md`.**

```markdown
# Simplify coverage gap tracking

`simplify_coverage_gaps.md` records simplify-pathway branches that the
`simplify_*` tests do not cover on their own.

## Regenerate

    source ~/pyenv/bin/activate      # gcovr 8.x
    ast_tester/coverage/run_simplify_coverage.sh

This builds `build-cov` with `--coverage -O0`, runs the simplify tests then
the full suite (capturing branch coverage after each via gcovr), and rewrites
the ledger.

## Reading the ledger

- **differential** — branch covered by the full suite but not by the simplify
  fixtures alone. Highest priority: author a fixture so the simplify set is
  self-contained.
- **absolute-only** — branch covered by nothing. Investigate: author a fixture
  if reachable, else set Status to `unreachable:<reason>`.

The **Status** column is the only field you edit by hand:
`open` | `fixture=<name>` | `unreachable:<reason>` | `wontfix:<reason>`.
Regeneration preserves your Status notes by `(file, function, line, branch)`.
```

- [ ] **Step 5: Generate the baseline ledger.**

```bash
source ~/pyenv/bin/activate 2>/dev/null || true
ast_tester/coverage/run_simplify_coverage.sh
sed -n '1,20p' ast_tester/coverage/simplify_coverage_gaps.md
```

Expected: the ledger is written and prints its open differential/absolute counts. This is the P1 baseline.

- [ ] **Step 6: Commit.**

```bash
git add .gitignore ast_tester/coverage/run_simplify_coverage.sh \
        ast_tester/coverage/README.md ast_tester/coverage/simplify_coverage_gaps.md
git commit -m "Add simplify coverage wrapper and baseline gap ledger"
```

---

## Task 4: Differential pass — close `differential` gaps (P2, repeatable)

**Goal.** Drive the ledger's open `differential` count to zero. This is the primary milestone. The task repeats the per-gap procedure below for each `differential` row; the gaps come from the ledger, so this is a recipe rather than a fixed step list. Batch commits per source class.

**Files (per gap):**
- Create: `ast_tester/simplify_fixtures/<name>.map` (+ `.simp` for positives)
- Modify: `ast_tester/simplify_tests.txt`, the inventory, `simplify_pathways.md`, and the ledger Status column.

**Per-gap procedure:**

- [ ] **Step 1: Pick the next open `differential` gap from the ledger.** Note its `file`, `line`, `function`.

- [ ] **Step 2: Read the branch in `src/<file>` around `:line`** to understand the condition and the input pipeline shape that triggers it.

```bash
sed -n '<line-12>,<line+12>p' src/<file>
```

- [ ] **Step 3: Author the smallest input `.map`.** Hand-write for simple shapes. For complex inputs use a throwaway generator (`ast_tester/gen_simplify_fixtures.c` is the precedent) compiled against the build, then discard it — do not commit the generator. Save as `ast_tester/simplify_fixtures/<name>.map`.

  Worked illustrative example (substitute the real gap). Suppose the gap is a `winmap.c MapMerge` branch that merges a WinMap with a following ZoomMap; author `win_zoom_absorb.map` as a 2-component series `CmpMap(WinMap, ZoomMap)` with matching `Nin`.

- [ ] **Step 4: Produce and verify the candidate output.**

```bash
./build/ast_tester/simplify ast_tester/simplify_fixtures/<name>.map /tmp/<name>.simp
cat /tmp/<name>.simp
```

Confirm: for a **positive** fixture the structure is the simpler form the branch produces; for a **negative** fixture the structure is unchanged from the input. If it does not match the branch's intent, the input shape is wrong — revise the `.map` and repeat.

- [ ] **Step 5: Install the fixture.** Positive: `cp /tmp/<name>.simp ast_tester/simplify_fixtures/<name>.simp`. Negative: no `.simp`; the row uses the `.map` in both columns with `yes`.

- [ ] **Step 6: Register and document.** Append the `simplify_tests.txt` row, an inventory row (next free `<class>-NN`, with the real `Lines` range), and a `simplify_pathways.md` row. Set the ledger Status for this gap's `(file,line,branch)` to `fixture=<name>`.

- [ ] **Step 7: Build and run the new test.**

```bash
cmake --build build >/dev/null
ctest --test-dir build -R "^simplify_<name>(_astequal)?$" --output-on-failure
```

Expected: pass (apply the Task 1 Step 5 ulp/skip logic if the string test fails but astequal passes).

- [ ] **Step 8: Regenerate the ledger and confirm the gap closed.**

```bash
source ~/pyenv/bin/activate 2>/dev/null || true
ast_tester/coverage/run_simplify_coverage.sh
```

Expected: the gap moves from "Open gaps" to "Closed"; open differential count drops by at least one. If the branch is still open, the fixture does not actually exercise it — return to Step 3.

- [ ] **Step 9: Commit (batched per class).** After finishing a class's differential gaps:

```bash
git add ast_tester/simplify_fixtures/ ast_tester/simplify_tests.txt \
        docs/superpowers/specs/2026-05-05-simplify-branch-coverage-inventory.md \
        ast_tester/simplify_pathways.md ast_tester/coverage/simplify_coverage_gaps.md
git commit -m "Close differential simplify gaps in <class>"
```

- [ ] **Step 10: Repeat Steps 1–9 until the ledger shows 0 open differential gaps.** This reaches the primary milestone; the effort may be paused here.

---

## Task 5: Absolute pass — close or annotate `absolute-only` gaps (P3, open-ended)

**Goal.** For each open `absolute-only` gap, either author a fixture (if the branch is reachable) or annotate it as unreachable with justification. Open-ended; stop after any commit.

**Per-gap procedure:**

- [ ] **Step 1: Pick the next open `absolute-only` gap.** Note `file`, `line`, `function`.

- [ ] **Step 2: Read the branch and judge reachability.**

```bash
sed -n '<line-15>,<line+15>p' src/<file>
```

Decide: can any Mapping pipeline drive `astSimplify` into this branch? Consider guard preconditions and what neighbor maps would be required.

- [ ] **Step 3a (reachable): author a fixture** following Task 4 Steps 3–9, set ledger Status to `fixture=<name>`, and confirm it closes via regeneration.

- [ ] **Step 3b (unreachable): annotate the ledger.** Set the Status to `unreachable:<one-line reason>` (e.g. `unreachable:guard above already returns for Nin!=Nout`). Do not modify `src/`. If you are not confident it is dead, prefer `wontfix:<reason>` and leave it for later analysis.

- [ ] **Step 4: Regenerate and commit (batched).**

```bash
source ~/pyenv/bin/activate 2>/dev/null || true
ast_tester/coverage/run_simplify_coverage.sh
git add ast_tester/ docs/superpowers/specs/2026-05-05-simplify-branch-coverage-inventory.md
git commit -m "Absolute-pass: close/annotate simplify gaps in <class>"
```

- [ ] **Step 5: Repeat until every open `absolute-only` gap is either `fixture=…`, `unreachable:…`, or `wontfix:…`.**

---

## Task 6: Milestone verification

**Goal.** Confirm a coherent, committable state at whichever milestone the effort pauses.

- [ ] **Step 1: Full simplify suite under sanitizers.**

```bash
cmake -B build-dev -DCMAKE_BUILD_TYPE=Debug -DAST_ENABLE_WARNINGS=ON -DAST_ENABLE_SANITIZERS=ON 2>&1 | tail -3
cmake --build build-dev 2>&1 | tail -5
ctest --test-dir build-dev -R '^simplify_' --output-on-failure 2>&1 | tail -15
```

Expected: all `simplify_*` tests pass under sanitizers and warnings.

- [ ] **Step 2: Full suite shows no new regressions.**

```bash
ctest --test-dir build-dev --output-on-failure 2>&1 | tail -20
```

Expected: same pass/fail set as `master` (no regressions introduced by the new fixtures).

- [ ] **Step 3: Confirm no source was modified.**

```bash
git diff --stat master..HEAD -- src/
```

Expected: empty output.

- [ ] **Step 4: Report milestone status.** State the open differential and absolute counts from the ledger header. If differential = 0, the self-contained milestone is reached. Hand back to the user for the push/PR decision (do not push).

---

## Self-review notes

- **Spec coverage:** import (Task 1), gcovr-collection + stdlib-diff tooling (Tasks 2–3), differential mechanism + classification (Task 2), gap ledger + annotation preservation (Tasks 2–3), differential pass / primary milestone (Task 4), absolute pass (Task 5), verification gate + no-src-change + done milestone (Task 6). All spec sections map to a task.
- **No-Rust:** the contributed fixtures are referenced only by their staging path; no characterization of their origin appears.
- **Type consistency:** `load_coverage`/`load_coverage_obj`, `function_for`, `classify_gaps`, `parse_prior`, `render_ledger` signatures match between the test (Task 2 Step 1), the implementation (Step 3), and the wrapper CLI (Task 3).
