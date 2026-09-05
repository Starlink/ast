# Test Tiering for Distribution — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Ship a ~5.5 MB default `make dist` whose tests adapt to the fixtures
present, with a `make fulldist` carrying the full 20 MB corpus.

**Architecture:** Tests are gated on whether their fixtures exist, not on a
configure switch, so one mechanism serves both a cut-down tarball and a git
checkout. A tracked `fixtures/DIST_MANIFEST` declares the default keep-set; a
generated `fixtures/DIST_TIER` marks a tarball and distinguishes a deliberate
omission from a manifest typo. CMake registers absent-fixture tests as
`DISABLED` so `ctest` still reports them; the autotools `.tap` producers emit
TAP `SKIP`.

**Tech Stack:** CMake 3.24+, GNU Autotools (automake 1.16, autoconf 2.71+),
POSIX shell, `awk`, `git`.

**Spec:** `docs/superpowers/specs/2026-09-05-test-tiering-design.md`

## Global Constraints

- CMake minimum is **3.24**. Do not use `cmake_language(EXIT)` (3.29+).
- Tests that serialize must build as C11: `-DAST_C_STANDARD=11`. `AST_DBL_DIG`
  is 18 under C99 and 17 under C11, so C99 changes reference output.
- No Starlink dependencies on the CMake path.
- `autoreconf` **cannot be run in this development environment**: system
  autoconf 2.69 with m4 1.4.19 loops and writes multi-gigabyte trace files into
  `autom4te.cache/`. Every autotools task in this plan is therefore verified in
  CI, not locally. Do not attempt `./bootstrap.local`; if you do and it hangs,
  kill the `m4` process and `rm -rf autom4te.cache`.
- Never push to a remote. Never post to GitHub.
- Use `git commit --fixup` for corrections to an earlier commit in this plan.
- Commit messages end with a blank line, then `Generated with AI`, then a blank
  line, then `Co-Authored-By: SLAC AI`.
- Prose documents use one sentence per line and American English spelling.
- Changes under `src/` require a prologue History entry. **No task in this plan
  touches `src/`.**

## Deviation from the spec — read before Task 5

The spec specifies an `awk` filter that **rewrites** a shipped oracle to contain
only sections whose fixture shipped.
Task 5 implements a **verification** instead: the dist hook fails if a shipped
oracle names a fixture that did not ship.

Reason: for the chosen keep-set the filter is a no-op, because all 277 sections
of `headers.oracle` reference fixtures that ship (259 in `wcsconv/inputs/`, 18 in
`programs/joye_car_headers/`).
A rewrite would also fail silently in the one case that matters — if someone
later trims `wcsconv/inputs/`, a filter quietly drops oracle sections and reduces
coverage with nothing to say so, whereas a check fails loudly at `make dist` and
forces a decision.
This matches the existing `dist-check-tracked` philosophy.

If the reviewer prefers the spec's rewrite, Task 5 Step 3 is the only step that
changes.

## File Structure

| File | Responsibility |
| --- | --- |
| `ast_tester/fixtures/DIST_MANIFEST` | **Create.** Tracked declaration of the default tier's keep-set. Single source of truth. |
| `ast_tester/fixtures/DIST_TIER` | **Generated into the tarball only**, never tracked. Contains `cutdown` or `full`. |
| `ast_tester/CMakeLists.txt` | Reads both files; `ast_fixture_status()` replaces `ast_require_fixture()`; registers `DISABLED` tests; prints the summary. |
| `ast_tester/test_simplify.tap` | Self-skips when its fixtures are absent. |
| `ast_tester/test_wcsconverter.tap` | Self-skips when its fixtures are absent. |
| `ast_tester/test_roundtrip.tap` | Tolerates an absent `simplify/cases.txt`; globs find only what is present. |
| `ast_tester/test_keymap_mapsz.sh` | Exits 77 when its fixture is absent. |
| `ast_tester/Makefile.am` | `dist-hook` honors `DIST_MANIFEST`, writes `DIST_TIER`, checks oracle consistency. |
| `Makefile.am` | Adds `fulldist`; splits `dist-check-tracked` into source and fixture checks. |
| `.github/workflows/autotools.yaml` | Asserts the tarball's `ctest` reports skips; adds a `fulldist` check. |

---

### Task 1: Delete the unused PostScript plot references

`fixtures/plot/expected/*.ps` is 1.7 MB referenced by no `cases.txt` row, no
CMake driver and no test source. They were a visual reference a human could open
and compare against by eye; the `.svg` references now fill that role and are
additionally compared automatically by the 20 `grid_*` tests. Deleting rather
than excluding them means they stop costing 1.7 MB in `fulldist` too.

**Files:**
- Delete: `ast_tester/fixtures/plot/expected/*.ps` (20 files)

**Interfaces:**
- Consumes: nothing.
- Produces: nothing. No later task depends on this; it is first because it is
  independent and shrinks the tree everything else measures.

- [ ] **Step 1: Prove nothing references them**

```bash
cd /sdf/home/t/timj/WORK/ast/ast_tester
grep -rn '\.ps' fixtures/plot/cases.txt CMakeLists.txt *.tap *.sh 2>/dev/null
grep -rln 'plot/expected' --include='*.c' --include='*.cmake' --include='*.txt' . ../cmake
```

Expected: the first command prints nothing. The second prints only
`fixtures/plot/cases.txt`.
If either prints anything else, **stop** — a reference exists and this task's
premise is wrong.

- [ ] **Step 2: Record the size before**

```bash
du -sh fixtures
du -ch fixtures/plot/expected/*.ps | tail -1
```

Expected: `20M` and `1.7M`.

- [ ] **Step 3: Delete them**

```bash
git rm -q fixtures/plot/expected/*.ps
du -sh fixtures
```

Expected: `18M`.

- [ ] **Step 4: Verify the grid tests still pass**

```bash
cd /sdf/home/t/timj/WORK/ast
cmake --build build-dev -j 16
ctest --test-dir build-dev -R '^grid_' --output-on-failure
```

Expected: 20 tests, all pass.

- [ ] **Step 5: Commit**

```bash
git add -A ast_tester/fixtures/plot/expected
git commit -F - <<'MSG'
test: drop the PostScript plot references

fixtures/plot/expected/*.ps was 1.7MB that no cases.txt row, CMake driver or
test source referenced.  The files were a visual reference for a human to
compare a plot against by eye; the .svg references now serve that purpose and,
unlike these, are compared automatically by the 20 grid tests, so the role is
filled by files a test keeps honest.

Generated with AI

Co-Authored-By: SLAC AI
MSG
```

---

### Task 2: The manifest, the tier marker, and fixture status

Introduces the two files that let a missing fixture be classified, and replaces
`ast_require_fixture()` with a function that returns a status instead of always
aborting. No test registration changes yet, so after this task the build behaves
exactly as before in a git checkout.

**Files:**
- Create: `ast_tester/fixtures/DIST_MANIFEST`
- Modify: `ast_tester/CMakeLists.txt` (replace the function at line 84; add
  manifest/tier reading near the fixture roots at lines 51-52)

**Interfaces:**
- Consumes: `AST_FIXTURE_SOURCE_ROOT` (already defined,
  `ast_tester/CMakeLists.txt:51`).
- Produces, for Task 3:
  - `AST_DIST_TIER` — cache-free variable, one of `cutdown`, `full`, or `""`
    (empty means git checkout).
  - `ast_manifest_covers(<relative_path> <out_var>)` — sets `out_var` to `TRUE`
    or `FALSE`.
  - `ast_fixture_status(<manifest> <relative_path> <out_var>)` — sets `out_var`
    to `TRUE` when the fixture exists, `FALSE` when it is a deliberate omission,
    and calls `message(FATAL_ERROR)` otherwise.

- [ ] **Step 1: Create the manifest**

A line ending in `/` includes everything beneath that directory. Any other line
is a glob matched against the fixture path relative to `fixtures/`, where `*`
does not cross a `/`.

Create `ast_tester/fixtures/DIST_MANIFEST`:

```
# Fixtures carried by the default `make dist`.  `make fulldist` ships
# everything and ignores this file.
#
# A line ending in "/" includes every file beneath that directory.  Any other
# line is a glob matched against the path relative to fixtures/, in which "*"
# does not cross a "/".
#
# Every entry must match at least one file and every shipped fixture must match
# an entry; `make dist` fails otherwise, so this list cannot drift from what is
# actually distributed.

# Inputs for the C and Fortran test programs.
programs/

# The 20 grid tests: their manifest and their SVG references.  Their head files
# come from wcsconv/inputs/ below.
plot/cases.txt
plot/expected/*.svg

# The header oracle corpus.  Also supplies the grid head files and the SIP
# headers that testfitschan.c and testfitschan.f read directly.
wcsconv/inputs/
oracle/headers.oracle
oracle/transform_oracle_overrides.txt
oracle/README.md

# Manifests for the families whose fixtures are not distributed.  They ship so
# that ctest can enumerate those tests and report them as skipped rather than
# omitting them silently; 84K buys an honest test count.
simplify/cases.txt
wcsconv/cases.txt
```

- [ ] **Step 2: Write the failing test**

The test is a shell script that configures a simulated cut-down tree and asserts
the classification. `tools/` does not exist yet, so create it:

```bash
mkdir -p /sdf/home/t/timj/WORK/ast/tools
```

Then create `tools/test-dist-tier.sh`:

```bash
#!/bin/sh
# Verify fixture-status classification: a deliberate omission must configure
# cleanly, a missing manifest-listed fixture must fail, and a git checkout must
# keep failing on any missing fixture.
set -e
cd "$(dirname "$0")/.."
FIX=ast_tester/fixtures
STASH=$(mktemp -d)
fail=0

cleanup() {
    rm -f "$FIX/DIST_TIER"
    test -d "$STASH/expected" && mv "$STASH/expected" "$FIX/wcsconv/expected"
    test -f "$STASH/headers.oracle" && \
        mv "$STASH/headers.oracle" "$FIX/oracle/headers.oracle"
    rm -rf "$STASH" _tiertest
}
trap cleanup EXIT

say() { printf '%s\n' "$*"; }
check() {
    if [ "$1" = "$2" ]; then say "  ok   - $3"; else
        say "  FAIL - $3 (expected $2, got $1)"; fail=1; fi
}

configure() {
    rm -rf _tiertest
    cmake -B _tiertest -DAST_C_STANDARD=11 > "$STASH/log" 2>&1 && echo pass || echo fail
}

say "case 1: git checkout, all fixtures present -> configure succeeds"
check "$(configure)" pass "baseline configures"

say "case 2: cutdown tarball, wcsconv/expected omitted -> configure succeeds"
mv "$FIX/wcsconv/expected" "$STASH/expected"
echo cutdown > "$FIX/DIST_TIER"
check "$(configure)" pass "deliberate omission is tolerated"

say "case 3: git checkout, wcsconv/expected missing -> configure FAILS"
rm -f "$FIX/DIST_TIER"
check "$(configure)" fail "missing fixture in a checkout is fatal"
mv "$STASH/expected" "$FIX/wcsconv/expected"

say "case 4: cutdown tarball, a manifest-listed fixture missing -> FAILS"
mv "$FIX/oracle/headers.oracle" "$STASH/headers.oracle"
echo cutdown > "$FIX/DIST_TIER"
check "$(configure)" fail "missing manifest-listed fixture is fatal"
mv "$STASH/headers.oracle" "$FIX/oracle/headers.oracle"

exit $fail
```

```bash
chmod +x tools/test-dist-tier.sh
```

- [ ] **Step 3: Run it to make sure it fails**

```bash
cd /sdf/home/t/timj/WORK/ast && ./tools/test-dist-tier.sh
```

Expected: cases 1, 3 and 4 report `ok`, and **case 2 FAILS**.

Cases 3 and 4 both expect a failed configure, which today's
`ast_require_fixture` delivers for any missing fixture, so they pass from the
start; they are here to stay green through the change rather than to drive it.
Case 2 is the one this task exists for: a fixture omitted on purpose must
configure cleanly, and today it aborts.

- [ ] **Step 4: Implement the manifest reader and status function**

In `ast_tester/CMakeLists.txt`, immediately after line 52
(`set(AST_FIXTURE_BINARY_ROOT ...)`), insert:

```cmake
# The distribution tier, written into a tarball by the dist hook and absent from
# a git checkout.  It is what lets a missing fixture be read as a deliberate
# omission rather than as a typo in a manifest.
set(AST_DIST_TIER "")
if(EXISTS "${AST_FIXTURE_SOURCE_ROOT}/DIST_TIER")
    file(READ "${AST_FIXTURE_SOURCE_ROOT}/DIST_TIER" AST_DIST_TIER)
    string(STRIP "${AST_DIST_TIER}" AST_DIST_TIER)
endif()

# The default tier's keep-set.  Comments and blank lines are dropped.
set(AST_DIST_MANIFEST_PATTERNS "")
if(EXISTS "${AST_FIXTURE_SOURCE_ROOT}/DIST_MANIFEST")
    file(STRINGS "${AST_FIXTURE_SOURCE_ROOT}/DIST_MANIFEST" _dm_lines)
    foreach(_line IN LISTS _dm_lines)
        string(STRIP "${_line}" _line)
        if(NOT _line STREQUAL "" AND NOT _line MATCHES "^#")
            list(APPEND AST_DIST_MANIFEST_PATTERNS "${_line}")
        endif()
    endforeach()
endif()
```

Then replace the whole `ast_require_fixture` function at line 84 with:

```cmake
# TRUE in out_var when "relative_path" matches a DIST_MANIFEST entry.  A pattern
# ending in "/" is a directory prefix; anything else is a glob in which "*" does
# not cross a "/".
function(ast_manifest_covers relative_path out_var)
    foreach(_pat IN LISTS AST_DIST_MANIFEST_PATTERNS)
        if(_pat MATCHES "/$")
            string(LENGTH "${_pat}" _n)
            string(LENGTH "${relative_path}" _rn)
            if(NOT _rn LESS _n)
                string(SUBSTRING "${relative_path}" 0 ${_n} _head)
                if(_head STREQUAL "${_pat}")
                    set(${out_var} TRUE PARENT_SCOPE)
                    return()
                endif()
            endif()
        else()
            string(REPLACE "." "\\." _re "${_pat}")
            string(REPLACE "*" "[^/]*" _re "${_re}")
            if(relative_path MATCHES "^${_re}$")
                set(${out_var} TRUE PARENT_SCOPE)
                return()
            endif()
        endif()
    endforeach()
    set(${out_var} FALSE PARENT_SCOPE)
endfunction()

# Classify a fixture named by a manifest.  Sets out_var TRUE when it is present,
# FALSE when it was deliberately left out of a cut-down distribution, and fails
# the configure otherwise.
#
# The FATAL_ERROR cases matter as much as the skip.  A typo in a cases.txt used
# to surface only as a harness failure part-way through a ctest run, with
# nothing to say the file was never there, and a git checkout -- where typos are
# actually introduced -- keeps exactly that strictness.  A fixture that the
# manifest says should have shipped but did not is a packaging defect and is
# likewise fatal, so the relaxation covers only files the manifest agrees are
# absent on purpose.
function(ast_fixture_status manifest relative_path out_var)
    if(EXISTS "${AST_FIXTURE_SOURCE_ROOT}/${relative_path}")
        set(${out_var} TRUE PARENT_SCOPE)
        return()
    endif()
    if(NOT AST_DIST_TIER STREQUAL "cutdown")
        message(FATAL_ERROR
            "${manifest}: no such fixture: fixtures/${relative_path}")
    endif()
    ast_manifest_covers("${relative_path}" _covered)
    if(_covered)
        message(FATAL_ERROR
            "${manifest}: fixtures/${relative_path} matches a DIST_MANIFEST "
            "entry but is missing from this distribution")
    endif()
    set(${out_var} FALSE PARENT_SCOPE)
endfunction()

# Retained so that call sites which cannot yet skip keep their old behavior.
function(ast_require_fixture manifest relative_path)
    ast_fixture_status("${manifest}" "${relative_path}" _ok)
    if(NOT _ok)
        message(FATAL_ERROR
            "${manifest}: fixtures/${relative_path} is not distributed, and "
            "this call site cannot skip")
    endif()
endfunction()
```

- [ ] **Step 5: Run the test to verify all four cases pass**

```bash
cd /sdf/home/t/timj/WORK/ast && ./tools/test-dist-tier.sh
```

Expected: cases 1, 3 and 4 `ok`, and **case 2 still FAILS**.

`ast_fixture_status` now classifies the omission correctly, but every call site
still goes through the `ast_require_fixture` wrapper, which aborts. Case 2 turns
green in Task 3 Step 9, once the call sites are converted. This is the one place
in the plan where a task ends with a known-red check; it is recorded here so
that it is not mistaken for a mistake.

- [ ] **Step 6: Verify the checkout build is unchanged**

```bash
cmake -B build-dev -DAST_C_STANDARD=11 && ctest --test-dir build-dev -j 16 | tail -3
```

Expected: `100% tests passed`, same total as before this task.

- [ ] **Step 7: Commit**

```bash
git add ast_tester/fixtures/DIST_MANIFEST ast_tester/CMakeLists.txt tools/test-dist-tier.sh
git commit -F - <<'MSG'
build: classify a missing fixture as omitted or as a defect

A cut-down distribution leaves fixtures out on purpose, which the fixture check
cannot tell from a typo in a cases.txt -- and catching that typo at configure
time, rather than part-way through a ctest run, is why the check exists.

DIST_MANIFEST states the default tier's keep-set and DIST_TIER marks a tarball,
so the two cases separate: a git checkout stays as strict as before, a tarball
tolerates only what the manifest agrees is absent, and a fixture the manifest
says should have shipped is a packaging defect either way.

Generated with AI

Co-Authored-By: SLAC AI
MSG
```

---

### Task 3: Register absent-fixture tests as skipped

Converts the call sites so a family whose fixtures are missing is registered and
reported by `ctest` rather than aborting the configure. Omitting the tests
entirely would print `113 tests` with nothing to say 2362 more exist, leaving
"did I run the full suite?" unanswerable.

**Files:**
- Modify: `ast_tester/CMakeLists.txt` — lines 315-316 (`add_wcsconv_test`),
  429 (`add_roundtrip_test`), 479-480 (`keymap_mapsz`), 551-552
  (`add_simplify_test`), 699 and 709 (`add_grid_test`), 836
  (`add_plotter_test`), and the ends of the `wcsconv`/`simplify` loops at 356
  and 602

**Interfaces:**
- Consumes: `ast_fixture_status()`, `AST_DIST_TIER` from Task 2.
- Produces: `ast_add_skipped_test(<name>)` — registers `name` as a disabled
  test and increments the global property `AST_SKIPPED_TEST_COUNT`.

- [ ] **Step 1: Add the skip helper**

In `ast_tester/CMakeLists.txt`, immediately after the `ast_fixture_status`
function from Task 2, insert:

```cmake
# Register a test that cannot run because its fixture is not distributed.
#
# DISABLED rather than SKIP_RETURN_CODE: the latter needs the test to run and
# return the code, which a test with no input cannot do without a stub program.
# ctest names a disabled test under "the following tests did not run", so the
# count stays visible instead of the suite silently shrinking.
function(ast_add_skipped_test name)
    add_test(NAME ${name} COMMAND "${CMAKE_COMMAND}" -E true)
    set_tests_properties(${name} PROPERTIES DISABLED TRUE)
    get_property(_n GLOBAL PROPERTY AST_SKIPPED_TEST_COUNT)
    if(NOT _n)
        set(_n 0)
    endif()
    math(EXPR _n "${_n} + 1")
    set_property(GLOBAL PROPERTY AST_SKIPPED_TEST_COUNT ${_n})
endfunction()
```

- [ ] **Step 2: Convert `add_wcsconv_test`**

Replace lines 315-316:

```cmake
    ast_require_fixture("wcsconv/cases.txt" "${W_INPUT}")
    ast_require_fixture("wcsconv/cases.txt" "${W_REFERENCE}")
```

with:

The skip branch must mirror the conditions the real registration uses, not assume
a fixed pair, or the reported count will not match the tests that would have
existed. `add_wcsconv_test` registers the bare test only when
`SKIP_STRING_COMPARE` is absent, and `_astequal` always:

```cmake
    ast_fixture_status("wcsconv/cases.txt" "${W_INPUT}" _in_ok)
    ast_fixture_status("wcsconv/cases.txt" "${W_REFERENCE}" _ref_ok)
    if(NOT _in_ok OR NOT _ref_ok)
        if(NOT W_SKIP_STRING_COMPARE)
            ast_add_skipped_test(wcsconv_${W_NAME})
        endif()
        ast_add_skipped_test(wcsconv_${W_NAME}_astequal)
        return()
    endif()
```

- [ ] **Step 3: Convert `add_simplify_test`**

Replace lines 551-552 the same way:

`add_simplify_test` registers three tests under three different conditions: the
bare test only when `SKIP_STRING_COMPARE` is absent, `_noop` only when input and
reference are the same file, and `_astequal` always. The skip branch mirrors all
three:

```cmake
    ast_fixture_status("simplify/cases.txt" "${S_INPUT}" _in_ok)
    ast_fixture_status("simplify/cases.txt" "${S_REFERENCE}" _ref_ok)
    if(NOT _in_ok OR NOT _ref_ok)
        if(NOT S_SKIP_STRING_COMPARE)
            ast_add_skipped_test(simplify_${S_NAME})
        endif()
        if(S_INPUT STREQUAL S_REFERENCE)
            ast_add_skipped_test(simplify_${S_NAME}_noop)
        endif()
        ast_add_skipped_test(simplify_${S_NAME}_astequal)
        return()
    endif()
```

Measured in a full checkout, these conditions yield 982 simplify tests: 490
`_astequal`, 157 `_noop`, and 335 bare. The skipped count must reproduce those
numbers exactly, which is why the conditions are mirrored rather than assumed.

- [ ] **Step 4: Convert `add_roundtrip_test`**

Replace line 429:

```cmake
    ast_require_fixture("roundtrip" "${R_FIXTURE}")
```

with:

```cmake
    ast_fixture_status("roundtrip" "${R_FIXTURE}" _rt_ok)
    if(NOT _rt_ok)
        ast_add_skipped_test(roundtrip_${R_NAME})
        return()
    endif()
```

- [ ] **Step 5: Convert `keymap_mapsz`**

Replace line 479 and wrap the `add_test` that follows it:

```cmake
ast_fixture_status("keymap_mapsz" "oracle/keymap_mapsz.txt" _km_ok)
if(_km_ok)
    add_test(NAME keymap_mapsz
        COMMAND "${CMAKE_COMMAND}"
                -DGENERATOR=$<TARGET_FILE:gen_keymap_mapsz>
                -DREF_FILE=${AST_FIXTURE_BINARY_ROOT}/oracle/keymap_mapsz.txt
                -DOUT_FILE=${CMAKE_CURRENT_BINARY_DIR}/keymap_mapsz.out
                -P "${CMAKE_SOURCE_DIR}/cmake/run_keymap_mapsz_test.cmake"
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
    if(AST_ENABLE_SANITIZERS)
        set_tests_properties(keymap_mapsz PROPERTIES
            ENVIRONMENT "ASAN_OPTIONS=detect_leaks=0")
    endif()
else()
    ast_add_skipped_test(keymap_mapsz)
endif()
add_dependencies(gen_keymap_mapsz ast_stage_fixtures)
```

- [ ] **Step 6: Convert the oracle tests whose oracle file may be absent**

Replace the three `add_test(NAME transform_oracle_*)` blocks at lines 188-202 and
the `set_tests_properties` call at lines 204-206 with the single loop below.
Leave `transform_oracle_selftest` at line 184 exactly as it is: it reads no
fixture and always runs.

Note the file name and the test name differ for one entry — the oracle file is
`simplify_fixtures.oracle` but the test is `transform_oracle_simplify` — so the
loop carries both:

```cmake
# Each oracle is gated on its own file: a cut-down distribution carries the
# header oracle and not the other two.
foreach(_orc simplify_fixtures headers framesets)
    if(_orc STREQUAL "simplify_fixtures")
        set(_tname "simplify")
    else()
        set(_tname "${_orc}")
    endif()
    if(EXISTS "${AST_FIXTURE_SOURCE_ROOT}/oracle/${_orc}.oracle")
        add_test(NAME transform_oracle_${_tname}
                 COMMAND check_transform_oracle
                         "${AST_FIXTURE_BINARY_ROOT}"
                         "${AST_FIXTURE_BINARY_ROOT}/oracle/${_orc}.oracle"
                         "${_oracle_overrides}")
        if(AST_ENABLE_SANITIZERS)
            set_tests_properties(transform_oracle_${_tname} PROPERTIES
                ENVIRONMENT "ASAN_OPTIONS=detect_leaks=0")
        endif()
    else()
        ast_add_skipped_test(transform_oracle_${_tname})
    endif()
endforeach()
```

Confirm the three test names are unchanged from before this step:

```bash
cd /sdf/home/t/timj/WORK/ast
cmake -B build-dev -DAST_C_STANDARD=11 >/dev/null
ctest --test-dir build-dev -N -R '^transform_oracle' | grep Test:
```

Expected, in some order: `transform_oracle_selftest`,
`transform_oracle_simplify`, `transform_oracle_headers`,
`transform_oracle_framesets`.

- [ ] **Step 7: Guard the manifest-reading generator**

`gen_frameset_fixtures` reads `wcsconv/cases.txt`, which always ships, so it
needs no guard. Confirm that is still true:

```bash
grep -n 'cases.txt' ast_tester/gen_frameset_fixtures.c ast_tester/gen_native_fixtures.c
```

If either names a manifest that `DIST_MANIFEST` does not ship, wrap its
`add_executable` in an `if(EXISTS ...)`. Otherwise make no change.

- [ ] **Step 8: Print the configure summary**

At the very end of `ast_tester/CMakeLists.txt`, append:

```cmake
get_property(_ast_skipped GLOBAL PROPERTY AST_SKIPPED_TEST_COUNT)
if(_ast_skipped)
    message(STATUS
        "AST test fixtures: cut-down distribution (fixtures/DIST_TIER)")
    message(STATUS
        "  ${_ast_skipped} tests disabled: fixtures not distributed")
    message(STATUS
        "  use `make fulldist` from a git checkout for the full corpus")
endif()
```

- [ ] **Step 9: Verify the checkout is unchanged and the cut-down case works**

```bash
cd /sdf/home/t/timj/WORK/ast
cmake -B build-dev -DAST_C_STANDARD=11 && ctest --test-dir build-dev -j 16 | tail -3
./tools/test-dist-tier.sh
```

Expected: `100% tests passed` with the same total as before (no summary printed,
because nothing was skipped), and all four cases of the tier test `ok`.

- [ ] **Step 10: Verify skip reporting end to end**

```bash
cd /sdf/home/t/timj/WORK/ast
mv ast_tester/fixtures/wcsconv/expected /tmp/wc-expected
echo cutdown > ast_tester/fixtures/DIST_TIER
cmake -B /tmp/cutdown-build -DAST_C_STANDARD=11 2>&1 | grep -A3 'cut-down'
ctest --test-dir /tmp/cutdown-build -N | tail -2
ctest --test-dir /tmp/cutdown-build -j 16 2>&1 | tail -6
rm -f ast_tester/fixtures/DIST_TIER
mv /tmp/wc-expected ast_tester/fixtures/wcsconv/expected
rm -rf /tmp/cutdown-build
```

Expected: the summary names 386 disabled tests; `ctest` lists them under "the
following tests did not run"; the remaining tests pass. The total test count
stays 2475.

- [ ] **Step 11: Commit**

```bash
git add ast_tester/CMakeLists.txt
git commit -F - <<'MSG'
build: report tests whose fixtures were not distributed

A family whose fixtures a cut-down distribution omits is now registered and
marked DISABLED rather than aborting the configure, so ctest still names it and
the count stays visible.  Omitting them would print a smaller suite with nothing
to say how much smaller, which is the one question an end user running the tests
from a tarball needs answered.

DISABLED rather than SKIP_RETURN_CODE because the latter needs the test to run
and return the code, which a test with no input cannot do without a stub.

Generated with AI

Co-Authored-By: SLAC AI
MSG
```

---

### Task 4: Make the autotools test scripts self-skip

The `TESTS` list stays static and the scripts decide at run time, so
`Makefile.am` needs no conditionals.

**Files:**
- Modify: `ast_tester/test_simplify.tap`
- Modify: `ast_tester/test_wcsconverter.tap`
- Modify: `ast_tester/test_roundtrip.tap`
- Modify: `ast_tester/test_keymap_mapsz.sh`

**Interfaces:**
- Consumes: nothing from earlier tasks. Independent of the CMake work.
- Produces: nothing later tasks depend on.

- [ ] **Step 1: Write the failing test**

Create `tools/test-tap-skip.sh`:

```bash
#!/bin/sh
# A .tap producer whose fixtures are absent must emit a valid skipped plan and
# exit 0; test_keymap_mapsz.sh must exit 77.
set -e
cd "$(dirname "$0")/.."
SRC=$PWD/ast_tester
BUILD=$PWD/build-dev/ast_tester
EMPTY=$(mktemp -d)
mkdir -p "$EMPTY/fixtures"
fail=0

check() {
    if [ "$1" = "$2" ]; then echo "  ok   - $3"; else
        echo "  FAIL - $3 (expected '$2', got '$1')"; fail=1; fi
}

cd "$BUILD"
for t in test_simplify test_wcsconverter test_roundtrip; do
    out=$(srcdir="$EMPTY" sh "$SRC/$t.tap" 2>&1); rc=$?
    check "$rc" 0 "$t.tap exits 0 with no fixtures"
    check "$(printf '%s' "$out" | head -1)" "1..0 # SKIP fixtures not distributed" \
          "$t.tap emits a skipped plan"
done

srcdir="$EMPTY" sh "$SRC/test_keymap_mapsz.sh" >/dev/null 2>&1; rc=$?
check "$rc" 77 "test_keymap_mapsz.sh exits 77 with no fixture"

rm -rf "$EMPTY"
exit $fail
```

```bash
chmod +x tools/test-tap-skip.sh
```

- [ ] **Step 2: Run it to verify it fails**

```bash
cd /sdf/home/t/timj/WORK/ast && ./tools/test-tap-skip.sh
```

Expected: all seven checks FAIL — the scripts currently crash or emit a nonzero
plan when the manifest is missing.

- [ ] **Step 3: Add the guard to `test_simplify.tap`**

After the line `data=${fixtures}/simplify/cases.txt`, insert:

```sh
# The cut-down distribution ships this manifest but not the fixtures it names,
# so report the family as skipped rather than failing every case.
if [ ! -f "$data" ] || [ ! -d "${fixtures}/simplify" ] || \
   ! ls "${fixtures}"/simplify/*.map >/dev/null 2>&1; then
    echo "1..0 # SKIP fixtures not distributed"
    exit 0
fi
```

- [ ] **Step 4: Add the guard to `test_wcsconverter.tap`**

After the line `data=${fixtures}/wcsconv/cases.txt`, insert:

```sh
# As test_simplify.tap: the manifest ships, the references may not.
if [ ! -f "$data" ] || [ ! -d "${fixtures}/wcsconv/expected" ]; then
    echo "1..0 # SKIP fixtures not distributed"
    exit 0
fi
```

- [ ] **Step 5: Add the guard to `test_roundtrip.tap`**

`test_roundtrip.tap` reads `simplify/cases.txt` and then globs. Replace the
`list=$(...)` assignment's `while` input so an absent manifest is not an error,
and skip only when the resulting list is empty. After the `cases=` line, insert:

```sh
# The simplify half of the corpus may not be distributed; the glob half then
# finds only the .ast fixtures that were.  Skip only if neither is present.
test -f "$cases" || cases=/dev/null
```

and replace the plan line

```sh
echo "1..$(printf '%s\n' "$list" | grep -c .)"
```

with

```sh
n=$(printf '%s\n' "$list" | grep -c . || true)
if [ "$n" -eq 0 ]; then
    echo "1..0 # SKIP fixtures not distributed"
    exit 0
fi
echo "1..$n"
```

- [ ] **Step 6: Add the guard to `test_keymap_mapsz.sh`**

After the `ref=` line, insert:

```sh
# Not distributed outside a git checkout: this checks a committed capture of
# KeyMap's growth policy, which a user's build cannot affect.
if [ ! -f "$ref" ]; then
    echo "keymap_mapsz: skipped (fixture not distributed)"
    exit 77
fi
```

- [ ] **Step 7: Run the test to verify it passes**

```bash
cd /sdf/home/t/timj/WORK/ast && ./tools/test-tap-skip.sh
```

Expected: all seven `ok`.

- [ ] **Step 8: Verify the normal path still works**

```bash
cd /sdf/home/t/timj/WORK/ast/build-dev/ast_tester
export srcdir=/sdf/home/t/timj/WORK/ast/ast_tester ASAN_OPTIONS=detect_leaks=0
sh $srcdir/test_wcsconverter.tap | tail -1
sh $srcdir/test_keymap_mapsz.sh; echo "keymap rc=$?"
```

Expected: the wcsconv producer's last line is `ok 193 - ...` and the keymap
script exits 0.

- [ ] **Step 9: Commit**

```bash
cd /sdf/home/t/timj/WORK/ast
git add ast_tester/test_simplify.tap ast_tester/test_wcsconverter.tap \
        ast_tester/test_roundtrip.tap ast_tester/test_keymap_mapsz.sh \
        tools/test-tap-skip.sh
git commit -F - <<'MSG'
test: skip the data-driven families when their fixtures are absent

The cut-down distribution ships the cases.txt manifests but not every fixture
they name, so each producer now emits a TAP skipped plan instead of failing
every row, and the keymap capture check exits 77 as testhuge.sh does.

Keeping the decision in the scripts rather than in Makefile.am means the TESTS
list stays static and the same suite adapts to whatever was distributed.

Generated with AI

Co-Authored-By: SLAC AI
MSG
```

---

### Task 5: Ship only the manifest's fixtures

**Files:**
- Modify: `ast_tester/Makefile.am` — the `dist-hook` at line 447

**Interfaces:**
- Consumes: `ast_tester/fixtures/DIST_MANIFEST` from Task 2.
- Produces: `fixtures/DIST_TIER` inside `$(distdir)`, read by Task 2's CMake
  code and by Task 7's guard. Honors `AST_DIST_FULL`, set by Task 6.

- [ ] **Step 1: Replace the fixture copy in the dist hook**

Replace the first `@list=` block of `dist-hook` (the `find fixtures -type f`
one) with:

```make
# Fixtures are shipped through fixtures/DIST_MANIFEST rather than wholesale.
# The tree is 18MB, larger than src/, and an end user verifying a build needs a
# fraction of it; `make fulldist' sets AST_DIST_FULL=1 to ship everything.
#
# DIST_TIER records which happened, so that the CMake build unpacked from this
# tarball can tell a fixture left out on purpose from one that went missing.
	@if test "x$(AST_DIST_FULL)" = x1; then \
	    list=`cd "$(srcdir)" && find fixtures -type f ! -name '.*'`; \
	    tier=full; \
	else \
	    list=`cd "$(srcdir)" && $(SHELL) ./dist-fixtures.sh fixtures/DIST_MANIFEST`; \
	    tier=cutdown; \
	fi; \
	for rel in $$list; do \
	    $(MKDIR_P) "$(distdir)/`dirname $$rel`"; \
	    cp -p "$(srcdir)/$$rel" "$(distdir)/$$rel" || exit 1; \
	done; \
	echo $$tier > "$(distdir)/fixtures/DIST_TIER"
```

- [ ] **Step 2: Create the manifest expander**

Create `ast_tester/dist-fixtures.sh` with exactly this content — it has been
run against the real fixture tree and against a deliberately broken manifest:

```sh
#!/bin/sh
# Print the fixture paths a manifest selects, relative to the ast_tester
# directory, one per line.
#
# A line ending in "/" includes everything beneath that directory; any other
# line is a glob matched against the path relative to fixtures/, in which "*"
# does not cross a "/".  Comments and blank lines are ignored.
#
# Every entry must match at least one file: a pattern that has stopped matching
# means the manifest no longer describes the tree, and shipping silently fewer
# fixtures than intended is the failure this exists to prevent.
set -e
manifest=${1:?usage: dist-fixtures.sh <manifest>}

pats=`sed -e 's/#.*//' -e 's/[[:space:]]*$//' "$manifest" | grep . || true`
out=`mktemp`
trap 'rm -f "$out"' EXIT

# A for-loop, not a "| while read" pipeline: the loop body must be able to fail
# the whole script, and a pipeline runs it in a subshell where exit cannot.
for pat in $pats; do
    case $pat in
        */) found=`find "fixtures/$pat" -type f ! -name '.*' 2>/dev/null || true` ;;
        *)  found=`for f in fixtures/$pat; do test -f "$f" && echo "$f"; done || true` ;;
    esac
    if test -z "$found"; then
        echo "dist-fixtures.sh: $manifest: pattern matches nothing: $pat" >&2
        exit 1
    fi
    printf '%s\n' "$found" >> "$out"
done
sort -u "$out"
```

```bash
chmod +x ast_tester/dist-fixtures.sh
```

Add it to `EXTRA_DIST` in `ast_tester/Makefile.am` by extending the line at 428:

```make
EXTRA_DIST = testhuge.sh test_keymap_mapsz.sh dist-fixtures.sh \
    test_wcsconverter.tap test_simplify.tap test_grid.tap test_roundtrip.tap
```

- [ ] **Step 3: Add the oracle consistency check**

Append to `dist-hook`:

```make
# A shipped oracle must not name a fixture that was not shipped.  The sections
# are keyed by fixture path, so this is a text check rather than a rewrite --
# and a check rather than a filter deliberately: silently dropping sections
# would reduce coverage with nothing to say so, whereas failing here forces the
# manifest and the oracle to be reconciled.
	@cd "$(distdir)/fixtures" && \
	find . -name '*.oracle' -type f | while read -r orc; do \
	    sed -n 's/^\[\([^ ]*\).*/\1/p' "$$orc" | sort -u | while read -r fx; do \
	        test -e "$$fx" || { \
	            echo "make dist: $$orc names undistributed fixture: $$fx" >&2; \
	            exit 1; }; \
	    done || exit 1; \
	done
```

- [ ] **Step 4: Verify the expander locally**

The dist hook itself cannot be run here (`autoreconf` loops, see Global
Constraints), but the expander is a plain script:

```bash
cd /sdf/home/t/timj/WORK/ast/ast_tester
./dist-fixtures.sh fixtures/DIST_MANIFEST > /tmp/shipped.txt
echo "files: $(wc -l < /tmp/shipped.txt)"
du -ch $(cat /tmp/shipped.txt) 2>/dev/null | tail -1
echo "--- every shipped path exists? ---"
while read -r f; do test -e "$f" || echo "MISSING $f"; done < /tmp/shipped.txt
echo "--- nothing from the dropped blocks? ---"
grep -E '^fixtures/(simplify/[^c]|wcsconv/(expected|framesets)|serialisation)' \
    /tmp/shipped.txt || echo "  none (correct)"
```

Expected: 209 files totaling 5.4M; no `MISSING` lines; no paths from the
dropped blocks.

- [ ] **Step 5: Verify the oracle check would pass**

```bash
cd /sdf/home/t/timj/WORK/ast/ast_tester
sed -n 's/^\[\([^ ]*\).*/\1/p' fixtures/oracle/headers.oracle | sort -u \
  | while read -r fx; do grep -qx "fixtures/$fx" /tmp/shipped.txt || echo "NOT SHIPPED $fx"; done
echo "(no output above means headers.oracle is internally consistent)"
```

Expected: no output. All 277 sections reference shipped fixtures.

- [ ] **Step 6: Commit**

```bash
cd /sdf/home/t/timj/WORK/ast
git add ast_tester/Makefile.am ast_tester/dist-fixtures.sh
git commit -F - <<'MSG'
build: ship the fixtures DIST_MANIFEST selects

The fixture tree is 18MB, larger than src/, and all of it shipped.  An end user
verifying a build needs a fraction: the binary tests' inputs, the SVG plot
references, and the FITS headers the transform oracle covers.

dist-fixtures.sh expands the manifest and fails on a pattern that matches
nothing, so the list cannot quietly stop describing the tree.  DIST_TIER records
which tier was built, and a shipped oracle is checked to name only shipped
fixtures -- a check rather than a filter, because dropping sections silently
would reduce coverage with nothing to say so.

Generated with AI

Co-Authored-By: SLAC AI
MSG
```

---

### Task 6: `make fulldist`

**Files:**
- Modify: `Makefile.am` — add the target near the existing `dist-hook`

**Interfaces:**
- Consumes: `AST_DIST_FULL` handling from Task 5.
- Produces: nothing later tasks depend on.

- [ ] **Step 1: Add the target**

Add to the top-level `Makefile.am`, immediately before `dist-check-tracked:` at
line 1093:

```make
# The full test corpus, for verifying a release rather than for shipping one.
# A distinct archive name so that a full tarball is never mistaken for the
# release artifact.
#
# Overriding distdir on the command line is what renames the archive; if a
# future automake computes the archive name independently of distdir, fall back
# to running dist and renaming the result.
fulldist:
	$(MAKE) AST_DIST_FULL=1 distdir=$(PACKAGE)-$(VERSION)-full dist

.PHONY: fulldist
```

- [ ] **Step 2: Verify it cannot be tested locally, and record that**

```bash
cd /sdf/home/t/timj/WORK/ast
ls configure 2>/dev/null || echo "no configure: autoreconf loops here, see Global Constraints"
```

Expected: the fallback message. This target is verified by CI in Task 8.

- [ ] **Step 3: Commit**

```bash
git add Makefile.am
git commit -F - <<'MSG'
build: add make fulldist for the complete fixture corpus

make dist now ships a cut-down fixture set, which is right for an end user
checking that a build works and wrong for verifying a release.  fulldist ships
everything under a distinct archive name so the two cannot be confused.

Generated with AI

Co-Authored-By: SLAC AI
MSG
```

---

### Task 7: Split the `dist-check-tracked` guard

The guard requires every tracked file under `cmake/`, `ast_tester/` and `src/`
to be in the tarball, which is what caught `testimmutable.c` and
`cmake/run_simplify_noop_test.cmake` before release. A cut-down dist omits
tracked fixtures on purpose, so the guard must be split rather than weakened.

**Files:**
- Modify: `Makefile.am` — `DIST_CHECK_EXCLUDE` at 1076, `dist-check-tracked` at
  1093

**Interfaces:**
- Consumes: `DIST_TIER` written by Task 5; `dist-fixtures.sh` from Task 5.
- Produces: nothing later tasks depend on.

- [ ] **Step 1: Exclude fixtures from the source check**

Add to `DIST_CHECK_EXCLUDE`, after the `':!ast_tester/simplify_pathways.md' \`
line:

```make
    ':!ast_tester/fixtures/*'
```

- [ ] **Step 2: Add the fixture check**

Replace the `dist-check-tracked:` recipe with:

```make
#  Two checks, because the two kinds of file fail differently.  A source file
#  must always ship: leaving one out makes the tarball configure and then fail
#  to build, which is how testimmutable.c and run_simplify_noop_test.cmake each
#  survived to a full distcheck.  A fixture may be left out on purpose, so what
#  matters there is that the tarball and DIST_MANIFEST agree.
dist-check-tracked:
	@test -d "$(top_srcdir)/.git" || exit 0; \
	tracked=`cd "$(top_srcdir)" && \
	    git ls-files cmake ast_tester src $(DIST_CHECK_EXCLUDE)`; \
	missing=; \
	for rel in $$tracked; do \
	    test -e "$(distdir)/$$rel" || missing="$$missing $$rel"; \
	done; \
	if test -n "$$missing"; then \
	    echo "make dist: tracked files missing from the distribution:" >&2; \
	    for rel in $$missing; do echo "    $$rel" >&2; done; \
	    echo "Add them to EXTRA_DIST, or widen the dist-hook patterns." >&2; \
	    exit 1; \
	fi
	@test -d "$(top_srcdir)/.git" || exit 0; \
	if test "x$(AST_DIST_FULL)" = x1; then \
	    expected=`cd "$(top_srcdir)" && git ls-files 'ast_tester/fixtures/*' \
	        | sed 's|^ast_tester/||'`; \
	else \
	    expected=`cd "$(top_srcdir)/ast_tester" && \
	        $(SHELL) ./dist-fixtures.sh fixtures/DIST_MANIFEST`; \
	fi; \
	missing=; \
	for rel in $$expected; do \
	    test -e "$(distdir)/ast_tester/$$rel" || missing="$$missing $$rel"; \
	done; \
	if test -n "$$missing"; then \
	    echo "make dist: fixtures missing from the distribution:" >&2; \
	    for rel in $$missing; do echo "    $$rel" >&2; done; \
	    exit 1; \
	fi; \
	got=`mktemp`; want=`mktemp`; \
	( cd "$(distdir)/ast_tester" && find fixtures -type f ! -name '.*' \
	    ! -name DIST_TIER ) | sort > "$$got"; \
	printf '%s\n' $$expected | sort > "$$want"; \
	extra=`comm -23 "$$got" "$$want"`; \
	rm -f "$$got" "$$want"; \
	if test -n "$$extra"; then \
	    echo "make dist: fixtures shipped but not in DIST_MANIFEST:" >&2; \
	    for rel in $$extra; do echo "    $$rel" >&2; done; \
	    exit 1; \
	fi
```

- [ ] **Step 3: Verify the expected-set computation locally**

The recipe needs automake, but its two halves are shell:

```bash
cd /sdf/home/t/timj/WORK/ast/ast_tester
./dist-fixtures.sh fixtures/DIST_MANIFEST | sort > /tmp/want.txt
git ls-files 'fixtures/*' | sort > /tmp/tracked.txt
echo "want: $(wc -l < /tmp/want.txt)  tracked: $(wc -l < /tmp/tracked.txt)"
echo "--- wanted but untracked (should be empty) ---"
comm -23 /tmp/want.txt /tmp/tracked.txt
```

Expected: `want` 209, `tracked` around 1450, and nothing in the third
section — every file the manifest selects is tracked in git.

- [ ] **Step 4: Commit**

```bash
cd /sdf/home/t/timj/WORK/ast
git add Makefile.am
git commit -F - <<'MSG'
build: check sources and fixtures separately in the dist guard

The guard required every tracked file under cmake/, ast_tester/ and src/ to be
in the tarball, and a cut-down dist omits tracked fixtures on purpose.
Weakening it would discard what it was built for -- it is what caught
testimmutable.c and run_simplify_noop_test.cmake.

Sources keep the existing check in both tiers, so a new .c, .tap or .cmake still
cannot silently fall out.  Fixtures are checked against DIST_MANIFEST in both
directions, so a file shipped without an entry and an entry matching nothing are
both errors, and the fixture set stays a stated decision rather than an
oversight.

Generated with AI

Co-Authored-By: SLAC AI
MSG
```

---

### Task 8: Prove it in CI

The autotools half of this plan cannot be verified locally. CI already unpacks
the tarball and runs `ctest` in it at `.github/workflows/autotools.yaml:117-124`,
which is the requirement driving the whole design; this task asserts the tier
behavior there and adds a `fulldist` check.

**Files:**
- Modify: `.github/workflows/autotools.yaml`

**Interfaces:**
- Consumes: everything from Tasks 1-7.
- Produces: nothing.

- [ ] **Step 1: Assert the tarball is cut down and reports its skips**

In the existing "Build the distribution tarball with CMake" step, after the
`tar xzf` line, insert:

```yaml
          test -f tarball-cmake/ast_tester/fixtures/DIST_TIER
          grep -qx cutdown tarball-cmake/ast_tester/fixtures/DIST_TIER
          size=$(du -sk tarball-cmake/ast_tester/fixtures | cut -f1)
          echo "fixture tree in tarball: ${size}K"
          test "$size" -lt 8000 \
            || { echo "cut-down fixture tree unexpectedly large" >&2; exit 1; }
          test ! -d tarball-cmake/ast_tester/fixtures/wcsconv/expected
```

and after the `ctest` line, insert:

```yaml
          ctest --test-dir tarball-cmake/_build -N | tail -2
```

- [ ] **Step 2: Assert `make check` in the tarball skips rather than fails**

After the CMake tarball step, add:

```yaml
      # The autotools suite must also adapt to the cut-down fixture set: the
      # .tap producers report a skipped plan rather than failing every row.
      - name: Run make check inside the distribution tarball
        run: |
          rm -rf tarball-am && mkdir tarball-am
          tar xzf ast-*.tar.gz -C tarball-am --strip-components=1
          cd tarball-am
          ./configure --without-starlink CC="${{ matrix.cc }}"
          make -j"$(nproc 2>/dev/null || sysctl -n hw.ncpu)"
          make check
          grep -l 'SKIP' ast_tester/*.log || true
```

Adjust the `configure` flags to match whatever the existing "Configure" step in
this workflow uses; read it before writing this step rather than assuming
`--without-starlink` exists.

- [ ] **Step 3: Add a `fulldist` check**

Add after the previous step:

```yaml
      # fulldist must produce a distinctly named archive carrying everything, so
      # that release verification has the full corpus and cannot be confused
      # with the release artifact.
      - name: Build the full distribution
        run: |
          make fulldist
          ls ast-*-full.tar.gz
          rm -rf tarball-full && mkdir tarball-full
          tar xzf ast-*-full.tar.gz -C tarball-full --strip-components=1
          grep -qx full tarball-full/ast_tester/fixtures/DIST_TIER
          test -d tarball-full/ast_tester/fixtures/wcsconv/expected
          test -d tarball-full/ast_tester/fixtures/simplify
          cmake -S tarball-full -B tarball-full/_build \
            -DCMAKE_BUILD_TYPE=Release -DAST_C_STANDARD=11
          cmake --build tarball-full/_build --parallel
          ctest --test-dir tarball-full/_build --output-on-failure
```

- [ ] **Step 4: Validate the workflow file parses**

```bash
cd /sdf/home/t/timj/WORK/ast
python3 -c "import yaml; yaml.safe_load(open('.github/workflows/autotools.yaml')); print('YAML OK')"
```

Expected: `YAML OK`.

- [ ] **Step 5: Commit**

```bash
git add .github/workflows/autotools.yaml
git commit -F - <<'MSG'
ci: check both tiers of the distribution

The tarball CMake build already existed and is what this design turns on:
ctest must work from an autoconf tarball that no longer carries every fixture.
It now also asserts the tarball is marked cutdown, that its fixture tree is
small, and that make check inside it skips rather than fails.

The fulldist job covers the other half, since a release is verified against the
full corpus and nothing else builds that archive.

Generated with AI

Co-Authored-By: SLAC AI
MSG
```

---

## Post-implementation verification

- [ ] `ctest --test-dir build-dev -j 16` in a git checkout: 2475 tests, all
      pass, no summary about a cut-down distribution.
- [ ] `./tools/test-dist-tier.sh`: all four cases `ok`.
- [ ] `./tools/test-tap-skip.sh`: all seven checks `ok`.
- [ ] `du -sh ast_tester/fixtures`: 18M (was 20M, after Task 1).
- [ ] `ast_tester/dist-fixtures.sh fixtures/DIST_MANIFEST | wc -l`: 209 files,
      5.4M.
- [ ] CI green on both `ubuntu-latest`/gcc and `macos-latest`/flang.
