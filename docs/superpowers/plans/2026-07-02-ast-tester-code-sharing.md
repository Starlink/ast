# ast_tester Shared Test Utilities Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Consolidate the duplicated helper code in `ast_tester/` (31 copies of `stopit`, ~13 `checkdump` round-trips, ~14 inline FITS-header readers, ~8 `.ast` readers, duplicated CMake tool blocks) into one shared `ast_test_util.{c,h}` unit plus a CMake `ast_add_tool()` helper.

**Architecture:** A new `ast_tester/ast_test_util.c` translation unit (public `ast.h` API only) is auto-appended to every `ast_add_test` target, following the existing `transform_oracle_util.c` precedent.
Migration is refactor-only — no test semantics change — and proceeds in mechanical waves with the full ctest suite run between waves.
A self-test executable `test_ast_test_util` grows alongside the util, mirroring `test_oracle_util`.

**Tech Stack:** C99 (C11 where native serialization is involved), CMake >= 3.24, ctest.

## Status / provenance

Written 2026-07-02 on branch `u/timj/transform-verification` (PR #63) as **future work**, capturing the duplication analysis done in review follow-up.
Execute on a fresh branch off `master` **after PR #63 merges** — do not mix this sweep into #63.
File/line references are current as of 2026-07-02 and may drift; re-verify with grep before editing.

## Global Constraints

- Refactor only: no test may change what it asserts. If a migration would change behavior, stop and flag it instead.
- C99 baseline; targets that round-trip native serialization keep/gain `C_STANDARD 11` (AST_DBL_DIG is 18 in C99 mode, 17 in C11 mode — see project CLAUDE.md).
- No new compiler warnings in changed files (check a build with `AST_ENABLE_WARNINGS=ON`).
- Full suite green after every task: `ctest --test-dir build-dev -j8` (1061+ tests at time of writing).
- Sanitizer suite for touched tests green at the end: `ctest --test-dir build-san -R '<names>'` (build-san must use homebrew clang, not conda — conda ASan hangs at startup on this machine).
- Edited test files keep their existing indentation style (mostly 3-space); the new util files use 4-space like `transform_oracle_util.c`.
- Each converted test's header comment documents differences from the Fortran original; update it only when a migration changes something an observer could notice (e.g. dropping the `fred.tmp` temp file).
- One commit per task, message style `refactor(ast_tester): ...` / `build: ...`.

## Duplication inventory (evidence base)

| Family | Copies | Where |
|---|---|---|
| `stopit(int *status, const char *text)` — byte-identical "variant A" | 19 | testcmpmap, testcmpframe, testflux, testfitstable, testframeset, testlutmap, testmapping, testplot3d, testregionmasking, testskyframe, testspecframe, testkeymap, testtable, testregions, testspecflux, testzoommap, testunitnormmap, testtime, teststc |
| `stopit` other signatures (leave alone) | ~12 | `(int,int*)`: testchannel, testchebymap, testnormmap, testpolymap, testyamlchan (non-static!); `(const char*,int*)`: testmoc, testmocchan, testhuge; `(int,double,int*)`: testrate, testratemap, testswitchmap; bespoke: testfitschan, test_bbox_values, testtrangrid |
| `checkdump` astToString/astFromString skeleton | ~7 | testflux + testspecflux (identical), testspecframe, testtime, testswitchmap, testregions, testframeset |
| `checkdump` via temp file (`fred.tmp`) + astEqual | 2 | testchebymap, testunitnormmap |
| fgets → strip CR/LF → astPutFits loop | ~14 | testgrid ≡ testplotter (identical), testtrangrid, testswitchmap, testhuge, test_bbox_values, tools |
| `readobj` (.ast via astChannel SourceFile) | ~8 | testcmpmap, testfitschan named helpers; rest inline |
| `alloc_cols`/`free_cols` | 2 | gen_transform_oracle.c, check_transform_oracle.c (despite transform_oracle_util.c existing) |
| raw `add_executable` 6-line blocks | 7 | gen_transform_oracle, check_transform_oracle, wcsconverter, ast_astequal, simplify, normalize_fixture, gen_simplify_fixtures |

Explicit **non-goals** (considered and rejected): macro-izing the `main()` prologue/epilogue (hides control flow, huge churn, zero behavior gain); renaming the non-variant-A `stopit` signatures; forcing old ad-hoc `fabs(a-b)<tol` assertions onto the shared tolerance helpers; consolidating the XML/STCS/MOC channel source/sink callbacks (each is genuinely channel-specific).

---

### Task 1: CMake `ast_add_tool()` helper

**Files:**
- Modify: `ast_tester/CMakeLists.txt` (helper functions near line 45; tool blocks near lines 113–125, 219–232, 335–352)

**Interfaces:**
- Produces: `ast_add_tool(name [C_STANDARD n] [SOURCES ...])` CMake function used by later tasks to register build-only executables.

- [ ] **Step 1: Add the function** directly below the existing `ast_add_test` function:

```cmake
# Register a build-only helper executable (no ctest entry).  Mirrors the
# compile/link setup of ast_add_test.
function(ast_add_tool name)
    cmake_parse_arguments(T "" "C_STANDARD" "SOURCES" ${ARGN})
    if(NOT T_SOURCES)
        set(T_SOURCES ${name}.c)
    endif()
    add_executable(${name} ${T_SOURCES})
    target_include_directories(${name} PRIVATE ${TEST_INCLUDES})
    target_link_libraries(${name} PRIVATE ${TEST_LIBS})
    target_compile_definitions(${name} PRIVATE HAVE_CONFIG_H)
    if(T_C_STANDARD)
        set_target_properties(${name} PROPERTIES C_STANDARD ${T_C_STANDARD})
    endif()
    ast_apply_dev_options(${name})
endfunction()
```

- [ ] **Step 2: Replace the seven raw blocks.** Each 5–6 line block collapses to one line:

```cmake
ast_add_tool(gen_transform_oracle SOURCES gen_transform_oracle.c transform_oracle_util.c C_STANDARD 11)
ast_add_tool(check_transform_oracle SOURCES check_transform_oracle.c transform_oracle_util.c C_STANDARD 11)
ast_add_tool(wcsconverter)
ast_add_tool(ast_astequal)
ast_add_tool(simplify)
ast_add_tool(normalize_fixture)
ast_add_tool(gen_simplify_fixtures)
```

Keep each call at its current location in the file (comments above the old blocks stay).

- [ ] **Step 3: Verify** — full reconfigure and build must be a no-op behaviorally:

Run: `cmake -B build-dev && cmake --build build-dev -j8 && ctest --test-dir build-dev -j8`
Expected: configures, builds, `100% tests passed`.

- [ ] **Step 4: Commit**

```bash
git add ast_tester/CMakeLists.txt
git commit -m "build: add ast_add_tool() and deduplicate tool add_executable blocks"
```

---

### Task 2: Create `ast_test_util` with `stopit`, self-test, and auto-append

**Files:**
- Create: `ast_tester/ast_test_util.h`, `ast_tester/ast_test_util.c`, `ast_tester/test_ast_test_util.c`
- Modify: `ast_tester/CMakeLists.txt` (`ast_add_test` body; register the self-test), `ast_tester/testyamlchan.c` (make its `stopit` static — pre-requisite, see Step 1)

**Interfaces:**
- Produces: `void stopit( int *status, const char *text )` — latch-first-failure reporter, linked into every test target from here on.

- [ ] **Step 1: Pre-requisite** — `testyamlchan.c:62` defines a NON-static `void stopit(int, int*)`. Once `ast_test_util.c` is linked into every test target its global `stopit` would collide at link time. Add `static` to testyamlchan's definition (and its forward declaration if present). All other duplicate `stopit`s are already `static` (file-scope) and cannot collide.

- [ ] **Step 2: Write the failing self-test** `ast_tester/test_ast_test_util.c`:

```c
/*
 * test_ast_test_util: self-tests for the shared ast_tester helpers
 * (ast_test_util.c).  New helper coverage is added here as the util grows.
 */
#include "ast_test_util.h"
#include <stdio.h>

int main( void ) {
    int fails = 0;

    /* stopit latches the first failure and ignores the rest. */
    {
        int st = 0;
        stopit( &st, "expected: first failure line" );
        if ( st != 1 ) { printf( "FAIL: stopit did not set status\n" ); fails++; }
        stopit( &st, "UNEXPECTED: must not print" );
        if ( st != 1 ) { printf( "FAIL: stopit changed set status\n" ); fails++; }
    }

    printf( fails ? "test_ast_test_util: FAILED\n"
                  : "test_ast_test_util: all passed\n" );
    return fails ? 1 : 0;
}
```

Register it in CMakeLists.txt next to `test_oracle_util`:

```cmake
ast_add_test(test_ast_test_util)
set_target_properties(test_ast_test_util PROPERTIES C_STANDARD 11)
```

- [ ] **Step 3: Verify it fails** — `cmake -B build-dev && cmake --build build-dev --target test_ast_test_util` must fail with `'ast_test_util.h' file not found`.

- [ ] **Step 4: Write the util.** `ast_tester/ast_test_util.h`:

```c
/* Shared helpers for the ast_tester test programs.  Public ast.h API only.
   This unit is compiled into every ast_add_test target; tests opt in by
   including this header instead of defining their own copy. */
#ifndef AST_TEST_UTIL_H
#define AST_TEST_UTIL_H

#include "ast.h"

/* Latch the first test failure: if *status is clear, print text and set
   *status to 1; otherwise do nothing, so the first failure is the one
   reported. */
void stopit( int *status, const char *text );

#endif
```

`ast_tester/ast_test_util.c`:

```c
/*
 * ast_test_util.c
 *
 * Shared helpers for the ast_tester test programs.  Kept dependency-free
 * (public ast.h API plus libc) so it can be linked into every test target.
 */
#include "ast_test_util.h"
#include <stdio.h>

void stopit( int *status, const char *text ) {
    if ( *status != 0 ) return;
    *status = 1;
    printf( "%s\n", text );
}
```

- [ ] **Step 5: Auto-append the util to every test target.** In `ast_add_test`, after the `if(NOT T_SOURCES)` default, add:

```cmake
    list(APPEND T_SOURCES ast_test_util.c)
```

Do not list `ast_test_util.c` manually in any `SOURCES` argument.

- [ ] **Step 6: Verify green** — `cmake -B build-dev && cmake --build build-dev -j8 && ctest --test-dir build-dev -j8`.
Expected: `test_ast_test_util` present and passing; total count grows by one; everything else unchanged.

- [ ] **Step 7: Commit**

```bash
git add ast_tester/ast_test_util.h ast_tester/ast_test_util.c \
        ast_tester/test_ast_test_util.c ast_tester/CMakeLists.txt \
        ast_tester/testyamlchan.c
git commit -m "refactor(ast_tester): introduce shared ast_test_util with stopit"
```

---

### Task 3: Migrate the 19 identical `stopit` copies

**Files:**
- Modify (delete local `static void stopit(...)` + any forward declaration; add `#include "ast_test_util.h"` after `ast.h`): `testcmpmap.c`, `testcmpframe.c`, `testflux.c`, `testfitstable.c`, `testframeset.c`, `testlutmap.c`, `testmapping.c`, `testplot3d.c`, `testregionmasking.c`, `testskyframe.c`, `testspecframe.c`, `testkeymap.c`, `testtable.c`, `testregions.c`, `testspecflux.c`, `testzoommap.c`, `testunitnormmap.c`, `testtime.c`, `teststc.c`

**Interfaces:**
- Consumes: `stopit` from Task 2. Call sites are untouched — the shared signature matches variant A exactly.

- [ ] **Step 1: Verify each candidate is byte-equivalent** before deleting (guard against drift since this plan was written):

Run: `grep -A5 'static void stopit' ast_tester/test<name>.c`
Expected body: `if( *status != 0 ) return; *status = 1; printf( "%s\n", text );` (whitespace variations fine; anything else → leave that file alone and note it).

- [ ] **Step 2: Per-file edit** (worked example, `testzoommap.c`; identical operation for every file in the list):

Remove:

```c
static void stopit( int *status, const char *text ) {
   if( *status != 0 ) return;
   *status = 1;
   printf( "%s\n", text );
}
```

Add below the existing `#include "ast.h"`:

```c
#include "ast_test_util.h"
```

Per-file checklist:
- [ ] testcmpmap.c
- [ ] testcmpframe.c
- [ ] testflux.c
- [ ] testfitstable.c
- [ ] testframeset.c
- [ ] testlutmap.c
- [ ] testmapping.c
- [ ] testplot3d.c
- [ ] testregionmasking.c
- [ ] testskyframe.c
- [ ] testspecframe.c
- [ ] testkeymap.c
- [ ] testtable.c
- [ ] testregions.c
- [ ] testspecflux.c
- [ ] testzoommap.c
- [ ] testunitnormmap.c
- [ ] testtime.c
- [ ] teststc.c

- [ ] **Step 3: Verify** — `cmake --build build-dev -j8 && ctest --test-dir build-dev -j8`: `100% tests passed`.

- [ ] **Step 4: Commit**

```bash
git add ast_tester/*.c
git commit -m "refactor(ast_tester): use shared stopit in the 19 identical copies"
```

---

### Task 4: Move `within_tol`/`within_tol_wrap` into the shared util

**Files:**
- Modify: `ast_tester/ast_test_util.{h,c}`, `ast_tester/transform_oracle.h`, `ast_tester/transform_oracle_util.c`, `ast_tester/check_transform_oracle.c`, `ast_tester/test_oracle_util.c`, `ast_tester/test_ast_test_util.c`, `ast_tester/CMakeLists.txt`

**Interfaces:**
- Produces: `int within_tol( double got, double ref, double rtol, double atol )` and `int within_tol_wrap( ... )` — same semantics as the current `oracle_within_tol`/`oracle_within_tol_wrap` (AST__BAD matches only AST__BAD; NaN never matches; wrap accepts a 2*pi offset).

- [ ] **Step 1: Add failing coverage** to `test_ast_test_util.c` by moving the `oracle_within_tol`/`oracle_within_tol_wrap` assertions out of `test_oracle_util.c` verbatim, renamed to `within_tol`/`within_tol_wrap`. Build must fail: undefined `within_tol`.

- [ ] **Step 2: Move the implementations** from `transform_oracle_util.c` into `ast_test_util.c` (rename `oracle_within_tol` → `within_tol`, `oracle_within_tol_wrap` → `within_tol_wrap`; bodies unchanged; declarations move from `transform_oracle.h` to `ast_test_util.h`; `transform_oracle.h` gains `#include "ast_test_util.h"`).
Update every call site: `grep -rn 'oracle_within_tol' ast_tester/` (expected: `check_transform_oracle.c` axis_match/compare_rtrip, `test_oracle_util.c`).
The two oracle **tools** are not `ast_add_test` targets, so append `ast_test_util.c` to their SOURCES:

```cmake
ast_add_tool(gen_transform_oracle SOURCES gen_transform_oracle.c transform_oracle_util.c ast_test_util.c C_STANDARD 11)
ast_add_tool(check_transform_oracle SOURCES check_transform_oracle.c transform_oracle_util.c ast_test_util.c C_STANDARD 11)
```

- [ ] **Step 3: Verify** — build all, run full suite, then regenerate nothing: the oracle files must NOT need regeneration (comparison code moved, not changed). Confirm: `ctest --test-dir build-dev -R 'oracle' --output-on-failure` passes against the checked-in oracles.

- [ ] **Step 4: Commit**

```bash
git add ast_tester/
git commit -m "refactor(ast_tester): move tolerance comparators into ast_test_util"
```

---

### Task 5: Deduplicate `alloc_cols`/`free_cols` between the oracle tools

**Files:**
- Modify: `ast_tester/transform_oracle.h`, `ast_tester/transform_oracle_util.c`, `ast_tester/gen_transform_oracle.c`, `ast_tester/check_transform_oracle.c`

**Interfaces:**
- Produces: `double **oracle_alloc_cols( int ncol, int nrow )` and `void oracle_free_cols( double **c, int ncol )` in transform_oracle_util (these are oracle-plumbing, not general test helpers — they stay in the oracle util).

- [ ] **Step 1:** Move one copy of `alloc_cols`/`free_cols` (they are identical; verify with diff) into `transform_oracle_util.c` as `oracle_alloc_cols`/`oracle_free_cols`, declare in `transform_oracle.h`, delete both local copies, update call sites in both tools (`grep -n 'alloc_cols\|free_cols' ast_tester/gen_transform_oracle.c ast_tester/check_transform_oracle.c`).

- [ ] **Step 2: Verify** — build; run `./build-dev/ast_tester/check_transform_oracle --selftest` (expect `selftest: ok`) and the full suite.

- [ ] **Step 3: Commit**

```bash
git add ast_tester/
git commit -m "refactor(ast_tester): share alloc_cols/free_cols via transform_oracle_util"
```

---

### Task 6: Shared `roundtrip()` and the string-based `checkdump` wave

**Files:**
- Modify: `ast_tester/ast_test_util.{h,c}`, `ast_tester/test_ast_test_util.c`, then `testflux.c`, `testspecflux.c`, `testspecframe.c`, `testtime.c`, `testswitchmap.c`, `testregions.c`, `testframeset.c`

**Interfaces:**
- Produces: `AstObject *roundtrip( AstObject *obj, int *status )` — astToString/astFromString round-trip returning the reloaded object (NULL + `stopit` on failure); caller annuls the result.

- [ ] **Step 1: Failing self-test** (add to `test_ast_test_util.c`; needs `astWatch` — add the standard status prologue to its `main` if not present):

```c
    /* roundtrip reproduces an object through astToString/astFromString. */
    {
        int st = 0;
        AstZoomMap *zm = astZoomMap( 2, 5.0, " " );
        AstObject *back = roundtrip( (AstObject *) zm, &st );
        if ( st || !back || !astEqual( zm, back ) ) {
            printf( "FAIL: roundtrip did not reproduce a ZoomMap\n" ); fails++;
        }
        if ( back ) back = astAnnul( back );
        zm = astAnnul( zm );
    }
```

Verify it fails to build (undefined `roundtrip`).

- [ ] **Step 2: Implement** in `ast_test_util.c` (declaration in the header with the comment from Interfaces):

```c
AstObject *roundtrip( AstObject *obj, int *status ) {
    char *pickle = astToString( obj );
    if ( !pickle ) {
        stopit( status, "roundtrip: astToString returned NULL" );
        return NULL;
    }
    AstObject *result = astFromString( pickle );
    pickle = astFree( pickle );
    if ( !result ) stopit( status, "roundtrip: astFromString returned NULL" );
    return result;
}
```

Verify the self-test passes.

- [ ] **Step 3: Migrate each file's local `checkdump`** to call `roundtrip` for the skeleton, keeping the type-specific comparison verbatim. Worked example — `testflux.c` (and `testspecflux.c`, which is identical):

Before:

```c
static void checkdump( AstObject *obj, const char *text, int *status ) {
   char *pickle;
   AstObject *result;
   if( *status != 0 ) return;
   pickle = astToString( obj );
   if( !pickle ) { stopit( status, text ); return; }
   result = astFromString( pickle );
   pickle = astFree( pickle );
   if( !result ) { stopit( status, text ); return; }
   if( astGetD( obj, "specval" ) != astGetD( result, "specval" ) ) {
      stopit( status, text );
   }
}
```

After:

```c
static void checkdump( AstObject *obj, const char *text, int *status ) {
   AstObject *result;
   if( *status != 0 ) return;
   result = roundtrip( obj, status );
   if( !result ) { stopit( status, text ); return; }
   if( astGetD( obj, "specval" ) != astGetD( result, "specval" ) ) {
      stopit( status, text );
   }
   result = astAnnul( result );
}
```

Note the added `astAnnul` — most current copies leak `result` inside `astBegin/astEnd` scopes (harmless but untidy); annul explicitly.
For `testspecframe.c`, `testtime.c`, `testswitchmap.c`, `testregions.c`, `testframeset.c`: apply the same replacement of the pickle block only; every line after `if( !result ) ...` stays byte-identical.
`testframeset.c`'s variant returns the restored object — have it `return roundtrip( obj, status );`-style but keep its existing signature and semantics.

Per-file checklist:
- [ ] testflux.c
- [ ] testspecflux.c
- [ ] testspecframe.c
- [ ] testtime.c
- [ ] testswitchmap.c
- [ ] testregions.c
- [ ] testframeset.c

- [ ] **Step 4: Verify** — full suite green.

- [ ] **Step 5: Commit**

```bash
git add ast_tester/
git commit -m "refactor(ast_tester): share the checkdump round-trip skeleton"
```

---

### Task 7: Convert the file-based checkdumps (`fred.tmp`) to `roundtrip()`

**Files:**
- Modify: `ast_tester/testchebymap.c`, `ast_tester/testunitnormmap.c`

**Interfaces:**
- Consumes: `roundtrip` from Task 6.

These two write the object to a hardcoded `fred.tmp` via `astChannel SinkFile=`, read it back with `SourceFile=`, and compare with `astEqual`.
That diverges from the project's stated conversion pattern (CLAUDE.md: astToString/astFromString replaces file/COMMON channels) and litters temp files in whatever CWD the binary runs from (see the stray `fred.txt` in the repo root).

- [ ] **Step 1:** Replace each file's `checkdump` body with:

```c
static void checkdump( AstObject *obj, int *status ) {
   AstObject *result;
   if( *status != 0 ) return;
   result = roundtrip( obj, status );
   if( !result ) return;
   if( !astEqual( obj, result ) ) stopit( status, "checkdump: astEqual failed" );
   result = astAnnul( result );
}
```

(Keep each file's actual signature and failure text if they differ; only the mechanism changes.)

- [ ] **Step 2:** Update each file's header comment: the round-trip now uses astToString/astFromString instead of a temporary channel file, matching the standard conversion pattern.

- [ ] **Step 3:** Check `C_STANDARD` — native-form serialization now flows through these tests; confirm their targets get C11 (add `set_target_properties(<name> PROPERTIES C_STANDARD 11)` in CMakeLists.txt if not already).

- [ ] **Step 4: Verify** — full suite green; confirm no `fred.tmp` appears in the build dir after `ctest --test-dir build-dev -R 'testchebymap|testunitnormmap'`.

- [ ] **Step 5: Commit**

```bash
git add ast_tester/
git commit -m "refactor(ast_tester): checkdump via roundtrip, dropping fred.tmp channel files"
```

---

### Task 8: Shared file readers (`read_fits_header`, `read_ast_object`)

**Files:**
- Modify: `ast_tester/ast_test_util.{h,c}`, `ast_tester/test_ast_test_util.c`, `ast_tester/CMakeLists.txt`, then `testgrid.c`, `testplotter.c`, `testtrangrid.c`, `testswitchmap.c` (FITS loops) and `testcmpmap.c`, `testfitschan.c` (`readobj`)

**Interfaces:**
- Produces:
  - `AstFitsChan *read_fits_header( const char *path, int *status )` — reads a text FITS header (one card per line, CR/LF stripped) into a new FitsChan, rewinds it (`astClear Card`), NULL + `stopit` if the file cannot be opened.
  - `AstObject *read_ast_object( const char *path, int *status )` — reads the first object from a native AST dump via `astChannel SourceFile=`, NULL + `stopit` on failure.

- [ ] **Step 1: Fixture copies for the self-test.** Next to the `test_ast_test_util` registration add (repeating a `configure_file` of the same source/dest used elsewhere is fine):

```cmake
configure_file(sip.head "${CMAKE_CURRENT_BINARY_DIR}/sip.head" COPYONLY)
configure_file(splittest1.ast "${CMAKE_CURRENT_BINARY_DIR}/splittest1.ast" COPYONLY)
```

- [ ] **Step 2: Failing self-test** additions:

```c
    /* read_fits_header loads a header and astRead yields a FrameSet. */
    {
        int st = 0;
        AstFitsChan *fc = read_fits_header( "sip.head", &st );
        AstObject *obj = fc ? astRead( fc ) : NULL;
        if ( st || !obj || !astIsAFrameSet( obj ) ) {
            printf( "FAIL: read_fits_header/sip.head\n" ); fails++;
        }
        if ( obj ) obj = astAnnul( obj );
        if ( fc ) fc = astAnnul( fc );
    }

    /* read_ast_object loads a native dump. */
    {
        int st = 0;
        AstObject *obj = read_ast_object( "splittest1.ast", &st );
        if ( st || !obj ) { printf( "FAIL: read_ast_object\n" ); fails++; }
        if ( obj ) obj = astAnnul( obj );
    }
```

Verify build failure (undefined functions).

- [ ] **Step 3: Implement** in `ast_test_util.c`:

```c
AstFitsChan *read_fits_header( const char *path, int *status ) {
    FILE *fp = fopen( path, "r" );
    if ( !fp ) {
        stopit( status, "read_fits_header: cannot open header file" );
        return NULL;
    }
    AstFitsChan *fc = astFitsChan( NULL, NULL, " " );
    char line[ 256 ];
    while ( fgets( line, (int) sizeof line, fp ) ) {
        size_t n = strlen( line );
        while ( n > 0 && ( line[ n - 1 ] == '\n' || line[ n - 1 ] == '\r' ) )
            line[ --n ] = '\0';
        astPutFits( fc, line, 0 );
    }
    fclose( fp );
    astClear( fc, "Card" );
    return fc;
}

AstObject *read_ast_object( const char *path, int *status ) {
    AstChannel *ch = astChannel( NULL, NULL, "SourceFile=%s", path );
    AstObject *obj = astRead( ch );
    ch = astAnnul( ch );
    if ( !obj ) stopit( status, "read_ast_object: astRead returned NULL" );
    return obj;
}
```

(`#include <string.h>` joins the util's includes.)
Verify self-test passes.

- [ ] **Step 4: Migrate call sites.** For each listed file, replace the inline fgets/astPutFits loop (or local `readobj`) with a call to the shared helper; delete the local helper.
Before deleting, diff the local copy against the shared one — if a local copy does anything extra (e.g. `testfitschan.c`'s `readobj` builds an options buffer), keep the extra behavior local and only replace the common core, or skip the file and note it.
Per-file checklist:
- [ ] testgrid.c
- [ ] testplotter.c
- [ ] testtrangrid.c
- [ ] testswitchmap.c
- [ ] testcmpmap.c
- [ ] testfitschan.c

- [ ] **Step 5: Verify** — full suite green.

- [ ] **Step 6: Commit**

```bash
git add ast_tester/
git commit -m "refactor(ast_tester): shared FITS-header and .ast dump readers"
```

---

### Task 9: Documentation and final verification

**Files:**
- Modify: `PLAN.md` (test inventory: add `test_ast_test_util`; note the shared util under a short "Shared test infrastructure" bullet), `docs/superpowers/plans/2026-07-02-ast-tester-code-sharing.md` (tick remaining checkboxes)

- [ ] **Step 1:** Update PLAN.md as above.

- [ ] **Step 2: Full verification battery:**

Run: `cmake --build build-dev -j8 && ctest --test-dir build-dev -j8`
Expected: `100% tests passed`.

Run: `cmake --build build-san -j8 && ctest --test-dir build-san -j8`
Expected: `100% tests passed` (build-san = homebrew clang + `AST_ENABLE_SANITIZERS=ON`).

Run: `git grep -c 'static void stopit' -- ast_tester/ | wc -l`
Expected: only the non-variant-A files remain (~12).

- [ ] **Step 3: Commit**

```bash
git add PLAN.md docs/superpowers/plans/2026-07-02-ast-tester-code-sharing.md
git commit -m "docs: record shared ast_tester test infrastructure in PLAN.md"
```
