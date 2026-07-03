# Transform Output Oracle Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a tolerance-based, architecture-independent regression oracle that records the current numerical output of AST transforms over the `simplify_fixtures` and FITS-header corpora, with a scan-only generator and a file-driven `ctest` checker.

**Architecture:** Two standalone C programs share a small utility translation unit. `gen_transform_oracle` scans the fixture tree, samples input points, transforms them, and writes two committed text oracle files. `check_transform_oracle` reads an oracle file as its sole manifest, reloads each named fixture, re-transforms the recorded inputs, and compares outputs within tolerance (golden checks A/B, round-trip, and `.map`/`.simp` equivalence check C).

**Tech Stack:** C (public `ast.h` API only), POSIX `dirent.h` for scanning, CMake/CTest, the existing `ast_tester` whole-archive link pattern.

## Global Constraints

- Use the **public C API** (`ast.h`) only — no internal headers, no `#define astCLASS`.
- All new code must add **no new compiler warnings** and pass the ASan/UBSan build (`-DAST_ENABLE_WARNINGS=ON -DAST_ENABLE_SANITIZERS=ON`).
- New test/program targets use the **whole-archive satellite-library linkage** via the existing `TEST_LIBS` list / `ast_add_test` — without it the executables segfault on null function pointers.
- All recorded numbers are written with **`printf("%.17g", x)`** (17 significant digits round-trips an IEEE-754 double exactly); oracle program targets are built with `C_STANDARD 11`.
- `AST__BAD` is serialized as the literal token **`BAD`**; `BAD` matches only `BAD`.
- Two oracle files: **`ast_tester/simplify_fixtures.oracle`** (native dumps) and **`ast_tester/headers.oracle`** (FITS headers).
- Comparison is **tolerance-based**: `|got - ref| <= atol + rtol*|ref|`. Defaults: golden/round-trip `rtol=atol=1e-12`; equivalence `equiv_rtol=equiv_atol=1e-9`. Final values are tuned in Task 9.
- These are new files under `ast_tester/`, not `src/`, so the "update the prologue history" rule does **not** apply; but each new `.c` gets a header comment in the style of `simplify.c`/`wcsconverter.c`.
- Update `PLAN.md` as work progresses (project convention).

---

## File Structure

- `ast_tester/transform_oracle.h` — shared declarations and default-tolerance constants.
- `ast_tester/transform_oracle_util.c` — shared implementations: tolerance compare, Halton sampling, fixture loader, double format/parse. Compiled into every program and unit test below.
- `ast_tester/gen_transform_oracle.c` — generator (scan → sample → transform → write). Not a `ctest` test.
- `ast_tester/check_transform_oracle.c` — checker (parse oracle → reload → transform → compare). Wired into `ctest`.
- `ast_tester/transform_oracle_rtrip_overrides.txt` — human-curated per-fixture round-trip tolerance overrides (e.g. PolyMap iterative inverses).
- `ast_tester/simplify_fixtures.oracle`, `ast_tester/headers.oracle` — generated, committed reference data.
- `ast_tester/test_oracle_util.c` — unit tests for the pure utility functions.
- `ast_tester/CMakeLists.txt` — new targets and tests (modify).
- `PLAN.md` — status (modify).

---

### Task 1: Shared header and tolerance-compare primitive

**Files:**
- Create: `ast_tester/transform_oracle.h`
- Create: `ast_tester/transform_oracle_util.c`
- Create: `ast_tester/test_oracle_util.c`
- Modify: `ast_tester/CMakeLists.txt`

**Interfaces:**
- Produces:
  - `#define ORACLE_BAD_TOKEN "BAD"`
  - default tolerance macros `ORACLE_DEF_RTOL`, `ORACLE_DEF_ATOL`, `ORACLE_DEF_EQUIV_RTOL`, `ORACLE_DEF_EQUIV_ATOL` (all `double` literals)
  - `int oracle_within_tol(double got, double ref, double rtol, double atol);` — returns 1 if within tolerance. `AST__BAD` matches only `AST__BAD`; any NaN never matches.

- [ ] **Step 1: Write the failing test**

Create `ast_tester/test_oracle_util.c`:

```c
/* Unit tests for transform_oracle_util.c pure helpers. */
#include "ast.h"
#include "transform_oracle.h"
#include <stdio.h>
#include <math.h>

static int failures = 0;
#define CHECK(cond) do { if(!(cond)) { \
    fprintf(stderr, "FAIL line %d: %s\n", __LINE__, #cond); failures++; } } while(0)

static void test_within_tol(void) {
    /* exact equality */
    CHECK(oracle_within_tol(1.0, 1.0, 0.0, 0.0) == 1);
    /* a few ULP apart, within default rtol */
    CHECK(oracle_within_tol(1.0 + 4e-16, 1.0, ORACLE_DEF_RTOL, ORACLE_DEF_ATOL) == 1);
    /* gross difference rejected */
    CHECK(oracle_within_tol(1.0001, 1.0, ORACLE_DEF_RTOL, ORACLE_DEF_ATOL) == 0);
    /* absolute floor near zero */
    CHECK(oracle_within_tol(5e-13, 0.0, ORACLE_DEF_RTOL, ORACLE_DEF_ATOL) == 1);
    /* BAD matches BAD only */
    CHECK(oracle_within_tol(AST__BAD, AST__BAD, ORACLE_DEF_RTOL, ORACLE_DEF_ATOL) == 1);
    CHECK(oracle_within_tol(AST__BAD, 1.0, ORACLE_DEF_RTOL, ORACLE_DEF_ATOL) == 0);
    CHECK(oracle_within_tol(1.0, AST__BAD, ORACLE_DEF_RTOL, ORACLE_DEF_ATOL) == 0);
    /* NaN never matches */
    CHECK(oracle_within_tol(NAN, 1.0, 1e9, 1e9) == 0);
}

int main(void) {
    test_within_tol();
    if (failures) { fprintf(stderr, "%d failure(s)\n", failures); return 1; }
    printf("test_oracle_util: all passed\n");
    return 0;
}
```

- [ ] **Step 2: Create the header with declarations and the util file with a deliberately wrong stub**

Create `ast_tester/transform_oracle.h`:

```c
/* Shared declarations for the transform-output oracle programs. */
#ifndef TRANSFORM_ORACLE_H
#define TRANSFORM_ORACLE_H

#include "ast.h"

#define ORACLE_BAD_TOKEN "BAD"

/* Default comparison tolerances (tuned in the tolerance-tuning task). */
#define ORACLE_DEF_RTOL        1e-12
#define ORACLE_DEF_ATOL        1e-12
#define ORACLE_DEF_EQUIV_RTOL  1e-9
#define ORACLE_DEF_EQUIV_ATOL  1e-9

/* Return 1 if got and ref agree within |got-ref| <= atol + rtol*|ref|.
   AST__BAD matches only AST__BAD; any NaN never matches. */
int oracle_within_tol( double got, double ref, double rtol, double atol );

#endif
```

Create `ast_tester/transform_oracle_util.c` with an intentionally wrong stub so the test fails first:

```c
/*
 * transform_oracle_util.c
 *
 * Shared helpers for gen_transform_oracle and check_transform_oracle:
 * tolerance comparison, low-discrepancy sampling, fixture loading, and
 * full-precision double formatting/parsing.  Public AST API only.
 */
#include "transform_oracle.h"
#include <math.h>

int oracle_within_tol( double got, double ref, double rtol, double atol ) {
    (void) rtol; (void) atol;
    return got == ref;   /* WRONG ON PURPOSE: tightened in step 4 */
}
```

- [ ] **Step 3: Register the unit-test target and build it, expecting failure**

Add to `ast_tester/CMakeLists.txt` (after the existing public-API test block, near the other `ast_add_test` calls):

```cmake
# --- Transform output oracle ---
ast_add_test(test_oracle_util SOURCES test_oracle_util.c transform_oracle_util.c)
set_target_properties(test_oracle_util PROPERTIES C_STANDARD 11)
```

Run:
```bash
cmake --build build --target test_oracle_util && ctest --test-dir build -R test_oracle_util --output-on-failure
```
Expected: build succeeds, test FAILS (the few-ULP and absolute-floor cases fail under exact `==`).

- [ ] **Step 4: Implement the real comparison**

Replace the stub body in `transform_oracle_util.c`:

```c
int oracle_within_tol( double got, double ref, double rtol, double atol ) {
    int got_bad = ( got == AST__BAD );
    int ref_bad = ( ref == AST__BAD );
    if ( got_bad || ref_bad ) return ( got_bad && ref_bad );
    if ( isnan( got ) || isnan( ref ) ) return 0;
    return fabs( got - ref ) <= atol + rtol * fabs( ref );
}
```

- [ ] **Step 5: Run the test, expect pass**

Run:
```bash
cmake --build build --target test_oracle_util && ctest --test-dir build -R test_oracle_util --output-on-failure
```
Expected: PASS (`test_oracle_util: all passed`).

- [ ] **Step 6: Commit**

```bash
git add ast_tester/transform_oracle.h ast_tester/transform_oracle_util.c \
        ast_tester/test_oracle_util.c ast_tester/CMakeLists.txt
git commit -m "test: add transform-oracle tolerance compare primitive"
```

---

### Task 2: Low-discrepancy sampling

**Files:**
- Modify: `ast_tester/transform_oracle.h`
- Modify: `ast_tester/transform_oracle_util.c`
- Modify: `ast_tester/test_oracle_util.c`

**Interfaces:**
- Consumes: nothing new.
- Produces:
  - `double oracle_halton( unsigned index, unsigned base );` — radical inverse, returns a value in `[0,1)`.
  - `int oracle_sample_axis_count(void);` — number of sampled points per fixture (returns the tunable constant).
  - `void oracle_sample_points( int naxis, const double *lo, const double *hi, int npoint, double **out );` — fills `out[axis][point]` with `npoint` samples per axis: a Halton spread (distinct prime base per axis) mapped into `[lo,hi]`, with the first few points forced to edge values (all-lo, all-hi, all-zero clamped into range, midpoint).

- [ ] **Step 1: Write the failing test**

Append to `test_oracle_util.c` (and call from `main`):

```c
static void test_halton(void) {
    /* Known radical-inverse base-2 values: 1->1/2, 2->1/4, 3->3/4. */
    CHECK(fabs(oracle_halton(1, 2) - 0.5)  < 1e-15);
    CHECK(fabs(oracle_halton(2, 2) - 0.25) < 1e-15);
    CHECK(fabs(oracle_halton(3, 2) - 0.75) < 1e-15);
    /* base-3: 1->1/3, 2->2/3. */
    CHECK(fabs(oracle_halton(1, 3) - 1.0/3.0) < 1e-15);
    CHECK(fabs(oracle_halton(2, 3) - 2.0/3.0) < 1e-15);
}

static void test_sample_points(void) {
    int n = oracle_sample_axis_count();
    CHECK(n >= 8);
    double lo[2] = {-10.0, -10.0}, hi[2] = {10.0, 10.0};
    double *col[2];
    col[0] = malloc(sizeof(double) * (size_t)n);
    col[1] = malloc(sizeof(double) * (size_t)n);
    oracle_sample_points(2, lo, hi, n, col);
    /* every sample inside the closed range */
    for (int a = 0; a < 2; a++)
        for (int i = 0; i < n; i++)
            CHECK(col[a][i] >= lo[a] - 1e-9 && col[a][i] <= hi[a] + 1e-9);
    /* edge points present: an all-lo and an all-hi row somewhere */
    int saw_lo = 0, saw_hi = 0;
    for (int i = 0; i < n; i++) {
        if (fabs(col[0][i]-lo[0])<1e-12 && fabs(col[1][i]-lo[1])<1e-12) saw_lo = 1;
        if (fabs(col[0][i]-hi[0])<1e-12 && fabs(col[1][i]-hi[1])<1e-12) saw_hi = 1;
    }
    CHECK(saw_lo && saw_hi);
    free(col[0]); free(col[1]);
}
```

Add `#include <stdlib.h>` at the top of `test_oracle_util.c`, and add `test_halton(); test_sample_points();` to `main` before the failure check.

- [ ] **Step 2: Run, expect failure**

Run:
```bash
cmake --build build --target test_oracle_util
```
Expected: FAIL to build/link — `oracle_halton`, `oracle_sample_axis_count`, `oracle_sample_points` undefined.

- [ ] **Step 3: Implement sampling**

Add declarations to `transform_oracle.h` (before `#endif`):

```c
double oracle_halton( unsigned index, unsigned base );
int    oracle_sample_axis_count( void );
void   oracle_sample_points( int naxis, const double *lo, const double *hi,
                             int npoint, double **out );
```

Add to `transform_oracle_util.c` (add `#include <stddef.h>`):

```c
/* Number of sampled points per fixture (tunable). */
#define ORACLE_NPOINT 24

/* Distinct prime bases per axis; extended cyclically beyond the table. */
static const unsigned oracle_primes[] = {
    2u, 3u, 5u, 7u, 11u, 13u, 17u, 19u, 23u, 29u, 31u, 37u
};

double oracle_halton( unsigned index, unsigned base ) {
    double f = 1.0, r = 0.0;
    while ( index > 0u ) {
        f /= (double) base;
        r += f * (double) ( index % base );
        index /= base;
    }
    return r;
}

int oracle_sample_axis_count( void ) { return ORACLE_NPOINT; }

static double clamp( double v, double lo, double hi ) {
    return v < lo ? lo : ( v > hi ? hi : v );
}

void oracle_sample_points( int naxis, const double *lo, const double *hi,
                           int npoint, double **out ) {
    /* First four rows are deterministic edge/structure points. */
    for ( int a = 0; a < naxis; a++ ) {
        if ( npoint > 0 ) out[a][0] = lo[a];                 /* all-lo  */
        if ( npoint > 1 ) out[a][1] = hi[a];                 /* all-hi  */
        if ( npoint > 2 ) out[a][2] = clamp(0.0, lo[a], hi[a]); /* zero */
        if ( npoint > 3 ) out[a][3] = 0.5 * ( lo[a] + hi[a] );  /* mid  */
    }
    /* Remaining rows: Halton spread, distinct prime base per axis. */
    int nedge = npoint < 4 ? npoint : 4;
    for ( int a = 0; a < naxis; a++ ) {
        unsigned base = oracle_primes[ (size_t) a % ( sizeof(oracle_primes)
                                       / sizeof(oracle_primes[0]) ) ];
        for ( int i = nedge; i < npoint; i++ ) {
            double u = oracle_halton( (unsigned)( i - nedge + 1 ), base );
            out[a][i] = lo[a] + u * ( hi[a] - lo[a] );
        }
    }
}
```

- [ ] **Step 4: Run, expect pass**

Run:
```bash
cmake --build build --target test_oracle_util && ctest --test-dir build -R test_oracle_util --output-on-failure
```
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add ast_tester/transform_oracle.h ast_tester/transform_oracle_util.c ast_tester/test_oracle_util.c
git commit -m "test: add Halton low-discrepancy sampling for the oracle"
```

---

### Task 3: Fixture loader and double format/parse helpers

**Files:**
- Modify: `ast_tester/transform_oracle.h`
- Modify: `ast_tester/transform_oracle_util.c`
- Modify: `ast_tester/test_oracle_util.c`

**Interfaces:**
- Produces:
  - `AstMapping *oracle_load_mapping( const char *root, const char *relpath );` — loads `<root>/<relpath>`. `.head` → `FitsChan` → `FrameSet` → `astGetMapping(fs, AST__BASE, AST__CURRENT)`; everything else → `Channel`/`astRead`. Returns NULL (after clearing status) if nothing usable could be read. Caller `astAnnul`s the result.
  - `void oracle_format_double( char *buf, size_t buflen, double v );` — writes `%.17g`, or `ORACLE_BAD_TOKEN` for `AST__BAD`.
  - `double oracle_parse_double( const char *tok, int *ok );` — parses a token; `BAD` → `AST__BAD`; sets `*ok` 0 on malformed input.

- [ ] **Step 1: Write the failing test**

Append to `test_oracle_util.c` and call from `main`:

```c
static void test_format_parse(void) {
    char buf[64];
    int ok = 0;
    oracle_format_double(buf, sizeof buf, AST__BAD);
    CHECK(strcmp(buf, "BAD") == 0);
    double back = oracle_parse_double("BAD", &ok);
    CHECK(ok && back == AST__BAD);
    /* round-trip a representative double exactly via %.17g */
    double v = 1.0/3.0;
    oracle_format_double(buf, sizeof buf, v);
    back = oracle_parse_double(buf, &ok);
    CHECK(ok && back == v);
    oracle_parse_double("not_a_number", &ok);
    CHECK(ok == 0);
}

static void test_load_mapping(void) {
    const char *root = getenv("ORACLE_TEST_ROOT");
    if (!root) { fprintf(stderr, "skip test_load_mapping (no ORACLE_TEST_ROOT)\n"); return; }
    AstMapping *m1 = oracle_load_mapping(root, "simplify_fixtures/matrix_diagonal_to_zoom.map");
    CHECK(m1 != NULL);
    if (m1) { CHECK(astGetI(m1, "Nin") == 2); m1 = astAnnul(m1); }
    AstMapping *m2 = oracle_load_mapping(root, "cobe.head");
    CHECK(m2 != NULL);
    if (m2) { CHECK(astGetI(m2, "Nin") == 2); m2 = astAnnul(m2); }
}
```

Add `#include <string.h>` to the test file; call `test_format_parse(); test_load_mapping();` in `main`. Register the test root by adding `set_tests_properties(test_oracle_util PROPERTIES ENVIRONMENT "ORACLE_TEST_ROOT=${CMAKE_CURRENT_SOURCE_DIR}")` right after the `ast_add_test(test_oracle_util ...)` line in `CMakeLists.txt`.

- [ ] **Step 2: Run, expect failure**

Run: `cmake --build build --target test_oracle_util`
Expected: FAIL to link — the three new functions are undefined.

- [ ] **Step 3: Implement the helpers**

Add declarations to `transform_oracle.h`:

```c
AstMapping *oracle_load_mapping( const char *root, const char *relpath );
void   oracle_format_double( char *buf, size_t buflen, double v );
double oracle_parse_double( const char *tok, int *ok );
```

Add to `transform_oracle_util.c` (`#include <stdio.h>`, `#include <string.h>`, `#include <stdlib.h>`):

```c
void oracle_format_double( char *buf, size_t buflen, double v ) {
    if ( v == AST__BAD ) { snprintf( buf, buflen, "%s", ORACLE_BAD_TOKEN ); }
    else                 { snprintf( buf, buflen, "%.17g", v ); }
}

double oracle_parse_double( const char *tok, int *ok ) {
    if ( strcmp( tok, ORACLE_BAD_TOKEN ) == 0 ) { *ok = 1; return AST__BAD; }
    char *end = NULL;
    double v = strtod( tok, &end );
    *ok = ( end != tok && *end == '\0' );
    return v;
}

AstMapping *oracle_load_mapping( const char *root, const char *relpath ) {
    char path[1024];
    snprintf( path, sizeof path, "%s/%s", root, relpath );

    const char *dot = strrchr( relpath, '.' );
    int is_head = ( dot && strcmp( dot, ".head" ) == 0 );

    AstMapping *result = NULL;

    if ( is_head ) {
        AstFitsChan *fc = astFitsChan( NULL, NULL, " " );
        FILE *fp = fopen( path, "r" );
        if ( fp ) {
            char line[256];
            while ( fgets( line, (int) sizeof line, fp ) ) {
                size_t n = strlen( line );
                while ( n > 0 && ( line[n-1] == '\n' || line[n-1] == '\r' ) )
                    line[--n] = '\0';
                astPutFits( fc, line, 0 );
            }
            fclose( fp );
            astClear( fc, "Card" );
            AstObject *obj = astRead( fc );
            if ( obj ) {
                if ( astIsAFrameSet( obj ) ) {
                    result = astGetMapping( (AstFrameSet *) obj,
                                            AST__BASE, AST__CURRENT );
                }
                obj = astAnnul( obj );
            }
        }
        fc = astAnnul( fc );
    } else {
        AstChannel *chan = astChannel( NULL, NULL, "SourceFile=%s", path );
        AstObject *obj = astRead( chan );
        chan = astAnnul( chan );
        if ( obj ) {
            if ( astIsAMapping( obj ) ) result = (AstMapping *) obj;
            else obj = astAnnul( obj );
        }
    }

    if ( !astOK ) astClearStatus;   /* a bad/unsupported fixture is not fatal */
    return result;
}
```

- [ ] **Step 4: Run, expect pass**

Run: `cmake --build build --target test_oracle_util && ctest --test-dir build -R test_oracle_util --output-on-failure`
Expected: PASS (the load cases run because `ORACLE_TEST_ROOT` is set by CMake).

- [ ] **Step 5: Commit**

```bash
git add ast_tester/transform_oracle.h ast_tester/transform_oracle_util.c \
        ast_tester/test_oracle_util.c ast_tester/CMakeLists.txt
git commit -m "test: add fixture loader and double format/parse helpers"
```

---

### Task 4: Generator program

**Files:**
- Create: `ast_tester/gen_transform_oracle.c`
- Modify: `ast_tester/CMakeLists.txt`

**Interfaces:**
- Consumes: everything from `transform_oracle_util.c`.
- Produces: executable `gen_transform_oracle <root> <simplify_oracle_out> <headers_oracle_out>`. Scans `<root>/simplify_fixtures/*.map`,`*.simp` and `<root>/*.head`; writes the two oracle files. Section format:
  `[<relpath>  nin=<N> nout=<M> dir=forward|inverse]` followed by one row per point: `N` input values then `M` output values, space-separated, `%.17g`/`BAD`.

- [ ] **Step 1: Write the generator**

Create `ast_tester/gen_transform_oracle.c`:

```c
/*
 * gen_transform_oracle: scan the transform fixture corpus, sample input
 * points, transform them, and write two text oracle files used by
 * check_transform_oracle as a tolerance-based regression reference.
 *
 * Usage:
 *   gen_transform_oracle <root> <simplify_oracle_out> <headers_oracle_out>
 *
 * <root> is the ast_tester source directory.  This program is the only
 * component that scans the filesystem; the checker is driven entirely by
 * the files written here.
 */
#include "ast.h"
#include "transform_oracle.h"

#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Pick the sampling range for one axis of one fixture.  FITS-header
   mappings sample the pixel grid 1..N (approximated by a fixed generous
   pixel range here, since NAXISn is not exposed on the bare Mapping);
   native-dump mappings use a fixed symmetric range. */
static void axis_range( int is_head, double *lo, double *hi ) {
    if ( is_head ) { *lo = 1.0;     *hi = 2000.0;  }    /* pixel-ish */
    else           { *lo = -1000.0; *hi =  1000.0; }
}

/* Transform npoint points through map in the given direction and write a
   section to fp.  in_n / out_n are the arities for that direction. */
static void write_section( FILE *fp, const char *relpath, AstMapping *map,
                           int forward, int in_n, int out_n, int is_head ) {
    int np = oracle_sample_axis_count();

    double **in  = malloc( sizeof(double *) * (size_t) in_n );
    double **out = malloc( sizeof(double *) * (size_t) out_n );
    double *lo = malloc( sizeof(double) * (size_t) in_n );
    double *hi = malloc( sizeof(double) * (size_t) in_n );
    for ( int a = 0; a < in_n;  a++ ) { in[a]  = malloc(sizeof(double)*(size_t)np);
                                        axis_range(is_head,&lo[a],&hi[a]); }
    for ( int a = 0; a < out_n; a++ )   out[a] = malloc(sizeof(double)*(size_t)np);

    oracle_sample_points( in_n, lo, hi, np, in );
    astTranP( map, np, in_n, (const double **) in, forward, out_n, out );

    fprintf( fp, "[%s  nin=%d nout=%d dir=%s]\n",
             relpath, in_n, out_n, forward ? "forward" : "inverse" );
    char buf[64];
    for ( int p = 0; p < np; p++ ) {
        for ( int a = 0; a < in_n;  a++ ) {
            oracle_format_double( buf, sizeof buf, in[a][p] );
            fprintf( fp, "%s%s", a ? " " : "  ", buf );
        }
        fprintf( fp, "  " );
        for ( int a = 0; a < out_n; a++ ) {
            oracle_format_double( buf, sizeof buf, out[a][p] );
            fprintf( fp, "%s%s", a ? " " : "", buf );
        }
        fprintf( fp, "\n" );
    }
    fprintf( fp, "\n" );

    for ( int a = 0; a < in_n;  a++ ) free( in[a] );
    for ( int a = 0; a < out_n; a++ ) free( out[a] );
    free( in ); free( out ); free( lo ); free( hi );
}

/* Emit one fixture (forward if defined, else inverse-only). Returns 1 if
   a section was written. */
static int emit_fixture( FILE *fp, const char *root, const char *relpath ) {
    int wrote = 0;
    astBegin;
    AstMapping *map = oracle_load_mapping( root, relpath );
    if ( map ) {
        int nin  = astGetI( map, "Nin" );
        int nout = astGetI( map, "Nout" );
        int is_head = ( strstr( relpath, ".head" ) != NULL );
        if ( astGetI( map, "TranForward" ) ) {
            write_section( fp, relpath, map, 1, nin, nout, is_head );
            wrote = 1;
        } else if ( astGetI( map, "TranInverse" ) ) {
            write_section( fp, relpath, map, 0, nout, nin, is_head );
            wrote = 1;
        }
    }
    if ( !astOK ) astClearStatus;
    astEnd;
    return wrote;
}

/* Comparator for qsort of strings, for deterministic file order. */
static int cmp_str( const void *a, const void *b ) {
    return strcmp( *(const char *const *) a, *(const char *const *) b );
}

/* Collect entries of dir matching suffix into a sorted relpath list.
   prefix is prepended to each name (e.g. "simplify_fixtures/" or ""). */
static char **scan_dir( const char *root, const char *subdir,
                        const char *prefix, const char *suffix, int *count ) {
    char path[1024];
    snprintf( path, sizeof path, "%s%s%s", root, subdir[0] ? "/" : "", subdir );
    DIR *d = opendir( path );
    *count = 0;
    if ( !d ) return NULL;
    size_t cap = 64, n = 0;
    char **list = malloc( sizeof(char *) * cap );
    struct dirent *e;
    size_t slen = strlen( suffix );
    while ( ( e = readdir( d ) ) ) {
        size_t nlen = strlen( e->d_name );
        if ( nlen <= slen || strcmp( e->d_name + nlen - slen, suffix ) != 0 )
            continue;
        if ( n == cap ) { cap *= 2; list = realloc( list, sizeof(char *) * cap ); }
        char rel[1024];
        snprintf( rel, sizeof rel, "%s%s", prefix, e->d_name );
        list[n++] = strdup( rel );
    }
    closedir( d );
    qsort( list, n, sizeof(char *), cmp_str );
    *count = (int) n;
    return list;
}

int main( int argc, char *argv[] ) {
    int status_value = 0;
    int *status = &status_value;
    if ( argc < 4 ) {
        fprintf( stderr,
            "Usage: gen_transform_oracle <root> <simplify_out> <headers_out>\n" );
        return 1;
    }
    const char *root         = argv[1];
    const char *simplify_out = argv[2];
    const char *headers_out  = argv[3];
    astWatch( status );

    const char *hdr =
        "# transform oracle - regenerate with gen_transform_oracle\n"
        "# tol: rtol=1e-12 atol=1e-12  equiv_rtol=1e-9 equiv_atol=1e-9\n\n";

    /* Native-dump corpus: .map and .simp interleaved by stem so paired
       sections sit together in the file. */
    FILE *fs = fopen( simplify_out, "w" );
    if ( !fs ) { fprintf( stderr, "cannot write %s\n", simplify_out ); return 1; }
    fputs( hdr, fs );
    int nmap = 0;
    char **maps = scan_dir( root, "simplify_fixtures", "simplify_fixtures/",
                            ".map", &nmap );
    int nwritten = 0;
    for ( int i = 0; i < nmap; i++ ) {
        char simp[1024];
        size_t L = strlen( maps[i] );
        snprintf( simp, sizeof simp, "%.*s.simp", (int)(L - 4), maps[i] );
        nwritten += emit_fixture( fs, root, maps[i] );
        nwritten += emit_fixture( fs, root, simp );
        free( maps[i] );
    }
    free( maps );
    fclose( fs );

    /* FITS-header corpus. */
    FILE *fh = fopen( headers_out, "w" );
    if ( !fh ) { fprintf( stderr, "cannot write %s\n", headers_out ); return 1; }
    fputs( hdr, fh );
    int nhead = 0;
    char **heads = scan_dir( root, "", "", ".head", &nhead );
    for ( int i = 0; i < nhead; i++ ) {
        nwritten += emit_fixture( fh, root, heads[i] );
        free( heads[i] );
    }
    free( heads );
    fclose( fh );

    printf( "gen_transform_oracle: wrote %d sections\n", nwritten );
    return astOK ? 0 : 1;
}
```

- [ ] **Step 2: Register the generator target (built, not a test)**

Add to `ast_tester/CMakeLists.txt` after the `test_oracle_util` block:

```cmake
add_executable(gen_transform_oracle gen_transform_oracle.c transform_oracle_util.c)
target_include_directories(gen_transform_oracle PRIVATE ${TEST_INCLUDES})
target_link_libraries(gen_transform_oracle PRIVATE ${TEST_LIBS})
target_compile_definitions(gen_transform_oracle PRIVATE HAVE_CONFIG_H)
set_target_properties(gen_transform_oracle PROPERTIES C_STANDARD 11)
ast_apply_dev_options(gen_transform_oracle)
```

- [ ] **Step 3: Build and run against the real corpus into a scratch file**

Run:
```bash
cmake --build build --target gen_transform_oracle
./build/ast_tester/gen_transform_oracle ast_tester /tmp/simp.oracle /tmp/head.oracle
head -20 /tmp/simp.oracle
grep -c '^\[' /tmp/simp.oracle /tmp/head.oracle
```
Expected: a non-zero "wrote N sections" line; `/tmp/simp.oracle` begins with the `#` header then `[simplify_fixtures/...map ...]` / `.simp` sections; section counts are in the hundreds for simp and ~100+ for head.

- [ ] **Step 4: Sanity-check the format**

Run:
```bash
awk '/^\[/{print; getline; print; exit}' /tmp/simp.oracle
```
Expected: a section header line `[simplify_fixtures/...  nin=.. nout=.. dir=forward]` followed by a data row whose column count equals `nin+nout`.

- [ ] **Step 5: Commit (program only; real oracle files are committed in Task 8)**

```bash
git add ast_tester/gen_transform_oracle.c ast_tester/CMakeLists.txt
git commit -m "feat: add transform-oracle generator"
```

---

### Task 5: Oracle file parser

**Files:**
- Create: `ast_tester/check_transform_oracle.c` (parser portion only this task)
- Modify: `ast_tester/CMakeLists.txt`

**Interfaces:**
- Produces (file-local to the checker, exercised via a `--selftest` mode this task):
  - a `Section` struct: `char relpath[1024]; int nin, nout, forward; int npoint; double **in; double **out;`
  - `int oracle_read_file(const char *path, Section **secs, int *nsec);` — parse all sections. Returns 0 on success.
  - `--selftest` CLI mode that parses a tiny inline fixture string written to a temp file and asserts the parsed shape.

- [ ] **Step 1: Write the checker skeleton with parser and a selftest**

Create `ast_tester/check_transform_oracle.c`:

```c
/*
 * check_transform_oracle: read a transform oracle file (written by
 * gen_transform_oracle), reload each named fixture, re-transform the
 * recorded inputs, and compare outputs within tolerance.  The oracle
 * file is the sole manifest; this program does no corpus scanning.
 *
 * Usage:
 *   check_transform_oracle <root> <oracle_file> [<rtrip_overrides>]
 *   check_transform_oracle --selftest
 */
#include "ast.h"
#include "transform_oracle.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    char relpath[1024];
    int nin, nout, forward;
    int npoint;
    double **in;
    double **out;
} Section;

static double **alloc_cols( int ncol, int nrow ) {
    double **c = malloc( sizeof(double *) * (size_t) ncol );
    for ( int a = 0; a < ncol; a++ ) c[a] = malloc( sizeof(double) * (size_t) nrow );
    return c;
}

static int oracle_read_file( const char *path, Section **out_secs, int *out_n ) {
    FILE *fp = fopen( path, "r" );
    if ( !fp ) { fprintf( stderr, "cannot open oracle %s\n", path ); return 1; }

    size_t cap = 64, n = 0;
    Section *secs = malloc( sizeof(Section) * cap );
    char line[8192];
    Section *cur = NULL;
    int row = 0, rowcap = 0;

    while ( fgets( line, (int) sizeof line, fp ) ) {
        if ( line[0] == '#' || line[0] == '\n' ) continue;
        if ( line[0] == '[' ) {
            if ( n == cap ) { cap *= 2; secs = realloc( secs, sizeof(Section)*cap ); }
            cur = &secs[n++];
            char dir[16] = "forward";
            /* header: [<relpath>  nin=N nout=M dir=word] */
            if ( sscanf( line, "[%1023s nin=%d nout=%d dir=%15[^]]",
                         cur->relpath, &cur->nin, &cur->nout, dir ) < 3 ) {
                fprintf( stderr, "bad section header: %s", line );
                fclose( fp ); return 1;
            }
            /* relpath may have captured nothing past whitespace; re-parse
               cleanly: the %s stops at first space which is what we want. */
            cur->forward = ( strncmp( dir, "inverse", 7 ) != 0 );
            cur->npoint = 0;
            rowcap = oracle_sample_axis_count() + 8;
            cur->in  = alloc_cols( cur->nin,  rowcap );
            cur->out = alloc_cols( cur->nout, rowcap );
            row = 0;
            continue;
        }
        if ( !cur ) continue;
        /* a data row: nin inputs then nout outputs */
        char *p = line, *tok;
        int col = 0, ok = 1;
        if ( row >= rowcap ) {
            rowcap *= 2;
            for ( int a = 0; a < cur->nin;  a++ )
                cur->in[a]  = realloc( cur->in[a],  sizeof(double)*(size_t)rowcap );
            for ( int a = 0; a < cur->nout; a++ )
                cur->out[a] = realloc( cur->out[a], sizeof(double)*(size_t)rowcap );
        }
        while ( ( tok = strtok( p, " \t\r\n" ) ) ) {
            p = NULL;
            int good = 0;
            double v = oracle_parse_double( tok, &good );
            if ( !good ) { ok = 0; break; }
            if ( col < cur->nin ) cur->in[col][row] = v;
            else if ( col < cur->nin + cur->nout ) cur->out[col - cur->nin][row] = v;
            col++;
        }
        if ( ok && col == cur->nin + cur->nout ) { row++; cur->npoint = row; }
    }
    fclose( fp );
    *out_secs = secs; *out_n = (int) n;
    return 0;
}

static int selftest( void ) {
    const char *tmp = "/tmp/oracle_selftest.txt";
    FILE *fp = fopen( tmp, "w" );
    fputs( "# header\n\n"
           "[foo.map  nin=2 nout=1 dir=forward]\n"
           "  1 2  3\n"
           "  4 5  BAD\n\n", fp );
    fclose( fp );
    Section *s = NULL; int n = 0;
    if ( oracle_read_file( tmp, &s, &n ) || n != 1 ) {
        fprintf( stderr, "selftest: parse failed\n" ); return 1;
    }
    int ok = ( strcmp( s[0].relpath, "foo.map" ) == 0 && s[0].nin == 2 &&
               s[0].nout == 1 && s[0].forward == 1 && s[0].npoint == 2 &&
               s[0].in[0][0] == 1 && s[0].in[1][1] == 5 &&
               s[0].out[0][0] == 3 && s[0].out[0][1] == AST__BAD );
    printf( "selftest: %s\n", ok ? "ok" : "FAIL" );
    return ok ? 0 : 1;
}

int main( int argc, char *argv[] ) {
    int status_value = 0; int *status = &status_value;
    astWatch( status );
    if ( argc == 2 && strcmp( argv[1], "--selftest" ) == 0 ) return selftest();
    fprintf( stderr, "checking not yet implemented\n" );
    return 2;   /* replaced in Task 6 */
}
```

> Note on the `relpath` header parse: `%1023s` reads up to the first whitespace, which is exactly the relpath (paths in this corpus contain no spaces). The trailing `nin=`/`nout=`/`dir=` fields then match regardless of the double-space separator.

- [ ] **Step 2: Register the checker target and a selftest test**

Add to `ast_tester/CMakeLists.txt`:

```cmake
add_executable(check_transform_oracle check_transform_oracle.c transform_oracle_util.c)
target_include_directories(check_transform_oracle PRIVATE ${TEST_INCLUDES})
target_link_libraries(check_transform_oracle PRIVATE ${TEST_LIBS})
target_compile_definitions(check_transform_oracle PRIVATE HAVE_CONFIG_H)
set_target_properties(check_transform_oracle PROPERTIES C_STANDARD 11)
ast_apply_dev_options(check_transform_oracle)

add_test(NAME transform_oracle_selftest
         COMMAND check_transform_oracle --selftest)
```

- [ ] **Step 3: Build and run the selftest, expect pass**

Run:
```bash
cmake --build build --target check_transform_oracle && ctest --test-dir build -R transform_oracle_selftest --output-on-failure
```
Expected: PASS (`selftest: ok`).

- [ ] **Step 4: Commit**

```bash
git add ast_tester/check_transform_oracle.c ast_tester/CMakeLists.txt
git commit -m "test: add transform-oracle file parser with selftest"
```

---

### Task 6: Checker comparisons (A, B, round-trip, C)

**Files:**
- Modify: `ast_tester/check_transform_oracle.c`

**Interfaces:**
- Consumes: `oracle_load_mapping`, `oracle_within_tol`, `astTranP`, parsed `Section`s.
- Produces: full `main` that runs golden checks A/B, round-trip, and `.map`/`.simp` equivalence check C; exits non-zero on any failure with per-failure diagnostics. CLI: `check_transform_oracle <root> <oracle_file> [<rtrip_overrides>]`.

- [ ] **Step 1: Add the comparison core**

In `check_transform_oracle.c`, add above `main` (after `oracle_read_file`):

```c
/* Re-transform a section's recorded inputs with the live mapping into
   freshly allocated columns the caller frees.  Returns NULL on arity
   mismatch (a structural change to the fixture). */
static double **transform_section( AstMapping *map, const Section *s ) {
    if ( astGetI( map, "Nin" ) != ( s->forward ? s->nin : s->nout ) ) return NULL;
    if ( astGetI( map, "Nout" ) != ( s->forward ? s->nout : s->nin ) ) return NULL;
    double **live = alloc_cols( s->nout, s->npoint );
    astTranP( map, s->npoint, s->nin, (const double **) s->in,
              s->forward, s->nout, live );
    return live;
}

static void free_cols( double **c, int ncol ) {
    if ( !c ) return;
    for ( int a = 0; a < ncol; a++ ) free( c[a] );
    free( c );
}

/* Compare two equal-shape output sets; report and count mismatches. */
static int compare_outputs( const char *label, const char *relpath,
                            double **got, double **ref, int ncol, int npoint,
                            double rtol, double atol ) {
    int fails = 0;
    for ( int p = 0; p < npoint; p++ ) {
        for ( int a = 0; a < ncol; a++ ) {
            if ( !oracle_within_tol( got[a][p], ref[a][p], rtol, atol ) ) {
                double d = fabs( got[a][p] - ref[a][p] );
                double rel = ref[a][p] != 0.0 ? d / fabs( ref[a][p] ) : d;
                fprintf( stderr,
                    "MISMATCH [%s] %s row=%d axis=%d ref=%.17g got=%.17g rel=%.3g\n",
                    label, relpath, p, a, ref[a][p], got[a][p], rel );
                fails++;
            }
        }
    }
    return fails;
}
```

Add `#include <math.h>` at the top of the file.

- [ ] **Step 2: Add round-trip overrides lookup**

Add above `main`:

```c
typedef struct { char relpath[1024]; double rtol, atol; int off; } Override;

static Override *g_over = NULL;
static int g_nover = 0;

static void load_overrides( const char *path ) {
    if ( !path ) return;
    FILE *fp = fopen( path, "r" );
    if ( !fp ) return;
    size_t cap = 16; g_over = malloc( sizeof(Override) * cap );
    char line[1024];
    while ( fgets( line, (int) sizeof line, fp ) ) {
        if ( line[0] == '#' || line[0] == '\n' ) continue;
        char rel[1024], what[64];
        /* "relpath off"  or  "relpath <rtol> <atol>" */
        if ( sscanf( line, "%1023s %63s", rel, what ) < 2 ) continue;
        if ( g_nover == (int) cap ) { cap *= 2; g_over = realloc( g_over, sizeof(Override)*cap ); }
        Override *o = &g_over[g_nover++];
        snprintf( o->relpath, sizeof o->relpath, "%s", rel );
        if ( strcmp( what, "off" ) == 0 ) { o->off = 1; o->rtol = o->atol = 0; }
        else { o->off = 0; sscanf( line, "%*s %lf %lf", &o->rtol, &o->atol ); }
    }
    fclose( fp );
}

/* Fill rtol/atol for a fixture's round-trip; return 0 if round-trip is
   switched off for it. */
static int rtrip_tol( const char *relpath, double *rtol, double *atol ) {
    *rtol = ORACLE_DEF_RTOL; *atol = ORACLE_DEF_ATOL;
    for ( int i = 0; i < g_nover; i++ ) {
        if ( strcmp( g_over[i].relpath, relpath ) == 0 ) {
            if ( g_over[i].off ) return 0;
            *rtol = g_over[i].rtol; *atol = g_over[i].atol; return 1;
        }
    }
    return 1;
}
```

- [ ] **Step 3: Replace `main` with the full checker**

Replace the body of `main` (keep the `--selftest` branch):

```c
int main( int argc, char *argv[] ) {
    int status_value = 0; int *status = &status_value;
    astWatch( status );
    if ( argc == 2 && strcmp( argv[1], "--selftest" ) == 0 ) return selftest();
    if ( argc < 3 ) {
        fprintf( stderr,
            "Usage: check_transform_oracle <root> <oracle_file> [<overrides>]\n" );
        return 2;
    }
    const char *root   = argv[1];
    const char *oracle = argv[2];
    if ( argc >= 4 ) load_overrides( argv[3] );

    Section *secs = NULL; int nsec = 0;
    if ( oracle_read_file( oracle, &secs, &nsec ) ) return 2;

    int total_fail = 0, checked = 0;

    /* Cache one stem's live .map outputs so check C can compare against the
       matching .simp without re-deriving. Keyed by stem (relpath minus
       extension). Since gen writes .map immediately before its .simp, a
       single one-slot cache suffices. */
    char prev_stem[1024] = "";
    double **prev_map_out = NULL; int prev_nout = 0, prev_np = 0;

    for ( int i = 0; i < nsec; i++ ) {
        Section *s = &secs[i];
        astBegin;
        AstMapping *map = oracle_load_mapping( root, s->relpath );
        if ( !map ) {
            fprintf( stderr, "LOADFAIL %s (fixture missing or unreadable)\n",
                     s->relpath );
            total_fail++;
            astEnd;
            continue;
        }
        double **live = transform_section( map, s );
        if ( !live ) {
            fprintf( stderr, "ARITY %s (Nin/Nout changed)\n", s->relpath );
            total_fail++;
            map = astAnnul( map ); astEnd; continue;
        }

        /* Golden check (A for .map, B for .simp, or the .head golden). */
        total_fail += compare_outputs( "golden", s->relpath, live, s->out,
                                       s->nout, s->npoint,
                                       ORACLE_DEF_RTOL, ORACLE_DEF_ATOL );
        checked++;

        /* Round-trip, where the live mapping defines the inverse and we
           recorded the forward direction. */
        if ( s->forward && astGetI( map, "TranInverse" ) ) {
            double rtol, atol;
            if ( rtrip_tol( s->relpath, &rtol, &atol ) ) {
                double **back = alloc_cols( s->nin, s->npoint );
                astTranP( map, s->npoint, s->nout, (const double **) live,
                          0, s->nin, back );
                total_fail += compare_outputs( "round-trip", s->relpath,
                                               back, s->in, s->nin, s->npoint,
                                               rtol, atol );
                free_cols( back, s->nin );
            }
        }

        /* Equivalence check C: when this section is a .simp whose stem
           matches the cached preceding .map, compare their live outputs. */
        size_t L = strlen( s->relpath );
        char stem[1024];
        const char *ext = strrchr( s->relpath, '.' );
        size_t stemlen = ext ? (size_t)( ext - s->relpath ) : L;
        snprintf( stem, sizeof stem, "%.*s", (int) stemlen, s->relpath );
        if ( ext && strcmp( ext, ".simp" ) == 0 && prev_map_out &&
             strcmp( stem, prev_stem ) == 0 &&
             prev_nout == s->nout && prev_np == s->npoint ) {
            total_fail += compare_outputs( "equiv", s->relpath, live,
                                           prev_map_out, s->nout, s->npoint,
                                           ORACLE_DEF_EQUIV_RTOL,
                                           ORACLE_DEF_EQUIV_ATOL );
        }

        /* Refresh the one-slot cache when this section is a .map. */
        free_cols( prev_map_out, prev_nout );
        prev_map_out = NULL;
        if ( ext && strcmp( ext, ".map" ) == 0 ) {
            prev_map_out = live;            /* hand ownership to the cache */
            prev_nout = s->nout; prev_np = s->npoint;
            snprintf( prev_stem, sizeof prev_stem, "%s", stem );
        } else {
            free_cols( live, s->nout );
        }

        map = astAnnul( map );
        astEnd;
    }
    free_cols( prev_map_out, prev_nout );

    printf( "check_transform_oracle: %d sections checked, %d mismatch(es)\n",
            checked, total_fail );
    return ( total_fail == 0 && astOK ) ? 0 : 1;
}
```

> The one-slot `.map` cache works because the generator emits each stem's `.map` immediately followed by its `.simp` with identical inputs. The cache holds the live `.map` outputs across exactly that gap.

- [ ] **Step 4: Build, then check against a freshly generated oracle (expect pass)**

Run:
```bash
cmake --build build --target gen_transform_oracle check_transform_oracle
./build/ast_tester/gen_transform_oracle ast_tester /tmp/simp.oracle /tmp/head.oracle
./build/ast_tester/check_transform_oracle ast_tester /tmp/simp.oracle
./build/ast_tester/check_transform_oracle ast_tester /tmp/head.oracle
```
Expected: both print `... 0 mismatch(es)` and exit 0. (Round-trip mismatches on iterative inverses such as PolyMap are expected here and are resolved in Task 7 via overrides; note which relpaths report `round-trip` mismatches.)

- [ ] **Step 5: Verify the checker actually fails on drift**

Run:
```bash
cp /tmp/simp.oracle /tmp/tampered.oracle
# perturb the first numeric output token of the first data row
awk 'BEGIN{done=0} /^\[/{print;next} (!done && NF>0 && $1 !~ /^#/){ $NF=$NF+1.0; done=1 } {print}' \
    /tmp/simp.oracle > /tmp/tampered.oracle
./build/ast_tester/check_transform_oracle ast_tester /tmp/tampered.oracle; echo "exit=$?"
```
Expected: a `MISMATCH [golden] ...` line and `exit=1`.

- [ ] **Step 6: Commit**

```bash
git add ast_tester/check_transform_oracle.c
git commit -m "feat: implement transform-oracle golden, round-trip, and equivalence checks"
```

---

### Task 7: Round-trip overrides file and wiring

**Files:**
- Create: `ast_tester/transform_oracle_rtrip_overrides.txt`
- Modify: `ast_tester/CMakeLists.txt`

**Interfaces:**
- Consumes: the `load_overrides`/`rtrip_tol` mechanism from Task 6.
- Produces: a committed, human-curated overrides file, applied by the registered `ctest` checks (added in Task 8).

- [ ] **Step 1: Identify fixtures whose round-trip exceeds the default tolerance**

Run (using the oracle generated in Task 6):
```bash
./build/ast_tester/check_transform_oracle ast_tester /tmp/simp.oracle 2>&1 \
  | awk '/round-trip/{print $3}' | sort -u
./build/ast_tester/check_transform_oracle ast_tester /tmp/head.oracle 2>&1 \
  | awk '/round-trip/{print $3}' | sort -u
```
Expected: a (possibly empty) list of relpaths. For each, decide per the spec: investigate whether it is a genuine iterative-inverse limitation (e.g. PolyMap) or a real bug. **Do not** blanket-loosen; set a per-fixture tolerance no looser than the observed residual needs, or `off` only when an inverse is legitimately not round-trippable.

- [ ] **Step 2: Write the overrides file**

Create `ast_tester/transform_oracle_rtrip_overrides.txt`. Format: `# comment`, blank lines ignored, otherwise `relpath off` or `relpath <rtol> <atol>`. Seed it with the fixtures found in Step 1, each with a one-line justification, for example:

```
# Per-fixture round-trip tolerance overrides for check_transform_oracle.
# Format:  <relpath> off            disable round-trip for this fixture
#          <relpath> <rtol> <atol>  use a looser round-trip tolerance
#
# PolyMap inverses are iterative (Newton), so round-trip accuracy is bounded
# by convergence, not ULP noise.  Values below are the measured residual
# rounded up, not a blanket loosening.
# (Fill in from Task 7 Step 1 output; leave empty if none exceed default.)
```

Add the concrete entries discovered in Step 1 beneath the header (one line each). If Step 1 produced no fixtures, leave only the comment header.

- [ ] **Step 3: Confirm overrides clear the round-trip mismatches**

Run:
```bash
./build/ast_tester/check_transform_oracle ast_tester /tmp/simp.oracle \
    ast_tester/transform_oracle_rtrip_overrides.txt; echo "exit=$?"
```
Expected: `0 mismatch(es)` and `exit=0` (golden and equivalence still enforced; only the listed round-trips are relaxed).

- [ ] **Step 4: Commit**

```bash
git add ast_tester/transform_oracle_rtrip_overrides.txt
git commit -m "feat: add round-trip tolerance overrides for iterative inverses"
```

---

### Task 8: Generate, commit reference data, and wire the ctest gate

**Files:**
- Create: `ast_tester/simplify_fixtures.oracle`
- Create: `ast_tester/headers.oracle`
- Modify: `ast_tester/CMakeLists.txt`
- Modify: `PLAN.md`

**Interfaces:**
- Consumes: `gen_transform_oracle`, `check_transform_oracle`, the overrides file.
- Produces: two committed oracle files and two registered `ctest` tests.

- [ ] **Step 1: Generate the committed oracle files into the source tree**

Run:
```bash
./build/ast_tester/gen_transform_oracle ast_tester \
    ast_tester/simplify_fixtures.oracle ast_tester/headers.oracle
wc -l ast_tester/simplify_fixtures.oracle ast_tester/headers.oracle
```
Expected: both files written; line counts in the thousands (hundreds of sections × tens of rows). Spot-check that the files are plain text and diff-friendly.

- [ ] **Step 2: Register the two gate tests**

Add to `ast_tester/CMakeLists.txt` after the checker target:

```cmake
add_test(NAME transform_oracle_simplify
         COMMAND check_transform_oracle
                 "${CMAKE_CURRENT_SOURCE_DIR}"
                 "${CMAKE_CURRENT_SOURCE_DIR}/simplify_fixtures.oracle"
                 "${CMAKE_CURRENT_SOURCE_DIR}/transform_oracle_rtrip_overrides.txt")
add_test(NAME transform_oracle_headers
         COMMAND check_transform_oracle
                 "${CMAKE_CURRENT_SOURCE_DIR}"
                 "${CMAKE_CURRENT_SOURCE_DIR}/headers.oracle"
                 "${CMAKE_CURRENT_SOURCE_DIR}/transform_oracle_rtrip_overrides.txt")
if(AST_ENABLE_SANITIZERS)
    set_tests_properties(transform_oracle_simplify transform_oracle_headers
        PROPERTIES ENVIRONMENT "ASAN_OPTIONS=detect_leaks=0")
endif()
```

- [ ] **Step 3: Run the gate tests, expect pass**

Run:
```bash
cmake --build build
ctest --test-dir build -R 'transform_oracle' --output-on-failure
```
Expected: `transform_oracle_selftest`, `transform_oracle_simplify`, `transform_oracle_headers` all pass.

- [ ] **Step 4: Update PLAN.md**

Add an entry to `PLAN.md` recording the new oracle: its two programs, the two committed oracle files, the three checks (golden/round-trip/equivalence), the regeneration command, and the overrides file. Match the surrounding `PLAN.md` formatting.

- [ ] **Step 5: Commit**

```bash
git add ast_tester/simplify_fixtures.oracle ast_tester/headers.oracle \
        ast_tester/CMakeLists.txt PLAN.md
git commit -m "feat: commit transform oracle reference data and ctest gate"
```

---

### Task 9: Tolerance tuning and full verification

**Files:**
- Modify: `ast_tester/transform_oracle.h` (final tolerance constants, if adjusted)
- Modify: `ast_tester/simplify_fixtures.oracle`, `ast_tester/headers.oracle` (regenerate if header text changes)
- Modify: `docs/superpowers/specs/2026-06-27-transform-oracle-design.md` (record final values)

**Interfaces:** none new.

- [ ] **Step 1: Measure the actual output spread to choose final tolerances**

Run a self-vs-self margin probe: regenerate the oracle, then re-run the checker but with tolerances temporarily tightened to 0 to print every nonzero residual, capturing the largest relative differences observed per check class.

```bash
./build/ast_tester/gen_transform_oracle ast_tester /tmp/a_simp.oracle /tmp/a_head.oracle
# golden self-consistency should be exact on the same machine/build:
./build/ast_tester/check_transform_oracle ast_tester /tmp/a_simp.oracle \
    ast_tester/transform_oracle_rtrip_overrides.txt
```
Expected: golden mismatches are zero on the same build (the checker recomputes identical numbers). The meaningful spreads are: round-trip residuals (from Task 7 Step 1) and equivalence residuals (`equiv` lines). Inspect those magnitudes.

- [ ] **Step 2: Set final tolerances with margin**

Choose final `ORACLE_DEF_*` values that sit comfortably above the largest legitimate residual observed for each class (round-trip and equivalence) and far below any plausible algorithmic change (≥ ~1e-6 relative). If the cross-architecture story is available, sanity-check on an ARM machine too; otherwise document the x86 basis. Edit the constants in `transform_oracle.h` and the header-comment defaults in `gen_transform_oracle.c` to match. If the comment line in the oracle header changes, regenerate both oracle files (Task 8 Step 1).

- [ ] **Step 3: Record final values in the spec**

Update the "Tolerances" and "Open implementation parameters" sections of `docs/superpowers/specs/2026-06-27-transform-oracle-design.md` to state the chosen final values and the basis (observed spread, machine(s) used).

- [ ] **Step 4: Full clean build + sanitizer run + whole suite**

Run:
```bash
cmake -B build-dev -DCMAKE_BUILD_TYPE=Debug -DAST_ENABLE_WARNINGS=ON -DAST_ENABLE_SANITIZERS=ON
cmake --build build-dev 2>&1 | tee /tmp/oracle_build.log
ctest --test-dir build-dev -R 'transform_oracle|test_oracle_util' --output-on-failure
```
Expected: no new warnings attributable to the new files in `/tmp/oracle_build.log`; all oracle tests pass under ASan/UBSan.

- [ ] **Step 5: Commit**

```bash
git add ast_tester/transform_oracle.h ast_tester/gen_transform_oracle.c \
        ast_tester/simplify_fixtures.oracle ast_tester/headers.oracle \
        docs/superpowers/specs/2026-06-27-transform-oracle-design.md
git commit -m "chore: finalize transform-oracle tolerances and document basis"
```

---

## Self-Review

**Spec coverage:**
- Tolerance-based, arch-independent guard → Task 1 (`oracle_within_tol`), Task 9 (tuning). ✓
- Two corpora (native dumps + FITS headers) → Task 3 loader, Task 4 scan, Task 8 two files. ✓
- `.map` and `.simp` both recorded as distinct entries → Task 4 emits both per stem. ✓
- Scan-only generator / file-driven checker → Task 4 scans; Task 5/6 checker only reads the oracle + named fixtures. ✓
- Self-describing sectioned text format, `%.17g`, `BAD` token → Tasks 3,4,5. ✓
- Inputs read from file, not regenerated → Task 5 parses inputs; Task 6 transforms them. ✓
- Forward stored; round-trip where both directions; inverse-only stored directly → Task 4 (`emit_fixture`), Task 6 (round-trip). ✓
- Checks A/B/C → Task 6 (`golden`, `equiv`). ✓
- Per-fixture round-trip tolerance for iterative inverses, treated as a signal → Task 7. ✓
- Two oracle files matching the two loader paths → Task 8. ✓
- Low-discrepancy + edge sampling, tunable → Task 2; adaptive broaden noted as soft/optional (see note below). 
- ctest integration, regeneration workflow, PLAN.md → Task 8. ✓
- Tolerances tuned empirically → Task 9. ✓

**Adaptive broadening:** The spec describes adaptive range broadening as a *soft* heuristic. This plan implements fixed per-corpus ranges (Task 4 `axis_range`) and does **not** implement auto-broadening in the first cut, to keep generation deterministic and machine-independent. If, during Task 4 Step 3/4 inspection, a class of fixtures is clearly never exercised (all outputs `BAD` or constant), widen the constants in `axis_range` and regenerate — a manual, reviewable adjustment rather than runtime adaptation. This is a deliberate YAGNI narrowing of the soft requirement; flagged here for the reviewer.

**Placeholder scan:** No `TBD`/`TODO`/"add error handling"-style placeholders; every code step contains complete code. The overrides file (Task 7) is intentionally seeded from measured output, with the procedure spelled out, not left vague.

**Type consistency:** `oracle_within_tol`, `oracle_halton`, `oracle_sample_points`, `oracle_sample_axis_count`, `oracle_load_mapping`, `oracle_format_double`, `oracle_parse_double` are declared in Task 1–3 and used with matching signatures in Tasks 4–6. The `Section` struct and `oracle_read_file` defined in Task 5 are consumed unchanged in Task 6. CLI shapes (`gen ... <root> <simp_out> <head_out>`, `check ... <root> <oracle> [overrides]`) are consistent between Tasks 4/6 and the CMake wiring in Task 8.
