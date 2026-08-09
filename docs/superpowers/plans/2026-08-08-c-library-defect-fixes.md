# C Library Defect Fixes Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Repair nineteen verified defects in the AST C library across `keymap.c`, `fitschan.c`, `cmpmap.c` and `mapping.c`, plus a tree-wide correction of a rounding idiom that is wrong for every negative value.

**Architecture:** Each defect lands as its own commit on the existing branch `u/timj/misc-bugs`, so a reviewer can accept or reject one without touching its neighbours. Where the repair is observable through the public API, a regression test is written first and confirmed failing against the unfixed library before the fix is applied. The rounding correction lands first because later tasks modify the same lines.

**Tech Stack:** C99 (C11 for serialization tests), CMake ≥ 3.24, CTest, clang/gcc with AddressSanitizer and UndefinedBehaviorSanitizer.

**Source spec:** `docs/superpowers/specs/2026-08-08-c-library-defect-fixes-design.md`

## Global Constraints

- Work on branch `u/timj/misc-bugs`. Never push to a remote.
- One commit per defect. Do not batch unrelated fixes.
- Every commit that touches a file under `src/` MUST add an entry to that file's prologue `History:` section, immediately before the closing `*class--` or `*--` line, in the established format:
  ```
  *     8-AUG-2026 (TIMJ):
  *        One-line description of the change.
  ```
- Commit messages must not reference any port of this library to another language, any plan or spec document, or the history of how the defect was found. Describe the defect and the fix only.
- Every commit message ends with:
  ```
  Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
  ```
- No change may introduce a new compiler warning. The AST source is not warning-clean at full warnings; your changes must not add to it.
- C99, no C++ features. Use the public API (`ast.h`) only in `ast_tester/` tests.
- Tests use `astWatch(&status)` + `astOK` for error checking, never Starlink EMS.
- Do NOT modify `c-bugs.md`. Do NOT commit it.

## Build and Test Commands

Two build trees are used throughout.

Plain build, for fast iteration:

```bash
cmake -B build
cmake --build build -j
ctest --test-dir build --output-on-failure
```

Sanitizer build, which every commit must pass before it is made:

```bash
cmake -B build-dev -DCMAKE_BUILD_TYPE=Debug \
  -DAST_ENABLE_WARNINGS=ON -DAST_ENABLE_SANITIZERS=ON
cmake --build build-dev -j
ctest --test-dir build-dev --output-on-failure
```

Run a single test:

```bash
ctest --test-dir build -R testkeymap --output-on-failure
```

**Note on the sanitizer build:** if AddressSanitizer hangs on startup, the compiler is probably a conda one. Use the homebrew clang explicitly:

```bash
cmake -B build-dev -DCMAKE_C_COMPILER=/opt/homebrew/opt/llvm/bin/clang \
  -DCMAKE_BUILD_TYPE=Debug -DAST_ENABLE_WARNINGS=ON -DAST_ENABLE_SANITIZERS=ON
```

## Writing a Scratch Probe

Several tasks call for a throwaway program that demonstrates a defect before it is fixed. Write these to a scratch directory, never into the repository, and never commit them.

Save as `/tmp/probe.c`, then build and run with:

```bash
cc -O0 -I build -o /tmp/probe /tmp/probe.c \
   -L build -last -Wl,-force_load,build/libast_err.a -lm
DYLD_LIBRARY_PATH=build /tmp/probe
```

On Linux use `LD_LIBRARY_PATH` in place of `DYLD_LIBRARY_PATH`, and `-Wl,--whole-archive build/libast_err.a -Wl,--no-whole-archive` in place of `-Wl,-force_load`.

The `-force_load` of `libast_err.a` is required: `libast` references `astPutErr` without defining it, and without whole-archive linkage the symbol resolves to a null pointer at runtime and the program crashes.

## Test File Conventions

Two different `stopit` helpers exist. Use the one belonging to the file you are editing.

`ast_tester/testkeymap.c:54`:

```c
static void stopit( int *status, const char *text );
/* called as: stopit( status, "Error 15" ); */
```

`ast_tester/testfitschan.c:232`:

```c
static void stopit( int errnum, const char *text, int *status );
/* called as: stopit( 12401, "grism order not preserved", status ); */
```

Error numbers already used in `testfitschan.c` go up to 12344. New checks in this plan use 12400 and above.

No new test executables are added. All new checks go into `ast_tester/testkeymap.c` and `ast_tester/testfitschan.c`, which are already registered in `ast_tester/CMakeLists.txt` (lines 192 and elsewhere). `PLAN.md` does not need updating, because that file tracks Fortran-to-C test conversions and no conversion is happening here.

## File Structure

| File | Responsibility in this plan |
|---|---|
| `src/keymap.c` | Tasks 2, 4-13. Ten defects plus the KeyMap share of the rounding sweep. |
| `src/fitschan.c` | Tasks 1, 3, 15, 17, 18. Grism order, `EncodeFloat`, `WcsNative`, `CDjjjiii`, duplicate consumption. |
| `src/mapping.c`, `src/cmpmap.c` | Task 16 only, and only if the investigation selects the flag-split fix. |
| 13 further `src/*.c` files | Task 3, the no-op share of the rounding sweep. |
| `ast_tester/testkeymap.c` | New regression checks for Tasks 2, 4, 5, 6, 7, 11. |
| `ast_tester/testfitschan.c` | New regression check for Task 1. |

---

## Task 1: Negative grism interference order

The rounding sweep's one behaviourally significant site outside `keymap.c`. Done first because it carries the clearest user-visible symptom.

**Files:**
- Modify: `src/fitschan.c:14322`
- Test: `ast_tester/testfitschan.c`

**Interfaces:**
- Consumes: nothing from earlier tasks.
- Produces: nothing later tasks depend on.

**Background.** `GrismSpecWcs` sets a GrismMap's interference order from FITS keyword `PVi_1`. `GrismM` is a plain signed `int` (`src/grismmap.c:1923-1927`) and negative diffraction orders are physically routine. The current cast truncates toward zero, so an order of −1 becomes 0, which is degenerate and makes `astRead` reject the header, and an order of −2 becomes −1, silently producing the wrong dispersion.

- [ ] **Step 1: Write the failing test**

Add this function to `ast_tester/testfitschan.c`, immediately above `int main( void )` at line 1218:

```c
/*
 * testgrismorder: a FITS grism header must preserve a negative
 * interference order. PVi_1 holds the order, which is commonly negative.
 */
static void testgrismorder( int errbase, double morder, int *status ) {
   AstFitsChan *fc;
   AstFitsChan *out;
   AstFrameSet *fs;
   char buf[ 81 ];
   double back;

   if( !astOK ) return;

   fc = astFitsChan( NULL, NULL, "" );
   astPutFits( fc, "CTYPE1  = 'WAVE-GRI'", 0 );
   astPutFits( fc, "CRVAL1  =            5.0E-7", 0 );
   astPutFits( fc, "CRPIX1  =                1.0", 0 );
   astPutFits( fc, "CDELT1  =            1.0E-9", 0 );
   astPutFits( fc, "CUNIT1  = 'm       '", 0 );
   astPutFits( fc, "PV1_0   =            1.0E6", 0 );
   sprintf( buf, "PV1_1   = %20.1f", morder );
   astPutFits( fc, buf, 0 );
   astPutFits( fc, "PV1_2   =                0.0", 0 );
   astPutFits( fc, "PV1_3   =                1.0", 0 );
   astPutFits( fc, "PV1_4   =                0.0", 0 );
   astPutFits( fc, "PV1_5   =                0.0", 0 );
   astPutFits( fc, "PV1_6   =                0.0", 0 );
   astClear( fc, "Card" );

   fs = astRead( fc );
   if( !astOK || !fs ) {
      stopit( errbase, "grism header could not be read", status );
      astClearStatus;
      return;
   }

   out = astFitsChan( NULL, NULL, "Encoding=FITS-WCS" );
   if( !astWrite( out, fs ) ) {
      stopit( errbase + 1, "grism FrameSet produced no header", status );
   } else {
      astClear( out, "Card" );
      if( !astGetFitsF( out, "PV1_1", &back ) ) {
         stopit( errbase + 2, "no PV1_1 written back", status );
      } else if( back != morder ) {
         printf( "PV1_1 in = %g, out = %g\n", morder, back );
         stopit( errbase + 3, "grism interference order not preserved",
                 status );
      }
   }
}
```

Then add these three calls inside `main`, immediately after `astBegin;` (line 1238):

```c
   /* A grism interference order is commonly negative and must survive a
      read/write round trip. The positive case is a control: it confirms
      the chosen grating parameters are physically usable, so a failure in
      the negative cases cannot be blamed on the parameter choice. */
   testgrismorder( 12400, 1.0, status );
   testgrismorder( 12410, -1.0, status );
   testgrismorder( 12420, -2.0, status );
```

- [ ] **Step 2: Run the test to verify it fails**

```bash
cmake --build build -j && ctest --test-dir build -R testfitschan --output-on-failure
```

Expected: FAIL. The `-1.0` case reports error 12410 ("grism header could not be read"), because the order rounds to 0 and `astRead` rejects it. The `-2.0` case reports error 12423 and prints `PV1_1 in = -2, out = -1`. The `1.0` control passes.

If the control at 12400 also fails, the grating parameters in the test are wrong, not the library. Stop and fix the test before going further.

- [ ] **Step 3: Apply the fix**

In `src/fitschan.c:14322`, change:

```c
      astSetGrismM( gmap, ( pv != AST__BAD )?(int) ( pv + 0.5 ):0);
```

to:

```c
      astSetGrismM( gmap, ( pv != AST__BAD )?(int) round( pv ):0);
```

`fitschan.c` already includes `<math.h>`; confirm with `grep -n '#include <math.h>' src/fitschan.c` before assuming it.

- [ ] **Step 4: Add the prologue history entry**

In `src/fitschan.c`, find the `History:` section and insert immediately before the closing `*class--` line:

```
*     8-AUG-2026 (TIMJ):
*        Use round() rather than truncation when reading the grism
*        interference order from PVi_1, so that negative orders are
*        preserved rather than being rounded toward zero.
```

- [ ] **Step 5: Run the test to verify it passes**

```bash
cmake --build build -j && ctest --test-dir build -R testfitschan --output-on-failure
```

Expected: PASS.

- [ ] **Step 6: Run the full suite and the sanitizer build**

```bash
ctest --test-dir build --output-on-failure
cmake --build build-dev -j && ctest --test-dir build-dev --output-on-failure
```

Expected: no failures, no new warnings.

- [ ] **Step 7: Commit**

```bash
git add src/fitschan.c ast_tester/testfitschan.c
git commit -F - <<'EOF'
fitschan: preserve negative grism interference orders

GrismSpecWcs read the interference order from FITS keyword PVi_1 with
(int)( pv + 0.5 ). A cast from floating point to integer truncates
toward zero, so this expression is correct only for non-negative values.

Negative diffraction orders are physically routine. An order of -1 was
rounded to 0, which is degenerate, so astRead rejected an otherwise
valid header with "unusable parameter values". An order of -2 was
stored as -1, producing a WCS with the wrong dispersion and no
diagnostic.

Use round() instead. Positive orders are unaffected.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
EOF
```

---

## Task 2: KeyMap value conversion rounds negatives toward zero

**Files:**
- Modify: `src/keymap.c:1953, 1957, 1961, 1965, 2012, 2016, 2020, 2024, 2069, 2085, 2101, 2117`
- Test: `ast_tester/testkeymap.c`

**Interfaces:**
- Consumes: nothing.
- Produces: the corrected `ConvertValue` arms that Task 6 wraps in range checks.

**Background.** `ConvertValue` converts a stored `Double`, `Float` or numeric `String` to a requested integer type using `(int)( dval + 0.5 )`. Every negative value is rounded toward zero by one, exact integers included: a stored `-2.0` is returned as `-1`.

The existing coverage at `ast_tester/testkeymap.c:559` uses `1999.9`, a positive value that rounds identically either way, which is why this survived.

- [ ] **Step 1: Write the failing test**

In `ast_tester/testkeymap.c`, find the block near line 377 that puts the `Fred*` entries:

```c
   astMapPut0D( map, "Fredd", 1999.9, "com2 " );
```

Add these puts immediately after it:

```c
   astMapPut0D( map, "Negd", -1999.9, "com2 " );
   astMapPut0F( map, "Negf", -1999.9f, "com2 " );
   astMapPut0D( map, "Negint", -2.0, "com2 " );
   astMapPut0D( map, "Neghalf", -2.5, "com2 " );
   astMapPut0D( map, "Negsmall", -0.6, "com2 " );
   astMapPut0D( map, "Negtiny", -0.4, "com2 " );
```

Then find the existing check near line 559:

```c
   /* double key read as int: 1999.9 -> 2000 */
   if( !astMapGet0I( map2, "Fredd", &ival ) ) {
      stopit( status, "Error 14" );
   } else if( ival != 2000 ) {
```

Add this block immediately after that check's closing brace:

```c
   /* Negative values must round to nearest, not toward zero. A cast of
      the form (int)( x + 0.5 ) truncates toward zero and so is wrong for
      every negative value, exact integers included. */
   {
      struct { const char *key; int expect; } negcases[] = {
         { "Negd",     -2000 },
         { "Negf",     -2000 },
         { "Negint",      -2 },
         { "Neghalf",     -3 },
         { "Negsmall",    -1 },
         { "Negtiny",      0 }
      };
      short int sval_n;
      int64_t kval_n;
      size_t ic;

      for( ic = 0; ic < sizeof( negcases ) / sizeof( negcases[ 0 ] ); ic++ ) {
         if( !astMapGet0I( map2, negcases[ ic ].key, &ival ) ) {
            printf( "%s\n", negcases[ ic ].key );
            stopit( status, "Error 14neg-get" );
         } else if( ival != negcases[ ic ].expect ) {
            printf( "%s: got %d want %d\n", negcases[ ic ].key, ival,
                    negcases[ ic ].expect );
            stopit( status, "Error 14neg-int" );
         }

         if( !astMapGet0S( map2, negcases[ ic ].key, &sval_n ) ) {
            stopit( status, "Error 14neg-gets" );
         } else if( sval_n != (short int) negcases[ ic ].expect ) {
            printf( "%s: got %d want %d\n", negcases[ ic ].key, (int) sval_n,
                    negcases[ ic ].expect );
            stopit( status, "Error 14neg-short" );
         }

         if( !astMapGet0K( map2, negcases[ ic ].key, &kval_n ) ) {
            stopit( status, "Error 14neg-getk" );
         } else if( kval_n != (int64_t) negcases[ ic ].expect ) {
            printf( "%s: got %" PRId64 " want %d\n", negcases[ ic ].key,
                    kval_n, negcases[ ic ].expect );
            stopit( status, "Error 14neg-int64" );
         }
      }
   }
```

`Negtiny` (−0.4 → 0) is a control: it rounds to the same value either way, so it would catch an over-correction that rounded away from zero unconditionally. The `Byte` arm is deliberately not tested here; a negative value converted to `unsigned char` is a range question and belongs to Task 6.

`testkeymap.c` already includes `<inttypes.h>` for `PRId64`; confirm with `grep -n inttypes ast_tester/testkeymap.c`.

- [ ] **Step 2: Run the test to verify it fails**

```bash
cmake --build build -j && ctest --test-dir build -R testkeymap --output-on-failure
```

Expected: FAIL, printing lines such as `Negint: got -1 want -2`, `Negd: got -1999 want -2000`, `Neghalf: got -2 want -3`, `Negsmall: got 0 want -1`. `Negtiny` must NOT appear.

- [ ] **Step 3: Apply the fix**

Replace the twelve casts in `src/keymap.c`. In the `Double` arm (lines 1953-1965):

```c
      if( out_type == AST__INTTYPE ) {
         if( out ) *( (int *) out ) = (int) round( dval );
      } else if( out_type == AST__SINTTYPE ) {
         if( out ) *( (short int *) out ) = (int) round( dval );
      } else if( out_type == AST__KINTTYPE ) {
         if( out ) *( (int64_t *) out ) = (int64_t) round( dval );
      } else if( out_type == AST__BYTETYPE ) {
         if( out ) *( (unsigned char *) out ) = (int) round( dval );
```

In the `Float` arm (lines 2012-2024), the operand is `fval`:

```c
         if( out ) *( (int *) out ) = (int) round( fval );
         if( out ) *( (short int *) out ) = (int) round( fval );
         if( out ) *( (int64_t *) out ) = (int64_t) round( fval );
         if( out ) *( (unsigned char *) out ) = (int) round( fval );
```

In the four `%lf` fallbacks inside the `String` arms (lines 2069, 2085, 2101, 2117), the operand is `dval`:

```c
               if( out ) *( (int *) out ) = (int) round( dval );
               if( out ) *( (short int *) out ) = (int) round( dval );
               if( out ) *( (int64_t *) out ) = (int64_t) round( dval );
               if( out ) *( (unsigned char *) out ) = (int) round( dval );
```

Confirm `keymap.c` includes `<math.h>`:

```bash
grep -n '#include <math.h>' src/keymap.c
```

If it does not, add it alongside the other system includes.

- [ ] **Step 4: Add the prologue history entry**

In `src/keymap.c`, immediately before the closing `*class--` line:

```
*     8-AUG-2026 (TIMJ):
*        Use round() rather than truncation in ConvertValue, so that
*        negative values convert to the nearest integer rather than
*        being rounded toward zero.
```

- [ ] **Step 5: Run the test to verify it passes**

```bash
cmake --build build -j && ctest --test-dir build -R testkeymap --output-on-failure
```

Expected: PASS.

- [ ] **Step 6: Run the full suite and the sanitizer build**

```bash
ctest --test-dir build --output-on-failure
cmake --build build-dev -j && ctest --test-dir build-dev --output-on-failure
```

- [ ] **Step 7: Commit**

```bash
git add src/keymap.c ast_tester/testkeymap.c
git commit -F - <<'EOF'
keymap: round negative values to nearest in ConvertValue

ConvertValue converted a stored Double, Float or numeric String to an
integer type with (int)( dval + 0.5 ). A cast from floating point to
integer truncates toward zero, so this expression is correct only for
non-negative values.

Every negative value was rounded toward zero by one, exact integers
included: astMapGet0I on an entry holding -2.0 returned -1.

Use round() instead. The existing coverage used 1999.9, a positive
value that rounds identically either way, so the defect was not
detected; the new checks cover an exact negative integer, values either
side of a tie, and a control that must not change.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
EOF
```

---

## Task 3: Remaining rounding sites

The 57 sites where the input cannot be negative. The substitution changes no behaviour; it is made so that one rounding idiom remains in the library instead of two, and so a reader need not re-derive the sign argument at each site.

**Files:** `src/axis.c` (1), `src/fitschan.c` (12 remaining), `src/grf3d_pgplot.c` (3), `src/grf3d_plplot.c` (3), `src/grf_log.c` (3), `src/grf_pgplot.c` (3), `src/mapping.c` (1), `src/mathmap.c` (4), `src/matrixmap.c` (2), `src/moc.c` (12), `src/pcdmap.c` (2), `src/polymap.c` (1), `src/wcsmap.c` (2), `src/winmap.c` (2), `src/yamlchan.c` (6)

**Interfaces:**
- Consumes: nothing.
- Produces: nothing.

- [ ] **Step 1: Enumerate the sites**

```bash
grep -nE '\((int|long|short int|unsigned char|int64_t|AstDim|size_t|dim_t)\)[[:space:]]*\([^;]*\+[[:space:]]*0\.5[[:space:]]*\)' src/*.c \
  | grep -v '^\S*:[0-9]*: *\*'
```

The second `grep -v` drops prologue comment lines, which match the pattern but are documentation. Expect 57 results at this point (70 total, minus the 1 fixed in Task 1 and the 12 fixed in Task 2).

- [ ] **Step 2: Convert each site**

For each result, rewrite `(T)( EXPR + 0.5 )` as `(T)round( EXPR )`, preserving the cast type `T` exactly. Two shapes need care:

`src/moc.c:2188` and its neighbours have trailing arithmetic outside the cast, which must be left alone:

```c
   glbnd_min[ 0 ] = (AstDim)( outa[ 0 ] + 0.5 ) - 1;   /* before */
   glbnd_min[ 0 ] = (AstDim)round( outa[ 0 ] ) - 1;    /* after  */
```

`src/mathmap.c:2940` likewise:

```c
   nstack -= (int) ( con[ icon++ ] + 0.5 ) - 1;   /* before */
   nstack -= (int) round( con[ icon++ ] ) - 1;    /* after  */
```

Note the `icon++` side effect: it must remain inside the cast and must not be duplicated.

Do NOT touch any other `+ 0.5` in these files. Those are genuine arithmetic — FITS pixel-centre conventions and midpoint calculations — not rounding.

- [ ] **Step 3: Confirm math.h is available in each file**

```bash
for f in src/axis.c src/fitschan.c src/grf3d_pgplot.c src/grf3d_plplot.c \
         src/grf_log.c src/grf_pgplot.c src/mapping.c src/mathmap.c \
         src/matrixmap.c src/moc.c src/pcdmap.c src/polymap.c \
         src/wcsmap.c src/winmap.c src/yamlchan.c; do
  grep -q '#include <math.h>' "$f" || echo "MISSING math.h: $f"
done
```

Add the include to any file the command names, alongside the other system includes.

- [ ] **Step 4: Verify no sites remain**

```bash
grep -nE '\((int|long|short int|unsigned char|int64_t|AstDim|size_t|dim_t)\)[[:space:]]*\([^;]*\+[[:space:]]*0\.5[[:space:]]*\)' src/*.c \
  | grep -v '^\S*:[0-9]*: *\*'
```

Expected: no output.

- [ ] **Step 5: Add prologue history entries**

Add an entry to each of the 15 modified files, immediately before the closing `*class--` or `*--` line:

```
*     8-AUG-2026 (TIMJ):
*        Use round() rather than (int)(x+0.5) for rounding, so that the
*        library uses a single rounding idiom that is correct for
*        negative values.
```

- [ ] **Step 6: Run the full suite and the sanitizer build**

```bash
cmake --build build -j && ctest --test-dir build --output-on-failure
cmake --build build-dev -j && ctest --test-dir build-dev --output-on-failure
```

Expected: no failures and NO fixture changes. The inputs at these sites cannot be negative, so the substitution is a no-op.

**If any test fails or any reference output changes, STOP.** That identifies a site that was rounding a negative value incorrectly — a finding to report to the user, not a diff to absorb. Report which site, what changed, and why the input can be negative there.

- [ ] **Step 7: Commit**

```bash
git add src/
git commit -F - <<'EOF'
Use round() consistently for floating point to integer rounding

The idiom (int)( x + 0.5 ) rounds toward zero for negative values,
because a cast from floating point to integer truncates. It is correct
only when the operand cannot be negative.

Convert the remaining occurrences to round(). The operands at these
sites are structurally non-negative -- axis indices, counts, polynomial
powers, colour indices and the like -- so behaviour is unchanged. The
value of the change is that one rounding idiom now remains in the
library instead of two, so a reader need not work out whether a given
+ 0.5 was chosen deliberately or merely inherited, and the remaining
occurrences of + 0.5 are all genuine arithmetic.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
EOF
```

---

## Task 4: astMapIterate never terminates when SortBy is set

**Files:**
- Modify: `src/keymap.c:8703-8706`
- Test: `ast_tester/testkeymap.c`

**Interfaces:**
- Consumes: nothing.
- Produces: nothing.

**Background.** `AddToSortedList` (`keymap.c:802`) builds a circular doubly linked list; `snext` is never NULL for an entry in the list. `MapKey` guards for the wrap at `keymap.c:3526-3533`. `MapIterate`'s sorted branch does not, and its only termination signal is `key == NULL`, so a caller's loop never ends.

- [ ] **Step 1: Write the failing test**

Add this function to `ast_tester/testkeymap.c`, immediately above `int main( void )`:

```c
/*
 * testiterate: astMapIterate must terminate for a KeyMap with SortBy set.
 * The sorted list is circular, so a walk that does not detect the wrap
 * runs forever. The loop below is bounded so a regression fails the test
 * rather than hanging the suite.
 */
static void testiterate( int *status ) {
   AstKeyMap *km;
   const char *key;
   int nseen;

   if( !astOK ) return;

   km = astKeyMap( "SortBy=KeyUp" );
   astMapPut0I( km, "AAA", 1, " " );
   astMapPut0I( km, "BBB", 2, " " );
   astMapPut0I( km, "CCC", 3, " " );

   nseen = 0;
   key = astMapIterate( km, 1 );
   while( key && astOK && nseen < 100 ) {
      nseen++;
      key = astMapIterate( km, 0 );
   }

   if( nseen != 3 ) {
      printf( "astMapIterate yielded %d keys, expected 3\n", nseen );
      stopit( status, "Error iterate-sorted" );
   }

   km = astAnnul( km );
}
```

Add the call inside `main`, immediately after `astBegin;`:

```c
   testiterate( status );
```

Check the exact signature before writing the call:

```bash
grep -n 'astMapIterate' build/ast.h
```

- [ ] **Step 2: Run the test to verify it fails**

```bash
cmake --build build -j && ctest --test-dir build -R testkeymap --output-on-failure
```

Expected: FAIL with `astMapIterate yielded 100 keys, expected 3`. The bound of 100 is what stops this hanging the suite.

- [ ] **Step 3: Apply the fix**

In `src/keymap.c`, the sorted branch of `MapIterate` currently reads:

```c
      if( entry ) {
         key = entry->key;
         this->iter_entry = entry->snext;
      }
```

Change it to:

```c
      if( entry ) {
         key = entry->key;

/* The sorted list is circular, so stop when the walk returns to the
   first entry rather than following the link back round. */
         this->iter_entry = ( entry->snext == this->first ) ? NULL :
                                                             entry->snext;
      }
```

Do NOT change the `SORTBY_NONE` branch above it. That branch walks the per-bucket `next` chains, where NULL genuinely marks the end, and is correct as written.

- [ ] **Step 4: Add the prologue history entry**

```
*     8-AUG-2026 (TIMJ):
*        Terminate astMapIterate correctly when SortBy is set. The
*        sorted list is circular, so the walk needs the same wrap guard
*        that astMapKey already uses.
```

- [ ] **Step 5: Run the test to verify it passes**

```bash
cmake --build build -j && ctest --test-dir build -R testkeymap --output-on-failure
```

Expected: PASS.

- [ ] **Step 6: Run the full suite and the sanitizer build**

```bash
ctest --test-dir build --output-on-failure
cmake --build build-dev -j && ctest --test-dir build-dev --output-on-failure
```

- [ ] **Step 7: Commit**

```bash
git add src/keymap.c ast_tester/testkeymap.c
git commit -F - <<'EOF'
keymap: terminate astMapIterate when SortBy is set

AddToSortedList builds a circular doubly linked list, so an entry's
snext pointer is never NULL. astMapKey accounts for this and stops when
the walk returns to the first entry; the sorted branch of MapIterate did
not.

MapIterate signals the end of iteration by returning a NULL key, and the
key was taken from an entry that is never NULL, so a caller's loop ran
forever, re-yielding the same keys in order.

Add the wrap guard astMapKey already uses. The unsorted branch is
unchanged: it walks the per-bucket chains, where NULL does mark the end.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
EOF
```

---

## Task 5: astMapRename on a locked KeyMap loses the entry

**Files:**
- Modify: `src/keymap.c:7688-7742`
- Test: `ast_tester/testkeymap.c`

**Interfaces:**
- Consumes: nothing.
- Produces: nothing.

**Background.** `MapRename` removes the entry from the table before it checks `MapLocked`. `AddTableEntry`, the only call that returns the entry to the table, is reached only `if( astOK )`. The `AST__BADKEY` branch sets the error status, so the entry falls through to `FreeMapEntry` and the KeyMap is left holding neither key.

- [ ] **Step 1: Write the failing test**

Add to `ast_tester/testkeymap.c`, immediately above `int main( void )`:

```c
/*
 * testrenamelocked: a rename refused because the KeyMap is locked must
 * leave the KeyMap unchanged, not consume the entry.
 */
static void testrenamelocked( int *status ) {
   AstKeyMap *km;
   int ival;
   int local_status;

   if( !astOK ) return;

   km = astKeyMap( " " );
   astMapPut0I( km, "KEEP", 42, " " );
   astSetI( km, "MapLocked", 1 );

   /* This rename must fail: NEWKEY is not already a known key. */
   local_status = 0;
   astWatch( &local_status );
   astMapRename( km, "KEEP", "NEWKEY" );
   if( local_status == 0 ) {
      stopit( status, "Error rename-locked-noerr" );
   }
   astClearStatus;
   astWatch( status );

   /* The entry must still be present under its original key. */
   if( !astMapHasKey( km, "KEEP" ) ) {
      stopit( status, "Error rename-locked-lost" );
   } else if( !astMapGet0I( km, "KEEP", &ival ) ) {
      stopit( status, "Error rename-locked-unreadable" );
   } else if( ival != 42 ) {
      printf( "got %d want 42\n", ival );
      stopit( status, "Error rename-locked-value" );
   }

   if( astMapHasKey( km, "NEWKEY" ) ) {
      stopit( status, "Error rename-locked-newkey" );
   }

   km = astAnnul( km );
}
```

Add the call inside `main`, immediately after `astBegin;`:

```c
   testrenamelocked( status );
```

- [ ] **Step 2: Run the test to verify it fails**

```bash
cmake --build build -j && ctest --test-dir build -R testkeymap --output-on-failure
```

Expected: FAIL with `Error rename-locked-lost` — the entry is gone.

- [ ] **Step 3: Apply the fix**

The current structure removes first and checks later. Restructure so the check happens before any mutation. In `src/keymap.c`, inside `MapRename`, replace the region that begins:

```c
/* Search the relevent table entry for the required MapEntry. Remove it
   from the list, but do not free it. */
      entry = RemoveTableEntry( this, itab, oldkey, status );
```

with:

```c
/* If the KeyMap is locked, renaming an entry to a key that is not
   already present would introduce a new key, which is not allowed.
   Check this before removing anything, so that a refused rename leaves
   the KeyMap unchanged. */
      if( astGetMapLocked( this ) ) {
         int newtab = HashFun( newkey, this->mapsize - 1, &hash, status );
         if( SearchTableEntry( this, itab, oldkey, status ) &&
             !SearchTableEntry( this, newtab, newkey, status ) ) {
            astError( AST__BADKEY, "astMapRename(%s): Failed to rename item "
                      "\"%s\" in a KeyMap to \"%s\": \"%s\" is not a known "
                      "item.", status, astGetClass( this ), oldkey, newkey,
                      newkey );
         }
      }

/* Search the relevent table entry for the required MapEntry. Remove it
   from the list, but do not free it. */
      entry = astOK ? RemoveTableEntry( this, itab, oldkey, status ) : NULL;
```

Then delete the now-redundant `MapLocked` check further down, the block reading:

```c
/* If the KeyMap is locked we report an error if an attempt is made to
   introduce a new key. */
         if( !there && astGetMapLocked( this ) ) {
            astError( AST__BADKEY, "astMapRename(%s): Failed to rename item "
                      "\"%s\" in a KeyMap to \"%s\": \"%s\" is not a known "
                      "item.", status, astGetClass( this ), oldkey, newkey,
                      newkey );
         }
```

Leave the `if( astOK ) ... else ... FreeMapEntry` structure below it untouched: it is still the right response to an error arising anywhere else.

Check `SearchTableEntry`'s signature before use:

```bash
grep -n 'static AstMapEntry \*SearchTableEntry' src/keymap.c
```

- [ ] **Step 4: Add the prologue history entry**

```
*     8-AUG-2026 (TIMJ):
*        Check MapLocked before removing the entry in astMapRename, so
*        that a refused rename leaves the KeyMap unchanged rather than
*        discarding the entry.
```

- [ ] **Step 5: Run the test to verify it passes**

```bash
cmake --build build -j && ctest --test-dir build -R testkeymap --output-on-failure
```

- [ ] **Step 6: Run the full suite and the sanitizer build**

```bash
ctest --test-dir build --output-on-failure
cmake --build build-dev -j && ctest --test-dir build-dev --output-on-failure
```

- [ ] **Step 7: Commit**

```bash
git add src/keymap.c ast_tester/testkeymap.c
git commit -F - <<'EOF'
keymap: do not discard the entry when astMapRename is refused

MapRename removed the entry from the hash table before testing whether
the KeyMap was locked. AddTableEntry, which puts the entry back, is
reached only when the status is still good, so the AST__BADKEY branch
fell through to FreeMapEntry instead. A refused rename left the KeyMap
holding neither the old key nor the new one, and the value was lost.

Move the MapLocked test ahead of the removal. Whether the rename would
introduce a new key can be determined without mutating anything, so the
error is now reported before the entry is touched.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
EOF
```

---

## Task 6: Out-of-range floating point to integer conversion

**Files:**
- Modify: `src/keymap.c` (helpers plus the twelve arms corrected in Task 2)
- Test: `ast_tester/testkeymap.c`

**Interfaces:**
- Consumes: the `round()`-corrected arms from Task 2.
- Produces: static helpers `DtoI`, `DtoS`, `DtoB`, `DtoK` used only within `keymap.c`.

**Background.** C11 6.3.1.4p1 makes conversion of a floating value to an integer type that cannot represent it undefined behaviour. On the toolchain this was observed on, `5000000100.0` converts to `INT_MIN` for `int` and to its low bytes for `short` and `unsigned char`, but nothing about that is guaranteed.

- [ ] **Step 1: Write the failing test**

Add to `ast_tester/testkeymap.c`, immediately above `int main( void )`:

```c
/*
 * testconvrange: converting an out-of-range Double to a narrow integer
 * type must saturate at the type's limit rather than producing an
 * undefined result.
 */
static void testconvrange( int *status ) {
   AstKeyMap *km;
   int ival;
   short int sval;
   unsigned char bval;

   if( !astOK ) return;

   km = astKeyMap( " " );
   astMapPut0D( km, "BIG", 5000000100.0, " " );
   astMapPut0D( km, "NEG", -5000000100.0, " " );

   if( !astMapGet0I( km, "BIG", &ival ) ) {
      stopit( status, "Error range-get-i" );
   } else if( ival != INT_MAX ) {
      printf( "int: got %d want %d\n", ival, INT_MAX );
      stopit( status, "Error range-int-max" );
   }

   if( !astMapGet0I( km, "NEG", &ival ) ) {
      stopit( status, "Error range-get-i2" );
   } else if( ival != INT_MIN ) {
      printf( "int: got %d want %d\n", ival, INT_MIN );
      stopit( status, "Error range-int-min" );
   }

   if( !astMapGet0S( km, "BIG", &sval ) ) {
      stopit( status, "Error range-get-s" );
   } else if( sval != SHRT_MAX ) {
      printf( "short: got %d want %d\n", (int) sval, SHRT_MAX );
      stopit( status, "Error range-short-max" );
   }

   if( !astMapGet0B( km, "BIG", &bval ) ) {
      stopit( status, "Error range-get-b" );
   } else if( bval != UCHAR_MAX ) {
      printf( "byte: got %d want %d\n", (int) bval, UCHAR_MAX );
      stopit( status, "Error range-byte-max" );
   }

   /* A negative value must saturate at zero for the unsigned type, not
      convert to a large unsigned value. */
   if( !astMapGet0B( km, "NEG", &bval ) ) {
      stopit( status, "Error range-get-b2" );
   } else if( bval != 0 ) {
      printf( "byte: got %d want 0\n", (int) bval );
      stopit( status, "Error range-byte-min" );
   }

   km = astAnnul( km );
}
```

Add the call inside `main`, immediately after `astBegin;`:

```c
   testconvrange( status );
```

`testkeymap.c` already includes `<limits.h>` (it defines `VAL__MAXI` from `INT_MAX` at line 51). Confirm `SHRT_MAX` and `UCHAR_MAX` are therefore available.

- [ ] **Step 2: Run the test to verify it fails**

```bash
cmake --build build -j && ctest --test-dir build -R testkeymap --output-on-failure
```

Expected: FAIL. The exact wrong values are toolchain-dependent — that is the point of the defect — so do not assert on them. Confirm only that the test fails before the fix and passes after.

The sanitizer build is the stronger signal here:

```bash
cmake --build build-dev -j && ctest --test-dir build-dev -R testkeymap --output-on-failure
```

Expected: UndefinedBehaviorSanitizer reports a runtime error at the `keymap.c` conversion line.

- [ ] **Step 3: Apply the fix**

Add these helpers to `src/keymap.c`, above `ConvertValue`:

```c
/* Convert a double to an integer type, saturating rather than relying on
   the undefined result of an out-of-range floating point to integer
   cast (C11 6.3.1.4p1). */

static int DtoI( double dval ) {
   double r = round( dval );
   return ( r >= (double) INT_MAX ) ? INT_MAX :
          ( r <= (double) INT_MIN ) ? INT_MIN : (int) r;
}

static short int DtoS( double dval ) {
   double r = round( dval );
   return ( r >= (double) SHRT_MAX ) ? SHRT_MAX :
          ( r <= (double) SHRT_MIN ) ? SHRT_MIN : (short int) r;
}

static unsigned char DtoB( double dval ) {
   double r = round( dval );
   return ( r >= (double) UCHAR_MAX ) ? UCHAR_MAX :
          ( r <= 0.0 ) ? 0 : (unsigned char) r;
}

static int64_t DtoK( double dval ) {
   double r = round( dval );

/* INT64_MAX is not exactly representable as a double, so compare against
   2^63 directly rather than against (double) INT64_MAX. */
   return ( r >= 9223372036854775808.0 ) ? INT64_MAX :
          ( r <= -9223372036854775808.0 ) ? INT64_MIN : (int64_t) r;
}
```

`INT_MAX`, `INT_MIN`, `SHRT_MAX`, `SHRT_MIN` and `UCHAR_MAX` need `<limits.h>`; `INT64_MAX` and `INT64_MIN` need `<stdint.h>`. Confirm both:

```bash
grep -nE '#include <(limits|stdint)\.h>' src/keymap.c
```

Now route the twelve arms through the helpers. `Double` arm:

```c
      if( out_type == AST__INTTYPE ) {
         if( out ) *( (int *) out ) = DtoI( dval );
      } else if( out_type == AST__SINTTYPE ) {
         if( out ) *( (short int *) out ) = DtoS( dval );
      } else if( out_type == AST__KINTTYPE ) {
         if( out ) *( (int64_t *) out ) = DtoK( dval );
      } else if( out_type == AST__BYTETYPE ) {
         if( out ) *( (unsigned char *) out ) = DtoB( dval );
```

`Float` arm — pass `(double) fval` so the same helpers serve:

```c
         if( out ) *( (int *) out ) = DtoI( (double) fval );
         if( out ) *( (short int *) out ) = DtoS( (double) fval );
         if( out ) *( (int64_t *) out ) = DtoK( (double) fval );
         if( out ) *( (unsigned char *) out ) = DtoB( (double) fval );
```

The four `%lf` fallbacks in the `String` arms:

```c
               if( out ) *( (int *) out ) = DtoI( dval );
               if( out ) *( (short int *) out ) = DtoS( dval );
               if( out ) *( (int64_t *) out ) = DtoK( dval );
               if( out ) *( (unsigned char *) out ) = DtoB( dval );
```

- [ ] **Step 4: Add the prologue history entry**

```
*     8-AUG-2026 (TIMJ):
*        Range check before converting a floating point value to an
*        integer type in ConvertValue, so that an out of range value
*        saturates rather than producing an undefined result.
```

- [ ] **Step 5: Run the test to verify it passes**

```bash
cmake --build build -j && ctest --test-dir build -R testkeymap --output-on-failure
cmake --build build-dev -j && ctest --test-dir build-dev -R testkeymap --output-on-failure
```

Expected: PASS in both, with no UndefinedBehaviorSanitizer report.

- [ ] **Step 6: Run the full suite and the sanitizer build**

```bash
ctest --test-dir build --output-on-failure
ctest --test-dir build-dev --output-on-failure
```

- [ ] **Step 7: Commit**

```bash
git add src/keymap.c ast_tester/testkeymap.c
git commit -F - <<'EOF'
keymap: range check floating point to integer conversion

ConvertValue cast a stored Double or Float straight to the requested
integer type. C11 6.3.1.4p1 makes that undefined when the value is
outside the destination's range, so astMapGet0I on an entry holding a
value too large for int returned whatever the toolchain happened to
produce, with nothing in the call signalling that it had happened.

Convert through helpers that saturate at the destination's limits. A
negative value now saturates at zero for the unsigned type rather than
converting to a large unsigned value. Every in-range value converts to
the same integer as before.

INT64_MAX is not exactly representable as a double, so the 64 bit helper
compares against 2^63 rather than against a cast of the limit.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
EOF
```

---

## Task 7: Numeric string overflow through %d

**Files:**
- Modify: `src/keymap.c:2060-2121`
- Test: `ast_tester/testkeymap.c`

**Interfaces:**
- Consumes: the `DtoI`/`DtoS`/`DtoB`/`DtoK` helpers from Task 6.
- Produces: nothing.

**Background.** C11 7.21.6.2p10 makes a `sscanf` numeric conversion undefined when the value is outside the destination's range. glibc does not fail on overflow; it wraps the accumulated value and reports success, so the acceptance test `nval == 1 && nc >= strlen( cval )` passes and the existing `%lf` fallback is never reached.

- [ ] **Step 1: Write the failing test**

Add to `ast_tester/testkeymap.c`, immediately above `int main( void )`:

```c
/*
 * testconvstring: a numeric string too large for int must not convert to
 * a wrapped value. The %lf fallback in ConvertValue handles it, and its
 * result saturates.
 */
static void testconvstring( int *status ) {
   AstKeyMap *km;
   int ival;

   if( !astOK ) return;

   km = astKeyMap( " " );
   astMapPut0C( km, "OVER", "9999999999", " " );
   astMapPut0C( km, "HUGE", "99999999999999999999999999999999999999", " " );
   astMapPut0C( km, "OK", "42", " " );
   astMapPut0C( km, "OKSPACE", " 42 ", " " );
   astMapPut0C( km, "BLANK", "   ", " " );
   astMapPut0C( km, "TRAIL", "12abc", " " );

   /* In range values, with and without surrounding space, are unchanged. */
   if( !astMapGet0I( km, "OK", &ival ) || ival != 42 ) {
      stopit( status, "Error convstr-ok" );
   }
   if( !astMapGet0I( km, "OKSPACE", &ival ) || ival != 42 ) {
      stopit( status, "Error convstr-okspace" );
   }

   /* Non numeric input is still refused. */
   if( astMapGet0I( km, "BLANK", &ival ) ) {
      stopit( status, "Error convstr-blank" );
   }
   if( astMapGet0I( km, "TRAIL", &ival ) ) {
      stopit( status, "Error convstr-trail" );
   }

   /* Magnitude overflow must saturate, not wrap. */
   if( !astMapGet0I( km, "OVER", &ival ) ) {
      stopit( status, "Error convstr-over-get" );
   } else if( ival != INT_MAX ) {
      printf( "over: got %d want %d\n", ival, INT_MAX );
      stopit( status, "Error convstr-over" );
   }

   if( !astMapGet0I( km, "HUGE", &ival ) ) {
      stopit( status, "Error convstr-huge-get" );
   } else if( ival != INT_MAX ) {
      printf( "huge: got %d want %d\n", ival, INT_MAX );
      stopit( status, "Error convstr-huge" );
   }

   km = astAnnul( km );
}
```

Add the call inside `main`, immediately after `astBegin;`:

```c
   testconvstring( status );
```

- [ ] **Step 2: Run the test to verify it fails**

```bash
cmake --build build -j && ctest --test-dir build -R testkeymap --output-on-failure
```

Expected: FAIL with `over: got 1410065407 want 2147483647` and a similar line for `huge`. The in-range and refusal cases must already pass — if they do not, the test is wrong, not the library.

- [ ] **Step 3: Apply the fix**

In the `String`-to-`Int` arm of `ConvertValue`, replace:

```c
      if( out_type == AST__INTTYPE ) {
         nc = 0;
         nval = astSscanf( cval, " %d %n", &ival, &nc );
         if( ( nval == 1 ) && ( nc >= (int) strlen( cval ) ) ) {
            if( out ) *( (int *) out ) = ival;
         } else {
```

with:

```c
      if( out_type == AST__INTTYPE ) {
         char *end;
         long lval;
         int consumed;

         errno = 0;
         lval = strtol( cval, &end, 10 );

/* Record whether any digits were consumed before skipping trailing
   space. Testing this after the skip would accept an all blank string,
   because the skip would advance "end" away from "cval" on its own. */
         consumed = ( end != cval );
         while( isspace( (int) *end ) ) end++;

         if( consumed && errno != ERANGE && *end == '\0' &&
             lval >= INT_MIN && lval <= INT_MAX ) {
            if( out ) *( (int *) out ) = (int) lval;
         } else {
```

Leave the `else` branch, which performs the `%lf` scan, exactly as it is. It now receives the overflow cases, and its `DtoI( dval )` call from Task 6 saturates them.

Apply the same shape to the three sibling arms:

- `AST__SINTTYPE`: `strtol`, bounds `SHRT_MIN`/`SHRT_MAX`, store `(short int) lval`.
- `AST__BYTETYPE`: `strtol`, bounds `0`/`UCHAR_MAX`, store `(unsigned char) lval`.
- `AST__KINTTYPE`: `strtoll` into a `long long`, bounds `INT64_MIN`/`INT64_MAX`, store `(int64_t) lval`.

Declare the block locals inside each arm's braces so the four arms do not collide.

`strtol`/`strtoll` need `<stdlib.h>`, `errno` needs `<errno.h>`, `isspace` needs `<ctype.h>`. Confirm all three:

```bash
grep -nE '#include <(stdlib|errno|ctype)\.h>' src/keymap.c
```

`errno`-based range detection is already used for this purpose in `axis.c`, `channel.c`, `fitschan.c`, `mathmap.c`, `memory.c` and `object.c`, so it is consistent with the codebase.

- [ ] **Step 4: Add the prologue history entry**

```
*     8-AUG-2026 (TIMJ):
*        Detect magnitude overflow when converting a numeric string to an
*        integer in ConvertValue, using strtol and errno rather than a
*        sscanf conversion whose overflow behaviour is undefined.
```

- [ ] **Step 5: Run the test to verify it passes**

```bash
cmake --build build -j && ctest --test-dir build -R testkeymap --output-on-failure
```

- [ ] **Step 6: Run the full suite and the sanitizer build**

```bash
ctest --test-dir build --output-on-failure
cmake --build build-dev -j && ctest --test-dir build-dev --output-on-failure
```

- [ ] **Step 7: Commit**

```bash
git add src/keymap.c ast_tester/testkeymap.c
git commit -F - <<'EOF'
keymap: detect numeric string overflow in ConvertValue

The String to integer arms of ConvertValue scanned with " %d %n" and
accepted the result when the whole string had been consumed. C11
7.21.6.2p10 makes a sscanf numeric conversion undefined when the value
is out of range, and glibc does not fail: it wraps the accumulated value
and reports success. The acceptance test therefore passed and the
existing "%lf" fallback was never reached, so "9999999999" converted to
1410065407.

Use strtol with an errno range check instead, which this codebase
already relies on in six other source files. Out of range strings now
fall through to the "%lf" attempt exactly as any other unparsable string
does.

The check for consumed digits is captured before the trailing whitespace
skip. Testing it afterwards would accept an all blank string, because
the skip advances the end pointer away from the start on its own.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
EOF
```

---

## Task 8: Zero-length vectors dump pointer bytes

**Files:**
- Modify: `src/keymap.c` (`MAKE_MAPPUT1` at 4894, `DumpEntry` String branch at 2880-2883)
- Test: `ast_tester/testkeymap.c`

**Interfaces:**
- Consumes: nothing.
- Produces: nothing.

**Background.** `AstMapEntry.nel` is documented as "0 => scalar, >0 => array" (`keymap.h:116`), so there is no encoding for a zero-length vector. `MAKE_MAPPUT1` applies no check to `size` because `CHECK_I` and its siblings are empty macros. A put with `size == 0` stores `nel = 0` and a NULL value pointer, and `DumpEntry` then takes the scalar branch and aliases the structure, handing that NULL to `astWriteString` — which dereferences it.

**This is the one change in the plan that alters a public API contract.** Rejecting the call is preferred to silently degrading it, because a caller passing zero has a bug that should be reported.

- [ ] **Step 1: Confirm the error code**

Use `AST__NELIN`, whose message text in `ast_err.msg` is "number of array elements invalid" — the exact condition. Confirm it exists:

```bash
grep -n 'NELIN' ast_err.msg
```

Do NOT use `AST__MPVIN`: its text is "invalid integer index supplied for a KeyMap vector element", which describes a bad index, not a bad length. Do NOT add a new error code — that changes the generated `ast_err.h` and is beyond this defect.

- [ ] **Step 2: Write the failing test**

Add to `ast_tester/testkeymap.c`, immediately above `int main( void )`:

```c
/*
 * testemptyvector: a zero length vector has no representation in a
 * KeyMap entry, where nel == 0 means scalar. The put must be refused
 * rather than creating an entry that cannot be dumped.
 */
static void testemptyvector( int *status ) {
   AstKeyMap *km;
   const char *cvec[ 1 ];
   int ivec[ 1 ];
   int local_status;
   char *dump;

   if( !astOK ) return;

   km = astKeyMap( " " );

   local_status = 0;
   astWatch( &local_status );
   astMapPut1I( km, "EMPTYI", 0, ivec, " " );
   if( local_status == 0 ) {
      stopit( status, "Error emptyvec-int-noerr" );
   }
   astClearStatus;

   local_status = 0;
   cvec[ 0 ] = "unused";
   astMapPut1C( km, "EMPTYC", 0, cvec, " " );
   if( local_status == 0 ) {
      stopit( status, "Error emptyvec-str-noerr" );
   }
   astClearStatus;
   astWatch( status );

   /* Neither key may have been created, and the KeyMap must still dump. */
   if( astMapHasKey( km, "EMPTYI" ) ) {
      stopit( status, "Error emptyvec-int-present" );
   }
   if( astMapHasKey( km, "EMPTYC" ) ) {
      stopit( status, "Error emptyvec-str-present" );
   }

   astMapPut0I( km, "REAL", 7, " " );
   dump = astToString( km );
   if( !dump ) {
      stopit( status, "Error emptyvec-dump" );
   } else {
      dump = astFree( dump );
   }

   km = astAnnul( km );
}
```

Add the call inside `main`, immediately after `astBegin;`:

```c
   testemptyvector( status );
```

- [ ] **Step 3: Run the test to verify it fails**

```bash
cmake --build build -j && ctest --test-dir build -R testkeymap --output-on-failure
```

Expected: FAIL. The `astMapPut1C` case will most likely **segfault** rather than report a failure, because the dump dereferences NULL. A crash here is the defect, confirmed. If it crashes, note that and continue to Step 4.

- [ ] **Step 4: Apply the fix**

In `MAKE_MAPPUT1` (`src/keymap.c:4894`), add a size check immediately after the `if ( !astOK ) return;` line and before the `CHECK_##X` invocation:

```c
/* A KeyMap entry has no representation for a zero length vector: nel of
   zero means the entry is a scalar. Reject the call rather than creating
   an entry that cannot be read back or dumped. */ \
   if( size < 1 ) { \
      astError( AST__NELIN, "astMapPut1" #X "(%s): Illegal vector length " \
                "%d supplied for KeyMap entry \"%s\" - must be at least " \
                "one.", status, astGetClass( this ), size, skey ); \
      return; \
   } \
```

Then add the guard to `DumpEntry`'s String scalar branch, matching the guard the Object branch already has at `keymap.c:2895`. Change:

```c
   } else if( type == AST__STRINGTYPE ) {
      if( entry->nel == 0 ) {
         (void) sprintf( buff, "Val%d", nentry );
         astWriteString( channel, buff, 1, 1, ((Entry0C *)entry)->value, "Item value" );
```

to:

```c
   } else if( type == AST__STRINGTYPE ) {
      if( entry->nel == 0 ) {
         if( ((Entry0C *)entry)->value ) {
            (void) sprintf( buff, "Val%d", nentry );
            astWriteString( channel, buff, 1, 1, ((Entry0C *)entry)->value, "Item value" );
         }
```

and add the matching closing brace. This is defence in depth: with the put refused, no such entry can be created, but the dump should not dereference a null pointer regardless.

- [ ] **Step 5: Add the prologue history entry**

```
*     8-AUG-2026 (TIMJ):
*        Reject a zero length vector in astMapPut1<X>, which has no
*        representation in a KeyMap entry, and guard against a null
*        string pointer when dumping a scalar entry.
```

- [ ] **Step 6: Run the test to verify it passes**

```bash
cmake --build build -j && ctest --test-dir build -R testkeymap --output-on-failure
```

- [ ] **Step 7: Run the full suite and the sanitizer build**

```bash
ctest --test-dir build --output-on-failure
cmake --build build-dev -j && ctest --test-dir build-dev --output-on-failure
```

Watch for any existing test that passes `size = 0` deliberately. If one exists, report it before committing — it changes the impact assessment of this API change.

- [ ] **Step 8: Commit**

```bash
git add src/keymap.c ast_tester/testkeymap.c
git commit -F - <<'EOF'
keymap: reject zero length vectors in astMapPut1<X>

A KeyMap entry records a vector length of zero to mean the entry is a
scalar, so there is no representation for a zero length vector. The put
macros applied no check to the supplied size, because the per type check
macros are empty for the numeric types, so a size of zero created an
entry claiming to be a scalar whose value pointer was null.

DumpEntry then took the scalar branch and read that pointer as the
payload. For the Object type an existing guard caught it; for the string
type there was no guard, and astWrite dereferenced the null pointer.

Reject the call instead, and add the missing null guard to the string
branch of DumpEntry so the dump path cannot dereference a null pointer
regardless of how an entry was created.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
EOF
```

---

## Task 9: astMapGetElem reports success on an undefined entry

**Files:**
- Modify: `src/keymap.c:6910` (`MAKE_MAPGETELEM`), `src/keymap.c:7116` (`MapGetElemC`)
- Test: `ast_tester/testkeymap.c`

**Interfaces:**
- Consumes: nothing.
- Produces: nothing.

**Background.** `MAKE_MAPGETELEM` sets `result = 1` as soon as the entry is found. The `AST__UNDEFTYPE` arm sets `raw = NULL`, and the conversion is guarded by `} else if( raw )`, so nothing writes the caller's output and nothing resets `result`. `MAKE_MAPGET0` (`keymap.c:5525`) and `MAKE_MAPGET1` (`keymap.c:6014`) both have the `if( !raw ) { result = 0; }` guard, so the three functions answer the same question differently.

- [ ] **Step 1: Write the failing test**

Add to `ast_tester/testkeymap.c`, immediately above `int main( void )`:

```c
/*
 * testundefelem: astMapGetElem<X> must report failure for an undefined
 * entry, as astMapGet0<X> and astMapGet1<X> already do.
 */
static void testundefelem( int *status ) {
   AstKeyMap *km;
   int ival;
   char cbuf[ 200 ];

   if( !astOK ) return;

   km = astKeyMap( " " );
   astMapPutU( km, "UNDEF", " " );

   /* The scalar and vector accessors already report failure. */
   if( astMapGet0I( km, "UNDEF", &ival ) ) {
      stopit( status, "Error undef-get0" );
   }

   /* The element accessor must agree with them. */
   if( astMapGetElemI( km, "UNDEF", 0, &ival ) ) {
      stopit( status, "Error undef-getelem-i" );
   }

   if( astMapGetElemC( km, "UNDEF", sizeof( cbuf ), 0, cbuf ) ) {
      stopit( status, "Error undef-getelem-c" );
   }

   km = astAnnul( km );
}
```

Add the call inside `main`, immediately after `astBegin;`:

```c
   testundefelem( status );
```

Confirm the exact signatures before writing the calls:

```bash
grep -n 'astMapGetElemI\|astMapGetElemC\|astMapPutU' build/ast.h
```

- [ ] **Step 2: Run the test to verify it fails**

```bash
cmake --build build -j && ctest --test-dir build -R testkeymap --output-on-failure
```

Expected: FAIL with `Error undef-getelem-i` and `Error undef-getelem-c`. The `astMapGet0I` check must already pass.

- [ ] **Step 3: Apply the fix**

In `MAKE_MAPGETELEM`, the bounds check and conversion currently read:

```c
      if( elem >= nel || elem < 0 ) { \
         ... astError( AST__MPVIN, ... ); \
\
/* Get a pointer to the requested raw value. */ \
      } else if( raw ) { \
```

Add an arm for the null raw pointer after the `else if( raw )` block's closing brace, so the ordering of the bounds check is preserved:

```c
/* An undefined entry has no value to return. Report failure, as \
   astMapGet0<X> and astMapGet1<X> do for the same entry. */ \
      } else { \
         result = 0; \
      } \
```

Apply the identical change to `MapGetElemC` at the corresponding point (`keymap.c:7116`).

Read both functions fully before editing. The bounds check must continue to fire first, so an out-of-range index on an undefined entry still reports `AST__MPVIN`.

- [ ] **Step 4: Add the prologue history entry**

```
*     8-AUG-2026 (TIMJ):
*        Report failure from astMapGetElem<X> for an undefined entry,
*        matching astMapGet0<X> and astMapGet1<X>.
```

- [ ] **Step 5: Run the test to verify it passes**

```bash
cmake --build build -j && ctest --test-dir build -R testkeymap --output-on-failure
```

- [ ] **Step 6: Run the full suite and the sanitizer build**

```bash
ctest --test-dir build --output-on-failure
cmake --build build-dev -j && ctest --test-dir build-dev --output-on-failure
```

- [ ] **Step 7: Commit**

```bash
git add src/keymap.c ast_tester/testkeymap.c
git commit -F - <<'EOF'
keymap: report failure from astMapGetElem<X> for an undefined entry

MAKE_MAPGETELEM set its result to one as soon as the entry was found.
For an undefined entry the raw value pointer is null and the conversion
step is skipped, so the function returned success without writing
anything to the caller's output buffer. A caller trusting the return
value read whatever was already there.

astMapGet0<X> and astMapGet1<X> both report failure for the same entry,
so three functions answering the same question gave two different
answers. Add the same guard to astMapGetElem<X> and astMapGetElemC.

An out of range index still reports AST__MPVIN, because the bounds check
continues to run first.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
EOF
```

---

## Task 10: LoadKeyMap injects a spurious KyCas card

**Files:**
- Modify: `src/keymap.c:11001`

**Interfaces:**
- Consumes: nothing.
- Produces: nothing.

**Background.** `KeyCase` is the one KeyMap attribute whose unset sentinel is `-1` (`keymap.c:10160-10161`); `KeyError`, `MapLocked` and `SortBy` all use `-INT_MAX`. Loading `KyCas` with a default of `-INT_MAX` makes `TestKeyCase` report an absent attribute as set, so the next dump writes it and no freshly written KeyMap matches its own re-dump.

No new test. The existing `checkdump` round-trips in `testkeymap.c` cover this, which is where the effect appears.

- [ ] **Step 1: Confirm the current behaviour**

Write this to `/tmp/probe.c` and build it per "Writing a Scratch Probe" above:

```c
#include <stdio.h>
#include <string.h>
#include "ast.h"

int main( void ) {
   int status = 0;
   AstKeyMap *km;
   AstKeyMap *km2;
   char *d1;
   char *d2;

   astWatch( &status );
   km = astKeyMap( " " );
   astMapPut0C( km, "ONLY", "x", " " );

   d1 = astToString( km );
   km2 = astFromString( d1 );
   d2 = astToString( km2 );

   printf( "dump1 has KyCas: %s\n", strstr( d1, "KyCas" ) ? "yes" : "no" );
   printf( "dump2 has KyCas: %s\n", strstr( d2, "KyCas" ) ? "yes" : "no" );
   return 0;
}
```

Expected before the fix: `dump1 has KyCas: no`, `dump2 has KyCas: yes`. A single entry is used so that no bucket collision can be the cause.

- [ ] **Step 2: Apply the fix**

In `src/keymap.c:11001`, change:

```c
      new->keycase = astReadInt( channel, "kycas", -INT_MAX );
```

to:

```c
      new->keycase = astReadInt( channel, "kycas", -1 );
```

Do NOT change the sentinel at `keymap.c:10160-10161`. That would touch the public-facing default value expression; the load path default is the smaller change with the same effect.

- [ ] **Step 3: Add the prologue history entry**

```
*     8-AUG-2026 (TIMJ):
*        Use the correct unset value when reading the KyCas card in
*        astLoadKeyMap, so that an absent card is not recorded as an
*        explicitly set KeyCase attribute.
```

- [ ] **Step 4: Verify the round trip is now stable**

Re-run the probe from Step 1. Expected: both dumps report `no`.

- [ ] **Step 5: Run the full suite and the sanitizer build**

```bash
cmake --build build -j && ctest --test-dir build --output-on-failure
cmake --build build-dev -j && ctest --test-dir build-dev --output-on-failure
```

If any reference `.ast` fixture containing a `KyCas` card now differs, review the diff and report it before committing.

- [ ] **Step 6: Commit**

```bash
git add src/keymap.c
git commit -F - <<'EOF'
keymap: do not record an absent KyCas card as a set attribute

KeyCase is the one KeyMap attribute whose unset value is -1; KeyError,
MapLocked and SortBy all use -INT_MAX. astLoadKeyMap read the KyCas card
with a default of -INT_MAX, which is not equal to -1, so TestKeyCase
reported an absent attribute as explicitly set and the next dump wrote
it out. No freshly written KeyMap matched its own re-dump.

Read the card with the correct unset value. The sentinel itself is
unchanged, since altering it would change the attribute's public default
value expression.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
EOF
```

---

## Task 11: MpLck is applied before any entry is loaded

**Files:**
- Modify: `src/keymap.c:11008-11012` and the region after the entry loop
- Test: `ast_tester/testkeymap.c`

**Interfaces:**
- Consumes: nothing.
- Produces: nothing.

**Background.** `astLoadKeyMap_` applies `MapLocked` before the entry loop. Every `MapPut0<X>`/`MapPut1<X>` then reports `AST__BADKEY` for a key not already present, and a KeyMap is empty when loading starts, so the first entry always fails. Dumping a locked, non-empty KeyMap and reading it back is impossible.

`SetMapLocked` (`keymap.c:9381`) does not error on a non-empty KeyMap, and it propagates the flag into nested KeyMaps held as entries (`keymap.c:9426`) — which the current placement cannot do, because those entries do not exist yet when it runs. Moving the call therefore repairs a second latent defect.

- [ ] **Step 1: Write the failing test**

Add to `ast_tester/testkeymap.c`, immediately above `int main( void )`:

```c
/*
 * testlockedroundtrip: a locked, non-empty KeyMap must survive a dump
 * and reload. MapLocked has to be applied after the entries are read,
 * or every put during the load is refused.
 */
static void testlockedroundtrip( int *status ) {
   AstKeyMap *km;
   AstKeyMap *km2;
   AstKeyMap *nested;
   AstObject *obj;
   char *dump;
   int ival;

   if( !astOK ) return;

   km = astKeyMap( " " );
   astMapPut0I( km, "AAA", 1, " " );
   astMapPut0I( km, "BBB", 2, " " );

   nested = astKeyMap( " " );
   astMapPut0I( nested, "INNER", 3, " " );
   astMapPut0A( km, "NEST", nested, " " );
   nested = astAnnul( nested );

   astSetI( km, "MapLocked", 1 );

   dump = astToString( km );
   if( !dump ) {
      stopit( status, "Error locked-dump" );
      km = astAnnul( km );
      return;
   }

   km2 = astFromString( dump );
   dump = astFree( dump );

   if( !astOK || !km2 ) {
      stopit( status, "Error locked-reload" );
      astClearStatus;
   } else {
      if( !astMapGet0I( km2, "AAA", &ival ) || ival != 1 ) {
         stopit( status, "Error locked-reload-aaa" );
      }
      if( !astGetI( km2, "MapLocked" ) ) {
         stopit( status, "Error locked-reload-flag" );
      }

      /* The flag must reach nested KeyMaps too. */
      if( !astMapGet0A( km2, "NEST", &obj ) ) {
         stopit( status, "Error locked-reload-nest" );
      } else {
         if( !astGetI( obj, "MapLocked" ) ) {
            stopit( status, "Error locked-reload-nestflag" );
         }
         obj = astAnnul( obj );
      }
      km2 = astAnnul( km2 );
   }

   km = astAnnul( km );
}
```

Add the call inside `main`, immediately after `astBegin;`:

```c
   testlockedroundtrip( status );
```

- [ ] **Step 2: Run the test to verify it fails**

```bash
cmake --build build -j && ctest --test-dir build -R testkeymap --output-on-failure
```

Expected: FAIL with `Error locked-reload` — `astFromString` returns nothing because the first put was refused.

- [ ] **Step 3: Apply the fix**

In `src/keymap.c`, cut these two lines from their position at 11011-11012:

```c
      new->maplocked = astReadInt( channel, "mplck", -INT_MAX );
      if ( TestMapLocked( new, status ) ) SetMapLocked( new, new->maplocked, status );
```

and paste them immediately after the entry-reading `while` loop terminates, keeping the surrounding comment structure:

```c
/* MapLocked. */
/* --------- */
/* This is read after the entries, because applying it beforehand would
   cause every put made by the entry loop to be refused: the KeyMap is
   empty at that point, so each key is a new key. Applying it here also
   propagates the flag into any nested KeyMaps that were read as
   entries. */
      new->maplocked = astReadInt( channel, "mplck", -INT_MAX );
      if ( TestMapLocked( new, status ) ) SetMapLocked( new, new->maplocked, status );
```

Locate the exact end of the entry loop first:

```bash
grep -n 'memcnt' src/keymap.c
```

Place the moved block adjacent to the `MemCnt` read, which is already after the loop.

- [ ] **Step 4: Add the prologue history entry**

```
*     8-AUG-2026 (TIMJ):
*        Apply the MpLck card after reading the entries in
*        astLoadKeyMap, so that a locked non-empty KeyMap can be read
*        back and so that nested KeyMaps inherit the flag.
```

- [ ] **Step 5: Run the test to verify it passes**

```bash
cmake --build build -j && ctest --test-dir build -R testkeymap --output-on-failure
```

- [ ] **Step 6: Run the full suite and the sanitizer build**

```bash
ctest --test-dir build --output-on-failure
cmake --build build-dev -j && ctest --test-dir build-dev --output-on-failure
```

- [ ] **Step 7: Commit**

```bash
git add src/keymap.c ast_tester/testkeymap.c
git commit -F - <<'EOF'
keymap: apply MpLck after reading entries in astLoadKeyMap

astLoadKeyMap applied every attribute, MapLocked included, before the
entry reading loop. A KeyMap is empty at that point, so every key the
loop then put was a new key, and the shared guard in the put macros
refused it with AST__BADKEY. Dumping a locked, non-empty KeyMap and
reading the result back was not possible: astWrite succeeded and
produced an MpLck card, but the matching astRead failed and returned a
null object.

Read and apply the card after the entries instead. This also propagates
the flag into nested KeyMaps held as entries, which the earlier
placement could not do because those entries did not exist yet.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
EOF
```

---

## Task 12: LoadKeyMap resets the member counter before every put

**Files:**
- Modify: `src/keymap.c:11065-11068`

**Interfaces:**
- Consumes: nothing.
- Produces: nothing.

**Background.** `astLoadKeyMap_` reads a per-entry `mem%d` card and assigns it to `new->member_count` before adding each entry. `DumpEntry` never writes such a card, so `astReadInt`'s default of 0 fires every iteration. `AddTableEntry` assigns `entry->member = (this->member_count)++` (`keymap.c:686`), so the counter is reset to 0 before each put and every loaded entry receives `member = 0`. This destroys the per-entry numbering that `SortBy=AgeUp` and `SortBy=AgeDown` order by.

No new test: the existing `testsorting` function in `testkeymap.c:131` covers sort ordering, and the `checkdump` round-trips cover the dump.

- [ ] **Step 1: Confirm no Mem card is written**

```bash
grep -n '"Mem%d"\|"mem%d"' src/keymap.c
```

Expected: exactly one hit, the read at line 11067. If a write exists, STOP — the analysis is wrong and the read should be kept.

- [ ] **Step 2: Apply the fix**

Delete these four lines from `src/keymap.c` (11065-11068):

```c
/* Get the entry member number. Set the KeyMap member count to this value
   so that the next entry added to the KeyMap will get this value as its
   member index. */
         (void) sprintf( buff, "mem%d", nentry );
         new->member_count = astReadInt( channel, buff, 0 );
```

Leave the `MemCnt` read after the entry loop untouched: it restores the total and is written by `Dump`.

- [ ] **Step 3: Add the prologue history entry**

```
*     8-AUG-2026 (TIMJ):
*        Remove the per entry member number read in astLoadKeyMap. No
*        such card is written, so the read reset the member counter
*        before every put and gave every loaded entry the same member
*        number.
```

- [ ] **Step 4: Run the full suite and the sanitizer build**

```bash
cmake --build build -j && ctest --test-dir build --output-on-failure
cmake --build build-dev -j && ctest --test-dir build-dev --output-on-failure
```

- [ ] **Step 5: Commit**

```bash
git add src/keymap.c
git commit -F - <<'EOF'
keymap: remove the per entry member number read in astLoadKeyMap

astLoadKeyMap read a mem%d card for each entry and used it to seed the
member counter before adding that entry. DumpEntry writes no such card,
so the default of zero was used on every iteration: the counter was
reset before each put and every loaded entry was given member number
zero.

The per entry member numbers determine SortBy=AgeUp and SortBy=AgeDown
ordering, so the read actively destroyed that ordering across a load
rather than merely failing to preserve it.

Remove the read. Entries now receive sequential member numbers in the
order they are read, which is the natural reconstruction. The MemCnt
total, which Dump does write, is still read after the loop.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
EOF
```

---

## Task 13: astMapCopyEntry uses the source KeyMap's KeyCase for both lookups

> **STATUS: attempted, reverted, needs re-planning before it is attempted again.**
>
> The fix below is correct as far as it goes, but it is not sufficient on its
> own, and applying it alone makes matters worse.
>
> `CopyMapEntry` (`src/keymap.c:2436`) copies the source entry's key text
> verbatim: `result->key = astStore( NULL, in->key, strlen( in->key ) + 1 )`.
> It has no knowledge of the destination KeyMap, so a copied entry is stored
> under the *source's* casing.
>
> With only the lookup fixed, the destination is searched under its own
> casing while the entry is stored under the source's. The two disagree, so
> a copied entry becomes unfindable in the map that now holds it. Before the
> fix both were consistently wrong, which is less harmful.
>
> The storage path is also shared: `CopyMapEntry` is called from
> `MapCopyEntry` (twice) and from the bulk `MapCopy` loop (twice), so
> `astMapCopy` between two KeyMaps with different `KeyCase` settings has the
> same defect. Any fix has to cover both callers or neither.
>
> The likely correct shape is to leave `CopyMapEntry` a pure copy helper and
> re-key in the callers, which know the destination — mirroring
> `astMapPut0<X>`, which stores the converted key rather than the caller's
> raw string. That is a larger change than this task specifies and needs its
> own design pass.
>
> Verified: the lookup fix alone leaves `astMapSize( dst ) == 1` (the
> duplicate-entry symptom is genuinely fixed) but the stored key is
> `MixedKey` where the destination's own rule requires `MIXEDKEY`.

**Files:**
- Modify: `src/keymap.c:4350`
- Test: `ast_tester/testkeymap.c`

**Interfaces:**
- Consumes: nothing.
- Produces: nothing.

**Background.** `MapCopyEntry` converts the caller's key using the source KeyMap's `KeyCase` and uses the same converted string for both the source lookup and the destination lookup. When the two KeyMaps differ in `KeyCase`, the destination lookup runs under the source's casing rule.

- [ ] **Step 1: Read the function fully**

```bash
sed -n '4264,4420p' src/keymap.c
```

Identify every use of the `key` variable and which KeyMap each belongs to. The destination is `this`; the source is `that`.

- [ ] **Step 2: Write the failing test**

Add to `ast_tester/testkeymap.c`, immediately above `int main( void )`:

```c
/*
 * testcopyentrycase: astMapCopyEntry must apply each KeyMap's own
 * KeyCase when looking up that KeyMap.
 */
static void testcopyentrycase( int *status ) {
   AstKeyMap *src;
   AstKeyMap *dst;
   int ival;

   if( !astOK ) return;

   /* Source is case sensitive; destination is not. */
   src = astKeyMap( "KeyCase=1" );
   dst = astKeyMap( "KeyCase=0" );

   astMapPut0I( src, "MixedKey", 11, " " );

   /* The destination already holds the key under its own folded form. */
   astMapPut0I( dst, "MIXEDKEY", 99, " " );

   /* Argument order is (destination, key, source, merge). The merge flag
      only affects what happens when the entry holds a KeyMap and the
      destination already has a KeyMap under that key, so it is
      irrelevant to this integer entry. */
   astMapCopyEntry( dst, "MixedKey", src, 0 );

   if( !astOK ) {
      stopit( status, "Error copyentry-status" );
      astClearStatus;
   } else if( !astMapGet0I( dst, "MIXEDKEY", &ival ) ) {
      stopit( status, "Error copyentry-missing" );
   } else if( ival != 11 ) {
      printf( "got %d want 11\n", ival );
      stopit( status, "Error copyentry-value" );
   }

   /* The destination must not have gained a second, differently cased
      copy of the same logical key. */
   if( astMapSize( dst ) != 1 ) {
      printf( "destination has %d entries, expected 1\n", astMapSize( dst ) );
      stopit( status, "Error copyentry-dup" );
   }

   src = astAnnul( src );
   dst = astAnnul( dst );
}
```

Add the call inside `main`, immediately after `astBegin;`:

```c
   testcopyentrycase( status );
```

The signature is `astMapCopyEntry( AstKeyMap *this, const char *key, AstKeyMap *that, int merge )`, where `this` is the destination and `that` is the source. The `merge` flag governs only the case where the entry holds a KeyMap and the destination already holds a KeyMap under that key, so it does not affect this test.

The value assertion (`ival == 11`) and the entry-count assertion (`astMapSize( dst ) == 1`) are both meaningful: a destination lookup performed under the source's casing rule misses the existing `MIXEDKEY` entry and adds a second one, so the count is the primary discriminator.

- [ ] **Step 3: Run the test to verify it fails**

```bash
cmake --build build -j && ctest --test-dir build -R testkeymap --output-on-failure
```

Expected: FAIL, most likely `copyentry-dup`, because the destination lookup used the source's unfolded key and so missed the existing entry.

**If the test passes before the fix, STOP and report.** The lookup paths may coincide for this input, and the test needs a different `KeyCase` combination to expose the defect.

- [ ] **Step 4: Apply the fix**

In `src/keymap.c:4350`, the single conversion currently reads:

```c
   key = ConvertKey( that, skey, keybuf, AST__MXKEYLEN + 1, "astMapCopyEntry",
                     status );
```

Convert once per KeyMap into separate buffers, and use each for its own lookup:

```c
/* Convert the supplied key using each KeyMap's own KeyCase, since the
   two KeyMaps may differ in that attribute. */
   inkey = ConvertKey( that, skey, inkeybuf, AST__MXKEYLEN + 1,
                       "astMapCopyEntry", status );
   outkey = ConvertKey( this, skey, outkeybuf, AST__MXKEYLEN + 1,
                        "astMapCopyEntry", status );
```

Declare `inkey`, `outkey`, `inkeybuf` and `outkeybuf` alongside the existing `key` and `keybuf` declarations, then replace each use of `key`: the `HashFun`/`SearchTableEntry` pair against `that` uses `inkey`; the pair against `this` uses `outkey`. Remove the now-unused `key` and `keybuf`.

- [ ] **Step 5: Add the prologue history entry**

```
*     8-AUG-2026 (TIMJ):
*        Apply each KeyMap's own KeyCase in astMapCopyEntry, rather than
*        using the source KeyMap's casing rule for the destination
*        lookup as well.
```

- [ ] **Step 6: Run the test to verify it passes**

```bash
cmake --build build -j && ctest --test-dir build -R testkeymap --output-on-failure
```

- [ ] **Step 7: Run the full suite and the sanitizer build**

```bash
ctest --test-dir build --output-on-failure
cmake --build build-dev -j && ctest --test-dir build-dev --output-on-failure
```

- [ ] **Step 8: Commit**

```bash
git add src/keymap.c ast_tester/testkeymap.c
git commit -F - <<'EOF'
keymap: apply each KeyMap's own KeyCase in astMapCopyEntry

MapCopyEntry converted the supplied key using the source KeyMap's
KeyCase and then used that one converted string for both the source
lookup and the destination lookup. When the two KeyMaps differ in
KeyCase, the destination was searched under the source's casing rule
and could miss a key its own rule would have matched.

Convert the key once per KeyMap and use each result for its own lookup.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
EOF
```

---

## Task 14: KeyMap dumps never converge

**Files:**
- Modify: `src/keymap.c:10515-10527` (`Dump`)

**Interfaces:**
- Consumes: nothing.
- Produces: nothing.

**Background.** `Dump` walks the hash table bucket by bucket and follows each bucket's chain from the head, ignoring `SortBy` entirely. `LoadKeyMap` re-puts entries in the order it reads them, and `AddTableEntry` head-inserts (`keymap.c:675`). Any two keys sharing a bucket swap places on every round trip, so the dump never stabilises. `ast_tester/stcschan-test1-doc3-props.ast` has `SIZE`/`FRAME` colliding and `TIME`/`REFPOS` colliding.

**This task changes serialized output.** Reference fixtures will change. Review every diff before committing.

- [ ] **Step 1: Read the entry structure and Dump's locals**

```bash
grep -n 'struct AstMapEntry' -A 12 src/keymap.h
sed -n '10470,10530p' src/keymap.c
```

Confirm the field holding the key is named `key` and is a `const char *`, and note which local variables `Dump` already declares so the new ones do not collide.

- [ ] **Step 2: Write a probe demonstrating the instability**

Write this to `/tmp/probe.c` and build it per "Writing a Scratch Probe" above. `SIZE` and `FRAME` collide in bucket 0 of a 16-entry table.

```c
#include <stdio.h>
#include <string.h>
#include "ast.h"

int main( void ) {
   int status = 0;
   AstKeyMap *km;
   AstKeyMap *k2;
   AstKeyMap *k3;
   char *d1;
   char *d2;
   char *d3;

   astWatch( &status );
   km = astKeyMap( " " );
   astMapPut0I( km, "SIZE", 1, " " );
   astMapPut0I( km, "FRAME", 2, " " );

   d1 = astToString( km );
   k2 = astFromString( d1 );
   d2 = astToString( k2 );
   k3 = astFromString( d2 );
   d3 = astToString( k3 );

   printf( "dump1 == dump2 ? %s\n", strcmp( d1, d2 ) ? "NO" : "yes" );
   printf( "dump1 == dump3 ? %s\n", strcmp( d1, d3 ) ? "NO" : "yes" );
   return 0;
}
```

Expected before the fix: `dump1 == dump2 ? NO` and `dump1 == dump3 ? yes` — the two keys swap places on every round trip, so the dump oscillates rather than converging.

- [ ] **Step 3: Apply the fix**

In `Dump`, replace the bucket walk:

```c
/* Loop round each entry in the hash table. */
   for( i = 0; i < this->mapsize; i++ ) {
      next = this->table[ i ];
      while( next && astOK ) {
         DumpEntry( next, channel, ++nentry, status );
         next = next->next;
      }
   }
```

with a pass that collects every entry, sorts by key, and dumps in that order:

```c
/* Gather pointers to every entry, so they can be dumped in key order
   rather than in hash table order. Bucket order depends on the KeyMap's
   insertion history, and since a load re-inserts entries in the order it
   reads them, dumping in bucket order means colliding keys swap places
   on every write and read cycle. */
   nent = 0;
   for( i = 0; i < this->mapsize; i++ ) {
      for( next = this->table[ i ]; next; next = next->next ) nent++;
   }

   entries = astMalloc( sizeof( AstMapEntry * )*(size_t) nent );
   if( entries ) {
      j = 0;
      for( i = 0; i < this->mapsize; i++ ) {
         for( next = this->table[ i ]; next; next = next->next ) {
            entries[ j++ ] = next;
         }
      }

      qsort( entries, (size_t) nent, sizeof( AstMapEntry * ), DumpCmp );

      for( j = 0; j < nent && astOK; j++ ) {
         DumpEntry( entries[ j ], channel, ++nentry, status );
      }

      entries = astFree( entries );
   }
```

Declare `nent`, `j` as `int` and `entries` as `AstMapEntry **` alongside the existing locals. `astMalloc` sets the status on failure, so the `if( entries )` guard is sufficient; do not add an error report.

Add the comparator above `Dump`:

```c
/* Compare two entries by key, for sorting the dump into a stable order. */
static int DumpCmp( const void *a, const void *b ) {
   const AstMapEntry *ea = *( (const AstMapEntry * const *) a );
   const AstMapEntry *eb = *( (const AstMapEntry * const *) b );
   return strcmp( ea->key, eb->key );
}
```

`qsort` needs `<stdlib.h>` and `strcmp` needs `<string.h>`; both are already included by `keymap.c`, but confirm with `grep`.

Use `strcmp` directly rather than calling `KeyCmp`: `KeyCmp` is written as a `SortBy` comparator and takes its arguments in the form that machinery supplies, so reusing it here would mean adapting its signature. A plain key comparison is what this needs, and it is what makes the output idempotent.

The `MapSz` card and the hash table itself are unaffected; only the order in which entries are written changes.

- [ ] **Step 4: Add the prologue history entry**

```
*     8-AUG-2026 (TIMJ):
*        Dump KeyMap entries in key order rather than hash table order,
*        so that a dump is stable under repeated write and read cycles.
```

- [ ] **Step 5: Verify the dump is now idempotent**

Re-run the probe from Step 2. Expected: dump 1, dump 2 and dump 3 all identical.

- [ ] **Step 6: Run the full suite and review fixture diffs**

```bash
cmake --build build -j && ctest --test-dir build --output-on-failure
git status --short
git diff -- ast_tester/
```

Expect `.ast` reference files containing KeyMaps to change. **Review each diff and confirm it is a reordering of entries and nothing else.** Report the list of changed fixtures to the user before committing.

- [ ] **Step 7: Run the sanitizer build**

```bash
cmake --build build-dev -j && ctest --test-dir build-dev --output-on-failure
```

- [ ] **Step 8: Commit**

```bash
git add src/keymap.c ast_tester/
git commit -F - <<'EOF'
keymap: dump entries in key order so that dumps are stable

Dump walked the hash table bucket by bucket and followed each bucket's
chain from the head. astLoadKeyMap re-puts the entries in the order it
reads them, and AddTableEntry inserts at the head of the bucket, so any
two keys sharing a bucket swapped places on every write and read cycle
and the dump never converged. Whether two logically identical KeyMaps
serialised identically depended on their insertion history and on how
many times each had been round tripped.

Emit the entries in key order instead. Dump already ignored the SortBy
attribute, so no documented behaviour changes, and ordered output is
idempotent under any insertion order. The hash table and the MapSz card
that records its size are unaffected.

Reference files containing KeyMaps are updated for the new entry order.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
EOF
```

---

## Task 15: EncodeFloat injects a stray space into wide values

**Files:**
- Modify: `src/fitschan.c:9866-9873`

**Interfaces:**
- Consumes: nothing.
- Produces: nothing.

**Background.** `EncodeFloat` formats into a right-justified 20-column field and then removes a redundant leading zero from the exponent with a length-preserving shift: it slides everything left of the zero one place right and writes a space into the vacated first character. That shift *is* the right-justification, so it is correct whenever padding exists to consume. When the two-digit-exponent form already exceeds the field width there is no padding, the space survives, and the value starts one column further right than every other value in the header.

The same function already uses the correct idiom twenty lines below, where the `".0"` insertion guards its shift on `buf[ 0 ] == ' ' && buf[ 1 ] == ' '` (`fitschan.c:9913`).

**This task changes serialized output.** Review fixture diffs before committing.

- [ ] **Step 1: Confirm the current behaviour**

Write this to `/tmp/probe.c` and build it per "Writing a Scratch Probe" above:

```c
#include <stdio.h>
#include "ast.h"

static void show( AstFitsChan *fc, const char *key, double value ) {
   char card[ 81 ];
   astSetFitsF( fc, key, value, "Scale factor", 1 );
   astClear( fc, "Card" );
   while( astFindFits( fc, key, card, 1 ) ) {
      printf( "[%s]\n", card );
   }
}

int main( void ) {
   int status = 0;
   AstFitsChan *fc;

   astWatch( &status );
   fc = astFitsChan( NULL, NULL, "FitsDigits=17" );

   /* Wide: the two digit exponent form is 23 characters, so there is no
      padding for the shift to consume. */
   show( fc, "SCL1_A", -1.7453292519943296E-04 );

   /* Narrow: the two digit exponent form fits, so padding absorbs the
      shift and the value must stay right justified. */
   show( fc, "SCL2_A", 1.0E-05 );
   return 0;
}
```

Expected before the fix: the `SCL1_A` card has **two** spaces after the `=` and its value runs one character past column 30, while `SCL2_A` has one space and ends in column 30.

Record both cards verbatim. Step 4 compares against them.

- [ ] **Step 2: Apply the fix**

The code currently reads:

```c
/* If a leading zero was found, shuffle everything down from the start of
   the string by one character, over-writing the redundant zero, and insert
   a space at the start of the string. */
      if( w ) {
         r = w - 1 ;
         while( w != buf ) *(w--) = *(r--);
         *w = ' ';
      }
```

Replace it with:

```c
/* If a leading zero was found, remove it. If there is a leading space to
   absorb the change, shuffle everything down from the start of the
   string by one character, over-writing the redundant zero, and insert a
   space at the start; this keeps the field right justified within the
   desired field width. Otherwise there is no padding to consume, so
   close the gap instead and let the string shorten. */
      if( w ) {
         if( buf[ 0 ] == ' ' ) {
            r = w - 1 ;
            while( w != buf ) *(w--) = *(r--);
            *w = ' ';
         } else {
            memmove( w, w + 1, strlen( w + 1 ) + 1 );
         }
      }
```

`memmove` and `strlen` need `<string.h>`; confirm with `grep -n '#include <string.h>' src/fitschan.c`.

- [ ] **Step 3: Add the prologue history entry**

```
*     8-AUG-2026 (TIMJ):
*        In EncodeFloat, close the gap when removing a redundant leading
*        zero from an exponent if there is no leading padding to absorb
*        the shift, rather than pushing the value one column right.
```

- [ ] **Step 4: Verify the fix**

Re-run the probe from Step 1. Expected: `SCL1_A` now has one space after the `=`, matching neighbouring cards, and `SCL2_A` is **byte-identical to what Step 1 recorded**. If `SCL2_A` changed, the guard on `buf[ 0 ] == ' ' ` is wrong and narrow values are being left-shifted out of their column.

- [ ] **Step 5: Run the full suite and review fixture diffs**

```bash
cmake --build build -j && ctest --test-dir build --output-on-failure
git diff -- ast_tester/
```

Expect FITS reference outputs containing wide float values to lose a leading space. Review the diffs and report them before committing.

- [ ] **Step 6: Run the sanitizer build**

```bash
cmake --build build-dev -j && ctest --test-dir build-dev --output-on-failure
```

- [ ] **Step 7: Commit**

```bash
git add src/fitschan.c ast_tester/
git commit -F - <<'EOF'
fitschan: keep wide float values right justified in EncodeFloat

EncodeFloat removes a redundant leading zero from an exponent by sliding
everything to the left of the zero one place right and writing a space
into the vacated first character. That shift is what keeps the value
right justified, and it works whenever there is padding to consume.

When the two digit exponent form already fills or exceeds the field
width there is no padding, so the space survived: the value was written
starting one column further right than every other value in the header,
and the field ran one character past its intended end.

Consume a leading pad space when one exists, and otherwise close the gap
and let the string shorten. This matches the idiom already used a few
lines below when inserting a decimal point. Values narrow enough to have
padding are unaffected.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
EOF
```

---

## Task 16: WcsNative stamps a vestigial Invert on a shared PermMap

**Files:**
- Modify: `src/fitschan.c:38396-38409`

**Interfaces:**
- Consumes: nothing.
- Produces: nothing.

**Background.** When the celestial axes are not already at slots 0 and 1, `WcsNative` builds one `AstPermMap` and uses it twice: forward for stage 1, then inverted for stage 3. `astCmpMap` clones rather than copies — `CombineMaps` (`cmpmap.c:616`) copies only when the same pointer is supplied for both components of one CmpMap, which is not the case here — so the stage-1 CmpMap holds the object that `astInvert` then flips. The stage-1 CmpMap records `invert1 = 0` and dumps no `InvA`, yet its child PermMap dumps `Invert = 1`, a state acquired after encapsulation.

The flag never reaches the mathematics, because `Transform` forces the parent's recorded `invert1` onto the child before transforming (`cmpmap.c:3073`), but it is serialized.

**This task changes serialized output.** Review fixture diffs before committing.

- [ ] **Step 1: Apply the fix**

Change:

```c
/* Now invert the PermMap, so that it re-arranges the axes back into
   their original order. This is the mapping described as stage 3 in
   the prologue. */
         astInvert( permmap );

/* And finally.... add this inverted PermMap onto the end of the CmpMap. */
         cmpmap = astCmpMap( new, permmap, 1, "", status );
         permmap = astAnnul( permmap );
```

to:

```c
/* Take a copy for stage 3 and invert that, so that it re-arranges the
   axes back into their original order. A copy is used because astCmpMap
   clones rather than copies its components, so inverting the original
   would alter the PermMap already encapsulated in the stage 1 CmpMap
   and leave it carrying a flag set after encapsulation. */
         permmap2 = astCopy( permmap );
         permmap = astAnnul( permmap );
         astInvert( permmap2 );

/* And finally.... add this inverted PermMap onto the end of the CmpMap. */
         cmpmap = astCmpMap( new, permmap2, 1, "", status );
         permmap2 = astAnnul( permmap2 );
```

Declare `permmap2` alongside the existing `permmap` declaration, with the same type.

- [ ] **Step 2: Add the prologue history entry**

```
*     8-AUG-2026 (TIMJ):
*        In WcsNative, use a copy of the PermMap for the axis
*        rearrangement stage, so that the PermMap encapsulated in the
*        earlier CmpMap does not acquire an Invert flag afterwards.
```

- [ ] **Step 3: Run the full suite and review fixture diffs**

```bash
cmake --build build -j && ctest --test-dir build --output-on-failure
git diff -- ast_tester/
```

Expect reference dumps for headers whose celestial axes are not at slots 0 and 1 to lose an `Invert` card under a `PermMap`, which shifts the suffix of later `Invert` cards. Review and report before committing.

- [ ] **Step 4: Run the sanitizer build**

```bash
cmake --build build-dev -j && ctest --test-dir build-dev --output-on-failure
```

Watch for leaks: `permmap` and `permmap2` must each be annulled exactly once on every path.

- [ ] **Step 5: Commit**

```bash
git add src/fitschan.c ast_tester/
git commit -F - <<'EOF'
fitschan: do not invert an encapsulated PermMap in WcsNative

WcsNative built one PermMap and used it twice: forward for the stage 1
CmpMap, then inverted for stage 3. astCmpMap clones rather than copies
its components, so the stage 1 CmpMap held the same object that was
subsequently inverted.

The stage 1 CmpMap recorded an invert flag of zero at construction and
so dumped no InvA card, yet its child PermMap dumped Invert = 1, a state
it acquired only after being encapsulated. The flag never reached the
arithmetic, because Transform forces the parent's recorded value onto
the child, but it was serialised, and a reader trusting it would need
the parent's InvA to cancel it back out.

Use a copy for stage 3, so that neither CmpMap's child carries a flag
set after encapsulation and the dump reflects the Mapping rather than
the construction sequence.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
EOF
```

---

## Task 17: CDjjjiii cards are consumed under encodings that discard them

**Files:**
- Modify: `src/fitschan.c:31284-31299`

**Interfaces:**
- Consumes: nothing.
- Produces: nothing.

**Background.** `SpecTrans` reads each element of the legacy six-digit `CDjjjiii` matrix with `GetValue2` unconditionally, but stores the result only when `encoding == FITSIRAF_ENCODING`. `GetValue2` marks the card used when it finds it in the second FitsChan (`fitschan.c:18376`). Under any other encoding the values are read, discarded, and the cards still fail to reach the output, so the rotation is silently lost.

- [ ] **Step 1: Apply the fix**

The loop currently reads:

```c
                  sprintf( keyname, "CD%.3d%.3d", j + 1, i + 1 );
                  if( GetValue2( ret, this, keyname, AST__FLOAT, (void *) &dval, 0,
                                method, class, status ) ){
                     if( encoding == FITSIRAF_ENCODING ){
                        SetValue( ret, FormatKey( "PC", j + 1, i + 1, ' ', status ),
                                  (void *) &dval, AST__FLOAT, NULL, status );
                        dval = 1.0;
                        SetValue( ret, FormatKey( "CDELT", j + 1, -1, s, status ),
                                  (void *) &dval, AST__FLOAT, NULL, status );
                        gotpcij = 1;
                     }
                  }
```

Move the encoding test outside, so the read happens only when its result is used:

```c
                  if( encoding == FITSIRAF_ENCODING ){
                     sprintf( keyname, "CD%.3d%.3d", j + 1, i + 1 );
                     if( GetValue2( ret, this, keyname, AST__FLOAT,
                                    (void *) &dval, 0, method, class,
                                    status ) ){
                        SetValue( ret, FormatKey( "PC", j + 1, i + 1, ' ', status ),
                                  (void *) &dval, AST__FLOAT, NULL, status );
                        dval = 1.0;
                        SetValue( ret, FormatKey( "CDELT", j + 1, -1, s, status ),
                                  (void *) &dval, AST__FLOAT, NULL, status );
                        gotpcij = 1;
                     }
                  }
```

Better still, hoist the encoding test outside the two nested `for` loops entirely, so the loops do not run at all under other encodings. Read the surrounding code and choose whichever is the smaller, clearer change.

- [ ] **Step 2: Add the prologue history entry**

```
*     8-AUG-2026 (TIMJ):
*        Only read the legacy CDjjjiii matrix in SpecTrans under the
*        FITS-IRAF encoding, which is the only encoding that uses the
*        values, so the cards are not consumed and discarded under other
*        encodings.
```

- [ ] **Step 3: Run the full suite and review fixture diffs**

```bash
cmake --build build -j && ctest --test-dir build --output-on-failure
git diff -- ast_tester/
```

Expect any reference output for a header carrying `CDjjjiii` cards under a non-IRAF encoding to now retain those cards. Review and report before committing.

- [ ] **Step 4: Run the sanitizer build**

```bash
cmake --build build-dev -j && ctest --test-dir build-dev --output-on-failure
```

- [ ] **Step 5: Commit**

```bash
git add src/fitschan.c ast_tester/
git commit -F - <<'EOF'
fitschan: do not consume CDjjjiii cards under encodings that ignore them

SpecTrans read every element of the legacy six digit CDjjjiii matrix
unconditionally, but stored the result, as PCj_i plus a unit CDELTj,
only under the FITS-IRAF encoding. Reading a card marks it as used, so a
header carrying CDjjjiii lost those cards under every encoding while
only FITS-IRAF gained the rotation matrix they describe. Under the other
encodings the values were read, discarded, and the cards still failed to
reach the output, so the rotation was silently lost rather than either
honoured or passed through.

Read the cards only under the encoding that uses them.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
EOF
```

---

## Task 18: Duplicate-keyword consumption differs between the two read paths

**Files:**
- Modify: `src/fitschan.c` — the `GetValue` marking path and `SpecTrans`'s MJD-OBS handling at 31703

**Interfaces:**
- Consumes: nothing.
- Produces: nothing.

**Background.** `WcsFcRead` (`fitschan.c:36233`) is a card sweep that marks every copy of a repeated keyword. `SpecTrans` (`fitschan.c:30880`) reads through `GetValue2` → `GetValue` → `SearchCard` → `FindKeyCard`, a first-match forward search (`fitschan.c:18185`), so it marks exactly one card however many copies exist.

For a header holding two HDUs' worth of WCS cards, one `astRead` strips both copies of `CRPIXi`, `CRVALi`, `CTYPEi` and `EQUINOX` but leaves the second copy of `CDi_j`, `RADECSYS` and `DATE-OBS` — a partial, self-inconsistent WCS description.

There is also an ordering asymmetry: `WcsFcRead`'s `MJD-OBS` branch deliberately sets `mark = 0` so the card survives a read, but `SpecTrans`'s translation-9 block reaches `MJD-OBS` through `GetValue2` with marking enabled and consumes it anyway.

`WcsFcRead`'s sweep is taken as the intended semantics, because it carries explicit per-keyword exceptions whereas `SpecTrans`'s single-copy marking is a consequence of the lookup helper it happens to use.

**This is the largest behavioural change in the plan.** It alters which cards survive `astRead` on production headers.

- [ ] **Step 1: Map the marking path**

```bash
sed -n '18185,18300p' src/fitschan.c   # GetValue / SearchCard / FindKeyCard
grep -n 'MarkCard' src/fitschan.c
```

Establish exactly where a card is marked used, and whether `MarkCard` is called once per `GetValue` or by a lower layer.

- [ ] **Step 2: Add a mark-all-copies helper**

Add a static function beside `MarkCard` that marks every card whose keyword matches a given name, modelled on `WcsFcRead`'s sweep: save the current card index, `astClearCard`, walk the FitsChan marking each match, then restore the saved index.

Preserving and restoring the caller's card position is essential — `SpecTrans` reads many keywords in sequence and relies on the current card.

- [ ] **Step 3: Route SpecTrans's consumption through it**

Where `GetValue` currently marks the single card it found, mark every copy instead. Confine the change to the marking step; the *value* returned must still come from the first match, so translations continue to read the same numbers they read today.

- [ ] **Step 4: Respect the MJD-OBS exception**

In `SpecTrans`'s translation-9 block at `fitschan.c:31703`, read `MJD-OBS` without marking it, matching `WcsFcRead`'s deliberate `mark = 0` for that keyword. Whether `MJD-OBS` survives a read must no longer depend on whether `DATE-OBS` sits alongside it.

- [ ] **Step 5: Add the prologue history entry**

```
*     8-AUG-2026 (TIMJ):
*        Consume every copy of a repeated keyword in SpecTrans, matching
*        WcsFcRead, and respect the MJD-OBS exception so that whether
*        that card survives a read no longer depends on whether DATE-OBS
*        is present.
```

- [ ] **Step 6: Run the full suite and review fixture diffs**

```bash
cmake --build build -j && ctest --test-dir build --output-on-failure
git diff -- ast_tester/
```

Expect reference outputs for multi-HDU headers to lose the surviving second copies of `CDi_j`, `RADECSYS` and `DATE-OBS`. **Review every diff and report the full list to the user before committing.** This is the change most likely to surprise, so do not absorb any diff silently.

- [ ] **Step 7: Run the sanitizer build**

```bash
cmake --build build-dev -j && ctest --test-dir build-dev --output-on-failure
```

- [ ] **Step 8: Commit**

```bash
git add src/fitschan.c ast_tester/
git commit -F - <<'EOF'
fitschan: consume every copy of a repeated keyword when reading WCS

Reading a FITS-WCS header consumes cards through two paths that treated
repeated keywords differently. WcsFcRead sweeps every card and marks
each match, so every copy of a repeated keyword is consumed. SpecTrans
reads through GetValue, a first match forward search, so it marked
exactly one card however many copies existed.

For a header holding two HDUs' worth of WCS cards, a single astRead
therefore stripped both copies of CRPIXi, CRVALi, CTYPEi and EQUINOX but
left the second copy of CDi_j, RADECSYS and DATE-OBS behind. What
survived was a partial, self inconsistent WCS description, and which
keywords survived was decided by which read path happened to handle
them.

Mark every copy in the SpecTrans path as well. The value used is still
that of the first match, so no translation reads a different number.

Also respect the MJD-OBS exception. WcsFcRead deliberately leaves that
card unmarked so it survives a read, but SpecTrans reached it with
marking enabled when DATE-OBS was also present, so whether it survived
depended on an unrelated keyword.

Co-Authored-By: Claude Opus 5 (1M context) <noreply@anthropic.com>
EOF
```

---

## Task 19: Investigate the IsSimple tag placement defect

**This task produces a report, not a fix.** Its output selects between two candidate repairs, and the choice is put to the user before any code is written.

**Files:**
- Create: a scratch instrumented build. Nothing under `src/` is committed by this task.

**Interfaces:**
- Consumes: nothing.
- Produces: a recommendation the user acts on.

**Background.** `IsSimple` is one bit of `Mapping.flags` (`mapping.h:410`) carrying two meanings. `astSimplify_` raises it on whatever it returns (`mapping.c:24759`) and short-circuits for any Mapping whose bit is already set (`mapping.c:24740`), so it means both "dump `IsSimp = 1`" and "do not simplify again". `astMAKE_SET(Mapping,Invert,...)` clears it on every set (`mapping.c:23256`), even when the value does not change; `astMAKE_CLEAR` does not.

Probes decompose a candidate with `astMapList`, which returns reference-counted pointers to the live components. They call `astSetInvert` to line the components up and restore the Invert *value* afterwards, but not the `IsSimple` bit that `astSetInvert` incidentally cleared. `cmpmap.c`'s `Equal` (`cmpmap.c:386-398`), `matrixmap.c`'s `CanSwap` (`416-421`, restored at `582-583`) and `winmap.c`'s `CanSwap` (`102-106`, restored at `201-202`) all follow this shape.

Which interior node ends up tagged in a dump is therefore decided by which nodes were aliased into a probe, not by the tree's shape.

- [ ] **Step 1: Build an instrumented library**

In a scratch build directory, not the working tree, add temporary instrumentation to `astMAKE_SET(Mapping,Invert,...)` in `src/mapping.c` that logs, for every set that clears an already-raised `IsSimple` bit:

- the Mapping's address and class;
- whether the stored Invert value actually changed;
- a backtrace, or at minimum the immediate caller resolved with `backtrace_symbols` or `addr2line`.

Disable object caching so addresses are not recycled and pointer identity is readable across a run.

- [ ] **Step 2: Run the suite against it**

```bash
ctest --test-dir <scratch-build> --output-on-failure
```

Collect the log.

- [ ] **Step 3: Analyse the distribution**

Produce a table of clears grouped by calling routine, distinguishing:

- clears where the stored Invert value genuinely changed, which are arguably correct, since a flipped Mapping is a different Mapping;
- clears where the value did not change, or was restored immediately afterwards — the save-flip-restore probes, which are the defect;
- how many distinct routines and files the second group spans.

- [ ] **Step 4: Report and recommend**

Report to the user:

1. the table from Step 3;
2. how many probe sites would need editing under the "make probes non-mutating" fix;
3. whether any probe restores Invert with `astClearInvert` rather than `astSetInvert`, since `astMAKE_CLEAR` is already tag-preserving and those sites need no change;
4. a recommendation between:
   - **splitting the flag** — keep `AST__ISSIMPLE_FLAG` as the simplify guard that `astSetInvert` clears, and add a separate bit that only `astSimplify_` sets and only `Dump` reads. Confined to `mapping.c` and `mapping.h`; every probe stays as written;
   - **making probes non-mutating** — save and restore `IsSimple` alongside `Invert` at each site. Truer to each routine's read-only contract, but spans a dozen or more files and each site needs individual judgement.

**Stop here and wait for the user's decision.** Do not implement either fix in this task.

---

## Self-Review Notes

Spec coverage check against `docs/superpowers/specs/2026-08-08-c-library-defect-fixes-design.md`:

| Spec item | Task |
|---|---|
| R1 rounding sweep | 2, 3 (and 1 for the grism site) |
| R1a grism order | 1 |
| A1 MapIterate | 4 |
| A2 MapRename | 5 |
| A3 zero-length vector | 8 |
| A4 range clamping | 6 |
| A5 string overflow | 7 |
| A6 KyCas | 10 |
| A7 MpLck | 11 |
| A8 MapGetElem | 9 |
| A9 mem%d | 12 |
| A10 MapCopyEntry | 13 |
| B1 Dump order | 14 |
| B2 EncodeFloat | 15 |
| B3 WcsNative | 16 |
| B4 MapSz | No task. The spec records a deliberate decision to make no change. |
| C1 IsSimple | 19 (investigation only) |
| C2 CDjjjiii | 17 |
| C3 duplicate consumption | 18 |
| CheckCircle | No task. Examined and correct. |

Ordering constraints:

- Task 2 must precede Task 6: the clamping helpers wrap the `round()`-corrected expressions.
- Task 6 must precede Task 7: the `%lf` fallback that Task 7 routes overflow into relies on Task 6's saturating helpers.
- Tasks 14, 15, 16, 17 and 18 change serialized output; each reviews its fixture diffs before committing.
- Task 19 is terminal and blocks on a user decision.

All other tasks are independent and may be reordered.
