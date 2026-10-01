# ChebyMap Iterative Inverse Review Fixes Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Resolve the ten review findings and four cleanup items on branch `u/timj/chebymap-inverse` by moving the bounded Newton inverse into ChebyMap and fixing the attribute, FitsChan, test and documentation defects around it.

**Architecture:** `IterInverse` becomes a protected virtual method of PolyMap. PolyMap keeps the unbounded algorithm from `master`. ChebyMap overrides it with the bounded algorithm, `IterSteps`, and a static domain helper, falling back to the parent when the forward series has no usable Chebyshev box. Attribute rules for `IterInverse`, `NiterInverse` and `TolInverse` live in `src/polymap.c`; ChebyMap overrides only the `NiterInverse` default. FitsChan excludes ChebyMaps from SIP output at the point where it picks the PolyMap.

**Tech Stack:** C99 library, C11 tests, CMake, ctest, AST object system (`astMAKE_*` macros, vtabs, `astBegin`/`astEnd`).

**Spec:** `docs/superpowers/specs/2026-09-10-chebymap-inverse-review-fixes-design.md`

## Global Constraints

- Build and test with the developer configuration from `CLAUDE.md`: `cmake -B build-dev -DCMAKE_BUILD_TYPE=Debug -DAST_ENABLE_WARNINGS=ON -DAST_ENABLE_SANITIZERS=ON`.
- No new compiler warnings and no sanitizer reports from changed code.
- Every change to a file under `src/` adds a History entry in that file's prologue.
- Public C API only in tests, except `ast_tester/testchebyinverse.c`, which is built with `INTERNAL_HEADERS` and may use protected functions.
- Prose in Markdown uses one sentence per line and American English spelling.
- Commit after each task with `git commit`; never push. End commit messages with the attribution lines given in the session.
- The regression test set for every task is:

```bash
cmake --build build-dev -j8 && ctest --test-dir build-dev --output-on-failure \
  -R 'testchebymap_c|testchebyinverse|testpolymap_c|testthreads|testfitschan_c|transform_oracle'
```

- Ordering: Task 1 before Task 2. Task 2 before Tasks 5, 6, 10. Task 9 before Task 12. Task 3 before Task 4. Everything before Task 13.

---

### Task 1: Make IterInverse a protected virtual method of PolyMap

**Files:**
- Modify: `src/polymap.h:204-206` (vtab), `src/polymap.h:279-281` (prototypes), `src/polymap.h:360-365` (macros), `src/polymap.h:114` (History)
- Modify: `src/polymap.c:348` (static prototype), `src/polymap.c:2495-2497` (vtab init), `src/polymap.c:6691` (Transform dispatch), `src/polymap.c:8116-8130` (wrappers)
- Test: `ast_tester/testchebyinverse.c`

**Interfaces:**
- Produces: vtab member `void (* IterInverse)( AstPolyMap *, AstPointSet *, AstPointSet *, int * );`, protected function `void astIterInverse_( AstPolyMap *, AstPointSet *, AstPointSet *, int * );`, macro `astIterInverse(this,out,result)`.
- The `GetIterDomain` vtab member, `astGetIterDomain_` and the `astGetIterDomain` macro are still present after this task. Task 2 removes them.

- [ ] **Step 1: Add the vtab member, prototype and macro in `src/polymap.h`**

After line 206 (`int (* GetIterDomain)...`) add:

```c
   void (* IterInverse)( AstPolyMap *, AstPointSet *, AstPointSet *, int * );
```

After line 281 (`int astGetIterDomain_...`) add:

```c
   void astIterInverse_( AstPolyMap *, AstPointSet *, AstPointSet *, int * );
```

After line 365 (the `astGetIterDomain` macro) add:

```c
#define astIterInverse(this,out,result) \
        astINVOKE(V,astIterInverse_(astCheckPolyMap(this),out,result,STATUS_PTR))
```

Add to the History block after `28-SEP-2003 (DSB): Original version.`:

```
*     10-SEP-2026 (TIMJ):
*        Add the protected virtual methods astGetJacobian, astLinearGuess
*        and astIterInverse, so that a subclass can supply its own
*        Jacobian, initial guess and iterative inverse algorithm.
```

- [ ] **Step 2: Wire the vtab in `src/polymap.c`**

At line 2497, after `vtab->GetIterDomain = GetIterDomain;` add:

```c
   vtab->IterInverse = IterInverse;
```

At line 6691 replace `IterInverse( map, in, result, status );` with:

```c
      astIterInverse( map, in, result );
```

After the `astGetIterDomain_` wrapper (line 8130) add:

```c
void astIterInverse_( AstPolyMap *this, AstPointSet *out, AstPointSet *result,
                      int *status ){
   if ( !astOK ) return;
   (**astMEMBER(this,PolyMap,IterInverse))( this, out, result, status );
}
```

Update the `IterInverse` prologue (`src/polymap.c:2810-2870`): change `Type: Private function.` to `Protected virtual function.` and the Class Membership line to `PolyMap method (implements the astIterInverse protected method)`.

- [ ] **Step 3: Build and run the regression set**

Run the Global Constraints test command. Expected: all listed tests pass, no new warnings.

- [ ] **Step 4: Commit**

```bash
git add src/polymap.h src/polymap.c
git commit -m "Make the PolyMap iterative inverse a protected virtual method"
```

---

### Task 2: Move the bounded algorithm into ChebyMap and restore the unbounded PolyMap algorithm

**Files:**
- Modify: `src/polymap.c:2562-2808` (delete `IterSteps`), `src/polymap.c:2810-3261` (`IterInverse`), `src/polymap.c:1983-2038` (delete `GetIterDomain`), `src/polymap.c:349`, `src/polymap.c:323`, `src/polymap.c:2497`, `src/polymap.c:8126-8130`
- Modify: `src/polymap.h:206`, `src/polymap.h:281`, `src/polymap.h:364-365`
- Modify: `src/chebymap.c` (new `IterInverse`, `IterSteps`, `MarkUnsolved`, `Usable`; `GetIterDomain` becomes static-only; vtab init at 1727-1736)
- Test: `ast_tester/testchebyinverse.c` (`seeds`, `bounded_solver`, `caches`)

**Interfaces:**
- Consumes: `astIterInverse` vtab slot from Task 1; existing `astGetJacobian(map)` and `astLinearGuess(map)` protected methods.
- Produces in `src/chebymap.c`: `static void IterInverse( AstPolyMap *, AstPointSet *, AstPointSet *, int * );`, `static void IterSteps( AstPolyMap *, int, int, double **, double **, const double *, const double *, const double *, const double *, const double *, const double *, int *, int *, int *, int * );`, `static int GetIterDomain( AstChebyMap *, double *, double *, int * );`, `static void MarkUnsolved( double **, int, int, int *, int * );`, `static int Usable( double );`, and file-static `parent_iterinverse`.
- Tasks 5, 6 and 10 edit the ChebyMap `IterInverse` written here.

- [ ] **Step 1: Rewrite the tests that used the removed hook**

In `ast_tester/testchebyinverse.c`:

Replace the two `astGetIterDomain` checks in `seeds` (lines 108-110) with:

```c
   astChebyDomain( cm, 0, &dlo, &dhi );
   near( dlo, lo, "Original forward lower bound" );
   near( dhi, hi, "Original forward upper bound" );
```

Delete the `testdomain` function (lines 160-167) and replace `bounded_solver` (lines 236-262) with a ChebyMap on [-1, 1]. The monomial -0.125 + x + 0.25 x^2 is 0 T0 + 1 T1 + 0.125 T2 because x^2 = (T0 + T2)/2:

```c
static void bounded_solver( int *status ) {
   double lo = -1, hi = 1;
   double coeffs[] = { 1, 1, 1, .125, 1, 2 };
   double target[] = { -.875, -.125, 0, .3, 1.125, -1, 2, AST__BAD, NAN, INFINITY };
   double mono[] = { -.125, 1, 0, 1, 1, 1, .25, 1, 2 };
   double got[10];
   AstChebyMap *cm = astChebyMap( 1, 1, 2, coeffs, 0, NULL, &lo, &hi,
                                  NULL, NULL, "", status );
   AstPolyMap *pm;
   int i;
   astSet( cm, "NiterInverse=20,TolInverse=1e-12", status );
   astTran1( cm, 10, target, 0, got );
   for( i = 0; i < 5; i++ ) {
      near( got[i], 2*(target[i]+.125)/(1+sqrt(1+target[i]+.125)),
            "Bounded quadratic inverse" );
   }
   for( i = 5; i < 10; i++ ) {
      check( got[i] == AST__BAD, "Unsolved bounded input must return BAD" );
   }
   astSet( cm, "NiterInverse=1", status );
   astTran1( cm, 1, target+3, 0, got );
   check( got[0] == AST__BAD, "Exhaustion must not return the last iterate" );
   cm = astAnnul( cm );

/* The same polynomial as an unbounded PolyMap retains its historical
   last-iterate behavior. */
   pm = astPolyMap( 1, 1, 3, mono, 0, NULL, "NiterInverse=1", status );
   astTran1( pm, 1, target+3, 0, got );
   check( got[0] != AST__BAD && isfinite(got[0]), "Legacy PolyMap exhaustion" );
   pm = astAnnul( pm );
}
```

In `caches`, replace the two `astGetIterDomain` blocks (lines 208-234) with `astChebyDomain` checks that expect `AST__BAD`:

```c
   if( loaded ) {
      double dlo = -99, dhi = 99;
      astChebyDomain( loaded, 1, &dlo, &dhi );
      check( dlo == AST__BAD && dhi == AST__BAD,
             "Unresolvable box has no evaluable domain" );
      loaded = astAnnul( loaded );
   }
```

and

```c
   if( loaded ) {
      double dlo = -99, dhi = 99;
      astChebyDomain( loaded, 1, &dlo, &dhi );
      check( dlo == AST__BAD && dhi == AST__BAD,
             "Zero scale has no evaluable domain" );
      loaded = astAnnul( loaded );
   }
```

Note: the "Unresolvable box" expectation only holds after Task 9 routes `ChebyDomain` through `AxisBounds`. Until then that one `check` fails; that is expected and is listed in Task 9.

- [ ] **Step 2: Run the test to see it fail to compile**

Run: `cmake --build build-dev -j8 2>&1 | grep -E 'error|warning' | head`
Expected: no errors yet (the hooks still exist). Proceed; the failure point is Step 5.

- [ ] **Step 3: Write ChebyMap's `IterInverse` and helpers in `src/chebymap.c`**

Add static prototypes next to line 272:

```c
static void IterInverse( AstPolyMap *, AstPointSet *, AstPointSet *, int * );
static void IterSteps( AstPolyMap *, int, int, double **, double **,
                       const double *, const double *, const double *,
                       const double *, const double *, const double *,
                       int *, int *, int *, int * );
static int GetIterDomain( AstChebyMap *, double *, double *, int * );
static void MarkUnsolved( double **, int, int, int *, int * );
static int Usable( double );
```

Add next to line 215:

```c
static void (*parent_iterinverse)( AstPolyMap *, AstPointSet *, AstPointSet *, int * );
```

Add the two small helpers after `IterBounds`:

```c
static int Usable( double value ) {
/*
*  Name:
*     Usable
*  Purpose:
*     Test whether a coordinate value can take part in the iteration.
*  Description:
*     Returns non-zero if "value" is neither AST__BAD nor a non-finite
*     floating point value.
*/
   return value != AST__BAD && isfinite( value );
}

static void MarkUnsolved( double **inputs, int ncoord, int ipoint,
                          int *flags, int *nconv ) {
/*
*  Name:
*     MarkUnsolved
*  Purpose:
*     Record that a batch position has no solution.
*  Description:
*     Sets every coordinate of position "ipoint" in "inputs" to AST__BAD,
*     sets its convergence flag and increments the count of resolved
*     positions, so the iteration does not touch it again.
*/
   int i;
   for( i = 0; i < ncoord; i++ ) inputs[ i ][ ipoint ] = AST__BAD;
   flags[ ipoint ] = 1;
   (*nconv)++;
}
```

Change `GetIterDomain` (chebymap.c:1117) to take `AstChebyMap *this` directly, drop the `AstChebyMap *this = (AstChebyMap *) map;` cast, keep its body, and change the prologue `Class Membership` to `ChebyMap member function.` and `Type` to `Private function.`

Move `IterSteps` from `src/polymap.c:2562-2808` into `src/chebymap.c` verbatim, changing every `for( i = 0; i < ncoord; i++ ) inputs[ i ][ ipoint ] = AST__BAD; flags[ ipoint ] = 1; (*nconv)++;` triple into `MarkUnsolved( inputs, ncoord, ipoint, flags, nconv );` and every `value == AST__BAD || !isfinite( value )` into `!Usable( value )`. Change the prologue `Parameters: this` text from "Pointer to the PolyMap" to "Pointer to the ChebyMap, supplied as a PolyMap pointer".

Write ChebyMap's `IterInverse` after `IterSteps`. Its body is the `bounded == 1` path of `src/polymap.c:2810-3261` with the unbounded branches removed. The full function:

```c
static void IterInverse( AstPolyMap *map, AstPointSet *out,
                         AstPointSet *result, int *status ){
/*
*  Name:
*     IterInverse
*  Purpose:
*     Evaluate the original inverse transformation of a ChebyMap by
*     bounded Newton iteration.
*  Type:
*     Private function.
*  Synopsis:
*     #include "polymap.h"
*     void IterInverse( AstPolyMap *map, AstPointSet *out,
*                       AstPointSet *result, int *status )
*  Class Membership:
*     ChebyMap member function (over-rides the astIterInverse protected
*     method inherited from the PolyMap class).
*  Description:
*     This function transforms a set of original output positions into
*     original input positions using Newton-Raphson iteration restricted
*     to the forward bounding box of the ChebyMap. Initial guesses come
*     from astLinearGuess and are clipped into the box. Each Newton
*     correction is normalised by the box half-widths, each residual by
*     its Jacobian row norm, and IterSteps backtracks corrections that
*     do not reduce the residual. A position converges when both the
*     normalised correction and the scaled residual are within
*     TolInverse; an exactly zero residual converges immediately.
*
*     After NiterInverse updates the final candidate is checked once
*     more. Any position that is still unsolved, was supplied with a bad
*     or non-finite coordinate, or met a singular Jacobian is returned
*     as AST__BAD without setting the inherited status.
*
*     If the forward series has no usable Chebyshev box on every axis,
*     the parent PolyMap algorithm is used instead.
*  Parameters:
*     map
*        Pointer to the ChebyMap, supplied as a PolyMap pointer.
*     out
*        PointSet holding the original output positions to invert.
*     result
*        PointSet to receive the original input positions.
*     status
*        Pointer to the inherited status variable.
*/
   AstChebyMap *this = (AstChebyMap *) map;
   AstMapping *lintrunc;
   AstPointSet *work;
   AstPointSet **ps_jac;
   AstPolyMap **jacob;
   double ***ptr_jac;
   double **ptr_in;
   double **ptr_out;
   double **ptr_work;
   double *lbnd;
   double *mat;
   double *norms;
   double *pa;
   double *pb;
   double *scale;
   double *scales;
   double *steps;
   double *ubnd;
   double *vec;
   double *width;
   double det;
   double norm;
   double stepnorm;
   double tol;
   double xx;
   int *flags;
   int *iw;
   int *stepping;
   int exact;
   int fwd;
   int icol;
   int icoord;
   int ipoint;
   int irow;
   int iter;
   int maxiter;
   int nconv;
   int ncoord;
   int npoint;
   int sing;
   int valid;

   if( !astOK ) return;

   ncoord = astGetNin( map );
   if( ncoord != astGetNout( map ) ) {
      astError( AST__INTER, "astTransform(%s): Supplied %s has unequal numbers"
                " of inputs and outputs and therefore an iterative inverse "
                "cannot be used (internal AST Programming error).", status,
                astGetClass(map), astGetClass(map) );
      return;
   }

/* Without a Chebyshev normalisation or a usable box there is nothing to
   restrict the iteration to. */
   lbnd = astMalloc( ncoord*sizeof( *lbnd ) );
   ubnd = astMalloc( ncoord*sizeof( *ubnd ) );
   if( !astOK || !this->scale_f || !GetIterDomain( this, lbnd, ubnd, status ) ) {
      lbnd = astFree( lbnd );
      ubnd = astFree( ubnd );
      (*parent_iterinverse)( map, out, result, status );
      return;
   }

   /* From here, paste src/polymap.c:2943-3260 and apply these edits:
      - delete the "bounded = astGetIterDomain(...)" block and the
        allocation of lbnd/ubnd (done above); allocate width and scale
        unconditionally;
      - delete "maxerr = ...; maxerr *= maxerr;" and keep "tol = astGetTolInverse( map );";
      - the for loop condition becomes "iter <= maxiter && nconv < npoint && astOK";
      - remove every "if( bounded )" and "else if( bounded )" wrapper,
        keeping the bounded body; delete the final "else { vlensq ... }"
        unbounded update block;
      - replace every "for( icoord ... ) ptr_in[icoord][ipoint] = AST__BAD; flags[ipoint] = 1; nconv++;"
        with "MarkUnsolved( ptr_in, ncoord, ipoint, flags, &nconv );";
      - replace "xx == AST__BAD || !isfinite(xx)" and the like with "!Usable( xx )";
      - "jacob = astGetJacobian( map );", "lintrunc = astLinearGuess( map );",
        "astTransform( map, ... )", "astGetInvert( map )", "astGetNiterInverse( map )"
        all use "map". */
}
```

- [ ] **Step 4: Register the override and remove the ChebyMap vtab use of GetIterDomain**

At `src/chebymap.c:1736` replace `polymap->GetIterDomain = GetIterDomain;` with:

```c
   parent_iterinverse = polymap->IterInverse;
   polymap->IterInverse = IterInverse;
```

Update the `src/chebymap.c` prologue History (the entry near line 131 that mentions astGetIterDomain) to say the class over-rides `astIterInverse` and keeps the bounded algorithm in this file.

- [ ] **Step 5: Restore the unbounded PolyMap algorithm and delete the hook**

In `src/polymap.c`:

- Replace the body of `IterInverse` (lines 2877-3261) with the `master` body from `git show master:src/polymap.c | sed -n 2497-2726p`, with two changes: use `astGetJacobian( this )` in place of `GetJacobian( this, status )` and `astLinearGuess( this )` in place of `LinearGuess( this, status )`. Keep the prologue text written in Task 1 but delete its paragraphs about bounded iteration and the two bounded Notes.
- Delete `IterSteps` (lines 2562-2808) and its prototype at line 349.
- Delete `GetIterDomain` (lines 1983-2038), its prototype at line 323, `vtab->GetIterDomain = GetIterDomain;` at 2497, and the `astGetIterDomain_` wrapper at 8126-8130.
- In the History entry near line 178, delete the mention of `astGetIterDomain` and the bounded iteration; state that `astIterInverse` is virtual and that PolyMap implements the unbounded algorithm.

In `src/polymap.h` delete line 206 (`GetIterDomain` member), line 281 (prototype) and lines 364-365 (macro). Adjust the History entry from Task 1 to omit nothing (it already lists only the three surviving methods).

- [ ] **Step 6: Build and run the regression set**

Run the Global Constraints test command.
Expected: every test passes except one `check` in `testchebyinverse` reporting "Unresolvable box has no evaluable domain", which Task 9 fixes. Confirm there is no other failure and no new warning.

- [ ] **Step 7: Commit**

```bash
git add src/polymap.c src/polymap.h src/chebymap.c ast_tester/testchebyinverse.c
git commit -m "Move the bounded iterative inverse into the ChebyMap class"
```

---

### Task 3: One rule for the IterInverse attribute

**Files:**
- Modify: `src/polymap.c:6929-6938` (getter and setter), `src/polymap.c:8075-8076` (loader)
- Modify: `src/chebymap.c:1483-1583` (delete `GetIterInverse` and `SetIterInverse` overrides), `src/chebymap.c:214-215`, `src/chebymap.c:267-268`, `src/chebymap.c:1727-1730`, `src/chebymap.c:3117-3122` (loader comment)
- Test: `ast_tester/testpolymap.c`, `ast_tester/testchebymap.c:556-570`

**Interfaces:**
- Produces: `astSetIterInverse` raises `AST__ATTIN` when a non-zero value is requested on a PolyMap or ChebyMap that has no forward coefficients or unequal Nin and Nout. `astGetIterInverse` returns an explicitly set value when set, and otherwise 1 only when `ncoeff_f` is non-NULL, `ncoeff_i` is NULL and Nin equals Nout.

- [ ] **Step 1: Write failing tests in `ast_tester/testpolymap.c`**

After the test that ends with `stopit( 8022, status )` (line 222) add:

```c
   /* An explicitly set IterInverse survives Get and a dump round trip when
      it is legal, and is rejected when there is no forward transformation
      to iterate on. */
   pm2 = astAnnul( pm2 );
   {
      double fwd3[] = { 2.0, 1, 1 };
      double inv3[] = { 0.5, 1, 1 };
      char *dump;
      AstPolyMap *back;

      pm2 = astPolyMap( 1, 1, 1, fwd3, 1, inv3, "IterInverse=1" );
      if( !astTest( pm2, "IterInverse" ) ) stopit( 8023, status );
      if( astGetI( pm2, "IterInverse" ) != 1 ) stopit( 8024, status );
      dump = astToString( pm2 );
      back = astFromString( dump );
      dump = astFree( dump );
      if( !back || !astTest( back, "IterInverse" ) ||
          astGetI( back, "IterInverse" ) != 1 ) stopit( 8025, status );
      back = astAnnul( back );
      pm2 = astAnnul( pm2 );

      if( *status == 0 ) {
         pm2 = astPolyMap( 1, 1, 0, NULL, 1, inv3, "IterInverse=1" );
         int expected = ( *status == AST__ATTIN && !pm2 );
         astClearStatus;
         if( !expected ) stopit( 8026, status );
         if( pm2 ) pm2 = astAnnul( pm2 );
      }

      if( *status == 0 ) {
         pm2 = astPolyMap( 1, 1, 0, NULL, 1, inv3, "" );
         astSetI( pm2, "IterInverse", 1 );
         int expected = ( *status == AST__ATTIN );
         astClearStatus;
         if( !expected ) stopit( 8027, status );
         if( astGetI( pm2, "IterInverse" ) != 0 ) stopit( 8028, status );
         pm2 = astAnnul( pm2 );
      }

      /* A dump recording IterInv = 1 with no forward coefficients still
         loads, and the unusable value is treated as unset. */
      back = astFromString( " Begin PolyMap\n Nin = 1\n IsA Mapping\n"
                            " MPI1 = 1\n NCI1 = 1\n CI1 = 0.5\n PI1 = 1\n"
                            " IterInv = 1\n End PolyMap\n" );
      if( !back ) stopit( 8029, status );
      if( back && ( astTest( back, "IterInverse" ) ||
                    astGetI( back, "IterInverse" ) ) ) stopit( 8030, status );
      if( back ) back = astAnnul( back );
      pm2 = NULL;
   }
```

Check the dump keyword names by running `astShow` on an existing inverse-only PolyMap in a scratch program if unsure; the loader reads `mpi<i>`, `nci<i>`, `ci<i>`, `pi<i>` and `iterinv` (see `src/polymap.c` `astLoadPolyMap`).

- [ ] **Step 2: Run the test to verify it fails**

Run: `cmake --build build-dev -j8 && ctest --test-dir build-dev -R testpolymap_c --output-on-failure`
Expected: FAIL at 8026 (constructor currently accepts `IterInverse=1` on a PolyMap with no forward coefficients).

- [ ] **Step 3: Implement the rule in `src/polymap.c`**

Replace lines 6929-6938 with:

```c
astMAKE_GET(PolyMap,IterInverse,int,0,( ( this->iterinverse == -INT_MAX ) ?
                                        ( this->ncoeff_f != NULL &&
                                          this->ncoeff_i == NULL &&
                                          astGetNin( this ) == astGetNout( this ) ) :
                                        this->iterinverse ))
astMAKE_SET1(PolyMap,IterInverse,int,iterinverse,
  ( !value || ( this->ncoeff_f && astGetNin(this) == astGetNout(this) ) ) ?
  (( (value?1:0) != this->iterinverse ) ? astClearIsSimple(this) : (void)0,(value?1:0)):
  ( !this->ncoeff_f ?
    astError(AST__ATTIN,"astSetIterInverse(%s): Cannot use an iterative "
             "inverse because the %s has no forward transformation.",
             status, astGetClass(this), astGetClass(this)) :
    astError(AST__ATTIN,"astSetIterInverse(%s): Cannot use an iterative "
             "inverse because the %s has unequal numbers of inputs and "
             "outputs.", status, astGetClass(this), astGetClass(this)),
    this->iterinverse ))
```

In `astLoadPolyMap` replace line 8076 with:

```c
/* A recorded non-zero value that cannot be honoured (no forward
   coefficients) is treated as unset rather than rejected, so that dumps
   written by earlier versions still load. */
      if( new->iterinverse != -INT_MAX && new->iterinverse && !new->ncoeff_f ) {
         new->iterinverse = -INT_MAX;
      }
      if ( TestIterInverse( new, status ) ) SetIterInverse( new, new->iterinverse, status );
```

Update the `IterInverse` attribute documentation (lines 6881-6903, the ChebyMap section) so it no longer says the ChebyMap has its own rule; add to the Notes: "An error is reported if IterInverse is set non-zero for a PolyMap that has no forward transformation."

- [ ] **Step 4: Delete the ChebyMap overrides**

In `src/chebymap.c` delete `GetIterInverse` and `SetIterInverse` (lines 1483-1583), their prototypes (267-268), the `parent_getiterinverse` and `parent_setiterinverse` statics (214-215), the four vtab lines (1727-1730), and the loader comment at 3117-3122. Remove the History sentence near line 135 that describes those overrides. Tests 781 and 783 in `ast_tester/testchebymap.c` continue to expect `AST__ATTIN`, now raised by PolyMap; leave them.

- [ ] **Step 5: Run the regression set**

Run the Global Constraints test command. Expected: pass (plus the known Task 9 `check`).

- [ ] **Step 6: Commit**

```bash
git add src/polymap.c src/chebymap.c ast_tester/testpolymap.c
git commit -m "Apply one IterInverse validity rule to PolyMap and ChebyMap"
```

---

### Task 4: Validate TolInverse and NiterInverse in the setters

**Files:**
- Modify: `src/polymap.c:6982-6984` (NiterInverse setter), `src/polymap.c:7033-7035` (TolInverse setter), attribute docs at 6960-6975 and 7000-7030
- Modify: `src/chebymap.c` ChebyMap `IterInverse` (from Task 2): the `valid = isfinite(tol) && tol > 0.0 && maxiter >= 0;` line
- Test: `ast_tester/testpolymap.c`

- [ ] **Step 1: Write failing tests in `ast_tester/testpolymap.c`**

After the block added in Task 3 add:

```c
   /* Iteration controls must be usable: a non-positive or non-finite
      tolerance and a negative iteration count are rejected. */
   {
      double fwd3[] = { 2.0, 1, 1 };
      const char *bad[] = { "TolInverse=0", "TolInverse=-1", "NiterInverse=-1" };
      int k;
      pm2 = astPolyMap( 1, 1, 1, fwd3, 0, NULL, "" );
      for( k = 0; k < 3 && *status == 0; k++ ) {
         astSet( pm2, bad[k] );
         int expected = ( *status == AST__ATTIN );
         astClearStatus;
         if( !expected ) stopit( 8031 + k, status );
      }
      if( *status == 0 ) {
         astSetD( pm2, "TolInverse", INFINITY );
         int expected = ( *status == AST__ATTIN );
         astClearStatus;
         if( !expected ) stopit( 8034, status );
      }
      if( astGetD( pm2, "TolInverse" ) != 1.0E-6 ) stopit( 8035, status );
      if( astGetI( pm2, "NiterInverse" ) != 4 ) stopit( 8036, status );
      astSet( pm2, "NiterInverse=0,TolInverse=1e-3" );
      if( astGetI( pm2, "NiterInverse" ) != 0 ) stopit( 8037, status );
      pm2 = astAnnul( pm2 );
   }
```

`INFINITY` needs `#include <math.h>`; confirm it is included at the top of `ast_tester/testpolymap.c` and add it if not.

- [ ] **Step 2: Run the test to verify it fails**

Run: `cmake --build build-dev -j8 && ctest --test-dir build-dev -R testpolymap_c --output-on-failure`
Expected: FAIL at 8031.

- [ ] **Step 3: Implement the setters**

Replace `src/polymap.c:6982-6984` with:

```c
astMAKE_SET1(PolyMap,NiterInverse,int,niterinverse,(
            ( value < 0 ) ?
            ( astError( AST__ATTIN, "astSetNiterInverse(%s): Invalid value %d "
                        "supplied for NiterInverse (must be zero or positive).",
                        status, astGetClass(this), value ), this->niterinverse ) :
            ( ( value != this->niterinverse ) ? astClearIsSimple(this) : (void)0,
              value ) ))
```

Replace `src/polymap.c:7033-7035` with:

```c
astMAKE_SET1(PolyMap,TolInverse,double,tolinverse,(
            ( !isfinite( value ) || value <= 0.0 ) ?
            ( astError( AST__ATTIN, "astSetTolInverse(%s): Invalid value %g "
                        "supplied for TolInverse (must be positive and finite).",
                        status, astGetClass(this), value ), this->tolinverse ) :
            ( ( value != this->tolinverse ) ? astClearIsSimple(this) : (void)0,
              value ) ))
```

Add one sentence to each attribute's Description: "An error is reported if a negative value is supplied." for NiterInverse and "An error is reported if the value is not positive and finite." for TolInverse.

In ChebyMap's `IterInverse` (`src/chebymap.c`), change `valid = isfinite(tol) && tol > 0.0 && maxiter >= 0;` to `valid = 1;` and leave the following per-axis bound checks.

- [ ] **Step 4: Run the regression set**

Run the Global Constraints test command. Expected: pass (plus the known Task 9 `check`). Confirm `testchebymap_c` and `testchebyinverse` did not depend on a non-positive tolerance.

- [ ] **Step 5: Commit**

```bash
git add src/polymap.c src/chebymap.c ast_tester/testpolymap.c
git commit -m "Reject unusable TolInverse and NiterInverse values"
```

---

### Task 5: Nudge a singular seed once before giving up

**Files:**
- Modify: `src/chebymap.c` ChebyMap `IterInverse` (from Task 2)
- Test: `ast_tester/testchebyinverse.c`

**Interfaces:**
- Produces: a per-position `nudged` flag array and a `stale` counter inside ChebyMap `IterInverse`. Task 10 reads `stale` to decide whether the whole-batch forward transform can be skipped.

- [ ] **Step 1: Write the failing test**

Add to `ast_tester/testchebyinverse.c` before `main`:

```c
/* A forward series with no linear term has a singular Jacobian at the
   domain midpoint, which is where the fallback seed lands. The solver must
   step off that point rather than report every target as unsolved. */
static void singular_seed( int *status ) {
   double lo = 0, hi = 10;
   double coeffs[] = { 1, 1, 2 };
   double target[] = { 0, .5, -.5, .9, -.99 };
   double got[5], back[5];
   double lo2[] = { 0, 0 }, hi2[] = { 10, 10 };
   double coeffs2[] = { 1, 1, 2, 0,  1, 2, 0, 1 };
   double tx[] = { 0, .5, -.5 }, ty[] = { .2, -.4, .9 };
   double gx[3], gy[3], bx[3], by[3];
   AstChebyMap *cm;
   int i;

   cm = astChebyMap( 1, 1, 1, coeffs, 0, NULL, &lo, &hi, NULL, NULL, "", status );
   check( astGetI( cm, "TranInverse" ) == 1, "Pure T2 offers an inverse" );
   astTran1( cm, 5, target, 0, got );
   astTran1( cm, 5, got, 1, back );
   for( i = 0; i < 5; i++ ) {
      check( got[i] != AST__BAD, "Pure T2 target solved" );
      near( back[i], target[i], "Pure T2 round trip" );
   }
   cm = astAnnul( cm );

   cm = astChebyMap( 2, 2, 1, coeffs2, 0, NULL, lo2, hi2, NULL, NULL, "", status );
   astTran2( cm, 3, tx, ty, 0, gx, gy );
   astTran2( cm, 3, gx, gy, 1, bx, by );
   for( i = 0; i < 3; i++ ) {
      check( gx[i] != AST__BAD && gy[i] != AST__BAD, "T2,T1 target solved" );
      near( bx[i], tx[i], "T2,T1 round trip x" );
      near( by[i], ty[i], "T2,T1 round trip y" );
   }
   cm = astAnnul( cm );
}
```

The 2-D coefficient array describes output 1 = T2(x') and output 2 = T1(y'): each term is `{ coefficient, output, power_x, power_y }`.

Insert `astSet( cm, "TolInverse=1e-12", status );` immediately after each `astChebyMap` call in this function, so `near` (1e-11 relative) is satisfiable.

Call `singular_seed( status );` from `main` after `bounded_solver( status );`.

- [ ] **Step 2: Run the test to verify it fails**

Run: `cmake --build build-dev -j8 && ctest --test-dir build-dev -R testchebyinverse --output-on-failure`
Expected: FAIL with "Pure T2 target solved" for every target.

- [ ] **Step 3: Implement the nudge**

In ChebyMap `IterInverse`:

- Declare `int *nudged;` and `int stale;`. Allocate `nudged = astCalloc( npoint, sizeof( int ) );` with the other per-point arrays and free it with them. Set `stale = 0;` before the iteration loop.
- Where `if( sing ) { MarkUnsolved(...) }` handles a singular matrix, replace with:

```c
               if( sing ) {
                  if( !nudged[ ipoint ] ) {
/* Move once off a stationary point, towards the side of the box with
   more room, and evaluate the forward transformation there next time. */
                     nudged[ ipoint ] = 1;
                     stale++;
                     for( icoord = 0; icoord < ncoord; icoord++ ) {
                        xx = ptr_in[ icoord ][ ipoint ];
                        if( xx - lbnd[ icoord ] >= ubnd[ icoord ] - xx ) {
                           xx -= 0.25*width[ icoord ];
                        } else {
                           xx += 0.25*width[ icoord ];
                        }
                        ptr_in[ icoord ][ ipoint ] = xx;
                     }
                  } else {
                     MarkUnsolved( ptr_in, ncoord, ipoint, flags, &nconv );
                  }
```

- Leave the loop bound alone. A nudge costs the position one iteration, which the ChebyMap default of 10 from Task 6 absorbs.

Update the function prologue Description with one sentence: "A position whose Jacobian is singular is moved once by a quarter of each half-width before being declared unsolved."

- [ ] **Step 4: Run the regression set**

Run the Global Constraints test command. Expected: `testchebyinverse` passes the new checks (the Task 9 `check` is still pending), other tests unchanged.

- [ ] **Step 5: Commit**

```bash
git add src/chebymap.c ast_tester/testchebyinverse.c
git commit -m "Step a ChebyMap inverse off a singular seed before giving up"
```

---

### Task 6: Raise the ChebyMap NiterInverse default and pin default behavior

**Files:**
- Modify: `src/chebymap.c` (new `GetNiterInverse` override, vtab init, prologue History)
- Modify: `src/polymap.c:6960-6975` (NiterInverse attribute docs, ChebyMap section)
- Test: `ast_tester/testchebyinverse.c`

- [ ] **Step 1: Write the failing test**

Add before `main` in `ast_tester/testchebyinverse.c`:

```c
/* A user who constructs a ChebyMap with no options must get a usable
   inverse across the whole box. x' = T1(u) + 0.2 T2(u),
   y' = T1(v) + 0.1 T2(u) T1(v) on [0,10]^2 has a triangular Jacobian with
   positive diagonal everywhere, so every in-box target has exactly one
   solution. */
static void default_attributes( int *status ) {
   double lo[] = { 0, 0 }, hi[] = { 10, 10 };
   double coeffs[] = { 1, 1, 1, 0,   .2, 1, 2, 0,
                       1, 2, 0, 1,   .1, 2, 2, 1 };
   double u[900], v[900], x[900], y[900], bu[900], bv[900];
   int i, j, n = 0, nbad = 0;
   AstChebyMap *cm = astChebyMap( 2, 2, 2, coeffs, 0, NULL, lo, hi,
                                  NULL, NULL, "", status );
   check( astGetI( cm, "NiterInverse" ) == 10, "ChebyMap NiterInverse default" );
   for( i = 0; i < 30; i++ ) {
      for( j = 0; j < 30; j++ ) {
         u[n] = 10.0*i/29.0;
         v[n] = 10.0*j/29.0;
         n++;
      }
   }
   astTran2( cm, n, u, v, 1, x, y );
   astTran2( cm, n, x, y, 0, bu, bv );
   for( i = 0; i < n; i++ ) {
      if( bu[i] == AST__BAD || bv[i] == AST__BAD ) {
         nbad++;
      } else if( fabs( bu[i] - u[i] ) > 1e-5 || fabs( bv[i] - v[i] ) > 1e-5 ) {
         check( 0, "Default-attribute round trip accuracy" );
      }
   }
   check( nbad == 0, "Default-attribute round trip leaves no BAD points" );
   cm = astAnnul( cm );
}
```

Call `default_attributes( status );` from `main`.

- [ ] **Step 2: Run the test to verify it fails**

Run: `cmake --build build-dev -j8 && ctest --test-dir build-dev -R testchebyinverse --output-on-failure`
Expected: FAIL with "ChebyMap NiterInverse default" and a non-zero BAD count.

- [ ] **Step 3: Implement the override**

In `src/chebymap.c` add prototype `static int GetNiterInverse( AstPolyMap *, int * );` and static `static int (*parent_getniterinverse)( AstPolyMap *, int * );`. Add after `GetObjSize`:

```c
static int GetNiterInverse( AstPolyMap *this, int *status ) {
/*
*  Name:
*     GetNiterInverse
*  Purpose:
*     Return the value of the NiterInverse attribute.
*  Type:
*     Private function.
*  Synopsis:
*     #include "polymap.h"
*     int GetNiterInverse( AstPolyMap *this, int *status )
*  Class Membership:
*     ChebyMap member function (over-rides the astGetNiterInverse
*     protected method inherited from the PolyMap class).
*  Description:
*     This function returns the NiterInverse value. An explicitly set
*     value is returned unchanged. The default for a ChebyMap is ten,
*     because the bounded algorithm checks the final candidate after the
*     last update and returns AST__BAD on exhaustion, so it needs more
*     headroom than the unbounded PolyMap algorithm, whose default is
*     four.
*/
   if( !astOK ) return 0;
   if( astTestNiterInverse( this ) ) return (*parent_getniterinverse)( this, status );
   return 10;
}
```

In the vtab initialiser add:

```c
   parent_getniterinverse = polymap->GetNiterInverse;
   polymap->GetNiterInverse = GetNiterInverse;
```

In `src/polymap.c` NiterInverse docs, change the ChebyMap paragraph sentence "The default is four, but stronger distortions may require more iterations." to "The default is ten. Stronger distortions may require more iterations." Add a History entry to `src/chebymap.c`.

- [ ] **Step 4: Run the regression set**

Run the Global Constraints test command. Expected: `default_attributes` passes. If the BAD count is not zero, print the failing grid positions and confirm they are not on the box edge before changing anything else; the spec expects zero.

- [ ] **Step 5: Commit**

```bash
git add src/chebymap.c src/polymap.c ast_tester/testchebyinverse.c
git commit -m "Give ChebyMap a NiterInverse default suited to bounded iteration"
```

---

### Task 7: Exclude ChebyMaps from the FitsChan SIP writer

**Files:**
- Modify: `src/fitschan.c` (include near the other class headers; the search loop in `SIPIntWorld`, around line 28925; prologue History)
- Test: `ast_tester/testfitschan.c`

- [ ] **Step 1: Write the failing test**

Add before `main` in `ast_tester/testfitschan.c`:

```c
/* A ChebyMap holds Chebyshev coefficients, which SIP cannot represent.
   FitsChan must decline to write a FITS-WCS header for a FrameSet whose
   pixel-to-projection Mapping contains one, even though the ChebyMap
   reports an (iterative) inverse. */
static void testsipcheby( int errbase, int *status ) {
   double coeffs[] = { 100, 1, 1, 0,  5, 1, 2, 0,  100, 2, 0, 1,  5, 2, 0, 2 };
   double lbnd[] = { -100, -100 }, ubnd[] = { 100, 100 };
   double shift[] = { -50, -50 };
   double matrix[] = { 1e-4, 0, 0, 1e-4 };
   AstFrame *grid = astFrame( 2, "Domain=GRID" );
   AstSkyFrame *sky = astSkyFrame( "" );
   AstChebyMap *cm = astChebyMap( 2, 2, 2, coeffs, 0, NULL, lbnd, ubnd,
                                  NULL, NULL, "" );
   AstShiftMap *sm = astShiftMap( 2, shift, "" );
   AstMatrixMap *mm = astMatrixMap( 2, 2, 0, matrix, "" );
   AstWcsMap *wm = astWcsMap( 2, AST__TAN, 1, 2, "" );
   AstCmpMap *c1, *c2, *c3;
   AstFrameSet *fs;
   AstFitsChan *fc;
   int nwrite;

   if( !astGetI( cm, "TranInverse" ) ) stopit( errbase, "ChebyMap has no inverse", status );
   astInvert( wm );
   c1 = astCmpMap( sm, cm, 1, "" );
   c2 = astCmpMap( c1, mm, 1, "" );
   c3 = astCmpMap( c2, wm, 1, "" );
   fs = astFrameSet( grid, "" );
   astAddFrame( fs, AST__BASE, c3, sky );
   fc = astFitsChan( NULL, NULL, "Encoding=FITS-WCS" );
   nwrite = astWrite( fc, fs );
   if( nwrite != 0 ) stopit( errbase + 1, "FITS-WCS written for a ChebyMap", status );
   fc = astAnnul( fc );
   fs = astAnnul( fs );
   c3 = astAnnul( c3 );
   c2 = astAnnul( c2 );
   c1 = astAnnul( c1 );
   wm = astAnnul( wm );
   mm = astAnnul( mm );
   sm = astAnnul( sm );
   cm = astAnnul( cm );
   sky = astAnnul( sky );
   grid = astAnnul( grid );
}
```

Call `testsipcheby( 12500, status );` from `main` after the `testgrismorder` calls.

- [ ] **Step 2: Run the test to verify it fails**

Run: `cmake --build build-dev -j8 && ctest --test-dir build-dev -R testfitschan_c --output-on-failure`
Expected: FAIL with "Error 12501: FITS-WCS written for a ChebyMap".

- [ ] **Step 3: Implement the exclusion**

In `src/fitschan.c` add `#include "chebymap.h"` next to `#include "polymap.h"`. In the `SIPIntWorld` search loop (the `if( astIsAPolyMap( map_list[ imap ] ) )` around line 28925) change the condition to:

```c
            if( astIsAPolyMap( map_list[ imap ] ) &&
                !astIsAChebyMap( map_list[ imap ] ) ) {
```

and add a comment above it:

```c
/* A ChebyMap is a PolyMap but its coefficients are Chebyshev, not
   monomial, so the SIP analysis below would misread them. Skip it. */
```

Add a History entry to the `src/fitschan.c` prologue.

- [ ] **Step 4: Run the regression set**

Run the Global Constraints test command. Expected: pass. Also run `ctest --test-dir build-dev -R 'fits|wcsconv|header' --output-on-failure` and expect no change.

- [ ] **Step 5: Commit**

```bash
git add src/fitschan.c ast_tester/testfitschan.c
git commit -m "Keep ChebyMaps out of FitsChan SIP output"
```

---

### Task 8: Make tests 610 to 612 able to fail

**Files:**
- Modify: `ast_tester/testchebymap.c:312-330`

- [ ] **Step 1: Prove the current checks cannot fail**

Temporarily change `AST__NOBOX` to `AST__ATTIN` in the 610 check, build, run `ctest --test-dir build-dev -R testchebymap_c`. Expected: still passes. Revert.

- [ ] **Step 2: Rewrite the three checks**

Replace lines 318-330 with:

```c
      cm = astChebyMap( 1, 1, 1, box_coeffs, 0, NULL, NULL, &bhi, NULL, NULL,
                        " " );
      expected = ( *status == AST__NOBOX && !cm );
      astClearStatus;
      if( !expected ) stopit( 610, status );
      if( cm ) cm = astAnnul( cm );

      cm = astChebyMap( 1, 1, 1, box_coeffs, 0, NULL, &blo, NULL, NULL, NULL,
                        " " );
      expected = ( *status == AST__NOBOX && !cm );
      astClearStatus;
      if( !expected ) stopit( 611, status );
      if( cm ) cm = astAnnul( cm );

      cm = astChebyMap( 1, 1, 0, NULL, 1, box_coeffs, NULL, NULL, NULL, &bhi,
                        " " );
      expected = ( *status == AST__NOBOX && !cm );
      astClearStatus;
      if( !expected ) stopit( 612, status );
      if( cm ) cm = astAnnul( cm );
```

and declare `int expected;` alongside `box_coeffs` at the top of the block.

- [ ] **Step 3: Verify the checks now fail when wrong**

Repeat Step 1's temporary edit. Expected: `Error 610` and a non-zero exit. Revert the temporary edit, rebuild, confirm pass.

- [ ] **Step 4: Commit**

```bash
git add ast_tester/testchebymap.c
git commit -m "Let the missing-box ChebyMap tests fail"
```

---

### Task 9: One evaluable-domain reconstruction in ChebyMap

**Files:**
- Modify: `src/chebymap.c:1253-1318` (`IterBounds` becomes `AxisBounds`), `src/chebymap.c:399-433` (`ChebyDomain`), `src/chebymap.c:2064-2091` (`PolyTran` NULL bounds), the static `GetIterDomain`
- Test: `ast_tester/testchebymap.c`

**Interfaces:**
- Produces: `static int AxisBounds( double scale, double offset, double *lbnd, double *ubnd );` returning 1 with bounds the forward evaluator accepts, or 0 when no usable interval exists.

- [ ] **Step 1: Write the failing test**

Add to `ast_tester/testchebymap.c` near the other domain tests (after the 612 block from Task 8):

```c
/* Bounds reported by astChebyDomain must be evaluable, and astPolyTran
   without user bounds must be able to sample at them. The box [0.1,1.3]
   reconstructs an upper bound one ulp outside the evaluator's range. */
   if( *status == 0 ) {
      double lin[] = { 1.0, 1, 1 };
      double blo = 0.1, bhi = 1.3, dlo, dhi, xin[2], xout[2];
      AstPolyMap *fit;

      cm = astChebyMap( 1, 1, 1, lin, 0, NULL, &blo, &bhi, NULL, NULL, " " );
      astChebyDomain( cm, 1, &dlo, &dhi );
      xin[0] = dlo;
      xin[1] = dhi;
      astTran1( cm, 2, xin, 1, xout );
      if( xout[0] == AST__BAD || xout[1] == AST__BAD ) stopit( 620, status );
      if( fabs( dlo - blo ) > 1e-12 || fabs( dhi - bhi ) > 1e-12 ) stopit( 621, status );
      fit = astPolyTran( cm, 0, 1e-8, 1e-6, 6, NULL, NULL );
      if( !fit ) stopit( 622, status );
      if( fit ) fit = astAnnul( fit );
      cm = astAnnul( cm );
   }
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `cmake --build build-dev -j8 && ctest --test-dir build-dev -R testchebymap_c --output-on-failure`
Expected: FAIL at 620 (upper bound evaluates to BAD).

- [ ] **Step 3: Implement**

Rename `IterBounds` to `AxisBounds` everywhere in `src/chebymap.c` (prototype, definition, prologue name, the two calls in `GetIterDomain`). Delete the line `if( !( lo <= hi ) ) return 0;`. Change the prologue Purpose to "Recover the evaluable interval on one axis from a Chebyshev normalization."

In `ChebyDomain`, replace both `scale[iax] != 0.0 ? (…)/scale : AST__BAD` blocks (lines 401-408 and 421-428) with:

```c
         if( !AxisBounds( scale[ iax ], offset[ iax ], lbnd + iax, ubnd + iax ) ) {
            lbnd[ iax ] = AST__BAD;
            ubnd[ iax ] = AST__BAD;
         }
```

and the equivalent for `scale_o`, `offset_o`, `lbnd_o`, `ubnd_o`.

In `PolyTran` (lines 2064-2091) replace the two `( -1.0 - offset[k] )/scale[k]` and `( 1.0 - offset[k] )/scale[k]` computations with:

```c
   } else if( scale && offset ) {
      for( k = 0; k < nax; k++ ) {
         if( !AxisBounds( scale[ k ], offset[ k ], this_lbnd + k, this_ubnd + k ) && astOK ) {
            astError( AST__NOBOX, "astPolyTran(%s): The %s transformation "
                      "has no usable bounding box on axis %d.", status,
                      astGetClass( this ), word, k + 1 );
         }
      }
```

Because both bounds now come from one call, restructure so `this_lbnd` and `this_ubnd` are filled together when neither user array is supplied, and the user-supplied array is copied when it is; keep the two existing `AST__NOBOX` errors for the case of no scale and no user bounds. Declare `int k;`.

- [ ] **Step 4: Run the regression set**

Run the Global Constraints test command. Expected: all pass, including the `testchebyinverse` check "Unresolvable box has no evaluable domain" pending since Task 2, and `transform_oracle_simplify` (bounds shift by ulps only; if a section differs, inspect it and regenerate in Task 12, do not edit the oracle by hand).

- [ ] **Step 5: Commit**

```bash
git add src/chebymap.c ast_tester/testchebymap.c
git commit -m "Reconstruct every ChebyMap bound through one evaluable-domain helper"
```

---

### Task 10: Remove redundant work from the bounded hot path

**Files:**
- Modify: `src/chebymap.c` `IterSteps` and ChebyMap `IterInverse` (from Tasks 2 and 5)
- Test: `ast_tester/testchebyinverse.c` (existing checks pin results)

**Interfaces:**
- Consumes: the `stale` counter from Task 5.
- Produces: `IterSteps` gains one parameter, `double **work`, the array of forward values for the whole batch (`ptr_work` in `IterInverse`), into which it copies the accepted trial's forward values.

- [ ] **Step 1: Record the baseline**

Run `ctest --test-dir build-dev -R 'testchebyinverse|testchebymap_c' --output-on-failure` and confirm pass. Save the output of `build-dev/ast_tester/testchebyinverse` for comparison.

- [ ] **Step 2: Allocate the trial PointSets once**

In `IterSteps` move `trial = astPointSet( nsearch, ncoord, "", status ); trial_out = astPointSet( nsearch, ncoord, "", status ); x = astGetPoints( trial );` to before the backtrack loop, using the initial `nsearch`. Inside the loop replace the per-level creation with `astSetNpoint( trial, nstart ); astSetNpoint( trial_out, nstart );` at the top of each level (the count never grows, and `astSetNpoint` only shrinks). Keep the `if( ntrial < nstart )` shrink. Move the two `astAnnul` calls after the loop. Keep the `if( !astOK ) break;` guard by testing `astOK` after the initial allocation and returning early after freeing `index`.

- [ ] **Step 3: Keep the accepted forward values**

Add the parameter `double **work` to `IterSteps` (prototype, definition, prologue Parameters, and the call in `IterInverse`, passing `ptr_work`). In the acceptance block, after copying `x[ i ][ jpoint ]` into `inputs`, add:

```c
            for( i = 0; i < ncoord; i++ ) {
               work[ i ][ ipoint ] = y[ i ][ jpoint ];
            }
```

In `IterInverse`, wrap the whole-batch forward transform at the top of the iteration:

```c
/* Every position still iterating was either just stepped by IterSteps,
   which recorded its forward value, or nudged, which did not. Only the
   first iteration and a nudge need a fresh evaluation of the batch. */
         if( iter == 0 || stale > 0 ) {
            (void) astTransform( map, result, fwd, work );
            stale = 0;
         }
```

The residual differencing loop that follows (`*pb = *pa - *pb`) already runs on `ptr_work`, so it now sees either the fresh transform or the copied trial values. Confirm that no code path between `IterSteps` and the top of the next iteration modifies `ptr_in` for an unflagged position other than the nudge.

- [ ] **Step 4: Compare against the baseline**

Run the regression set. Expected: pass, and `testchebyinverse` prints the same output as in Step 1. Run `ctest --test-dir build-dev -R transform_oracle_simplify` and expect no mismatch.

- [ ] **Step 5: Commit**

```bash
git add src/chebymap.c
git commit -m "Reuse trial PointSets and accepted forward values in the bounded inverse"
```

---

### Task 11: Sample from the pristine original in ReplaceTransformation

**Files:**
- Modify: `src/polymap.c:337` (prototype), `src/polymap.c:5525-5526` (call), `src/polymap.c:5543-5765` (`ReplaceTransformation`)

- [ ] **Step 1: Change the signature**

`PolyTran` already holds the caller's unmodified PolyMap in `this` and passes a copy in `result`. Add a `AstPolyMap *source` parameter after `this`:

```c
static int ReplaceTransformation( AstPolyMap *, AstPolyMap *, int, double, double, int, const double *, const double *, int * );
```

Call site:

```c
   ok = ReplaceTransformation( result, this, forward, acc, maxacc, maxorder,
                               lbnd, ubnd, status );
```

In the function, delete `source = astCopy( this );` and `source = astAnnul( source );`, delete the local `AstPolyMap *source;`, and keep the two `SamplePoly*D( source, ... )` calls. Rewrite the comment above the former copy to say the sampling uses the caller's original PolyMap, which is never initialised by the fitting code, so every order samples the transformation as supplied. Update the prologue Parameters to describe `source`.

- [ ] **Step 2: Run the regression set**

Run the Global Constraints test command plus `ctest --test-dir build-dev -R 'polytran|testpolymap|testchebymap' --output-on-failure`. Expected: pass. The `caches` test in `testchebyinverse` ("Forward fit invalidates the copied seed") covers the case the copy protected.

- [ ] **Step 3: Commit**

```bash
git add src/polymap.c
git commit -m "Sample the caller's PolyMap directly when fitting a replacement"
```

---

### Task 12: Sample ChebyMap fixtures over their box and regenerate the oracle

**Files:**
- Modify: `ast_tester/gen_transform_oracle.c:138-177` (`axis_bounds`)
- Modify: `ast_tester/fixtures/oracle/README.md` (domain rules)
- Regenerate: `ast_tester/fixtures/oracle/simplify_fixtures.oracle`

- [ ] **Step 1: Extend the native-domain rule**

`axis_bounds` receives no object. Add a parameter `AstObject *obj` (pass the loaded fixture object from the caller) and in the `DOM_NATIVE` branch:

```c
    } else {
        for ( int a = 0; a < naxis; a++ ) { lo[a] = -1000.0; hi[a] = 1000.0; }
/* A ChebyMap is defined only inside its forward box. Sampling the box gives
   its inverse real coverage; the symmetric default lands almost entirely on
   BAD. Compound Mappings that contain a ChebyMap keep the default. */
        if ( obj && astIsAChebyMap( obj ) ) {
            double *blo = astMalloc( naxis*sizeof( double ) );
            double *bhi = astMalloc( naxis*sizeof( double ) );
            int usable = astOK;
            astChebyDomain( (AstChebyMap *) obj, 1, blo, bhi );
            for ( int a = 0; a < naxis && usable; a++ ) {
                if ( blo[a] == AST__BAD || bhi[a] == AST__BAD ) usable = 0;
            }
            if ( usable ) {
                for ( int a = 0; a < naxis; a++ ) { lo[a] = blo[a]; hi[a] = bhi[a]; }
            }
            blo = astFree( blo );
            bhi = astFree( bhi );
            if ( !astOK ) astClearStatus;
        }
    }
```

Include `ast.h` already provides `astIsAChebyMap` and `astChebyDomain`. Update every caller of `axis_bounds` to pass the object it is about to sample. `check_transform_oracle.c` regenerates rather than reads, so confirm whether it shares `axis_bounds` through `transform_oracle_util.c`; if it does, the change lives there and both binaries pick it up.

- [ ] **Step 2: Regenerate to temporary files and inspect**

```bash
cmake --build build-dev -j8 --target gen_transform_oracle
build-dev/ast_tester/gen_transform_oracle ast_tester/fixtures \
    /tmp/simplify.oracle /tmp/headers.oracle /tmp/framesets.oracle
diff ast_tester/fixtures/oracle/headers.oracle /tmp/headers.oracle && echo headers-unchanged
diff ast_tester/fixtures/oracle/framesets.oracle /tmp/framesets.oracle && echo framesets-unchanged
diff ast_tester/fixtures/oracle/simplify_fixtures.oracle /tmp/simplify.oracle | grep '^[<>] #' | sort | uniq -c
```

Expected: headers and framesets unchanged. The simplify diff touches only the sections for `cheby_linear_lone_reduce.map`, `cap_cheby_consolidate.simp`, `neg_cheby_quadratic_no_reduce.map`, `neg_chebymap_standalone.map` (new `dir=inverse` sections; the three top-level ChebyMaps also gain in-box forward samples). Any other changed section is a regression to investigate before continuing.

- [ ] **Step 3: Install and test**

```bash
cp /tmp/simplify.oracle ast_tester/fixtures/oracle/simplify_fixtures.oracle
ctest --test-dir build-dev -R transform_oracle --output-on-failure
```

Expected: pass.

- [ ] **Step 4: Document the rule**

In `ast_tester/fixtures/oracle/README.md`, in the paragraph describing native-dump sampling, add: "A native dump whose top-level object is a ChebyMap is sampled over its forward bounding box from `astChebyDomain`, because a Chebyshev series is undefined outside it. Compound Mappings containing a ChebyMap keep the symmetric default."

- [ ] **Step 5: Commit**

```bash
git add ast_tester/gen_transform_oracle.c ast_tester/transform_oracle_util.c \
        ast_tester/fixtures/oracle/README.md ast_tester/fixtures/oracle/simplify_fixtures.oracle
git commit -m "Pin the ChebyMap iterative inverse in the simplify oracle"
```

(Omit `transform_oracle_util.c` from `git add` if it was not modified.)

---

### Task 13: Documentation and conventions

**Files:**
- Modify: `src/polymap.h` History (verify the Task 1 entry is present and accurate after Task 2)
- Modify: `src/polymap.c`, `src/chebymap.c`, `src/fitschan.c` prologue History (verify one entry per task above)
- Modify: `docs/plans/PLAN-chebymap-iterative-inverse.md`
- Modify: `PLAN.md:87`
- Modify: `ast.news`

- [ ] **Step 1: Reflow the plan document**

Rewrite `docs/plans/PLAN-chebymap-iterative-inverse.md` so every sentence is on its own line. Use this check and require zero output:

```bash
grep -nP '(?<!e\.g|i\.e|etc|vs|\d)\. [A-Z]' docs/plans/PLAN-chebymap-iterative-inverse.md
```

Replace `centre` with `center`, `normalisation` with `normalization`, `behaviour` with `behavior`, `initialise` with `initialize` in that file. Add a closing section "Follow-up (2026-09-10)" that names the spec `docs/superpowers/specs/2026-09-10-chebymap-inverse-review-fixes-design.md` and states that the bounded algorithm now lives in `src/chebymap.c`.

- [ ] **Step 2: Fix `PLAN.md`**

Replace `behaviour` with `behavior` on line 87. Run `grep -n 'behaviour\|centre\|normalis' PLAN.md` and expect no output.

- [ ] **Step 3: Update `ast.news`**

In the entry added by this branch for ChebyMap iterative inverses, add: "Setting IterInverse to a non-zero value on a PolyMap or ChebyMap with no forward transformation, a non-positive or non-finite TolInverse, or a negative NiterInverse now reports an error. The ChebyMap default for NiterInverse is 10. FitsChan no longer writes SIP headers for FrameSets containing a ChebyMap."

- [ ] **Step 4: Verify the prologue histories**

```bash
for f in src/polymap.c src/polymap.h src/chebymap.c src/fitschan.c; do
  echo "== $f"; grep -n '10-SEP-2026' "$f" | head -3
done
```

Expected: at least one match per file. Add any missing entry.

- [ ] **Step 5: Full build and test**

```bash
cmake --build build-dev -j8 2>&1 | grep -c warning
ctest --test-dir build-dev --output-on-failure
```

Expected: the warning count is no higher than on `master` built the same way, and every test passes.

- [ ] **Step 6: Commit**

```bash
git add docs/plans/PLAN-chebymap-iterative-inverse.md PLAN.md ast.news src/polymap.h src/polymap.c src/chebymap.c src/fitschan.c
git commit -m "Record the ChebyMap inverse review fixes in the histories and plans"
```

---

## Deferred

- The `tnx-cheb` round-trip inaccuracy recorded in `ast_tester/fixtures/oracle/transform_oracle_overrides.txt` concerns the TPN projection, not ChebyMap, and is not part of this plan.
