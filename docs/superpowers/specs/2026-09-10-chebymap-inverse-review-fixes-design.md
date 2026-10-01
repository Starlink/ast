# ChebyMap Iterative Inverse: Review Fixes Design

Date: 2026-09-10.
Branch: `u/timj/chebymap-inverse`.
Scope: the ten findings and four cleanup items from the code review of this branch, plus the class-boundary question they raise.

## Problem

The branch gives ChebyMap an iterative inverse by extending PolyMap's `IterInverse` with a `bounded` flag and three new protected virtual methods (`astGetJacobian`, `astLinearGuess`, `astGetIterDomain`).
The bounded and unbounded Newton algorithms now share one function in `src/polymap.c`, with `if( bounded )` branches interleaved through it, and the bounded-only helper `IterSteps` also lives in `src/polymap.c`.

Several review findings are direct consequences of that arrangement.
The two algorithms have different stopping conventions, different exhaustion behavior, and different needs for `NiterInverse`, yet they share one default and one set of attribute checks.
The remaining findings are independent defects: a FitsChan SIP writer regression, an all-BAD result at a singular seed, inconsistent `IterInverse` attribute semantics, silent all-BAD output for unusable tolerances, tests that cannot fail, an oracle that was not regenerated, two definitions of the evaluable domain in one class, redundant work on the hot path, and documentation conventions not followed.

## Class-boundary decision

The bounded algorithm belongs to ChebyMap, and it is moved to `src/chebymap.c`.

`IterInverse` becomes the single protected virtual method of PolyMap for evaluating the original inverse transformation by iteration.
PolyMap's implementation returns to the unbounded algorithm as it stood on `master`, with the branch's two genuine bug fixes retained (discarding a partly built Jacobian after an error, and invalidating cached Jacobian and initial-guess Mappings when coefficients change).
ChebyMap overrides `IterInverse` with the bounded algorithm.
When the forward series has no Chebyshev normalization (a legacy ChebyMap holding an ordinary polynomial), the override calls the parent implementation.

The protected hooks `astGetJacobian` and `astLinearGuess` stay.
They are clean per-class overrides with no guards, `ast_tester/testchebyinverse.c` exercises them directly, and both algorithms consume them the same way.
The protected hook `astGetIterDomain` is removed from the PolyMap vtab and header, because only the bounded algorithm asks the question it answers.
ChebyMap's `GetIterDomain` becomes a static helper in `src/chebymap.c`.
When it reports no usable domain on some axis, ChebyMap's `IterInverse` calls the parent implementation, which is what the branch does today.
The Jacobian and initial-guess caches stay in the PolyMap structure (`jacobian`, `lintrunc`), because PolyMap's copy, delete, lock and size code already manages them and subclasses access parent members directly throughout AST.
`IterSteps` moves to `src/chebymap.c` unchanged in behavior.
The tests in `ast_tester/testchebyinverse.c` that drove the bounded solver through a fake PolyMap subclass with a `GetIterDomain` hook are rewritten to drive it through a ChebyMap on the box [-1, 1] with the same polynomial expressed in the Chebyshev basis.
The tests that called `astGetIterDomain` use the public `astChebyDomain`, which after F8 reports the same evaluable bounds.

Rationale.
Facts about the algorithm (seeding, backtracking, exhaustion, singular handling, work-space reuse) live with the class whose transformation they invert.
Facts about the attributes (`IterInverse`, `NiterInverse`, `TolInverse` validity and defaults) live in PolyMap, where the attributes are defined, and ChebyMap overrides only the one default that differs.
The cost is roughly one hundred lines of shared Newton skeleton (Jacobian evaluation and `palDmat` solve) appearing in both files.
That duplication is accepted because the two loops have different stopping rules and different per-point outcomes, and a shared skeleton with callbacks would reintroduce the guards this change removes.

## Fixes

Each item names the review finding it resolves.

### F1. FitsChan SIP writer must not accept a ChebyMap

`SIPIntWorld` in `src/fitschan.c` searches the simplified Mapping list for a PolyMap using `astIsAPolyMap`.
A ChebyMap passes that test, and now that forward-only square ChebyMaps report `TranInverse=1`, `AnalysePoly` proceeds and `astMergeShift` and `ScalePolyInputs` treat Chebyshev coefficients as monomials.
The search skips any Mapping for which `astIsAChebyMap` is true, so no SIP header is written and `astWrite` returns zero, as on `master`.
The exclusion is made at the search, the single point where the class is chosen, rather than in three downstream routines.

### F2. Singular Jacobian at the seed

When the affine and linear seeds are both unusable, ChebyMap seeds every target at the domain midpoint.
A forward series with no linear terms has a singular Jacobian there, so every position is marked BAD on the first iteration.
ChebyMap's `IterInverse` gives each position one nudge: when `palDmat` reports a singular matrix and the position has not been nudged, it is moved on every axis by one quarter of the axis half-width toward the side of the box with more room, the position is marked as needing a fresh forward evaluation, and iteration continues.
A second singularity marks the position BAD as before.

### F3. IterInverse attribute semantics

One rule, in PolyMap's `SetIterInverse`: a non-zero value is rejected with `AST__ATTIN` when the numbers of inputs and outputs differ, or when the original forward transformation is undefined.
The getter lets an explicitly set value win, as on `master`.
The default selects iteration only for a square Mapping that has forward coefficients and no inverse coefficients, so a PolyMap with no coefficients at all (which the constructor allows) defaults to zero.
With the setter enforcing the same condition, `astTest` and `astGet` agree and a dump round trip preserves an explicitly set value.
The loader in `astLoadPolyMap` treats a recorded non-zero `IterInv` on an object with no forward coefficients as unset, without error, so dumps written by earlier versions still load.
ChebyMap's `GetIterInverse` and `SetIterInverse` overrides, their parent pointers, and the loader comment in `src/chebymap.c` are deleted.

### F4. NiterInverse default for ChebyMap

ChebyMap overrides `GetNiterInverse` to return 10 when the attribute is unset.
PolyMap keeps 4.
The bounded algorithm checks the final candidate after the last update and returns BAD on exhaustion, so it needs more headroom than the unbounded algorithm, which returns the last iterate.
A new test round-trips a 30 by 30 grid through a two-dimensional map with default attributes and requires zero BAD results.
The attribute documentation states both defaults.

### F5. TolInverse and NiterInverse validation

PolyMap's setters reject `TolInverse` values that are not finite or not positive, and `NiterInverse` values that are negative, with `AST__ATTIN`.
This applies to both classes.
ChebyMap's `IterInverse` no longer needs its own checks of these two values; its remaining validity checks concern the bounds only.

### F6. Tests 610 to 612 must be able to fail

`stopit` returns without action when status is already set, and the following `astClearStatus` erases the failure.
The three checks adopt the idiom already used by tests 781 and 783: capture the expectation, clear status, then call `stopit` if the expectation failed.

### F7. Regenerate the simplify oracle and sample ChebyMap fixtures over their box

The four fixtures whose ChebyMaps gained an inverse have no `dir=inverse` sections in the committed oracle.
The oracle is regenerated and committed.
The generator's native-domain rule is extended: when a native fixture's top-level object is a ChebyMap with a defined forward box, the sampling domain is that box instead of plus or minus 1000, so the inverse is pinned over points where it is defined.
The README under `ast_tester/fixtures/oracle/` records the new rule.

### F8. One definition of the evaluable domain

`IterBounds` is renamed `AxisBounds` and becomes the only way `src/chebymap.c` reconstructs a bound from a scale and offset.
`ChebyDomain` (behind the public `astChebyDomain`), ChebyMap's `PolyTran` with NULL bounds, and the static `GetIterDomain` use it.
Bounds move inward by at most a few ulps, so a bound returned by `astChebyDomain` always evaluates without BAD.
An axis whose normalization describes no usable interval (zero or non-finite scale, or an interval too narrow to hold a representable evaluable value) is reported by `astChebyDomain` as `AST__BAD` on that axis, which the zero-scale case already did.
The unreachable `!( lo <= hi )` test is removed; the two `fabs` checks stay.

### F9. Hot path

In `IterSteps` the two trial PointSets are created once at the initial search size and shrunk with `astSetNpoint` at each backtrack level.
The forward values of an accepted trial are copied into the work array, and the whole-batch forward transform at the top of the next iteration is skipped when no position is stale.
A position becomes stale when it is nudged (F2).

### F10. Documentation conventions

`src/polymap.h` gains a History entry for the vtab change.
`src/polymap.c`, `src/chebymap.c` and `src/fitschan.c` prologue histories record this work.
`docs/plans/PLAN-chebymap-iterative-inverse.md` is reflowed to one sentence per line and uses American spelling.
`PLAN.md` uses American spelling.

### Cleanup items

`ReplaceTransformation` samples from the pristine original passed in by `PolyTran` instead of taking a second deep copy.
The five copies of the mark-unsolved loop in the bounded algorithm become one static helper in `src/chebymap.c`.
The hand-rolled `x == AST__BAD || !isfinite(x)` predicates in the moved code use one static inline helper.

## Testing

Every task runs the existing ChebyMap, PolyMap, threads, FitsChan and oracle tests in the developer build with warnings and sanitizers enabled.
New behavior is pinned in `ast_tester/testchebyinverse.c`, `ast_tester/testchebymap.c`, `ast_tester/testpolymap.c` and `ast_tester/testfitschan.c`.

## Out of scope

The TPN round-trip inaccuracy noted in the oracle overrides for `tnx-cheb` remains a separate investigation.
No change to the Fortran interface or to `ast.h`.
