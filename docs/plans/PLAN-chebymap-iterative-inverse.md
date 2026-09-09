Implement ChebyMap's iterative inverse by retaining the Chebyshev first-kind
basis and reusing PolyMap's Newton solver. Do not introduce negative orders or
second-kind coefficient storage. The derivative representation is a small part
of the work; initial guesses, domain constraints and convergence handling need
equal attention.

This is a design plan, not an implementation. Source inspection used local
commit `68b452f6` and the discussion and changed-file list of PR #80. Numerical
probes used the existing local build and NumPy 2.5.3.

**Correction to the motivating discussion.** The two `tnx-cheb` fixtures do not
contain a ChebyMap. `FitsChan::WATCoeffs` calls `Cheb2Poly` at
`src/fitschan.c:35618` and represents the result as TPN projection parameters.
The stored FrameSet contains a WcsMap with `Type = "TPN"` at
`ast_tester/fixtures/wcsconv/framesets/tnx-cheb-head.ast:106`; a freshly loaded
FITS header produces the same projection. Both have `TranInverse=1`.

A live round trip of pixel `(1024, 2048)` through either fixture returns roughly
`(1659.92538, -2813.10529)`, with finite coordinates. Their override comments
incorrectly attribute the failures to an absent ChebyMap inverse. Investigate
the TPN inverse, polynomial branches and fixture validity separately. Do not
remove those overrides as an acceptance condition for this feature.

**What PolyMap already provides.**

- `IterInverse` (`src/polymap.c:2457`) implements batched Newton iteration:
  evaluate the original forward transformation, form `target - f(x)`, evaluate
  the Jacobian, solve `J dx = residual` with `palDmat`, and update each point.
  It handles original versus inverted Mapping direction and tracks completed
  points separately. The forward evaluation already dispatches to ChebyMap's
  `PolyPowers` implementation.
- `GetJacobian` (`src/polymap.c:1949`) lazily caches one PolyMap per Jacobian
  column. Each takes all original inputs and returns all output derivatives
  with respect to one input. ChebyMap inherits from PolyMap, so ordinary
  first-kind derivative ChebyMaps fit the existing `AstPolyMap **` cache.
- `LinearGuess` (`src/polymap.c:2727`) caches an affine seed Mapping. Its present
  coefficient interpretation is specific to ordinary polynomials.
- Attribute dispatch, transformation availability, serialization and cleanup
  machinery already exist. Defaults are four iterations and relative tolerance
  `1e-6`; explicit inverse coefficients take precedence by default, while an
  explicit `IterInverse=1` selects iteration instead.

`GetJacobian` and `LinearGuess` are private static functions, not virtual
methods. Simply removing ChebyMap's always-zero getter would run the wrong
derivative and seed calculations. The MINPACK/Levenberg–Marquardt callbacks
elsewhere in PolyMap fit polynomial coefficients; their Jacobian is with
respect to coefficients, not input coordinates. They are not the inverse
solver to reuse.

**Derivative options and recommendation.**

| Option | Benefits | Costs and recommendation |
| --- | --- | --- |
| Differentiate into first-kind coefficients and cache derivative ChebyMaps | Exact polynomial algebra; existing storage, evaluator and cache lifecycle | Derivatives can have more terms. Recommended first implementation. |
| Evaluate values and coordinate derivatives together using recurrences | Shares basis tables; avoids derivative coefficient expansion | Requires a new numerical Jacobian hook and more restructuring of evaluation. Best alternative if high-order sparse maps make cache size or evaluation expensive. |
| Convert to an ordinary PolyMap and reuse its inverse | An existing low-dimensional conversion exists inside FitsChan | Monomial expansion can grow substantially and suffer cancellation; would also need domain and inverse-policy handling. Useful as a low-order test oracle, not the general implementation. |
| Store mixed first/second-kind terms, e.g. negative orders | Compact symbolic derivative of each input term | Changes order validation, indexing, evaluation, interchange assumptions and simplification. Unnecessary for this feature. |
| Finite-difference Jacobian via Mapping rate evaluation | Possible prototype with little symbolic mathematics | Extra evaluations, step-size sensitivity and boundary problems. Use as an independent interior test, not the production default. |

For a stored term `c * product(T[n_j](z_j))`, where
`z_j = scale_f[j] * x_j + offset_f[j]`, its derivative in input axis `j` is

`c * scale_f[j] * T'[n_j](z_j) * product(T[n_k](z_k), k != j)`.

There is no need to store `U`: expand the derivative in `T` instead. For degree
`n > 0`, visit `k = n-1, n-3, ...` down to zero or one. The coefficient of
`T[k]` in `T'[n]` is `2*n` when `k > 0`, and `n` when `k == 0`.
Thus `T1' = T0`, `T2' = 4*T1`, `T3' = 6*T2 + 3*T0`, and
`T4' = 8*T3 + 8*T1`. Only the differentiated axis changes order.

Accumulate equal output/order tuples and preserve sparsity. For densely
populated order sequences, grouping by output and all other axis orders allows
a backward coefficient recurrence instead of repeated term expansion. Avoid
allocating a full multidimensional coefficient tensor. A single sparse term
produces `ceil(n_j/2)` derivative terms along axis `j`, not a Cartesian product
of expansions on every axis.

This is the approach embodied by
[NumPy chebder](https://numpy.org/doc/stable/reference/generated/numpy.polynomial.chebyshev.chebder.html),
which supports differentiation along a coefficient-array axis and a scale
factor for a linear change of variable. GSL likewise produces a new Chebyshev
series with
[gsl_cheb_calc_deriv](https://www.gnu.org/software/gsl/doc/html/cheb.html).
No new runtime dependency is needed. Be careful when borrowing formulas:
some libraries use a half-weighted constant coefficient; AST uses `c0*T0`.

For the direct-evaluation alternative, differentiate the existing recurrence:
`D0=0`, `D1=1`, `D[n+1]=2*T[n]+2*z*D[n]-D[n-1]`, then multiply by the physical
axis scale. Form tensor-product derivatives without dividing by `T[n]`, which
can be zero. This also avoids the apparent endpoint singularities of formulas
that divide by `1-z*z`. Boost provides
[Chebyshev derivative and Clenshaw evaluation routines](https://www.boost.org/doc/libs/latest/libs/math/doc/html/math_toolkit/sf_poly/chebyshev.html)
as another implementation reference; Clenshaw evaluation is an optional later
performance/numerical refinement, not a prerequisite.

**Implementation sequence.**

1. Expose protected virtual hooks for Jacobian construction and the affine
   initial-guess Mapping in `polymap.h`/`polymap.c`, using AST's existing vtable
   conventions. Keep the parent implementations and behavior intact in this
   refactoring commit. Route the shared solver through the hooks. Add an
   optional protected iteration-domain hook: PolyMap has no finite bounds;
   ChebyMap supplies its original forward bounds. This keeps ChebyMap knowledge
   out of the parent class and permits one shared iteration loop.

2. Implement ChebyMap's Jacobian hook using the expansion above. Preserve the
   original forward normalization in every derivative map, and include the
   chain-rule scale in each coefficient. Prefer copying the normalization
   exactly to needlessly reconstructing it through rounded bounds. Disable
   iterative inversion on these internal derivative maps. Delegate to the
   parent for a forward transformation with `scale_f == NULL`: ChebyMap can
   contain an ordinary polynomial in either direction.

   A zero derivative column must be a defined zero Mapping, not a constructor
   with zero supplied coefficients (which means undefined transformation).
   An explicit zero constant record can establish the forward arrays, even
   though `AddCoeff` drops the zero term. Test this deliberately.

3. Implement a seed in physical input coordinates. Prefer the affine
   linearization at the forward domain midpoint: evaluate `f(c)` and `J(c)`,
   then solve `J(c) * (x0-c) = target-f(c)`. It accounts for higher Chebyshev
   terms that contribute to value or slope at the centre. When that matrix is
   singular, try the normalized constant/linear truncation; if unusable, use
   the midpoint as a documented fallback. An out-of-domain proposed seed must
   be projected into the box or replaced by an interior seed before evaluation.

   The truncation itself requires `A[i,j] += coeff*scale_f[j]` and
   `b[i] += coeff*offset_f[j]`, plus constant coefficients. It shares this
   arithmetic with the affine simplification fix discussed in PR #80, but
   seeding should not invoke `astSimplify`: simplification has other semantics,
   including treatment of explicit inverse coefficients and domain coverage.

4. Add safeguarded steps for finite-domain iteration. Work with normalized
   input steps, or equivalently scale Jacobian columns by domain half-widths.
   Limit trial steps to the box, backtrack until the scaled residual improves,
   and keep the current valid point when a trial fails. A tiny clipped step
   alone must not count as convergence. Preserve forward out-of-domain behavior;
   do not enable general extrapolation as a side effect of inversion. Specify
   a small roundoff allowance at endpoints and test it without accepting
   materially out-of-domain roots.

   Check input, residual, Jacobian and proposed steps for `AST__BAD` and
   non-finite values before linear algebra. Recognize an already satisfactory
   residual before rejecting a singular Jacobian. Verify the final forward
   residual as well as the correction size. For ChebyMap, measure input
   accuracy against domain widths so zero and very large translated coordinates
   do not distort the stopping rule; scale output residuals consistently with
   those input tolerances and the Jacobian. Document this ChebyMap-specific
   meaning of `TolInverse` without silently changing PolyMap's tolerance units.

   Return `AST__BAD` for individual unsolved ChebyMap points after exhaustion,
   invalid arithmetic or failure to find an acceptable step. A local numerical
   failure must not set the global AST status or contaminate other points.
   Configuration/allocation errors continue to use normal AST errors.

5. Enable the attribute once the preceding pieces work. Match PolyMap's default
   policy for square, forward-defined ChebyMaps without explicit inverse
   coefficients. Retain explicit inverse coefficients by default and respect
   `IterInverse=0`, explicit selection of iteration, `astClear`, and `Invert`.
   Reject an iterative inverse for unequal dimensions or a missing original
   forward transformation. An inverse-only map does not acquire the opposite
   iterative direction merely by setting `Invert`.

   Audit load/copy/dump and effective attributes: old dumps may contain an
   explicitly set IterInverse value that the old ChebyMap getter ignored.
   Document that enabling the new default changes inverse availability for
   existing forward-only square ChebyMaps. `TranInverse=1` indicates an
   available algorithm, not convergence or uniqueness at every target.

6. Audit derived caches after coefficient or normalization replacement,
   including `astPolyTran`, copying, merging and serialization. PolyMap already
   drops its Jacobian on copy but copies `lintrunc`; a subsequent fit can make
   that affine cache stale. Invalidate both caches when the forward definition
   changes. Review lock propagation and memory accounting for cached children,
   not just their destruction. Ensure new protected methods and generated
   interfaces are included in both supported build workflows.

7. Update ChebyMap's class/constructor text, IterInverse applicability,
   NiterInverse and TolInverse descriptions, and remove the obsolete storage
   TODO. Update the existing C and Fortran tests that require IterInverse to
   be zero. Add release notes explaining the default and failure semantics.
   Keep the TPN override diagnosis correction separate from the feature.

**PolyMap behavior that must be treated explicitly.** Its existing solver
uses only relative correction length for convergence, has no damping, passes
BAD Jacobian/residual values to `palDmat` without an explicit precheck, and
returns its last iterate when the iteration budget runs out. A probe of
`f(x)=x*x`, target `4`, `NiterInverse=1` returns `2.5` with residual `2.25`.
Consequently, reusing it unchanged would not provide the failure contract
proposed above. Keep legacy PolyMap exhaustion/tolerance behavior initially;
apply the finite-domain safeguards through the new domain-aware path. Any
global change to PolyMap's contract should be a separate tested change.

Four iterations may suffice for weak distortions but is not a general
convergence guarantee. Start with inherited defaults for API consistency and
measure representative ChebyMaps before considering a documented class-specific
default. A local solver also cannot promise the globally closest inverse branch.
For multiple roots, document dependence on the initial guess and require any
returned root to satisfy the forward equation and domain. Global root finding,
automatic inverse fitting and multistart search are outside this first scope.

**Independent analytic inverse fixtures.**

Use direct C arithmetic as the reference inverse, with neither an inverse fit
nor AST's derivative machinery involved. Construct forward-only ChebyMaps and
select iteration explicitly once supported. Test their inverse directly, before
any simplification, and also compare their forward evaluations against the
independent formulas. This prevents a matching forward/inverse error from
passing a round-trip-only test.

1. A monotonic quadratic with a square-root inverse. Let `x` be in `[0,10]`,
   `z=x/5-1`, and

   `y = T1(z) + T2(z)/8 = z + z*z/4 - 1/8`.

   Its derivative with respect to `z` is `1+z/2`, between `1/2` and `3/2`,
   so the inverse is unique throughout the domain. The image is
   `[-7/8,9/8]`. Use the rationalized expression

   `z = 2*(y+1/8)/(1+sqrt(1+y+1/8))`, then `x=5*(z+1)`.

   This avoids subtracting nearly equal quantities near `z=0`. Forward
   constructor records are `{1,1,1, 0.125,1,2}`. Generate both known input
   samples and independent targets across the image, including endpoints and
   values near `y=-1/8`. Targets outside the image must fail even where the
   square-root formula still has a real algebraic solution outside the domain.

2. A coupled 2D map made from two polynomial shears, giving an explicit
   polynomial inverse and higher-order/mixed terms. Normalize physical axes
   independently, for example `z=x/5-1` for `[0,10]` and `w=(y+3)/4-1` for
   `[-3,5]`. Choose `a=1/8`, `b=1/16` and define

   `u = z + a*T3(w)`

   `v = w + b*T2(u)`.

   Reverse the shears analytically:

   `w = v - b*(2*u*u-1)`

   `z = u - a*(4*w*w*w-3*w)`

   `x = 5*(z+1)`, `y = 4*(w+1)-3`.

   Each shear is globally invertible; their composition has normalized
   Jacobian determinant one. Thus this is a nonlinear, coupled, nonsingular
   reference with no branch-selection ambiguity. A single ChebyMap represents
   the composition exactly, using

   `u = T1(z) + a*T3(w)`

   `v = T1(w) + b*T2(z) + 4*a*b*T1(z)*T3(w)`
   `    + a*a*b*T6(w) + a*a*b*T0`.

   These explicit coefficients exercise derivatives through degree six and a
   mixed `T1*T3` term. Compute the reference forward map from the two shear
   formulas rather than from the expanded coefficient table. Require the
   analytic inverse to lie in the original input box; the valid output region
   is its curved image, not merely an axis-aligned output bounding box. Test
   valid interior/edge/corner targets and targets whose unique analytic inverse
   is outside the box. Check independent output samples as well as round trips.

   A simpler single-shear variant (`b=0`) isolates cross-axis derivatives and
   permits testing arbitrary Chebyshev orders with the same explicit inverse.

The reference formulas were checked over 101 quadratic samples and a 101-by-101
normalized 2D grid. These are planning checks of the formulas and the expanded
forward coefficients, not a test of the as-yet unimplemented iterative inverse.

**Acceptance tests.**

- Analytical derivatives for degrees 0 through 4, asymmetric/shifted bounds,
  mixed cross terms, duplicate terms, cancellations, constant outputs and zero
  columns. Compare higher-degree 1D/2D/3D derivatives with an independent
  coefficient recurrence; use finite differences only in the domain interior.
- The motivating `2*T1` on `[0,10]`: without simplification or fitted inverse,
  forward values of `0,1,10` invert to the same inputs. Also test `[-1,1]`,
  coupled affine maps and singular affine maps.
- Nonlinear monotonic 1D and mildly distorted coupled 2D/3D maps with known
  inputs. Include the midpoint, endpoints/faces/corners, differing axis widths,
  large domain offsets and near-zero coordinates. Assert both input recovery
  for injective examples and final forward residual.
- Overshooting Newton steps: `T1(z)+0.25*T3(z)` is strictly increasing on the
  interval, yet a step from zero toward the output at `z=0.9` proposes `3.816`.
  It must converge using safeguards rather than become BAD at the first trial.
- Unreachable outputs, singular points, multiple roots, too-small iteration
  budgets, BAD inputs, non-finite values and mixed valid/invalid batches.
  Reject false convergence at a boundary. Require BAD for unsolved points.
- Attribute defaults and set/clear behavior, explicit/fitted inverse priority,
  missing forward coefficients, unequal dimensions, inverted objects and
  ordinary-polynomial forward transformations inside ChebyMap.
- Copy and Channel round trips before/after warming caches; fit a new forward
  transformation after warming the affine cache; verify independent copies,
  cleanup and thread-lock transfer where the test infrastructure supports it.
- Run C/Fortran ChebyMap and PolyMap tests and the transform/simplification
  oracle suite. Review changed outputs individually. Benchmark cold and warm
  batches, derivative-cache size and iteration counts before claiming a
  performance improvement.

The planning probe checked the sparse first-kind derivative expansion against
NumPy for six differentiated axes across 1D, 2D and 3D coefficient arrays with
different physical scales; all checks passed. This validates the algebraic
choice, not an AST iterative-inverse implementation. The remaining engineering
risk is bounded nonlinear convergence and compatibility, rather than the
ability to represent a Chebyshev derivative.
