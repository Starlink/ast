# `astSimplify` Branch Coverage Pathways

**Branch**: `u/timj/simplify-testing`
**Maintained by**: AST developers
**Last inventory audit**: 2026-05-05

---

## What is `astSimplify`?

`astSimplify` is the core Mapping optimization routine in the AST library.
Given a Mapping (typically a compound `CmpMap` representing a chain of
coordinate transformations), it reduces it to a simpler equivalent by
repeatedly invoking each component's `MapMerge` virtual method until a
fixed point is reached.

The pipeline works as follows:

1. The input Mapping is decomposed into a flat list of atomic Mappings
   (series or parallel, depending on the CmpMap combination mode).
2. For each Mapping in the list, its `MapMerge` method is called. This
   method examines its neighbours in the list and attempts to merge,
   cancel, swap, or self-simplify.
3. If any `MapMerge` call succeeds (returns a modified list), the loop
   restarts from the beginning of the new list.
4. The process repeats until no `MapMerge` call modifies the list (fixed
   point).
5. The resulting list is reassembled into a single Mapping (or returned
   as-is if it contains only one element).

Each class's `MapMerge` implements class-specific simplification rules:
inverse-pair cancellation, absorption into neighbours, self-reduction to
simpler classes, or swapping past intervening Mappings to reach a merge
target.

### References

- [SUN/210](https://www.starlink.co.uk/docs/sun210.htx/sun210.html) --
  AST programmer's guide
- [SUN/211](https://www.starlink.co.uk/docs/sun211.htx/sun211.html) --
  AST class reference
- Fixture directory: `ast_tester/simplify_fixtures/`
- Design document:
  `docs/superpowers/specs/2026-05-05-simplify-branch-coverage-design.md`

Each table row below represents a distinct branch (code path) within a
class's `MapMerge` implementation that `astSimplify` can exercise.

---

## Coverage Status Legend

| Symbol | Meaning |
|--------|---------|
| `+` | Covered -- a fixture in `simplify_tests.txt` exercises this branch |
| `-` | NOT covered |

Reasons a branch may be uncovered:

- **no fixture** -- branch is reachable via public API but no test
  fixture has been written yet.
- **infeasible (protected constructor)** -- the class cannot be
  instantiated via the public C API (DssMap requires a DSS plate
  solution in a FitsChan; XphMap and IntraMap require internal
  registration that is not exposed).
- **class lacks astEqual** -- the output Mapping type does not implement
  `astEqual`, so the test harness cannot verify correctness of the
  simplified result.
- **requires deep internal nesting** -- the branch depends on internal
  list-position state or multi-level decomposition that cannot be
  triggered by a single focused fixture.
- **requires compound FrameSet** -- Region parallel merging requires
  Regions embedded in compound FrameSets with compatible domains.

The `Lines` column from the source inventory is intentionally omitted
from this document because line numbers are volatile across edits.

---

## Scenario Fixtures

These are pre-existing integration regressions that exercise multiple
branches across many classes simultaneously. They serve as integration
safety nets rather than targeted branch coverage.

| ID | Fixture | Type | Polarity | Status | Description | Trigger |
|---|---|---|---|---|---|---|
| scenario-01 | brad.map | scenario | positive | `+` | Multi-class integration regression from Brad Warren's FITS-WCS pipeline | Complex FrameSet with multiple WCS conversions |
| scenario-02 | lsst1.map | scenario | positive | `+` | LSST camera-geometry pipeline regression testing CmpMap decomposition | Nested CmpMap chain from LSST focal-plane transform |
| scenario-03 | rigby.map | scenario | positive | `+` | Deep nested Mapping chain from Jane Rigby's spectroscopic pipeline | Large multi-level CmpMap with WcsMaps, MatrixMaps, ShiftMaps |

---

## cmpmap.c

CmpMap's MapMerge works in three stages: (1) self-simplification via
astSimplify or decomposition into components if the combination mode
matches the list mode; (2) merging with a neighbouring CmpMap -- series
CmpMaps in a parallel list restructure into parallel-then-series,
parallel CmpMaps in a series list pair components by dimension;
(3) swapping a PermMap past a parallel CmpMap by restructuring
permutation arrays.

### Positive branches

| ID | Fixture | Type | Status | Description | Trigger |
|---|---|---|---|---|---|
| cmpmap-01 | cmpmap_self_simplify.map | focused | `+` | CmpMap simplifies on its own (astSimplify returns a different/simpler mapping) | CmpMap whose internal components simplify when composed |
| cmpmap-03 | cmpmap_nested_parallel_flatten.map | focused | `+` | CmpMap decomposed into components when combination mode matches list mode | Series CmpMap in a series list, or parallel CmpMap in a parallel list |
| cmpmap-07 | cmpmap_parallel_series_components.map | focused | `+` | Two series CmpMaps in parallel list restructured into parallel-then-series and at least one simplifies | Two series CmpMaps combined in parallel with simplifiable pairings |
| cmpmap-09 | cmpmap_parallel_in_series_merge.map | focused | `+` | Two parallel CmpMaps in series list paired by dimension, at least one pair simplifies | Two parallel CmpMaps in series with simplifiable component pairs |
| cmpmap-17 | cmpmap_perm_parallel_swap.map | focused | `+` | PermMap and parallel CmpMap swapped (no constants), producing reordered CmpMap + new PermMap | PermMap swapping two contiguous blocks feeding a parallel CmpMap |
| cmpmap-18 | cmpmap_perm_swap_aconstants.map | focused | `+` | PermMap swap with aconstants: first component gets all-constant outputs | PermMap with first block all constants before parallel CmpMap |
| cmpmap-19 | cmpmap_perm_swap_bconstants.map | focused | `+` | PermMap swap with bconstants: second component gets all-constant outputs | PermMap with second block all constants before parallel CmpMap |
| cmpmap-20 | reconstruct_right_assoc.map | cascade | `+` | Right-nested CmpMap chain reassociated so neighbouring components can merge | Right-associative CmpMap chain of WcsMap/ZoomMap/UnitMap/WcsMap |
| cmpmap-21 | series_parallel_absorb.map | cascade | `+` | UnitMap absorbed from nested series/parallel CmpMap, collapsing to WcsMap+ShiftMap | Nested CmpMap containing WcsMap, UnitMap and ShiftMap |
| cmpmap-22 | spec_cel_invert.map | scenario | `+` | Large celestial FITS-WCS pipeline collapses via repeated merge to CmpMap(WinMap, WcsMap) | Multi-class FITS-WCS conversion pipeline (>3 components) |
| cmpmap-23 | wcsconv_gappt_iwc_residue.map | scenario | `+` | GAPPT/IWC conversion residue: large pipeline reassociated and partly merged | Multi-class FITS-WCS conversion pipeline (>3 components) |

### Negative branches

| ID | Fixture | Type | Status | Description | Trigger |
|---|---|---|---|---|---|
| cmpmap-02 | -- | focused | `- (requires deep internal nesting)` | CmpMap does not simplify on its own (astSimplify returns same pointer unchanged) | CmpMap with irreducible components |
| cmpmap-04 | neg_cmpmap_mode_mismatch.map | focused | `+` | Guard rejects decomposition: CmpMap mode does not match list mode | Series CmpMap in a parallel list |
| cmpmap-05 | cmpmap_solo_after_unit.map | focused | `-` | Guard rejects merging: only one mapping in list (nmap <= 1) | Series CmpMap(UnitMap, parallel CmpMap); UnitMap elided leaves the parallel CmpMap alone |
| cmpmap-06 | neg_cmpmap_neighbour_nonexcmpmap.map | focused | `- (requires deep internal nesting)` | Guard rejects merging: neighbour is not a CmpMap | CmpMap adjacent to a non-CmpMap in the list |
| cmpmap-08 | -- | focused | `- (no fixture)` | Guard: re-arranged parallel CmpMaps do not simplify | Two series CmpMaps in parallel whose rearranged pairings remain irreducible |
| cmpmap-10 | -- | cascade | `- (no fixture)` | Guard: two CmpMaps are not both parallel, or list is not series | Two adjacent CmpMaps where at least one is series |
| cmpmap-11 | -- | cascade | `- (no fixture)` | Parallel-in-series pairing produces no simplification | Two parallel CmpMaps in series with all irreducible sub-mappings |
| cmpmap-12 | -- | focused | `- (requires deep internal nesting)` | Guard: earlier branch already succeeded (result != -1) | Any case where a prior branch already simplified |
| cmpmap-13 | -- | focused | `- (requires deep internal nesting)` | Guard: CmpMap at position 0 (no preceding neighbour) | CmpMap first in mapping list |
| cmpmap-14 | -- | focused | `- (no fixture)` | Guard: predecessor is not a PermMap | Parallel CmpMap preceded by non-PermMap |
| cmpmap-15 | -- | focused | `- (no fixture)` | Guard: CmpMap is series (not parallel) | Series CmpMap preceded by PermMap |
| cmpmap-16 | -- | focused | `- (no fixture)` | Guard: PermMap axes don't satisfy contiguous-block condition (canswap=0) | PermMap with non-contiguous axis routing before parallel CmpMap |

---

## unitmap.c

UnitMap's MapMerge: in series, removes UnitMap from the list. In parallel,
merges adjacent UnitMaps into one with combined dimensionality. Also handles
Invert flag normalization.

### Positive branches

| ID | Fixture | Type | Status | Description | Trigger |
|---|---|---|---|---|---|
| unitmap-01 | unit_invert_clear.map | focused | `+` | Single UnitMap with Invert flag set: flag is cleared | Lone UnitMap with invert_list=1 |
| unitmap-03 | unit_series_elision.map | focused | `+` | UnitMap removed from series composition | CmpMap(ShiftMap(2), UnitMap(2)), Series=1 |
| unitmap-04 | unit_parallel_merge.map | focused | `+` | Adjacent UnitMaps in parallel merged into one wider UnitMap | CmpMap(UnitMap(1), UnitMap(2)), Series=0 |
| unitmap-05 | unit_parallel_invert_clear.map | focused | `+` | Parallel UnitMap with no adjacent UnitMaps but Invert set: flag cleared | Single UnitMap in parallel with invert=1, flanked by non-UnitMaps |

### Negative branches

| ID | Fixture | Type | Status | Description | Trigger |
|---|---|---|---|---|---|
| unitmap-02 | neg_unit_lone.map | focused | `- (no fixture)` | Single UnitMap with Invert=0: no change needed | Lone UnitMap already canonical |
| unitmap-06 | -- | focused | `- (no fixture)` | Parallel UnitMap with no adjacent UnitMaps and Invert=0: no simplification | Single UnitMap in parallel already canonical |

---

## zoommap.c

ZoomMap's MapMerge: in series, accumulates zoom factors of adjacent
ZoomMaps/UnitMaps (product to ZoomMap, product=1 to UnitMap). If no
same-class merge, absorbs into neighbouring MatrixMap, WinMap, or ShiftMap.
In parallel, collects adjacent ZoomMaps into one ZoomMap (same factor) or
diagonal MatrixMap (different factors).

### Positive branches

| ID | Fixture | Type | Status | Description | Trigger |
|---|---|---|---|---|---|
| zoommap-01 | zoom_series_merge.map | focused | `+` | Adjacent series ZoomMaps multiplied into single ZoomMap | CmpMap(ZoomMap(2), ZoomMap(3)), Series=1 |
| zoommap-02 | zoom_series_cancel.map | focused | `+` | Adjacent ZoomMaps with product=1 collapse to UnitMap | CmpMap(ZoomMap(2), ZoomMap(2,Invert=1)), Series=1 |
| zoommap-03 | zoom_invert_normalize.map | focused | `+` | Single inverted ZoomMap normalized to forward-only (1/zoom) | ZoomMap with invert_list=1, no adjacent ZoomMaps |
| zoommap-05 | zoom_parallel_all_unit.map | focused | `+` | All parallel ZoomMaps/UnitMaps have factor 1 to UnitMap | CmpMap(UnitMap(2), ZoomMap(1,Nin=3)), Series=0 |
| zoommap-06 | zoom_parallel_same_factor.map | focused | `+` | All parallel ZoomMaps have same non-unity factor to single ZoomMap | CmpMap(ZoomMap(2,Nin=1), ZoomMap(2,Nin=2)), Series=0 |
| zoommap-07 | zoom_parallel_to_matrix.map | focused | `+` | Parallel ZoomMaps with different factors to diagonal MatrixMap | CmpMap(ZoomMap(2,Nin=1), ZoomMap(3,Nin=1)), Series=0 |
| zoommap-09 | zoom_absorb_prev_matrix.map | focused | `+` | ZoomMap absorbed into previous MatrixMap (elements scaled) | CmpMap(MatrixMap, ZoomMap), Series=1 |
| zoommap-10 | zoom_absorb_prev_win.map | focused | `+` | ZoomMap absorbed into previous WinMap (shifts and scales multiplied) | CmpMap(WinMap, ZoomMap), Series=1 |
| zoommap-11 | zoom_absorb_prev_shift.map | focused | `+` | ZoomMap absorbed into previous ShiftMap to WinMap | CmpMap(ShiftMap, ZoomMap), Series=1 |
| zoommap-12 | zoom_absorb_next_matrix.map | focused | `+` | ZoomMap absorbed into next MatrixMap (elements scaled) | CmpMap(ZoomMap, MatrixMap), Series=1 |
| zoommap-13 | zoom_absorb_next_win.map | focused | `+` | ZoomMap absorbed into next WinMap (scales only) | CmpMap(ZoomMap, WinMap), Series=1 |
| zoommap-14 | zoom_absorb_next_shift.map | focused | `+` | ZoomMap absorbed into next ShiftMap to WinMap | CmpMap(ZoomMap, ShiftMap), Series=1 |

### Negative branches

| ID | Fixture | Type | Status | Description | Trigger |
|---|---|---|---|---|---|
| zoommap-04 | neg_zoom_lone.map | focused | `- (no fixture)` | Single forward ZoomMap with no adjacent ZoomMaps: nothing to accumulate | Lone forward ZoomMap in series (falls through to absorb) |
| zoommap-08 | neg_zoom_parallel_lone.map | focused | `-` | Single ZoomMap in parallel with no adjacent ZoomMaps: no simplification | Lone ZoomMap in parallel |
| zoommap-15 | neg_zoom_no_absorb.map | focused | `+` | ZoomMap cannot be absorbed: neither neighbour is MatrixMap/WinMap/ShiftMap | CmpMap(SphMap, ZoomMap, MathMap), Series=1 |
| zoommap-16 | neg_zoom_lower_nonzoom.map | focused | `-` | Backward neighbour search stops at a non-ZoomMap/UnitMap lower neighbour; no merge | CmpMap(PermMap, ZoomMap), Series=1 |

---

## shiftmap.c

ShiftMap's MapMerge delegates to WinMap. It converts itself to WinMap with
unit scale, calls WinMap's MapMerge. If that only converts back to ShiftMap,
checks for Invert flag normalization.

### Positive branches

| ID | Fixture | Type | Status | Description | Trigger |
|---|---|---|---|---|---|
| shiftmap-01 | shift_invert_normalize.map | focused | `+` | Normalizes an inverted ShiftMap by negating shifts and clearing Invert | ShiftMap with Invert=1 |

---

## winmap.c

WinMap's MapMerge: (1) self-simplification -- all-zero-shift to MatrixMap,
all-unit-scale to ShiftMap; (2) series direct merge with WinMap/ZoomMap/
ShiftMap/MatrixMap(diag)/UnitMap, or merge with parallel CmpMap neighbour,
or swap past PermMap/MatrixMap/WcsMap to reach target; (3) parallel merge
with same classes.

### Positive branches

| ID | Fixture | Type | Status | Description | Trigger |
|---|---|---|---|---|---|
| winmap-01 | win_to_matrix.map | focused | `+` | WinMap with all shift=0 replaced by diagonal MatrixMap | Standalone WinMap with Scl only |
| winmap-02 | win_to_shift.map | focused | `+` | WinMap with all scale=1 replaced by ShiftMap | Standalone WinMap with Sft only |
| winmap-05 | win_win_series_merge.map | focused | `+` | WinMap + WinMap in series merged | CmpMap(WinMap, WinMap), Series=1 |
| winmap-06 | win_zoom_series_merge.map | focused | `+` | WinMap + ZoomMap in series merged (WinMap first) | CmpMap(WinMap, ZoomMap), Series=1 |
| winmap-07 | win_zoom_series_merge_rev.map | focused | `+` | ZoomMap + WinMap in series merged (ZoomMap first) | CmpMap(ZoomMap, WinMap), Series=1 |
| winmap-08 | win_shift_series_merge.map | focused | `+` | WinMap + ShiftMap in series merged (WinMap first) | CmpMap(WinMap, ShiftMap), Series=1 |
| winmap-09 | win_shift_series_merge_rev.map | focused | `+` | ShiftMap + WinMap in series merged (ShiftMap first) | CmpMap(ShiftMap, WinMap), Series=1 |
| winmap-10 | win_matrix_series_merge.map | focused | `+` | WinMap + diagonal MatrixMap series merged (WinMap first) | CmpMap(WinMap, MatrixMap[Diagonal]), Series=1 |
| winmap-11 | win_matrix_series_merge_rev.map | focused | `+` | Diagonal MatrixMap + WinMap series merged (MatrixMap first) | CmpMap(MatrixMap[Diagonal], WinMap), Series=1 |
| winmap-12 | win_unit_series_merge.map | focused | `+` | WinMap + UnitMap in series: UnitMap removed | CmpMap(WinMap, UnitMap), Series=1 |
| winmap-14 | win_cmpmap_parallel_merge.map | cascade | `+` | WinMap merges with neighbouring parallel CmpMap (lower neighbour) | CmpMap(parallel CmpMap, WinMap), Series=1 |
| winmap-15 | win_upper_cmpmap_parallel_merge.map | cascade | `+` | WinMap merges with neighbouring parallel CmpMap (upper neighbour) | CmpMap(WinMap, parallel CmpMap), Series=1 |
| winmap-18 | win_swap_past_matrix.map | cascade | `+` | WinMap swaps past MatrixMap to reach merge target | CmpMap(WinMap, MatrixMap, WinMap), Series=1 |
| winmap-19 | win_perm_swap_merge.map | cascade | `+` | WinMap swaps past PermMap to reach merge target | CmpMap(WinMap, PermMap, WinMap), Series=1 |
| winmap-20 | win_swap_past_wcsmap.map | cascade | `+` | WinMap swaps past WcsMap to reach merge target | CmpMap(WinMap, WcsMap, WinMap), Series=1 |
| winmap-26 | win_swap_simplifies.map | cascade | `+` | Swap accepted because swapped Mapping simplifies | PermMap that strips axes swapped past WinMap |
| winmap-27 | win_swap_outer_merge.map | cascade | `+` | Swap accepted because outer neighbours can merge after swap | Three-map series where outer pair merges once WinMap moves |
| winmap-29 | win_win_parallel_merge.map | focused | `+` | WinMap + WinMap in parallel merged | CmpMap(WinMap, WinMap), Series=0 |
| winmap-30 | win_zoom_parallel_merge.map | focused | `+` | WinMap + ZoomMap in parallel merged (WinMap first) | CmpMap(WinMap, ZoomMap), Series=0 |
| winmap-31 | win_zoom_parallel_merge_rev.map | focused | `+` | ZoomMap + WinMap in parallel merged (ZoomMap first) | CmpMap(ZoomMap, WinMap), Series=0 |
| winmap-32 | win_parallel_merge.map | focused | `+` | WinMap + ShiftMap in parallel merged (WinMap first) | CmpMap(WinMap, ShiftMap), Series=0 |
| winmap-33 | win_shift_parallel_merge_rev.map | focused | `+` | ShiftMap + WinMap in parallel merged (ShiftMap first) | CmpMap(ShiftMap, WinMap), Series=0 |
| winmap-34 | win_diagmatrix_parallel_merge.map | focused | `+` | WinMap + diagonal MatrixMap in parallel merged (WinMap first) | CmpMap(WinMap, DiagMatrixMap), Series=0 |
| winmap-35 | win_diagmatrix_parallel_merge_rev.map | focused | `+` | Diagonal MatrixMap + WinMap in parallel merged (MatrixMap first) | CmpMap(DiagMatrixMap, WinMap), Series=0 |
| winmap-36 | win_unit_parallel_merge.map | focused | `+` | WinMap + UnitMap in parallel merged (WinMap first) | CmpMap(WinMap, UnitMap), Series=0 |
| winmap-37 | win_unit_parallel_merge_rev.map | focused | `+` | UnitMap + WinMap in parallel merged (UnitMap first) | CmpMap(UnitMap, WinMap), Series=0 |

### Negative branches

| ID | Fixture | Type | Status | Description | Trigger |
|---|---|---|---|---|---|
| winmap-03 | neg_win_mixed_scale_shift.map | focused | `- (no fixture)` | Not all shifts are zero: MatrixMap replacement refused | WinMap with mixed shifts |
| winmap-04 | neg_win_mixed_scale_shift.map | focused | `- (no fixture)` | Not all scales are 1: ShiftMap replacement refused | WinMap with mixed scales, alone |
| winmap-13 | neg_win_nonmergeable_series.map | focused | `+` | Neither neighbour is a directly-mergeable class | CmpMap(WinMap, FullMatrixMap), Series=1 |
| winmap-16 | -- | cascade | `- (no fixture)` | CmpMap neighbour is series (not parallel): no merge | CmpMap(series CmpMap, WinMap), Series=1 |
| winmap-17 | -- | cascade | `- (no fixture)` | Parallel CmpMap split doesn't simplify: refused | Parallel CmpMap with non-simplifiable components next to WinMap |
| winmap-21 | -- | cascade | `- (no fixture)` | No higher neighbour exists (WinMap last in list) | WinMap at end of series list |
| winmap-22 | -- | cascade | `- (no fixture)` | No lower neighbour exists (WinMap first in list) | WinMap at start of series list |
| winmap-23 | -- | cascade | `- (no fixture)` | Forward scan hits non-swappable class before target | CmpMap(WinMap, SpecMap, WinMap), Series=1 |
| winmap-24 | -- | cascade | `- (no fixture)` | Backward scan hits non-swappable class before target | Same from other direction |
| winmap-25 | -- | cascade | `- (no fixture)` | Both swap directions find no reachable target | WinMap adjacent to non-swappable non-mergeable class |
| winmap-28 | neg_win_swap_no_simplify.map | cascade | `+` | Swap refused: neither swapped Mapping simplifies and no outer merge | WinMap + MatrixMap where swap produces equivalent complexity |
| winmap-38 | neg_win_nonmergeable_parallel.map | focused | `+` | Neither parallel neighbour is a mergeable class | CmpMap(WinMap, FullMatrixMap), Series=0 |

---

## matrixmap.c

MatrixMap's MapMerge: (1) self-simplification -- unit to UnitMap,
equal-diagonal to ZoomMap, zero-off-diagonal to Diagonal; (2) series merge
with MatrixMap/ZoomMap/PermMap/WinMap(diag)/UnitMap; (3) cascade swap past
WinMap/PermMap toward merge target or for local simplification.

### Positive branches

| ID | Fixture | Type | Status | Description | Trigger |
|---|---|---|---|---|---|
| matrixmap-01 | matrix_unit_to_unit.map | focused | `+` | Unit-form square MatrixMap replaced by UnitMap | MatrixMap(Form="Unit", Nin=Nout) |
| matrixmap-03 | matrix_diagonal_to_zoom.map | focused | `+` | Diagonal MatrixMap with equal elements replaced by ZoomMap | MatrixMap(diag, [4,4]) |
| matrixmap-05 | matrix_full_to_diagonal.map | focused | `+` | Full MatrixMap with zero off-diagonals replaced by Diagonal | MatrixMap(full, [2,0,0,3]) |
| matrixmap-07 | matrix_matrix_series_merge.map | focused | `+` | Two MatrixMaps in series merged via matrix multiplication | CmpMap(MatrixMap, MatrixMap), Series=1 |
| matrixmap-08 | matrix_zoom_series_merge.map | focused | `+` | MatrixMap + ZoomMap in series: elements scaled | CmpMap(MatrixMap, ZoomMap), Series=1 |
| matrixmap-09 | matrix_perm_series_merge.map | focused | `+` | MatrixMap + bidirectional PermMap in series merged via MatPerm | CmpMap(MatrixMap, PermMap[bidirectional]), Series=1 |
| matrixmap-11 | matrix_diagwin_series_merge.map | focused | `+` | Diagonal MatrixMap + WinMap in series merged via MatWin2 | CmpMap(MatrixMap[diag], WinMap), Series=1 |
| matrixmap-12 | matrix_unit_series_merge.map | focused | `+` | MatrixMap + UnitMap in series: UnitMap eliminated | CmpMap(MatrixMap, UnitMap), Series=1 |
| matrixmap-13 | matrix_swap_past_win.map | cascade | `+` | MatrixMap swaps past WinMap to reach merge target | CmpMap(MatrixMap, WinMap, MatrixMap), Series=1 |
| matrixmap-14 | matrix_swap_past_perm.map | cascade | `+` | MatrixMap swaps past PermMap to reach merge target | CmpMap(MatrixMap, PermMap, MatrixMap), Series=1 |
| matrixmap-15 | matrix_swap_win_simplifies.map | cascade | `+` | Swap produces simpler Mapping (local simplification) | CmpMap(PermMap[dropping axes], MatrixMap), Series=1 |

### Negative branches

| ID | Fixture | Type | Status | Description | Trigger |
|---|---|---|---|---|---|
| matrixmap-02 | neg_matrix_singular_diag.map | focused | `- (no fixture)` | Diagonal MatrixMap with NULL i_matrix (singular) cannot become ZoomMap | MatrixMap(diag, singular) |
| matrixmap-04 | neg_matrix_diag_unequal.map | focused | `+` | Diagonal MatrixMap with unequal elements cannot become ZoomMap | MatrixMap(diag, [2,3]) with non-mergeable neighbours |
| matrixmap-06 | neg_matrix_full_offdiag.map | focused | `+` | Full MatrixMap with non-zero off-diagonal cannot self-simplify | MatrixMap(full, [1,2,3,4]) alone |
| matrixmap-10 | neg_matrix_perm_not_bidirectional.map | focused | `+` | Adjacent PermMap has inconsistent fwd/inv axes: merge blocked | CmpMap(MatrixMap, PermMap[non-bidirectional]), Series=1 |
| matrixmap-16 | neg_matrix_swap_refused.map | cascade | `+` | Swap refused: neither swapped Mapping simplifies | CmpMap(WinMap, MatrixMap[full]), Series=1 |
| matrixmap-17 | neg_matrix_parallel_no_merge.map | focused | `+` | Parallel mode with non-self-simplifiable MatrixMap: returns -1 | MatrixMap(full) in parallel |

---

## permmap.c

PermMap's MapMerge: accumulates permutation effects of adjacent PermMaps/
UnitMaps. Series composes sequentially; parallel concatenates side-by-side.
Result may be UnitMap (null permutation), simplified PermMap (cleared invert,
simplified arrays), or unchanged.

### Positive branches

| ID | Fixture | Type | Status | Description | Trigger |
|---|---|---|---|---|---|
| permmap-01 | perm_series_merge.map | cascade | `+` | Two or more adjacent PermMaps in series compose into one | CmpMap(PermMap, PermMap), Series=1 |
| permmap-02 | perm_parallel_merge.map | cascade | `+` | Adjacent PermMaps/UnitMaps in parallel compose into one wider PermMap | CmpMap(PermMap, UnitMap), Series=0 |
| permmap-03 | perm_cancel_to_unit.map | focused | `+` | Composed PermMap reduces to UnitMap (both permutations null, nin==nout) | Two inverse PermMaps in series producing identity |
| permmap-04 | perm_invert_normalize.map | focused | `+` | Single PermMap with Invert flag normalized (flag cleared, arrays swapped) | Lone PermMap with invert_list=1 |
| permmap-05 | perm_array_simplify.map | focused | `+` | PermMap simplified: previously-stored array now null after composition | PermMap + UnitMap in series where array becomes identity |
| permmap-06 | perm_inperm_constant_fold.map | focused | `+` | PermMap simplified: inperm array differs after constant folding | PermMap with constants composed with routing PermMap |
| permmap-07 | perm_outperm_constant_fold.map | focused | `+` | PermMap simplified: outperm array differs after re-computation | Similar to permmap-06 affecting outperm |
| permmap-09 | perm_constant_propagation.map | cascade | `+` | Series composition propagates constants through merged PermMap | PermMap with constant outputs + PermMap routing those outputs |
| permmap-10 | perm_bad_propagation.map | cascade | `+` | Series composition propagates AST__BAD through (negative perm index) | PermMap with out-of-range indices in series |

### Negative branches

| ID | Fixture | Type | Status | Description | Trigger |
|---|---|---|---|---|---|
| permmap-08 | neg_perm_no_merge.map | focused | `+` | No simplification: single canonical PermMap with no mergeable neighbours | Lone forward PermMap, invert=0 |

---

## lutmap.c

LutMap's MapMerge: (1) linear LutMap (detected by GetLinear) replaced by
WinMap; (2) inverse pair cancellation to UnitMap.

### Positive branches

| ID | Fixture | Type | Status | Description | Trigger |
|---|---|---|---|---|---|
| lutmap-01 | lut_linear_to_win.map | focused | `+` | Linear non-inverted LutMap replaced by WinMap | LutMap with linear table values |
| lutmap-02 | lut_linear_to_win.map | focused | `+` | Linear inverted LutMap replaced by reversed WinMap | LutMap(linear, Invert=1) |
| lutmap-03 | lut_inverse_cancel.map | focused | `+` | LutMap cancels with equal upper-neighbour in opposite direction to UnitMap | CmpMap(LutMap, Inverse(LutMap)), Series=1 |
| lutmap-04 | lut_inverse_cancel.map | focused | `+` | LutMap cancels with equal lower-neighbour in opposite direction to UnitMap | CmpMap(Inverse(LutMap), LutMap), Series=1 |

### Negative branches

| ID | Fixture | Type | Status | Description | Trigger |
|---|---|---|---|---|---|
| lutmap-05 | neg_lut_nonlinear.map | focused | `+` | LutMap is not linear: WinMap replacement refused | LutMap with non-linear table |
| lutmap-06 | neg_lut_constant.map | focused | `-` | LutMap is linear but constant (b1==b2): WinMap impossible | LutMap([5,5,5,...]) |
| lutmap-07 | neg_lut_parallel.map | focused | `+` | Not in series: cancellation skipped | Two LutMaps in parallel |
| lutmap-08 | neg_lut_nonlut_neighbour.map | focused | `+` | Neither neighbour is a LutMap | CmpMap(LutMap, ZoomMap), Series=1 |
| lutmap-09 | neg_lut_different_tables.map | focused | `+` | Neighbouring LutMap not inverse-equal (different tables) | Two different LutMaps in opposite directions |

---

## polymap.c

PolyMap's MapMerge: (1) combine duplicate terms; (2) linearize if only
constant and first-power terms to ShiftMap+MatrixMap; (3) inverse-pair
cancellation.

### Positive branches

| ID | Fixture | Type | Status | Description | Trigger |
|---|---|---|---|---|---|
| polymap-01 | poly_duplicate_terms.map | focused | `+` | Duplicate-power terms combined, reducing coefficient count | PolyMap with coefficients sharing identical power vectors |
| polymap-02 | poly_duplicate_terms.map | focused | `+` | Linear PolyMap (after dedup) replaced by ShiftMap+MatrixMap | PolyMap with nin==nout, only constant and x^1 terms |
| polymap-06 | poly_inverse_cancel.map | focused | `+` | PolyMap + Inverse(PolyMap) cancel to UnitMap | CmpMap(PolyMap, Inverse(PolyMap)), Series=1 |

### Negative branches

| ID | Fixture | Type | Status | Description | Trigger |
|---|---|---|---|---|---|
| polymap-03 | neg_poly_no_forward.map | focused | `+` | Linearization refused: forward transform undefined (ncoeff_f NULL) | PolyMap defined only with inverse coefficients |
| polymap-04 | neg_poly_nin_ne_nout.map | focused | `+` | Linearization refused: nin != nout | PolyMap(2-in, 3-out) with linear terms |
| polymap-05 | neg_poly_nonlinear.map | focused | `+` | Linearization refused: term has power > 1 or multiple inputs | PolyMap with quadratic term |
| polymap-07 | neg_poly_parallel_nonlinear.map | focused | `+` | Inverse-cancel refused: combination is parallel | Two PolyMaps in parallel |
| polymap-08 | neg_poly_nonpoly_neighbour.map | focused | `+` | Inverse-cancel refused: neighbour is not PolyMap | CmpMap(PolyMap, ZoomMap), Series=1 |
| polymap-09 | neg_poly_same_direction.map | focused | `+` | Inverse-cancel refused: neighbour has same invert direction | Two forward PolyMaps in series |
| polymap-10 | neg_poly_different_coeffs.map | focused | `+` | Inverse-cancel refused: astEqual fails (different coefficients) | Two different PolyMaps in opposite directions |

---

## mathmap.c

MathMap's MapMerge: checks SimpFI/SimpIF permission, then compares function
text of adjacent MathMaps. If forward(first)==inverse(second) and vice versa,
pair cancels to UnitMap.

### Positive branches

| ID | Fixture | Type | Status | Description | Trigger |
|---|---|---|---|---|---|
| mathmap-01 | math_inverse_cancel.map | focused | `+` | Two MathMaps with matching fwd/inv text cancel to UnitMap | CmpMap(MathMap[SimpFI=1], Inverse(MathMap[SimpIF=1])), Series=1 |

### Negative branches

| ID | Fixture | Type | Status | Description | Trigger |
|---|---|---|---|---|---|
| mathmap-02 | -- | focused | `- (no fixture)` | Parallel mode: refuses | MathMaps in parallel |
| mathmap-03 | neg_math_lone.map | focused | `- (no fixture)` | No following Mapping (last in list) | Single MathMap |
| mathmap-04 | neg_math_nonmath_neighbour.map | focused | `+` | Neighbour is not a MathMap | CmpMap(MathMap, ZoomMap), Series=1 |
| mathmap-05 | neg_math_no_simpfi.map | focused | `+` | SimpFI/SimpIF not set: simplification refused | CmpMap(MathMap[SimpFI=0], Inverse(MathMap)), Series=1 |
| mathmap-06 | -- | focused | `- (no fixture)` | Dimension mismatch: nin(first) != nout(second) | MathMaps of different dimensionality |
| mathmap-07 | -- | focused | `- (no fixture)` | Forward function count mismatch | MathMaps with different output counts |
| mathmap-08 | -- | focused | `- (no fixture)` | Forward function text of first != inverse text of second | Different MathMaps with SimpFI set |
| mathmap-09 | -- | focused | `- (no fixture)` | Inverse function count mismatch | MathMaps whose inv/fwd counts differ |
| mathmap-10 | neg_math_inv_text_mismatch.map | focused | `-` | Inverse function text of first != forward text of second | MathMaps whose inv(first) != fwd(second) |

---

## ratemap.c

RateMap's MapMerge: simplifies encapsulated Mapping, then checks both
neighbours for inverse pair cancellation.

### Positive branches

| ID | Fixture | Type | Status | Description | Trigger |
|---|---|---|---|---|---|
| ratemap-01 | ratemap_simplify_interior.map | focused | `+` | Encapsulated Mapping simplifies to new RateMap with simplified interior | RateMap(CmpMap(ZoomMap,ZoomMap)) |
| ratemap-02 | ratemap_inverse_cancel.map | focused | `+` | RateMap cancels with equal lower-neighbour in opposite direction | CmpMap(Inverse(RateMap), RateMap), Series=1 |
| ratemap-03 | ratemap_inverse_cancel.map | focused | `+` | RateMap cancels with equal upper-neighbour in opposite direction | CmpMap(RateMap, Inverse(RateMap)), Series=1 |

### Negative branches

| ID | Fixture | Type | Status | Description | Trigger |
|---|---|---|---|---|---|
| ratemap-04 | neg_ratemap_parallel.map | focused | `+` | Not in series: skipped | RateMaps in parallel |
| ratemap-05 | neg_ratemap_lower_nonrate.map | focused | `-` | Lower neighbour not a RateMap | CmpMap(ZoomMap, RateMap), Series=1 |
| ratemap-06 | neg_ratemap_same_invert.map | focused | `-` | Lower neighbour same invert flag | Two forward RateMaps in series |
| ratemap-07 | neg_ratemap_diff_indices.map | focused | `-` | Lower neighbour different iin/iout indices | RateMaps with different axis indices |
| ratemap-08 | neg_ratemap_different_inner.map | focused | `+` | Lower neighbour non-equal encapsulated Mapping | RateMaps wrapping different Mappings |
| ratemap-09 | neg_ratemap_nonratemap_neighbour.map | focused | `+` | Upper neighbour not a RateMap | CmpMap(RateMap, ZoomMap), Series=1 |
| ratemap-10 | neg_ratemap_same_invert.map | focused | `-` | Upper neighbour same invert flag | Two forward RateMaps |
| ratemap-11 | neg_ratemap_diff_indices.map | focused | `-` | Upper neighbour different iin/iout | Different axis RateMaps |
| ratemap-12 | neg_ratemap_diff_inner.map | focused | `-` | Upper neighbour non-equal encapsulated Mapping | Different inner Mappings |

---

## intramap.c

IntraMap's MapMerge: checks following neighbour for same function/IntraFlag
and opposite direction, with AST__SIMPFI/SIMPIF permission.

### Positive branches

| ID | Fixture | Type | Status | Description | Trigger |
|---|---|---|---|---|---|
| intramap-01 | intramap_inverse_cancel.map | focused | `+` | Forward + Inverse IntraMap with SIMPFI cancels to UnitMap | CmpMap(IntraMap, Inverse(IntraMap)), Series=1, SIMPFI set |
| intramap-02 | -- | focused | `- (infeasible: protected constructor)` | Inverse + Forward IntraMap with SIMPIF cancels to UnitMap | CmpMap(Inverse(IntraMap), IntraMap), Series=1, SIMPIF set |

### Negative branches

| ID | Fixture | Type | Status | Description | Trigger |
|---|---|---|---|---|---|
| intramap-03 | -- | focused | `- (infeasible: protected constructor)` | Not in series or no following Mapping | IntraMaps in parallel |
| intramap-04 | -- | focused | `- (infeasible: protected constructor)` | Following Mapping is not IntraMap | CmpMap(IntraMap, ZoomMap), Series=1 |
| intramap-05 | -- | focused | `- (infeasible: protected constructor)` | Different transformation functions (ifun differs) | IntraMaps with different registered functions |
| intramap-06 | -- | focused | `- (infeasible: protected constructor)` | IntraFlag strings differ | Same function but different flags |
| intramap-07 | -- | focused | `- (infeasible: protected constructor)` | Dimension mismatch | Asymmetric IntraMaps |
| intramap-08 | -- | focused | `- (infeasible: protected constructor)` | Same direction (both forward or both inverse) | Two forward IntraMaps in series |
| intramap-09 | -- | focused | `- (infeasible: protected constructor)` | SIMPFI flag not set on forward-then-inverse pair | IntraMaps without SIMPFI permission |
| intramap-10 | -- | focused | `- (infeasible: protected constructor)` | SIMPIF flag not set on inverse-then-forward pair | IntraMaps without SIMPIF permission |

---

## splinemap.c

SplineMap's MapMerge: checks both neighbours for equal SplineMap in opposite
direction via astEqual.

### Positive branches

| ID | Fixture | Type | Status | Description | Trigger |
|---|---|---|---|---|---|
| splinemap-01 | spline_inverse_cancel.map | focused | `+` | Lower-neighbour SplineMap equal and opposite to UnitMaps | CmpMap(Inverse(SplineMap), SplineMap), Series=1 |
| splinemap-02 | spline_inverse_cancel.map | focused | `+` | Upper-neighbour SplineMap equal and opposite to UnitMaps | CmpMap(SplineMap, Inverse(SplineMap)), Series=1 |

### Negative branches

| ID | Fixture | Type | Status | Description | Trigger |
|---|---|---|---|---|---|
| splinemap-03 | neg_spline_parallel.map | focused | `+` | Not in series: skipped | SplineMaps in parallel |
| splinemap-04 | -- | focused | `- (no fixture)` | No neighbour (boundary) | Single SplineMap |
| splinemap-05 | neg_spline_nonspline_neighbour.map | focused | `+` | Neighbour is not a SplineMap | CmpMap(SplineMap, ZoomMap), Series=1 |
| splinemap-06 | neg_spline_same_direction.map | focused | `+` | Same invert flag (same direction) | Two forward SplineMaps in series |
| splinemap-07 | neg_spline_different_coeffs.map | focused | `+` | astEqual fails (different coefficients) | Different SplineMaps in opposite directions |

---

## normmap.c

NormMap's MapMerge: (1) simplify encapsulated Frame; (2) basic Frame to
UnitMap; (3) inverse-pair cancellation; (4) duplicate-NormMap elimination.

### Positive branches

| ID | Fixture | Type | Status | Description | Trigger |
|---|---|---|---|---|---|
| normmap-01 | -- | focused | `- (unreachable: API pre-simplifies)` | Encapsulated Frame simplifies to new NormMap with simplified Frame | NormMap wrapping compound Frame that simplifies |
| normmap-02 | normmap_basic_frame_to_unit.map | focused | `+` | NormMap encapsulating basic Frame replaced by UnitMap (astNorm is no-op) | NormMap wrapping plain Frame |
| normmap-03 | normmap_inverse_cancel.map | focused | `+` | NormMap cancels with inverse lower-neighbour NormMap | NormMap preceded by Inverse(NormMap) with same Frame |
| normmap-04 | normmap_inverse_cancel_upper.map | focused | `+` | NormMap cancels with inverse upper-neighbour NormMap | NormMap followed by Inverse(NormMap) with same Frame |
| normmap-05 | normmap_duplicate_elim.map | focused | `+` | Duplicate adjacent NormMaps (same Frame, same direction) to UnitMaps | Two identical NormMaps in series |

### Negative branches

| ID | Fixture | Type | Status | Description | Trigger |
|---|---|---|---|---|---|
| normmap-06 | -- | focused | `- (no fixture)` | Lower neighbour NormMap: invert flags not opposite | NormMap preceded by same-direction NormMap with different Frame |
| normmap-07 | neg_normmap_different_frames.map | focused | `+` | Lower inverse NormMap: Frames not equal | NormMap preceded by Inverse(NormMap) with different Frame |
| normmap-08 | neg_normmap_nonnorm_neighbour.map | focused | `+` | Upper neighbour is not a NormMap | NormMap followed by non-NormMap |
| normmap-09 | -- | focused | `- (no fixture)` | Upper inverse NormMap: Frames differ | NormMap followed by Inverse(NormMap) with different Frame |
| normmap-10 | -- | focused | `- (no fixture)` | Adjacent same-direction NormMap: Frames differ | Two NormMaps same direction, different Frames |
| normmap-11 | neg_normmap_parallel.map | focused | `+` | Parallel mode: no simplification beyond Frame-level | NormMap in parallel |
| normmap-12 | -- | focused | `- (no fixture)` | Non-basic Frame, doesn't simplify, not in series | NormMap(SkyFrame) in parallel |

---

## unitnormmap.c

UnitNormMap's MapMerge: merges with adjacent ShiftMap/WinMap (unit-scale) by
adjusting centre, or cancels with inverse UnitNormMap (same centre to UnitMap,
different to ShiftMap).

### Positive branches

| ID | Fixture | Type | Status | Description | Trigger |
|---|---|---|---|---|---|
| unitnormmap-01 | unitnormmap_shift_fwd_merge.map | focused | `+` | ShiftMap + forward UnitNormMap to new UnitNormMap with adjusted centre | ShiftMap followed by UnitNormMap(fwd) |
| unitnormmap-02 | unitnormmap_winmap_fwd_merge.map | focused | `+` | WinMap(unit scale) + forward UnitNormMap to new UnitNormMap with adjusted centre | WinMap(scale=1) followed by UnitNormMap(fwd) |
| unitnormmap-04 | unitnormmap_inv_shift_merge.map | focused | `+` | Inverse UnitNormMap + ShiftMap to new inverse UnitNormMap with adjusted centre | UnitNormMap(inv) followed by ShiftMap |
| unitnormmap-05 | unitnormmap_inv_winmap_merge.map | focused | `+` | Inverse UnitNormMap + WinMap(unit scale) to new inverse UnitNormMap | UnitNormMap(inv) followed by WinMap(scale=1) |
| unitnormmap-07 | unitnormmap_inverse_cancel.map | focused | `+` | Forward + Inverse UnitNormMap with same centre to UnitMap | UnitNormMap(fwd) + Inverse(UnitNormMap), same centre |
| unitnormmap-08 | unitnormmap_inv_fwd_cancel.map | focused | `+` | Inverse + Forward UnitNormMap with same centre to UnitMap | Inverse(UnitNormMap) + UnitNormMap(fwd), same centre |
| unitnormmap-09 | unitnormmap_diff_centre_to_shift.map | focused | `+` | Forward + Inverse UnitNormMap with different centres to ShiftMap | UnitNormMap(fwd,c1) + Inverse(UnitNormMap,c2) |

### Negative branches

| ID | Fixture | Type | Status | Description | Trigger |
|---|---|---|---|---|---|
| unitnormmap-03 | neg_unitnormmap_nonunit_scale.map | focused | `+` | WinMap(non-unit scale) + UnitNormMap: refused | WinMap(scale!=1) followed by UnitNormMap(fwd) |
| unitnormmap-06 | -- | focused | `- (no fixture)` | UnitNormMap(inv) + WinMap(non-unit scale): refused | UnitNormMap(inv) followed by WinMap(scale!=1) |
| unitnormmap-10 | -- | focused | `- (no fixture)` | Inverse + Forward with different centres: no merge (asymmetric) | Inverse(UnitNormMap,c1) + UnitNormMap(fwd,c2) |
| unitnormmap-11 | -- | focused | `- (no fixture)` | ShiftMap + UnitNormMap(inv): refused | ShiftMap followed by Inverse(UnitNormMap) |
| unitnormmap-12 | -- | focused | `- (no fixture)` | WinMap + UnitNormMap(inv): refused | WinMap followed by Inverse(UnitNormMap) |
| unitnormmap-13 | -- | focused | `- (no fixture)` | UnitNormMap(fwd) + ShiftMap: refused | Forward UnitNormMap followed by ShiftMap |
| unitnormmap-14 | -- | focused | `- (no fixture)` | UnitNormMap(fwd) + WinMap: refused | Forward UnitNormMap followed by WinMap |
| unitnormmap-15 | -- | focused | `- (no fixture)` | Two UnitNormMaps in same direction: refused | Two forward UnitNormMaps in series |
| unitnormmap-16 | -- | focused | `- (no fixture)` | Neighbour is not ShiftMap/WinMap/UnitNormMap: refused | UnitNormMap with non-mergeable neighbour |
| unitnormmap-17 | -- | focused | `- (no fixture)` | Parallel mode: never simplifies | UnitNormMap in parallel |
| unitnormmap-18 | -- | focused | `- (no fixture)` | No neighbour pair produces valid merge | UnitNormMap flanked by non-mergeable classes |

---

## selectormap.c

SelectorMap's MapMerge: (1) simplify internal regions; (2) inverse-pair
cancellation.

### Positive branches

| ID | Fixture | Type | Status | Description | Trigger |
|---|---|---|---|---|---|
| selectormap-01 | -- | focused | `- (unreachable: API pre-simplifies)` | Internal regions simplify to new SelectorMap with simplified regions | SelectorMap containing simplifiable Regions |
| selectormap-05 | selectormap_inverse_cancel.map | focused | `+` | Inverse-pair cancellation to UnitMap | SelectorMap + Inverse(SelectorMap), identical regions |

### Negative branches

| ID | Fixture | Type | Status | Description | Trigger |
|---|---|---|---|---|---|
| selectormap-02 | -- | focused | `- (no fixture)` | No region simplification and no adjacent SelectorMap | SelectorMap with already-simple regions, alone |
| selectormap-03 | -- | focused | `- (no fixture)` | No adjacent SelectorMap found in series | SelectorMap flanked by non-SelectorMaps |
| selectormap-04 | -- | focused | `- (no fixture)` | Adjacent SelectorMap not equal-and-opposite | Two SelectorMaps with different regions |
| selectormap-06 | -- | focused | `- (no fixture)` | Parallel mode: Phase 2 skipped | SelectorMap in parallel |

---

## switchmap.c

SwitchMap's MapMerge: (1) inverse-pair cancellation; (2) invert-flag
normalization (swap selectors, invert routes); (3) simplify internal
selectors/routes.

### Positive branches

| ID | Fixture | Type | Status | Description | Trigger |
|---|---|---|---|---|---|
| switchmap-01 | switchmap_inverse_cancel.map | focused | `+` | Inverse-pair cancellation to UnitMap | SwitchMap + Inverse(SwitchMap), identical selectors/routes |
| switchmap-03 | switchmap_invert_normalize.map | focused | `+` | Inverted SwitchMap normalized to non-inverted equivalent | SwitchMap with Invert=1 |
| switchmap-04 | switchmap_internal_simplify.map | focused | `+` | Internal selectors/routes simplify to new SwitchMap | SwitchMap with simplifiable route maps |

### Negative branches

| ID | Fixture | Type | Status | Description | Trigger |
|---|---|---|---|---|---|
| switchmap-02 | -- | focused | `- (no fixture)` | Adjacent SwitchMap not equal-and-opposite | Two SwitchMaps with different routes |
| switchmap-05 | -- | focused | `- (no fixture)` | Series: no adjacent match, non-inverted, internals don't simplify | Lone simple SwitchMap in series |
| switchmap-06 | -- | focused | `- (no fixture)` | Parallel: same as above | Lone simple SwitchMap in parallel |

---

## tranmap.c

TranMap's MapMerge: (1) invert normalization (swap+invert components);
(2) simplify each component; (3) equal-component detection (to single
Mapping); (4) adjacent TranMap series merge.

### Positive branches

| ID | Fixture | Type | Status | Description | Trigger |
|---|---|---|---|---|---|
| tranmap-01 | tranmap_invert_normalize.map | focused | `+` | Inverted TranMap normalized by swapping and inverting components | TranMap with Invert=1 |
| tranmap-02 | tranmap_component_simplify.map | focused | `+` | Component Mappings individually simplified, TranMap rebuilt | TranMap(CmpMap(Z,Z), UnitMap) |
| tranmap-03 | tranmap_equal_components.map | focused | `+` | Both components bidirectional and equal to single component Mapping | TranMap(ZoomMap[2], ZoomMap[2]) |
| tranmap-04 | tranmap_adjacent_merge.map | cascade | `+` | Two adjacent TranMaps in series merge by combining fwd/inv legs | CmpMap(TranMap(A,B), TranMap(C,D)), Series=1 |

### Negative branches

| ID | Fixture | Type | Status | Description | Trigger |
|---|---|---|---|---|---|
| tranmap-05 | neg_tranmap_parallel.map | focused | `+` | Parallel mode: adjacent merge skipped | TranMaps in parallel |
| tranmap-06 | neg_tranmap_oneway.map | focused | `-` | Equal-component check skipped: components lack bidirectional transforms | TranMap(OneWayFwd, OneWayInv) |
| tranmap-07 | neg_tranmap_unequal_components.map | focused | `+` | CmpMap(fwd,inv(inv)) doesn't simplify to UnitMap: components not equal | TranMap(ShiftMap(1), ShiftMap(2)) |
| tranmap-08 | neg_tranmap_nontranmap_neighbour.map | focused | `+` | Higher neighbour is not a TranMap | CmpMap(TranMap, ZoomMap), Series=1 |
| tranmap-09 | neg_tranmap_adjacent_no_merge.map | focused | `-` | Neither fwd nor inv series combination simplified | CmpMap(TranMap(A,B), TranMap(C,D)) with irreducible legs |

---

## timemap.c

TimeMap's MapMerge: gather adjacent TimeMap steps, eliminate inverse pairs by
argument count (0,1,2,3,5-arg), eliminate no-op steps.

### Positive branches

| ID | Fixture | Type | Status | Description | Trigger |
|---|---|---|---|---|---|
| timemap-01 | time_inverse_cancel.map | focused | `+` | Full cancellation: all steps cancel to UnitMap | Two TimeMaps with inverse steps |
| timemap-02 | time_partial_cancel.map | focused | `+` | Partial cancellation: some steps cancel, result is simplified TimeMap | TimeMap with 3+ steps, one pair cancels |
| timemap-03 | time_merge_no_cancel.map | cascade | `+` | Multi-map merge without step reduction: adjacent TimeMaps merged | Two TimeMaps with non-cancelling steps |
| timemap-04 | time_invert_normalize.map | focused | `+` | Invert-flag clearing: single inverted TimeMap rebuilt with invert=0 | Single TimeMap with Invert=1 |
| timemap-07 | time_noop_eliminate.map | focused | `+` | No-op step elimination: MJDTOMJD with zero offset removed | TimeMap with MJDTOMJD(0,0) step |
| timemap-08 | time_inverse_cancel.map | focused | `+` | 1-arg pair cancellation (TAITOTT+TTTOTAI etc.) | Adjacent 1-arg inverse steps with matching arg |
| timemap-09 | time_2arg_swapped_cancel.map | focused | `+` | 2-arg pair cancellation (swapped args: MJDTOJD+JDTOMJD) | Adjacent 2-arg steps with swapped matching args |
| timemap-10 | time_2arg_same_cancel.map | focused | `+` | 2-arg pair cancellation (same order: TAITOUTC+UTCTOTAI) | Adjacent 2-arg steps with same-order args |
| timemap-11 | time_3arg_cancel.map | focused | `+` | 3-arg pair cancellation (GMSTTOLMST+LMSTTOGMST) | Adjacent 3-arg steps with matching args |
| timemap-12 | time_5arg_cancel.map | focused | `+` | 5-arg pair cancellation (TTTOTDB+TDBTOTT) | Adjacent 5-arg steps with matching args |
| timemap-17 | time_run_identity.map | cascade | `+` | Run of three TimeMaps whose combined steps cancel to a UnitMap | Three adjacent TimeMaps forming an identity conversion |

### Negative branches

| ID | Fixture | Type | Status | Description | Trigger |
|---|---|---|---|---|---|
| timemap-05 | neg_time_parallel.map | focused | `+` | Parallel-mode guard | TimeMap in parallel |
| timemap-06 | neg_time_lone_forward.map | focused | `+` | No simplification: single forward TimeMap, no neighbours | Lone forward TimeMap |
| timemap-13 | neg_time_arg_mismatch.map | focused | `+` | 1-arg pair with mismatched argument | Steps with different DUT1 values |
| timemap-14 | -- | focused | `- (no fixture)` | 2-arg swapped pair with mismatched arguments | MJDTOJD + JDTOMJD with different offsets |
| timemap-15 | neg_time_3arg_mismatch.map | focused | `+` | 3-arg pair with mismatched arguments | GMSTTOLMST + LMSTTOGMST with different lon |
| timemap-16 | neg_time_5arg_mismatch.map | focused | `+` | 5-arg pair with mismatched arguments | TTTOTDB + TDBTOTT with one arg different |

---

## slamap.c

SlaMap's MapMerge: gather adjacent SlaMap steps, eliminate inverse pairs by
argument count (0,1,2,4-arg), eliminate redundant precession, merge adjacent
precession.

### Positive branches

| ID | Fixture | Type | Status | Description | Trigger |
|---|---|---|---|---|---|
| slamap-01 | sla_inverse_cancel.map | focused | `+` | Full cancellation: all steps cancel to UnitMap | Two SlaMaps with inverse steps |
| slamap-02 | sla_partial_cancel.map | focused | `+` | Partial cancellation: some steps cancel, result is simplified SlaMap | SlaMap with 3+ steps, one pair cancels |
| slamap-03 | sla_merge_no_cancel.map | cascade | `+` | Multi-map merge without step reduction: adjacent SlaMaps merged | Two SlaMaps with non-cancelling steps |
| slamap-04 | sla_invert_normalize.map | focused | `+` | Invert-flag clearing: single inverted SlaMap rebuilt with invert=0 | Single SlaMap with Invert=1 |
| slamap-07 | sla_inverse_cancel.map | focused | `+` | 0-arg pair: galactic (EQGAL+GALEQ) | Adjacent galactic steps |
| slamap-08 | sla_supergalactic_cancel.map | focused | `+` | 0-arg pair: supergalactic (GALSUP+SUPGAL) | Adjacent supergalactic steps |
| slamap-09 | sla_j2000_cancel.map | focused | `+` | 0-arg pair: dynamical J2000 (J2000H+HJ2000) | Adjacent J2000 steps |
| slamap-10 | sla_eterms_cancel.map | focused | `+` | 1-arg pair: E-terms (ADDET+SUBET) | Adjacent E-term steps, same epoch |
| slamap-11 | sla_fk45_cancel.map | focused | `+` | 1-arg pair: FK4/FK5 (FK45Z+FK54Z) | Adjacent FK conversion steps |
| slamap-12 | sla_icrs_cancel.map | focused | `+` | 1-arg pair: ICRS/FK5 (HFK5Z+FK5HZ) | Adjacent ICRS steps |
| slamap-13 | sla_ecliptic_cancel.map | focused | `+` | 1-arg pair: ecliptic (ECLEQ+EQECL) | Adjacent ecliptic steps |
| slamap-14 | sla_helioecl_cancel.map | focused | `+` | 1-arg pair: helio-ecliptic (EQHE+HEEQ) | Adjacent helio-ecliptic steps |
| slamap-15 | sla_ha_cancel.map | focused | `+` | 1-arg pair: HA (R2H+H2R) | Adjacent HA steps |
| slamap-16 | sla_geocentric_cancel.map | focused | `+` | 2-arg pair (cross-matched): geocentric (AMP+MAP) | Adjacent AMP/MAP with crossed args |
| slamap-17 | sla_azel_cancel.map | focused | `+` | 2-arg pair (same order): AzEl (DH2E+DE2H) | Adjacent AzEl steps |
| slamap-18 | sla_hpc_cancel.map | focused | `+` | 4-arg pair: helioprojective-Cartesian (HPCEQ+EQHPC) | Adjacent HPC steps |
| slamap-19 | sla_hpr_cancel.map | focused | `+` | 4-arg pair: helioprojective-Radial (HPREQ+EQHPR) | Adjacent HPR steps |
| slamap-23 | sla_prec_redundant.map | focused | `+` | Redundant precession: PREC/PREBN with start==end eliminated | PREC(2000,2000) |
| slamap-24 | sla_prec_merge.map | focused | `+` | Adjacent precession merge: PREC(a,b)+PREC(b,c) to PREC(a,c) | Adjacent PREC steps with common equinox |
| slamap-26 | sla_run_identity.map | cascade | `+` | Run of three SlaMaps whose combined steps cancel to a UnitMap | Three adjacent SlaMaps forming an identity conversion |
| slamap-27 | sla_run_partial.map | cascade | `+` | Run of three SlaMaps partially cancels to a single SlaMap | Three adjacent SlaMaps, a subset of steps cancel |

### Negative branches

| ID | Fixture | Type | Status | Description | Trigger |
|---|---|---|---|---|---|
| slamap-05 | neg_sla_parallel.map | focused | `+` | Parallel-mode guard | SlaMap in parallel |
| slamap-06 | neg_sla_lone_forward.map | focused | `+` | No simplification: single forward SlaMap, no neighbours | Lone forward SlaMap |
| slamap-20 | neg_sla_arg_mismatch.map | focused | `+` | 1-arg pair: mismatched argument | Steps with different epochs |
| slamap-21 | neg_sla_2arg_mismatch.map | focused | `+` | 2-arg pair: mismatched arguments | AMP+MAP with non-matching args |
| slamap-22 | neg_sla_4arg_mismatch.map | focused | `+` | 4-arg pair: mismatched arguments | HPC steps with one arg different |
| slamap-25 | neg_sla_prec_no_common.map | focused | `+` | Adjacent precession: non-common equinox prevents merge | PREC(1950,1975) + PREC(2000,2025) |

---

## specmap.c

SpecMap's MapMerge: gather adjacent SpecMap steps, eliminate inverse pairs by
argument count (0,1,2,3,6-arg).

### Positive branches

| ID | Fixture | Type | Status | Description | Trigger |
|---|---|---|---|---|---|
| specmap-01 | spec_inverse_cancel.map | focused | `+` | Full cancellation: all steps cancel to UnitMap | Two SpecMaps with inverse steps |
| specmap-02 | spec_partial_cancel.map | focused | `+` | Partial cancellation: some steps cancel, simplified SpecMap | SpecMap with 3+ steps, one pair cancels |
| specmap-03 | spec_merge_no_cancel.map | cascade | `+` | Multi-map merge without step reduction | Two adjacent SpecMaps with non-inverse steps |
| specmap-04 | spec_invert_normalize.map | focused | `+` | Invert-flag clearing | Single SpecMap with Invert=1 |
| specmap-08 | spec_unit_cancel.map | focused | `+` | 0-arg pair: unit conversions (ENTOFR+FRTOEN etc.) | Adjacent 0-arg inverse steps |
| specmap-09 | spec_inverse_cancel.map | focused | `+` | 1-arg pair (FRTOVL+VLTOFR) | Adjacent 1-arg steps, matching arg |
| specmap-10 | spec_lsr_cancel.map | focused | `+` | 2-arg pair: local-standard (LKF2HL+HLF2LK) | Adjacent 2-arg steps, matching args |
| specmap-11 | spec_geocentric_cancel.map | focused | `+` | 3-arg pair: geocentric/barycentric (GEF2HL+HLF2GE) | Adjacent 3-arg steps |
| specmap-12 | spec_topocentric_cancel.map | focused | `+` | 6-arg pair: topocentric (TPF2HL+HLF2TP) | Adjacent 6-arg steps |
| specmap-17 | spec_run4_identity.map | scenario | `+` | Run of four SpecMaps whose combined steps cancel to a UnitMap | Four adjacent SpecMaps forming an identity conversion |
| specmap-18 | spec_run_inverted_mid.map | cascade | `+` | Run of three SpecMaps with an inverted middle step partially cancels to one SpecMap | Three adjacent SpecMaps, middle one inverted |

### Negative branches

| ID | Fixture | Type | Status | Description | Trigger |
|---|---|---|---|---|---|
| specmap-05 | neg_spec_parallel.map | focused | `+` | Parallel-mode guard | SpecMap in parallel |
| specmap-06 | neg_spec_lone_forward.map | focused | `+` | No simplification: single forward SpecMap, no neighbours | Lone forward SpecMap |
| specmap-07 | -- | focused | `- (no fixture)` | Nin-mismatch guard: adjacent SpecMap with different nin | SpecMap(nin=1) adjacent to SpecMap(nin=3) |
| specmap-13 | neg_spec_arg_mismatch.map | focused | `+` | 1-arg pair: mismatched argument | Steps with different rest frequencies |
| specmap-14 | neg_spec_2arg_mismatch.map | focused | `+` | 2-arg pair: mismatched arguments | Steps with different RA/dec |
| specmap-15 | neg_spec_3arg_mismatch.map | focused | `+` | 3-arg pair: mismatched arguments | Steps with different epoch |
| specmap-16 | neg_spec_6arg_mismatch.map | focused | `+` | 6-arg pair: mismatched arguments | Steps with one arg different |

---

## wcsmap.c

WcsMap's MapMerge: (1) AST__WCSBAD to UnitMap; (2) inverse-pair cancellation
via CanMerge; (3) swap past PermMap toward merge target or for local
simplification.

### Positive branches

| ID | Fixture | Type | Status | Description | Trigger |
|---|---|---|---|---|---|
| wcsmap-01 | wcsmap_bad_to_unit.map | focused | `+` | WcsMap with AST__WCSBAD type replaced by UnitMap | WcsMap(AST__WCSBAD) |
| wcsmap-02 | wcs_inverse_cancel.map | focused | `+` | Adjacent inverse WcsMap pair cancels to UnitMap | CmpMap(WcsMap[TAN], Inverse(WcsMap[TAN])), Series=1 |
| wcsmap-03 | wcsmap_perm_swap_cancel.map | cascade | `+` | WcsMap swaps past PermMap to reach inverse merge target | WcsMap + PermMap + Inverse(WcsMap), Series=1 |
| wcsmap-04 | wcsmap_perm_swap_simplify.map | cascade | `+` | Swap with PermMap for local simplification (no merge target) | WcsMap + PermMap(adds axes) where swap simplifies |

### Negative branches

| ID | Fixture | Type | Status | Description | Trigger |
|---|---|---|---|---|---|
| wcsmap-05 | -- | cascade | `- (no fixture)` | Speculative swap refused: neither Mapping simplifies | WcsMap + PermMap(identity) with no merge target |
| wcsmap-06 | neg_wcs_parallel.map | focused | `+` | Parallel mode or nmap==1: refused | WcsMap in parallel |
| wcsmap-07 | neg_wcs_different_projection.map | focused | `+` | Neighbour not a WcsMap or different projection type | WcsMap[TAN] adjacent to WcsMap[SIN] |
| wcsmap-08 | neg_wcs_same_direction.map | focused | `+` | Two WcsMaps same invert direction | Two forward WcsMap[TAN] in series |
| wcsmap-09 | -- | focused | `- (no fixture)` | Lon/lat axis indices differ | WcsMap pair with different axis assignments |
| wcsmap-10 | neg_wcs_different_params.map | focused | `+` | Projection parameters differ | WcsMap pair with different PV values |
| wcsmap-11 | neg_wcs_nonperm_between.map | cascade | `+` | Intervening Mapping not a PermMap: swap blocked | WcsMap + MatrixMap + Inverse(WcsMap) |
| wcsmap-12 | -- | cascade | `- (no fixture)` | PermMap has non-bidirectional links: swap refused | WcsMap + PermMap(one-way) |
| wcsmap-13 | -- | cascade | `- (no fixture)` | PermMap doesn't pass through both lon/lat: swap refused | WcsMap + PermMap(disconnects one WCS axis) |

---

## sphmap.c

SphMap's MapMerge: (1) inverse-pair cancellation (UnitRadius or matching
PolarLong); (2) matrix sandwich -- Inv(SphMap)+DiagMatrix/ZoomMap+SphMap to
WinMap (coordinate reflection).

### Positive branches

| ID | Fixture | Type | Status | Description | Trigger |
|---|---|---|---|---|---|
| sphmap-01 | sph_inverse_cancel.map | focused | `+` | Inverse(SphMap) + SphMap with matching PolarLong cancels to UnitMap | CmpMap(Inverse(SphMap[UntRd=1]), SphMap), Series=1 |
| sphmap-02 | sph_fwd_inv_unitradius_cancel.map | focused | `+` | SphMap(UnitRadius) + Inverse(SphMap) cancels to UnitMap | CmpMap(SphMap(UntRd=1), Inverse(SphMap)), Series=1 |
| sphmap-08 | sph_matrix_sandwich.map | cascade | `+` | Inv(SphMap) + DiagMatrixMap + SphMap to WinMap | Inv(SphMap) + MatrixMap(diag,equal mag) + SphMap |
| sphmap-09 | sph_zoom_sandwich.map | cascade | `+` | Sandwich with ZoomMap instead of MatrixMap | Inv(SphMap) + ZoomMap + SphMap |

### Negative branches

| ID | Fixture | Type | Status | Description | Trigger |
|---|---|---|---|---|---|
| sphmap-03 | neg_sph_polarlong_mismatch.map | focused | `+` | Inverse+Forward but PolarLong values differ | Inverse(SphMap(PolarLong=0)) + SphMap(PolarLong=pi) |
| sphmap-04 | neg_sph_no_unitradius.map | focused | `+` | Forward+Inverse but UnitRadius not set | SphMap(UntRd=0) + Inverse(SphMap) |
| sphmap-05 | -- | focused | `- (no fixture)` | Same direction (both forward or both inverse) | Two forward SphMaps in series |
| sphmap-06 | neg_sph_parallel.map | focused | `+` | Parallel mode or last in list | SphMap in parallel |
| sphmap-07 | neg_sph_non_sphmap_neighbour.map | focused | `+` | Following Mapping is not a SphMap | SphMap + ZoomMap in series |
| sphmap-10 | -- | cascade | `- (no fixture)` | Third Mapping not a non-inverted SphMap | Inv(SphMap) + MatrixMap + ZoomMap |
| sphmap-11 | neg_sph_sandwich_wrong_middle.map | cascade | `+` | Middle not ZoomMap or diagonal MatrixMap | Inv(SphMap) + ShiftMap + SphMap |
| sphmap-12 | -- | cascade | `- (no fixture)` | ZoomMap has zero factor | Inv(SphMap) + ZoomMap(0) + SphMap |
| sphmap-13 | neg_sph_sandwich_full_matrix.map | cascade | `+` | MatrixMap not diagonal or null | Inv(SphMap) + FullMatrixMap + SphMap |
| sphmap-14 | neg_sph_sandwich_unequal_diag.map | cascade | `+` | MatrixMap diagonal: unequal magnitude | Inv(SphMap) + MatrixMap(diag=[1,2,3]) + SphMap |
| sphmap-15 | -- | cascade | `- (no fixture)` | MatrixMap first diagonal is zero | Inv(SphMap) + MatrixMap(diag=[0,1,1]) + SphMap |
| sphmap-16 | -- | cascade | `- (no fixture)` | Adjusted PolarLong doesn't match third SphMap | PolarLong mismatch after sign adjustment |
| sphmap-17 | -- | cascade | `- (no fixture)` | Nominated SphMap not inverted | Forward SphMap + MatrixMap + SphMap |
| sphmap-18 | -- | cascade | `- (no fixture)` | Fewer than 3 Mappings remain | Inv(SphMap) at end of list |

---

## pcdmap.c

PcdMap's MapMerge: (1) Disco=0 to UnitMap; (2) inverse-pair cancellation;
(3) swap past ZoomMap/PermMap toward merge target or for local simplification.

### Positive branches

| ID | Fixture | Type | Status | Description | Trigger |
|---|---|---|---|---|---|
| pcdmap-01 | pcd_zero_to_unit.map | focused | `+` | PcdMap with Disco=0 replaced by UnitMap | PcdMap(Disco=0) |
| pcdmap-02 | pcd_inverse_cancel.map | focused | `+` | PcdMap + Inverse(PcdMap) cancel to UnitMap | CmpMap(PcdMap, Inverse(PcdMap)), Series=1 |
| pcdmap-03 | pcd_unit_series_merge.map | focused | `+` | PcdMap + UnitMap neighbour: UnitMap eliminated | CmpMap(PcdMap, UnitMap), Series=1 |
| pcdmap-04 | pcd_zoom_swap_cancel.map | cascade | `+` | PcdMap swaps with ZoomMap toward merge target | PcdMap + ZoomMap + Inverse(PcdMap) |
| pcdmap-05 | pcd_perm_swap_cancel.map | cascade | `+` | PcdMap swaps with axis-swapping PermMap toward merge target | PcdMap + PermMap(swap) + Inverse(PcdMap) |
| pcdmap-06 | pcd_swap_zoom_simplifies.map | cascade | `+` | Swap without target if it simplifies one Mapping | PcdMap + ZoomMap(1) where ZoomMap becomes UnitMap after swap |

### Negative branches

| ID | Fixture | Type | Status | Description | Trigger |
|---|---|---|---|---|---|
| pcdmap-07 | neg_pcd_zoom_no_target.map | cascade | `-` | Speculative swap refused: neither simplifies | PcdMap + ZoomMap(2) with no target |
| pcdmap-08 | neg_pcd_parallel.map | focused | `+` | Parallel mode: refused | PcdMaps in parallel |
| pcdmap-09 | neg_pcd_nonpcd_neighbour.map | focused | `+` | Neighbour not PcdMap/UnitMap/inverse-PcdMap | PcdMap + ShiftMap in series |
| pcdmap-10 | neg_pcd_nonswappable_between.map | cascade | `+` | Intervening Mapping not ZoomMap or PermMap: swap blocked | PcdMap + ShiftMap + Inverse(PcdMap) |
| pcdmap-11 | pcd_search_blocked.map | cascade | `-` | Non-swappable class blocks the forward swap search | PcdMap + ZoomMap + MatrixMap (Zoom swappable, Matrix blocks) |
| pcdmap-12 | -- | focused | `- (no fixture)` | CanSwap false: PermMap doesn't simply swap axes | PcdMap + PermMap(identity) |
| pcdmap-13 | neg_pcd_zoom_no_target.map | cascade | `-` | Backward (swaplo) swap search: swappable lower ZoomMap but no merge target below | CmpMap(ZoomMap, PcdMap), Series=1 |

---

## dssmap.c

DssMap's MapMerge: absorbs a neighbouring WinMap into modified pixel
parameters (CNPIX, XPIXELSZ, YPIXELSZ) within its FitsChan. Has many guards
(integer CNPIX, non-zero scale, required keywords present).

### Positive branches

| ID | Fixture | Type | Status | Description | Trigger |
|---|---|---|---|---|---|
| dssmap-07 | dssmap_winmap_absorb.map | focused | `+` | Non-inverted DssMap absorbs preceding WinMap | WinMap + DssMap(Invert=0) with valid integer CNPIX |
| dssmap-08 | dssmap_inv_winmap_absorb.map | focused | `+` | Inverted DssMap absorbs following WinMap | DssMap(Invert=1) + WinMap with valid integer CNPIX |

### Negative branches

| ID | Fixture | Type | Status | Description | Trigger |
|---|---|---|---|---|---|
| dssmap-01 | -- | focused | `- (infeasible: protected constructor)` | Parallel mode: refused | DssMap in parallel |
| dssmap-02 | -- | focused | `- (infeasible: protected constructor)` | No adjacent mapping at expected index | DssMap alone or at boundary |
| dssmap-03 | dssmap_zoom_no_merge.map | focused | `-` | Adjacent mapping is not a WinMap | DssMap + ZoomMap in series |
| dssmap-04 | -- | focused | `- (infeasible: protected constructor)` | WinMap scale/shift unusable (AST__BAD or zero scale) | WinMap with zero scale preceding DssMap |
| dssmap-05 | -- | focused | `- (infeasible: protected constructor)` | Computed CNPIX values non-integer beyond tolerance | WinMap producing fractional CNPIX |
| dssmap-06 | -- | focused | `- (infeasible: protected constructor)` | Required FITS keywords missing from DssMap FitsChan | DssMap with incomplete FitsChan |

---

## grismmap.c

GrismMap's MapMerge: (1) inverse-pair cancellation; (2) ZoomMap absorption
into wavelength parameters (forward GrismMap+Zoom or Zoom+inverse GrismMap).

### Positive branches

| ID | Fixture | Type | Status | Description | Trigger |
|---|---|---|---|---|---|
| grismmap-01 | grism_inverse_cancel.map | focused | `+` | Two GrismMaps with matching attributes and opposite invert cancel to UnitMap | CmpMap(GrismMap, Inverse(GrismMap)), Series=1 |
| grismmap-02 | grism_zoom_merge.map | focused | `+` | Forward GrismMap + ZoomMap: zoom absorbed into wavelength params | CmpMap(GrismMap(fwd), ZoomMap), Series=1 |
| grismmap-03 | grism_zoom_inv_merge.map | focused | `+` | ZoomMap + Inverse(GrismMap): zoom absorbed into wavelength params | CmpMap(ZoomMap, GrismMap(inv)), Series=1 |

### Negative branches

| ID | Fixture | Type | Status | Description | Trigger |
|---|---|---|---|---|---|
| grismmap-04 | neg_grism_parallel.map | focused | `+` | Parallel mode: refused | GrismMaps in parallel |
| grismmap-05 | -- | focused | `- (no fixture)` | No mergeable neighbour found | Single GrismMap with non-mergeable neighbours |
| grismmap-06 | neg_grism_different_attrs.map | focused | `+` | Two GrismMaps with differing attributes (NR, NRP, WaveR, etc.) | GrismMaps with different parameters |
| grismmap-07 | -- | focused | `- (no fixture)` | GrismM values equal: merge blocked | Two GrismMaps with same M |
| grismmap-08 | neg_grism_same_direction.map | focused | `+` | Same invert flag (same direction) | Two forward GrismMaps in series |
| grismmap-09 | neg_grism_inv_then_zoom.map | focused | `+` | First GrismMap is inverted: ZoomMap merge N/A | Inverse(GrismMap) + ZoomMap |
| grismmap-10 | neg_grism_fwd_then_nonzoom.map | focused | `+` | First is forward GrismMap but second not ZoomMap | GrismMap(fwd) + ShiftMap |
| grismmap-11 | neg_grism_zoom_before_fwd.map | focused | `+` | Second GrismMap not inverted: pre-ZoomMap merge N/A | ZoomMap + GrismMap(fwd) |
| grismmap-12 | neg_grism_nonzoom_before_inv.map | focused | `+` | Second is inverted GrismMap but first not ZoomMap | ShiftMap + Inverse(GrismMap) |
| grismmap-13 | -- | focused | `- (no fixture)` | ZoomMap zoom factor is zero | GrismMap(fwd) + ZoomMap(0) |

---

## xphmap.c

XphMap's MapMerge: inverse-pair cancellation only (same order, opposite
direction via astEqual).

### Positive branches

| ID | Fixture | Type | Status | Description | Trigger |
|---|---|---|---|---|---|
| xphmap-04 | -- | focused | `- (infeasible: protected constructor)` | Inverse-pair cancellation to UnitMap | XphMap + Inverse(XphMap), same order |

### Negative branches

| ID | Fixture | Type | Status | Description | Trigger |
|---|---|---|---|---|---|
| xphmap-01 | -- | focused | `- (infeasible: protected constructor)` | Parallel mode: refused | XphMap in parallel |
| xphmap-02 | -- | focused | `- (infeasible: protected constructor)` | No adjacent XphMap found | XphMap flanked by non-XphMaps |
| xphmap-03 | -- | focused | `- (infeasible: protected constructor)` | Adjacent XphMap not equal-and-opposite | XphMaps with different order or same direction |

---

## Region-as-Mapping Classes

All four Region classes below share structurally identical MapMerge logic:
(1) self-simplification via astSimplify; (2) parallel merge with adjacent
Region via class-specific MergeXxx helper.

### box.c

#### Positive branches

| ID | Fixture | Type | Status | Description | Trigger |
|---|---|---|---|---|---|
| box-01 | -- | focused | `- (unreachable: API pre-simplifies)` | Self-simplification succeeds (astSimplify returns different pointer) | Box with non-trivial base-to-current FrameSet |
| box-03 | box_parallel_merge.map | cascade | `+` | Parallel merge with lower Region via MergeBox | Box in parallel with compatible Region (lower) |
| box-04 | box_parallel_merge.map | cascade | `+` | Parallel merge with upper Region via MergeBox | Box in parallel with compatible Region (upper) |

#### Negative branches

| ID | Fixture | Type | Status | Description | Trigger |
|---|---|---|---|---|---|
| box-02 | -- | focused | `- (requires compound FrameSet)` | No self-simplification and series mode | Already-simple Box in series |
| box-05 | -- | focused | `- (requires compound FrameSet)` | Parallel but no Region neighbour or MergeBox returns NULL | Box in parallel with incompatible Region |
| box-06 | neg_box_asymmetric_2d.map | focused | `-` | Standalone 2-D Box with asymmetric axis intervals does not self-simplify | Lone Box, differing intervals per axis, no neighbour |

### interval.c

#### Positive branches

| ID | Fixture | Type | Status | Description | Trigger |
|---|---|---|---|---|---|
| interval-01 | -- | focused | `- (unreachable: API pre-simplifies)` | Self-simplification succeeds | Interval with non-trivial FrameSet |
| interval-03 | interval_parallel_merge.map | cascade | `+` | Parallel merge with lower Region | Interval + compatible Region (lower) |
| interval-04 | interval_parallel_merge.map | cascade | `+` | Parallel merge with upper Region | Interval + compatible Region (upper) |

#### Negative branches

| ID | Fixture | Type | Status | Description | Trigger |
|---|---|---|---|---|---|
| interval-02 | -- | focused | `- (requires compound FrameSet)` | No self-simplification, series mode | Simple Interval in series |
| interval-05 | -- | focused | `- (requires compound FrameSet)` | Parallel but no compatible Region | Interval + non-Region in parallel |

### nullregion.c

#### Positive branches

| ID | Fixture | Type | Status | Description | Trigger |
|---|---|---|---|---|---|
| nullregion-01 | -- | focused | `- (unreachable: API pre-simplifies)` | Self-simplification succeeds | NullRegion with non-trivial FrameSet |
| nullregion-03 | -- | cascade | `- (no structural change)` | Parallel merge with lower Region | NullRegion + compatible Region (lower) |
| nullregion-04 | -- | cascade | `- (no structural change)` | Parallel merge with upper Region | NullRegion + compatible Region (upper) |

#### Negative branches

| ID | Fixture | Type | Status | Description | Trigger |
|---|---|---|---|---|---|
| nullregion-02 | -- | focused | `- (requires compound FrameSet)` | No self-simplification, series mode | Simple NullRegion in series |
| nullregion-05 | -- | focused | `- (requires compound FrameSet)` | Parallel but no compatible Region | NullRegion + non-Region in parallel |

### pointlist.c

#### Positive branches

| ID | Fixture | Type | Status | Description | Trigger |
|---|---|---|---|---|---|
| pointlist-01 | -- | focused | `- (unreachable: API pre-simplifies)` | Self-simplification succeeds | PointList with non-trivial FrameSet |
| pointlist-03 | -- | cascade | `- (no structural change)` | Parallel merge with lower Region | PointList + compatible Region (lower) |
| pointlist-04 | -- | cascade | `- (no structural change)` | Parallel merge with upper Region | PointList + compatible Region (upper) |

#### Negative branches

| ID | Fixture | Type | Status | Description | Trigger |
|---|---|---|---|---|---|
| pointlist-02 | -- | focused | `- (requires compound FrameSet)` | No self-simplification, series mode | Simple PointList in series |
| pointlist-05 | -- | focused | `- (requires compound FrameSet)` | Parallel but no compatible Region | PointList + non-Region in parallel |

---

## Maintenance Contract

When modifying a class's `MapMerge` method:

1. If you add a new branch, add a row to the appropriate table in this
   document with a unique ID (next sequential number for that class).
2. Create a corresponding `.map` fixture in `ast_tester/simplify_fixtures/`
   and add it to `ast_tester/simplify_tests.txt`.
3. For negative (guard) branches, the `.map` file is its own reference
   (input == expected output).
4. For positive (simplification fires) branches, generate a `.simp`
   reference file showing the expected simplified output.
5. Update the Status column to `+` once the fixture passes.

---

## Cross-class Index

Alphabetical list of all inventory IDs with one-line descriptions.

| ID | Description |
|---|---|
| box-01 | Box self-simplification succeeds |
| box-02 | No self-simplification, series mode |
| box-03 | Parallel merge with lower Region via MergeBox |
| box-04 | Parallel merge with upper Region via MergeBox |
| box-05 | Parallel but no compatible Region neighbour |
| box-06 | Standalone asymmetric 2-D Box does not self-simplify |
| cmpmap-01 | CmpMap self-simplifies via astSimplify |
| cmpmap-02 | CmpMap does not self-simplify |
| cmpmap-03 | CmpMap decomposed when mode matches list mode |
| cmpmap-04 | Decomposition refused: mode mismatch |
| cmpmap-05 | Merging refused: only one mapping in list |
| cmpmap-06 | Merging refused: neighbour is not CmpMap |
| cmpmap-07 | Series CmpMaps in parallel restructured to parallel-then-series |
| cmpmap-08 | Rearranged parallel CmpMaps do not simplify |
| cmpmap-09 | Parallel CmpMaps in series paired by dimension |
| cmpmap-10 | Not both parallel or not series list |
| cmpmap-11 | Parallel-in-series pairing produces no simplification |
| cmpmap-12 | Earlier branch already succeeded |
| cmpmap-13 | CmpMap at position 0 (no preceding neighbour) |
| cmpmap-14 | Predecessor is not PermMap |
| cmpmap-15 | CmpMap is series (not parallel) |
| cmpmap-16 | PermMap axes not contiguous-block |
| cmpmap-17 | PermMap and parallel CmpMap swapped (no constants) |
| cmpmap-18 | PermMap swap with aconstants |
| cmpmap-19 | PermMap swap with bconstants |
| cmpmap-20 | Right-associative CmpMap chain reassociated |
| cmpmap-21 | UnitMap absorbed from nested CmpMap |
| cmpmap-22 | Large celestial FITS-WCS pipeline collapses |
| cmpmap-23 | GAPPT/IWC conversion residue reassociated |
| dssmap-01 | Parallel mode refused |
| dssmap-02 | No adjacent mapping at expected index |
| dssmap-03 | Adjacent mapping not a WinMap |
| dssmap-04 | WinMap scale/shift unusable |
| dssmap-05 | Computed CNPIX non-integer |
| dssmap-06 | Required FITS keywords missing |
| dssmap-07 | Non-inverted DssMap absorbs preceding WinMap |
| dssmap-08 | Inverted DssMap absorbs following WinMap |
| grismmap-01 | Inverse-pair cancellation to UnitMap |
| grismmap-02 | Forward GrismMap + ZoomMap absorbed |
| grismmap-03 | ZoomMap + Inverse(GrismMap) absorbed |
| grismmap-04 | Parallel mode refused |
| grismmap-05 | No mergeable neighbour |
| grismmap-06 | Differing GrismMap attributes |
| grismmap-07 | GrismM values equal: merge blocked |
| grismmap-08 | Same invert flag |
| grismmap-09 | First GrismMap inverted: ZoomMap merge N/A |
| grismmap-10 | Second not ZoomMap after forward GrismMap |
| grismmap-11 | Second GrismMap not inverted |
| grismmap-12 | First not ZoomMap before inverse GrismMap |
| grismmap-13 | ZoomMap factor is zero |
| interval-01 | Interval self-simplification succeeds |
| interval-02 | No self-simplification, series mode |
| interval-03 | Parallel merge with lower Region |
| interval-04 | Parallel merge with upper Region |
| interval-05 | Parallel but no compatible Region |
| intramap-01 | Forward+Inverse IntraMap with SIMPFI cancels |
| intramap-02 | Inverse+Forward IntraMap with SIMPIF cancels |
| intramap-03 | Not in series or no following Mapping |
| intramap-04 | Following Mapping is not IntraMap |
| intramap-05 | Different transformation functions |
| intramap-06 | IntraFlag strings differ |
| intramap-07 | Dimension mismatch |
| intramap-08 | Same direction |
| intramap-09 | SIMPFI flag not set |
| intramap-10 | SIMPIF flag not set |
| lutmap-01 | Linear non-inverted LutMap to WinMap |
| lutmap-02 | Linear inverted LutMap to reversed WinMap |
| lutmap-03 | LutMap cancels with upper-neighbour |
| lutmap-04 | LutMap cancels with lower-neighbour |
| lutmap-05 | LutMap not linear |
| lutmap-06 | LutMap linear but constant |
| lutmap-07 | Not in series: cancellation skipped |
| lutmap-08 | Neither neighbour is a LutMap |
| lutmap-09 | Neighbouring LutMap not inverse-equal |
| mathmap-01 | Matching fwd/inv text cancel to UnitMap |
| mathmap-02 | Parallel mode refuses |
| mathmap-03 | No following Mapping |
| mathmap-04 | Neighbour is not a MathMap |
| mathmap-05 | SimpFI/SimpIF not set |
| mathmap-06 | Dimension mismatch |
| mathmap-07 | Forward function count mismatch |
| mathmap-08 | Forward function text mismatch |
| mathmap-09 | Inverse function count mismatch |
| mathmap-10 | Inverse function text mismatch |
| matrixmap-01 | Unit-form MatrixMap to UnitMap |
| matrixmap-02 | Singular diagonal cannot become ZoomMap |
| matrixmap-03 | Equal-diagonal MatrixMap to ZoomMap |
| matrixmap-04 | Unequal diagonal cannot become ZoomMap |
| matrixmap-05 | Zero off-diagonals to Diagonal form |
| matrixmap-06 | Non-zero off-diagonal: no self-simplify |
| matrixmap-07 | Two MatrixMaps series merged |
| matrixmap-08 | MatrixMap + ZoomMap series merged |
| matrixmap-09 | MatrixMap + bidirectional PermMap merged |
| matrixmap-10 | PermMap non-bidirectional: merge blocked |
| matrixmap-11 | Diagonal MatrixMap + WinMap merged |
| matrixmap-12 | MatrixMap + UnitMap: UnitMap eliminated |
| matrixmap-13 | MatrixMap swaps past WinMap |
| matrixmap-14 | MatrixMap swaps past PermMap |
| matrixmap-15 | Swap produces simpler Mapping |
| matrixmap-16 | Swap refused: neither simplifies |
| matrixmap-17 | Parallel: non-simplifiable, returns -1 |
| normmap-01 | Encapsulated Frame simplifies |
| normmap-02 | Basic Frame to UnitMap |
| normmap-03 | Inverse lower-neighbour cancellation |
| normmap-04 | Inverse upper-neighbour cancellation |
| normmap-05 | Duplicate adjacent NormMaps eliminated |
| normmap-06 | Lower neighbour: invert flags not opposite |
| normmap-07 | Lower inverse NormMap: Frames not equal |
| normmap-08 | Upper neighbour not a NormMap |
| normmap-09 | Upper inverse NormMap: Frames differ |
| normmap-10 | Same-direction NormMap: Frames differ |
| normmap-11 | Parallel mode: no simplification |
| normmap-12 | Non-basic Frame, parallel, no simplification |
| nullregion-01 | NullRegion self-simplification succeeds |
| nullregion-02 | No self-simplification, series mode |
| nullregion-03 | Parallel merge with lower Region |
| nullregion-04 | Parallel merge with upper Region |
| nullregion-05 | Parallel but no compatible Region |
| pcdmap-01 | Disco=0 to UnitMap |
| pcdmap-02 | Inverse-pair cancellation to UnitMap |
| pcdmap-03 | UnitMap neighbour eliminated |
| pcdmap-04 | Swaps with ZoomMap toward target |
| pcdmap-05 | Swaps with PermMap toward target |
| pcdmap-06 | Swap simplifies one Mapping |
| pcdmap-07 | Speculative swap refused |
| pcdmap-13 | Backward swap search, no merge target below |
| pcdmap-08 | Parallel mode refused |
| pcdmap-09 | Neighbour not PcdMap/UnitMap |
| pcdmap-10 | Intervening not ZoomMap/PermMap |
| pcdmap-11 | Non-swappable class blocks search |
| pcdmap-12 | PermMap doesn't simply swap axes |
| permmap-01 | Adjacent PermMaps series composed |
| permmap-02 | Adjacent PermMaps parallel composed |
| permmap-03 | Composed PermMap reduces to UnitMap |
| permmap-04 | Invert flag normalized |
| permmap-05 | Array becomes null after composition |
| permmap-06 | Inperm array differs after folding |
| permmap-07 | Outperm array differs after re-computation |
| permmap-08 | No simplification: canonical lone PermMap |
| permmap-09 | Constants propagated through series merge |
| permmap-10 | AST__BAD propagated through merge |
| pointlist-01 | PointList self-simplification succeeds |
| pointlist-02 | No self-simplification, series mode |
| pointlist-03 | Parallel merge with lower Region |
| pointlist-04 | Parallel merge with upper Region |
| pointlist-05 | Parallel but no compatible Region |
| polymap-01 | Duplicate terms combined |
| polymap-02 | Linear PolyMap to ShiftMap+MatrixMap |
| polymap-03 | No forward transform defined |
| polymap-04 | nin != nout: linearization refused |
| polymap-05 | Power > 1: linearization refused |
| polymap-06 | Inverse-pair cancellation to UnitMap |
| polymap-07 | Parallel: inverse-cancel refused |
| polymap-08 | Neighbour not PolyMap |
| polymap-09 | Same invert direction |
| polymap-10 | astEqual fails (different coefficients) |
| ratemap-01 | Encapsulated Mapping simplifies |
| ratemap-02 | Lower-neighbour inverse cancellation |
| ratemap-03 | Upper-neighbour inverse cancellation |
| ratemap-04 | Not in series: skipped |
| ratemap-05 | Lower neighbour not RateMap |
| ratemap-06 | Lower neighbour same invert flag |
| ratemap-07 | Lower neighbour different indices |
| ratemap-08 | Lower neighbour non-equal Mapping |
| ratemap-09 | Upper neighbour not RateMap |
| ratemap-10 | Upper neighbour same invert flag |
| ratemap-11 | Upper neighbour different indices |
| ratemap-12 | Upper neighbour non-equal Mapping |
| scenario-01 | Brad Warren FITS-WCS integration |
| scenario-02 | LSST camera-geometry integration |
| scenario-03 | Rigby spectroscopic pipeline integration |
| selectormap-01 | Internal regions simplify |
| selectormap-02 | No region simplification, no adjacent SelectorMap |
| selectormap-03 | No adjacent SelectorMap in series |
| selectormap-04 | Adjacent SelectorMap not equal-and-opposite |
| selectormap-05 | Inverse-pair cancellation to UnitMap |
| selectormap-06 | Parallel mode: skipped |
| shiftmap-01 | Inverted ShiftMap normalized |
| slamap-01 | Full cancellation to UnitMap |
| slamap-02 | Partial cancellation |
| slamap-03 | Adjacent SlaMaps merged without reduction |
| slamap-04 | Invert-flag clearing |
| slamap-05 | Parallel-mode guard |
| slamap-06 | Lone forward SlaMap: no simplification |
| slamap-07 | 0-arg: galactic pair |
| slamap-08 | 0-arg: supergalactic pair |
| slamap-09 | 0-arg: dynamical J2000 pair |
| slamap-10 | 1-arg: E-terms pair |
| slamap-11 | 1-arg: FK4/FK5 pair |
| slamap-12 | 1-arg: ICRS/FK5 pair |
| slamap-13 | 1-arg: ecliptic pair |
| slamap-14 | 1-arg: helio-ecliptic pair |
| slamap-15 | 1-arg: HA pair |
| slamap-16 | 2-arg: geocentric pair (crossed) |
| slamap-17 | 2-arg: AzEl pair (same order) |
| slamap-18 | 4-arg: HPC pair |
| slamap-19 | 4-arg: HPR pair |
| slamap-20 | 1-arg: mismatched argument |
| slamap-21 | 2-arg: mismatched arguments |
| slamap-22 | 4-arg: mismatched arguments |
| slamap-23 | Redundant precession eliminated |
| slamap-24 | Adjacent precession merged |
| slamap-25 | Precession: non-common equinox |
| slamap-26 | Run of three SlaMaps cancels to UnitMap |
| slamap-27 | Run of three SlaMaps partially cancels |
| specmap-01 | Full cancellation to UnitMap |
| specmap-02 | Partial cancellation |
| specmap-03 | Adjacent SpecMaps merged without reduction |
| specmap-04 | Invert-flag clearing |
| specmap-05 | Parallel-mode guard |
| specmap-06 | Lone forward SpecMap: no simplification |
| specmap-07 | Nin-mismatch guard |
| specmap-08 | 0-arg: unit conversion pair |
| specmap-09 | 1-arg: frequency/velocity pair |
| specmap-10 | 2-arg: local-standard pair |
| specmap-11 | 3-arg: geocentric/barycentric pair |
| specmap-12 | 6-arg: topocentric pair |
| specmap-13 | 1-arg: mismatched argument |
| specmap-14 | 2-arg: mismatched arguments |
| specmap-15 | 3-arg: mismatched arguments |
| specmap-16 | 6-arg: mismatched arguments |
| specmap-17 | Run of four SpecMaps cancels to UnitMap |
| specmap-18 | Run of three SpecMaps partially cancels |
| sphmap-01 | Inverse+Forward with matching PolarLong cancels |
| sphmap-02 | Forward(UnitRadius)+Inverse cancels |
| sphmap-03 | PolarLong values differ |
| sphmap-04 | UnitRadius not set |
| sphmap-05 | Same direction |
| sphmap-06 | Parallel mode or last in list |
| sphmap-07 | Following Mapping not SphMap |
| sphmap-08 | Matrix sandwich to WinMap |
| sphmap-09 | ZoomMap sandwich |
| sphmap-10 | Third Mapping not non-inverted SphMap |
| sphmap-11 | Middle not ZoomMap/diagonal MatrixMap |
| sphmap-12 | ZoomMap has zero factor |
| sphmap-13 | MatrixMap not diagonal |
| sphmap-14 | Diagonal unequal magnitude |
| sphmap-15 | First diagonal is zero |
| sphmap-16 | Adjusted PolarLong mismatch |
| sphmap-17 | Nominated SphMap not inverted |
| sphmap-18 | Fewer than 3 Mappings remain |
| splinemap-01 | Lower-neighbour inverse cancellation |
| splinemap-02 | Upper-neighbour inverse cancellation |
| splinemap-03 | Not in series |
| splinemap-04 | No neighbour (boundary) |
| splinemap-05 | Neighbour not SplineMap |
| splinemap-06 | Same invert flag |
| splinemap-07 | astEqual fails (different coefficients) |
| switchmap-01 | Inverse-pair cancellation to UnitMap |
| switchmap-02 | Adjacent SwitchMap not equal-and-opposite |
| switchmap-03 | Inverted SwitchMap normalized |
| switchmap-04 | Internal selectors/routes simplify |
| switchmap-05 | Series: no simplification possible |
| switchmap-06 | Parallel: no simplification possible |
| timemap-01 | Full cancellation to UnitMap |
| timemap-02 | Partial cancellation |
| timemap-03 | Adjacent TimeMaps merged without reduction |
| timemap-04 | Invert-flag clearing |
| timemap-05 | Parallel-mode guard |
| timemap-06 | Lone forward TimeMap: no simplification |
| timemap-07 | No-op step eliminated |
| timemap-08 | 1-arg pair cancellation |
| timemap-09 | 2-arg swapped pair cancellation |
| timemap-10 | 2-arg same-order pair cancellation |
| timemap-11 | 3-arg pair cancellation |
| timemap-12 | 5-arg pair cancellation |
| timemap-13 | 1-arg: mismatched argument |
| timemap-14 | 2-arg swapped: mismatched arguments |
| timemap-15 | 3-arg: mismatched arguments |
| timemap-16 | 5-arg: mismatched arguments |
| timemap-17 | Run of three TimeMaps cancels to UnitMap |
| tranmap-01 | Inverted TranMap normalized |
| tranmap-02 | Components individually simplified |
| tranmap-03 | Equal components to single Mapping |
| tranmap-04 | Adjacent TranMaps series merged |
| tranmap-05 | Parallel mode: merge skipped |
| tranmap-06 | Components lack bidirectional transforms |
| tranmap-07 | Components not equal |
| tranmap-08 | Higher neighbour not TranMap |
| tranmap-09 | Neither leg simplified |
| unitmap-01 | Invert flag cleared |
| unitmap-02 | Invert=0: no change |
| unitmap-03 | UnitMap removed from series |
| unitmap-04 | Adjacent UnitMaps parallel merged |
| unitmap-05 | Parallel Invert cleared |
| unitmap-06 | Parallel Invert=0: no simplification |
| unitnormmap-01 | ShiftMap + forward UnitNormMap merged |
| unitnormmap-02 | WinMap(unit) + forward UnitNormMap merged |
| unitnormmap-03 | WinMap(non-unit) + UnitNormMap refused |
| unitnormmap-04 | Inverse UnitNormMap + ShiftMap merged |
| unitnormmap-05 | Inverse UnitNormMap + WinMap(unit) merged |
| unitnormmap-06 | Inverse UnitNormMap + WinMap(non-unit) refused |
| unitnormmap-07 | Same-centre inverse pair to UnitMap |
| unitnormmap-08 | Same-centre inverse pair (reversed) to UnitMap |
| unitnormmap-09 | Different-centre inverse pair to ShiftMap |
| unitnormmap-10 | Different-centre reversed: no merge |
| unitnormmap-11 | ShiftMap + UnitNormMap(inv) refused |
| unitnormmap-12 | WinMap + UnitNormMap(inv) refused |
| unitnormmap-13 | UnitNormMap(fwd) + ShiftMap refused |
| unitnormmap-14 | UnitNormMap(fwd) + WinMap refused |
| unitnormmap-15 | Two same-direction UnitNormMaps refused |
| unitnormmap-16 | Non-mergeable neighbour class |
| unitnormmap-17 | Parallel mode: never simplifies |
| unitnormmap-18 | No valid merge produced |
| wcsmap-01 | AST__WCSBAD to UnitMap |
| wcsmap-02 | Inverse-pair cancellation to UnitMap |
| wcsmap-03 | Swaps past PermMap to merge target |
| wcsmap-04 | Swap for local simplification |
| wcsmap-05 | Speculative swap refused |
| wcsmap-06 | Parallel mode refused |
| wcsmap-07 | Different projection type |
| wcsmap-08 | Same invert direction |
| wcsmap-09 | Lon/lat axis indices differ |
| wcsmap-10 | Projection parameters differ |
| wcsmap-11 | Intervening Mapping not PermMap |
| wcsmap-12 | PermMap non-bidirectional |
| wcsmap-13 | PermMap disconnects lon/lat |
| winmap-01 | All shifts zero: WinMap to MatrixMap |
| winmap-02 | All scales unity: WinMap to ShiftMap |
| winmap-03 | Not all shifts zero: MatrixMap refused |
| winmap-04 | Not all scales unity: ShiftMap refused |
| winmap-05 | WinMap + WinMap series merged |
| winmap-06 | WinMap + ZoomMap series merged |
| winmap-07 | ZoomMap + WinMap series merged |
| winmap-08 | WinMap + ShiftMap series merged |
| winmap-09 | ShiftMap + WinMap series merged |
| winmap-10 | WinMap + diagonal MatrixMap series merged |
| winmap-11 | Diagonal MatrixMap + WinMap series merged |
| winmap-12 | WinMap + UnitMap: UnitMap removed |
| winmap-13 | Neither neighbour is directly mergeable |
| winmap-14 | WinMap merges with parallel CmpMap (lower) |
| winmap-15 | WinMap merges with parallel CmpMap (upper) |
| winmap-16 | CmpMap neighbour is series: no merge |
| winmap-17 | Parallel CmpMap split doesn't simplify |
| winmap-18 | WinMap swaps past MatrixMap |
| winmap-19 | WinMap swaps past PermMap |
| winmap-20 | WinMap swaps past WcsMap |
| winmap-21 | No higher neighbour exists |
| winmap-22 | No lower neighbour exists |
| winmap-23 | Forward scan hits non-swappable class |
| winmap-24 | Backward scan hits non-swappable class |
| winmap-25 | Both swap directions find no target |
| winmap-26 | Swap accepted: swapped Mapping simplifies |
| winmap-27 | Swap accepted: outer neighbours merge |
| winmap-28 | Swap refused: neither simplifies |
| winmap-29 | WinMap + WinMap parallel merged |
| winmap-30 | WinMap + ZoomMap parallel merged |
| winmap-31 | ZoomMap + WinMap parallel merged |
| winmap-32 | WinMap + ShiftMap parallel merged |
| winmap-33 | ShiftMap + WinMap parallel merged |
| winmap-34 | WinMap + diagonal MatrixMap parallel merged |
| winmap-35 | Diagonal MatrixMap + WinMap parallel merged |
| winmap-36 | WinMap + UnitMap parallel merged |
| winmap-37 | UnitMap + WinMap parallel merged |
| winmap-38 | Neither parallel neighbour mergeable |
| xphmap-01 | Parallel mode refused |
| xphmap-02 | No adjacent XphMap |
| xphmap-03 | Not equal-and-opposite |
| xphmap-04 | Inverse-pair cancellation to UnitMap |
| zoommap-01 | Adjacent ZoomMaps multiplied |
| zoommap-02 | Product=1 collapse to UnitMap |
| zoommap-03 | Inverted ZoomMap normalized |
| zoommap-04 | Lone forward ZoomMap: nothing to accumulate |
| zoommap-05 | All parallel factors=1 to UnitMap |
| zoommap-06 | Same parallel factor to single ZoomMap |
| zoommap-07 | Different parallel factors to MatrixMap |
| zoommap-08 | Lone ZoomMap in parallel: no simplification |
| zoommap-09 | Absorbed into previous MatrixMap |
| zoommap-10 | Absorbed into previous WinMap |
| zoommap-11 | Absorbed into previous ShiftMap |
| zoommap-12 | Absorbed into next MatrixMap |
| zoommap-13 | Absorbed into next WinMap |
| zoommap-14 | Absorbed into next ShiftMap |
| zoommap-15 | Cannot absorb: no compatible neighbour |
| zoommap-16 | Backward neighbour search stops at non-zoom lower neighbour |

## Captured differential-coverage fixtures

In addition to the hand-authored rule fixtures catalogued above, the
`cap_*` fixtures are real Mappings captured from the wider AST test suite as
it drives `astSimplify`, then replayed as self-contained fixtures to close
coverage that previously depended on running those other tests. Because each
capture exercises several branches at once, they are catalogued separately in
`ast_tester/coverage/captured_fixtures.md` (with the capture method in
`ast_tester/coverage/README.md`) rather than as per-branch rows here.
