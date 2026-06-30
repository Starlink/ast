# `astSimplify` Branch Coverage Inventory

**Campaign**: `astSimplify` branch coverage  
**Branch**: `u/timj/simplify-testing`  
**Phase**: 2 — Inventory expansion (full audit)  
**Design document**: `2026-05-05-simplify-branch-coverage-design.md` (same directory)

This inventory maps every `MapMerge` branch in scope to either an existing
fixture or a gap marker (`—`). It is a working artifact for the coverage
campaign; the permanent reference document will be
`ast_tester/simplify_pathways.md` (Phase 5).

## Row schema

| Column | Purpose |
| --- | --- |
| `ID` | Stable identifier `<class>-NN` |
| `Fixture` | `.map` filename, or `—` if uncovered |
| `Type` | `focused` \| `cascade` \| `scenario` |
| `Polarity` | `positive` \| `negative` |
| `Lines` | Source line range in MapMerge |
| `Description` | One sentence |
| `Trigger` | Minimum input shape |

---

## Scenario fixtures

These are pre-existing integration regressions that exercise multiple branches
across many classes. They are kept as integration safety nets, not targeted
branch coverage.

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
|---|---|---|---|---|---|---|
| scenario-01 | brad.map | scenario | positive | N/A | Multi-class integration regression from Brad Warren's FITS-WCS pipeline | Complex FrameSet with multiple WCS conversions |
| scenario-02 | lsst1.map | scenario | positive | N/A | LSST camera-geometry pipeline regression testing CmpMap decomposition | Nested CmpMap chain from LSST focal-plane transform |
| scenario-03 | rigby.map | scenario | positive | N/A | Deep nested Mapping chain from Jane Rigby's spectroscopic pipeline | Large multi-level CmpMap with WcsMaps, MatrixMaps, ShiftMaps |

---

## cmpmap.c

CmpMap's MapMerge works in three stages: (1) self-simplification via
astSimplify or decomposition into components if the combination mode matches
the list mode; (2) merging with a neighbouring CmpMap — series CmpMaps in a
parallel list restructure into parallel-then-series, parallel CmpMaps in a
series list pair components by dimension; (3) swapping a PermMap past a
parallel CmpMap by restructuring permutation arrays.

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
|---|---|---|---|---|---|---|
| cmpmap-01 | cmpmap_self_simplify.map | focused | positive | cmpmap.c:1456-1467 | CmpMap simplifies on its own (astSimplify returns a different/simpler mapping) | CmpMap whose internal components simplify when composed |
| cmpmap-02 | — | focused | negative | cmpmap.c:1456-1471 | CmpMap does not simplify on its own (astSimplify returns same pointer unchanged) | CmpMap with irreducible components |
| cmpmap-03 | cmpmap_nested_parallel_flatten.map | focused | positive | cmpmap.c:1477-1524 | CmpMap decomposed into components when combination mode matches list mode | Series CmpMap in a series list, or parallel CmpMap in a parallel list |
| cmpmap-04 | neg_cmpmap_mode_mismatch.map | focused | negative | cmpmap.c:1477-1478 | Guard rejects decomposition: CmpMap mode does not match list mode | Series CmpMap in a parallel list |
| cmpmap-05 | cmpmap_solo_after_unit.map | focused | negative | cmpmap.c:1534-1535 | Guard rejects merging: only one mapping in list (nmap <= 1) | Series CmpMap(UnitMap, parallel CmpMap); UnitMap elided leaves the parallel CmpMap alone |
| cmpmap-06 | neg_cmpmap_neighbour_nonexcmpmap.map | focused | negative | cmpmap.c:1551 | Guard rejects merging: neighbour is not a CmpMap | CmpMap adjacent to a non-CmpMap in the list |
| cmpmap-07 | cmpmap_parallel_series_components.map | focused | positive | cmpmap.c:1582-1622 | Two series CmpMaps in parallel list restructured into parallel-then-series and at least one simplifies | Two series CmpMaps combined in parallel with simplifiable pairings |
| cmpmap-08 | — | focused | negative | cmpmap.c:1610-1611 | Guard: re-arranged parallel CmpMaps do not simplify | Two series CmpMaps in parallel whose rearranged pairings remain irreducible |
| cmpmap-09 | cmpmap_parallel_in_series_merge.map | focused | positive | cmpmap.c:1629-1780 | Two parallel CmpMaps in series list paired by dimension, at least one pair simplifies | Two parallel CmpMaps in series with simplifiable component pairs |
| cmpmap-10 | — | cascade | negative | cmpmap.c:1629-1653 | Guard: two CmpMaps are not both parallel, or list is not series | Two adjacent CmpMaps where at least one is series |
| cmpmap-11 | — | cascade | negative | cmpmap.c:1726-1728 | Parallel-in-series pairing produces no simplification | Two parallel CmpMaps in series with all irreducible sub-mappings |
| cmpmap-12 | — | focused | negative | cmpmap.c:1825 | Guard: earlier branch already succeeded (result != -1) | Any case where a prior branch already simplified |
| cmpmap-13 | — | focused | negative | cmpmap.c:1825 | Guard: CmpMap at position 0 (no preceding neighbour) | CmpMap first in mapping list |
| cmpmap-14 | — | focused | negative | cmpmap.c:1838 | Guard: predecessor is not a PermMap | Parallel CmpMap preceded by non-PermMap |
| cmpmap-15 | — | focused | negative | cmpmap.c:1838 | Guard: CmpMap is series (not parallel) | Series CmpMap preceded by PermMap |
| cmpmap-16 | — | focused | negative | cmpmap.c:1882-1911 | Guard: PermMap axes don't satisfy contiguous-block condition (canswap=0) | PermMap with non-contiguous axis routing before parallel CmpMap |
| cmpmap-17 | cmpmap_perm_parallel_swap.map | focused | positive | cmpmap.c:1918-2066 | PermMap and parallel CmpMap swapped (no constants), producing reordered CmpMap + new PermMap | PermMap swapping two contiguous blocks feeding a parallel CmpMap |
| cmpmap-18 | cmpmap_perm_swap_aconstants.map | focused | positive | cmpmap.c:1955-1957 | PermMap swap with aconstants: first component gets all-constant outputs | PermMap with first block all constants before parallel CmpMap |
| cmpmap-19 | cmpmap_perm_swap_bconstants.map | focused | positive | cmpmap.c:1958-1960 | PermMap swap with bconstants: second component gets all-constant outputs | PermMap with second block all constants before parallel CmpMap |
| cmpmap-20 | reconstruct_right_assoc.map | cascade | positive | cmpmap.c:1205-2125 | Right-nested CmpMap chain reassociated so neighbouring components can merge | Right-associative CmpMap chain of WcsMap/ZoomMap/UnitMap/WcsMap |
| cmpmap-21 | series_parallel_absorb.map | cascade | positive | cmpmap.c:1205-2125 | UnitMap absorbed from nested series/parallel CmpMap, collapsing to WcsMap+ShiftMap | Nested CmpMap containing WcsMap, UnitMap and ShiftMap |
| cmpmap-22 | spec_cel_invert.map | scenario | positive | cmpmap.c:1205-2125 | Large celestial FITS-WCS pipeline collapses via repeated merge to CmpMap(WinMap, WcsMap) | Multi-class FITS-WCS conversion pipeline (>3 components) |
| cmpmap-23 | wcsconv_gappt_iwc_residue.map | scenario | positive | cmpmap.c:1205-2125 | GAPPT/IWC conversion residue: large pipeline reassociated and partly merged | Multi-class FITS-WCS conversion pipeline (>3 components) |

## unitmap.c

UnitMap's MapMerge: in series, removes UnitMap from the list. In parallel,
merges adjacent UnitMaps into one with combined dimensionality. Also handles
Invert flag normalization.

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
|---|---|---|---|---|---|---|
| unitmap-01 | unit_invert_clear.map | focused | positive | unitmap.c:508-514 | Single UnitMap with Invert flag set: flag is cleared | Lone UnitMap with invert_list=1 |
| unitmap-02 | neg_unit_lone.map | focused | negative | unitmap.c:508-514 | Single UnitMap with Invert=0: no change needed | Lone UnitMap already canonical |
| unitmap-03 | unit_series_elision.map | focused | positive | unitmap.c:521-537 | UnitMap removed from series composition | CmpMap(ShiftMap(2), UnitMap(2)), Series=1 |
| unitmap-04 | unit_parallel_merge.map | focused | positive | unitmap.c:589-618 | Adjacent UnitMaps in parallel merged into one wider UnitMap | CmpMap(UnitMap(1), UnitMap(2)), Series=0 |
| unitmap-05 | unit_parallel_invert_clear.map | focused | positive | unitmap.c:581-585 | Parallel UnitMap with no adjacent UnitMaps but Invert set: flag cleared | Single UnitMap in parallel with invert=1, flanked by non-UnitMaps |
| unitmap-06 | — | focused | negative | unitmap.c:577-581 | Parallel UnitMap with no adjacent UnitMaps and Invert=0: no simplification | Single UnitMap in parallel already canonical |

## zoommap.c

ZoomMap's MapMerge: in series, accumulates zoom factors of adjacent
ZoomMaps/UnitMaps (product→ZoomMap, product=1→UnitMap). If no same-class
merge, absorbs into neighbouring MatrixMap, WinMap, or ShiftMap. In parallel,
collects adjacent ZoomMaps into one ZoomMap (same factor) or diagonal
MatrixMap (different factors).

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
|---|---|---|---|---|---|---|
| zoommap-01 | zoom_series_merge.map | focused | positive | zoommap.c:806-817 | Adjacent series ZoomMaps multiplied into single ZoomMap | CmpMap(ZoomMap(2), ZoomMap(3)), Series=1 |
| zoommap-02 | zoom_series_cancel.map | focused | positive | zoommap.c:799-815 | Adjacent ZoomMaps with product=1 collapse to UnitMap | CmpMap(ZoomMap(2), ZoomMap(2,Invert=1)), Series=1 |
| zoommap-03 | zoom_invert_normalize.map | focused | positive | zoommap.c:806-817 | Single inverted ZoomMap normalized to forward-only (1/zoom) | ZoomMap with invert_list=1, no adjacent ZoomMaps |
| zoommap-04 | neg_zoom_lone.map | focused | negative | zoommap.c:806-809 | Single forward ZoomMap with no adjacent ZoomMaps: nothing to accumulate | Lone forward ZoomMap in series (falls through to absorb) |
| zoommap-05 | zoom_parallel_all_unit.map | focused | positive | zoommap.c:916-930 | All parallel ZoomMaps/UnitMaps have factor 1 → UnitMap | CmpMap(UnitMap(2), ZoomMap(1,Nin=3)), Series=0 |
| zoommap-06 | zoom_parallel_same_factor.map | focused | positive | zoommap.c:912-932 | All parallel ZoomMaps have same non-unity factor → single ZoomMap | CmpMap(ZoomMap(2,Nin=1), ZoomMap(2,Nin=2)), Series=0 |
| zoommap-07 | zoom_parallel_to_matrix.map | focused | positive | zoommap.c:937-938 | Parallel ZoomMaps with different factors → diagonal MatrixMap | CmpMap(ZoomMap(2,Nin=1), ZoomMap(3,Nin=1)), Series=0 |
| zoommap-08 | neg_zoom_parallel_lone.map | focused | negative | zoommap.c:923-926 | Single ZoomMap in parallel with no adjacent ZoomMaps: no simplification | Lone ZoomMap in parallel |
| zoommap-09 | zoom_absorb_prev_matrix.map | focused | positive | zoommap.c:1001-1011 | ZoomMap absorbed into previous MatrixMap (elements scaled) | CmpMap(MatrixMap, ZoomMap), Series=1 |
| zoommap-10 | zoom_absorb_prev_win.map | focused | positive | zoommap.c:1016-1028 | ZoomMap absorbed into previous WinMap (shifts and scales multiplied) | CmpMap(WinMap, ZoomMap), Series=1 |
| zoommap-11 | zoom_absorb_prev_shift.map | focused | positive | zoommap.c:1033-1057 | ZoomMap absorbed into previous ShiftMap → WinMap | CmpMap(ShiftMap, ZoomMap), Series=1 |
| zoommap-12 | zoom_absorb_next_matrix.map | focused | positive | zoommap.c:1068-1078 | ZoomMap absorbed into next MatrixMap (elements scaled) | CmpMap(ZoomMap, MatrixMap), Series=1 |
| zoommap-13 | zoom_absorb_next_win.map | focused | positive | zoommap.c:1083-1094 | ZoomMap absorbed into next WinMap (scales only) | CmpMap(ZoomMap, WinMap), Series=1 |
| zoommap-14 | zoom_absorb_next_shift.map | focused | positive | zoommap.c:1098-1121 | ZoomMap absorbed into next ShiftMap → WinMap | CmpMap(ZoomMap, ShiftMap), Series=1 |
| zoommap-15 | neg_zoom_no_absorb.map | focused | negative | zoommap.c:989-1155 | ZoomMap cannot be absorbed: neither neighbour is MatrixMap/WinMap/ShiftMap | CmpMap(SphMap, ZoomMap, MathMap), Series=1 |
| zoommap-16 | neg_zoom_lower_nonzoom.map | focused | negative | zoommap.c:747-766 | Backward neighbour search: lower neighbour is not a ZoomMap/UnitMap, so it stops with no merge | ZoomMap nominated at index>0 with a non-zoom lower neighbour, e.g. CmpMap(PermMap, ZoomMap) Series=1 |

## shiftmap.c

ShiftMap's MapMerge delegates to WinMap. It converts itself to WinMap with
unit scale, calls WinMap's MapMerge. If that only converts back to ShiftMap,
checks for Invert flag normalization.

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
|---|---|---|---|---|---|---|
| shiftmap-01 | shift_invert_normalize.map | focused | positive | shiftmap.c:721-738 | Normalizes an inverted ShiftMap by negating shifts and clearing Invert | ShiftMap with Invert=1 |

## winmap.c

WinMap's MapMerge: (1) self-simplification — all-zero-shift→MatrixMap,
all-unit-scale→ShiftMap; (2) series direct merge with WinMap/ZoomMap/ShiftMap/
MatrixMap(diag)/UnitMap, or merge with parallel CmpMap neighbour, or swap past
PermMap/MatrixMap/WcsMap to reach target; (3) parallel merge with same classes.

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
|---|---|---|---|---|---|---|
| winmap-01 | win_to_matrix.map | focused | positive | winmap.c:1086-1123 | WinMap with all shift=0 replaced by diagonal MatrixMap | Standalone WinMap with Scl only |
| winmap-02 | win_to_shift.map | focused | positive | winmap.c:1819-1857 | WinMap with all scale=1 replaced by ShiftMap | Standalone WinMap with Sft only |
| winmap-03 | neg_win_mixed_scale_shift.map | focused | negative | winmap.c:1088-1092 | Not all shifts are zero: MatrixMap replacement refused | WinMap with mixed shifts |
| winmap-04 | neg_win_mixed_scale_shift.map | focused | negative | winmap.c:1821-1825 | Not all scales are 1: ShiftMap replacement refused (and nothing else simplified) | WinMap with mixed scales, alone |
| winmap-05 | win_win_series_merge.map | focused | positive | winmap.c:1190-1194 | WinMap + WinMap in series merged | CmpMap(WinMap, WinMap), Series=1 |
| winmap-06 | win_zoom_series_merge.map | focused | positive | winmap.c:1196-1206 | WinMap + ZoomMap in series merged (WinMap first) | CmpMap(WinMap, ZoomMap), Series=1 |
| winmap-07 | win_zoom_series_merge_rev.map | focused | positive | winmap.c:1196-1206 | ZoomMap + WinMap in series merged (ZoomMap first) | CmpMap(ZoomMap, WinMap), Series=1 |
| winmap-08 | win_shift_series_merge.map | focused | positive | winmap.c:1208-1217 | WinMap + ShiftMap in series merged (WinMap first) | CmpMap(WinMap, ShiftMap), Series=1 |
| winmap-09 | win_shift_series_merge_rev.map | focused | positive | winmap.c:1208-1217 | ShiftMap + WinMap in series merged (ShiftMap first) | CmpMap(ShiftMap, WinMap), Series=1 |
| winmap-10 | win_matrix_series_merge.map | focused | positive | winmap.c:1220-1230 | WinMap + diagonal MatrixMap series merged (WinMap first) | CmpMap(WinMap, MatrixMap[Diagonal]), Series=1 |
| winmap-11 | win_matrix_series_merge_rev.map | focused | positive | winmap.c:1220-1230 | Diagonal MatrixMap + WinMap series merged (MatrixMap first) | CmpMap(MatrixMap[Diagonal], WinMap), Series=1 |
| winmap-12 | win_unit_series_merge.map | focused | positive | winmap.c:1232-1234 | WinMap + UnitMap in series: UnitMap removed | CmpMap(WinMap, UnitMap), Series=1 |
| winmap-13 | neg_win_nonmergeable_series.map | focused | negative | winmap.c:1164-1184 | Neither neighbour is a directly-mergeable class | CmpMap(WinMap, FullMatrixMap), Series=1 |
| winmap-14 | win_cmpmap_parallel_merge.map | cascade | positive | winmap.c:1268-1416 | WinMap merges with neighbouring parallel CmpMap (lower neighbour) | CmpMap(parallel CmpMap, WinMap), Series=1 |
| winmap-15 | win_upper_cmpmap_parallel_merge.map | cascade | positive | winmap.c:1268-1416 | WinMap merges with neighbouring parallel CmpMap (upper neighbour) | CmpMap(WinMap, parallel CmpMap), Series=1 |
| winmap-16 | — | cascade | negative | winmap.c:1303 | CmpMap neighbour is series (not parallel): no merge | CmpMap(series CmpMap, WinMap), Series=1 |
| winmap-17 | — | cascade | negative | winmap.c:1369-1373 | Parallel CmpMap split doesn't simplify (simpler==0): refused | Parallel CmpMap with non-simplifiable components next to WinMap |
| winmap-18 | win_swap_past_matrix.map | cascade | positive | winmap.c:1536-1537 | WinMap swaps past MatrixMap to reach merge target | CmpMap(WinMap, MatrixMap, WinMap), Series=1 |
| winmap-19 | win_perm_swap_merge.map | cascade | positive | winmap.c:1539-1540 | WinMap swaps past PermMap to reach merge target | CmpMap(WinMap, PermMap, WinMap), Series=1 |
| winmap-20 | win_swap_past_wcsmap.map | cascade | positive | winmap.c:1542-1543 | WinMap swaps past WcsMap to reach merge target | CmpMap(WinMap, WcsMap, WinMap), Series=1 |
| winmap-21 | — | cascade | negative | winmap.c:1432-1439 | No higher neighbour exists (WinMap last in list) | WinMap at end of series list |
| winmap-22 | — | cascade | negative | winmap.c:1480-1488 | No lower neighbour exists (WinMap first in list) | WinMap at start of series list |
| winmap-23 | — | cascade | negative | winmap.c:1468-1472 | Forward scan hits non-swappable class before target | CmpMap(WinMap, SpecMap, WinMap), Series=1 |
| winmap-24 | — | cascade | negative | winmap.c:1505-1509 | Backward scan hits non-swappable class before target | Same from other direction |
| winmap-25 | — | cascade | negative | winmap.c:1517-1528 | Both swap directions find no reachable target | WinMap adjacent to non-swappable non-mergeable class |
| winmap-26 | win_swap_simplifies.map | cascade | positive | winmap.c:1575-1682 | Swap accepted because swapped Mapping simplifies | PermMap that strips axes swapped past WinMap |
| winmap-27 | win_swap_outer_merge.map | cascade | positive | winmap.c:1619-1651 | Swap accepted because outer neighbours can merge after swap | Three-map series where outer pair merges once WinMap moves |
| winmap-28 | neg_win_swap_no_simplify.map | cascade | negative | winmap.c:1660-1669 | Swap refused: neither swapped Mapping simplifies and no outer merge | WinMap + MatrixMap where swap produces equivalent complexity |
| winmap-29 | win_win_parallel_merge.map | focused | positive | winmap.c:1729-1732 | WinMap + WinMap in parallel merged | CmpMap(WinMap, WinMap), Series=0 |
| winmap-30 | win_zoom_parallel_merge.map | focused | positive | winmap.c:1735-1745 | WinMap + ZoomMap in parallel merged (WinMap first) | CmpMap(WinMap, ZoomMap), Series=0 |
| winmap-31 | win_zoom_parallel_merge_rev.map | focused | positive | winmap.c:1735-1745 | ZoomMap + WinMap in parallel merged (ZoomMap first) | CmpMap(ZoomMap, WinMap), Series=0 |
| winmap-32 | win_parallel_merge.map | focused | positive | winmap.c:1747-1757 | WinMap + ShiftMap in parallel merged (WinMap first) | CmpMap(WinMap, ShiftMap), Series=0 |
| winmap-33 | win_shift_parallel_merge_rev.map | focused | positive | winmap.c:1747-1757 | ShiftMap + WinMap in parallel merged (ShiftMap first) | CmpMap(ShiftMap, WinMap), Series=0 |
| winmap-34 | win_diagmatrix_parallel_merge.map | focused | positive | winmap.c:1760-1770 | WinMap + diagonal MatrixMap in parallel merged (WinMap first) | CmpMap(WinMap, DiagMatrixMap), Series=0 |
| winmap-35 | win_diagmatrix_parallel_merge_rev.map | focused | positive | winmap.c:1760-1770 | Diagonal MatrixMap + WinMap in parallel merged (MatrixMap first) | CmpMap(DiagMatrixMap, WinMap), Series=0 |
| winmap-36 | win_unit_parallel_merge.map | focused | positive | winmap.c:1772-1783 | WinMap + UnitMap in parallel merged (WinMap first) | CmpMap(WinMap, UnitMap), Series=0 |
| winmap-37 | win_unit_parallel_merge_rev.map | focused | positive | winmap.c:1772-1783 | UnitMap + WinMap in parallel merged (UnitMap first) | CmpMap(UnitMap, WinMap), Series=0 |
| winmap-38 | neg_win_nonmergeable_parallel.map | focused | negative | winmap.c:1703-1723 | Neither parallel neighbour is a mergeable class | CmpMap(WinMap, FullMatrixMap), Series=0 |

## matrixmap.c

MatrixMap's MapMerge: (1) self-simplification — unit→UnitMap,
equal-diagonal→ZoomMap, zero-off-diagonal→Diagonal; (2) series merge with
MatrixMap/ZoomMap/PermMap/WinMap(diag)/UnitMap; (3) cascade swap past
WinMap/PermMap toward merge target or for local simplification.

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
|---|---|---|---|---|---|---|
| matrixmap-01 | matrix_unit_to_unit.map | focused | positive | matrixmap.c:1854-1855 | Unit-form square MatrixMap replaced by UnitMap | MatrixMap(Form="Unit", Nin=Nout) |
| matrixmap-02 | neg_matrix_singular_diag.map | focused | negative | matrixmap.c:1860-1862 | Diagonal MatrixMap with NULL i_matrix (singular) cannot become ZoomMap | MatrixMap(diag, singular) |
| matrixmap-03 | matrix_diagonal_to_zoom.map | focused | positive | matrixmap.c:1873-1881 | Diagonal MatrixMap with equal elements replaced by ZoomMap | MatrixMap(diag, [4,4]) |
| matrixmap-04 | neg_matrix_diag_unequal.map | focused | negative | matrixmap.c:1865-1868 | Diagonal MatrixMap with unequal elements cannot become ZoomMap | MatrixMap(diag, [2,3]) with non-mergeable neighbours |
| matrixmap-05 | matrix_full_to_diagonal.map | focused | positive | matrixmap.c:1887-1904 | Full MatrixMap with zero off-diagonals replaced by Diagonal | MatrixMap(full, [2,0,0,3]) |
| matrixmap-06 | neg_matrix_full_offdiag.map | focused | negative | matrixmap.c:1894-1896 | Full MatrixMap with non-zero off-diagonal cannot self-simplify | MatrixMap(full, [1,2,3,4]) alone |
| matrixmap-07 | matrix_matrix_series_merge.map | focused | positive | matrixmap.c:1981-1984 | Two MatrixMaps in series merged via matrix multiplication | CmpMap(MatrixMap, MatrixMap), Series=1 |
| matrixmap-08 | matrix_zoom_series_merge.map | focused | positive | matrixmap.c:1986-1996 | MatrixMap + ZoomMap in series: elements scaled | CmpMap(MatrixMap, ZoomMap), Series=1 |
| matrixmap-09 | matrix_perm_series_merge.map | focused | positive | matrixmap.c:1998-2007 | MatrixMap + bidirectional PermMap in series merged via MatPerm | CmpMap(MatrixMap, PermMap[bidirectional]), Series=1 |
| matrixmap-10 | neg_matrix_perm_not_bidirectional.map | focused | negative | matrixmap.c:1959-1960 | Adjacent PermMap has inconsistent fwd/inv axes: merge blocked | CmpMap(MatrixMap, PermMap[non-bidirectional]), Series=1 |
| matrixmap-11 | matrix_diagwin_series_merge.map | focused | positive | matrixmap.c:2010-2019 | Diagonal MatrixMap + WinMap in series merged via MatWin2 | CmpMap(MatrixMap[diag], WinMap), Series=1 |
| matrixmap-12 | matrix_unit_series_merge.map | focused | positive | matrixmap.c:2022-2025 | MatrixMap + UnitMap in series: UnitMap eliminated | CmpMap(MatrixMap, UnitMap), Series=1 |
| matrixmap-13 | matrix_swap_past_win.map | cascade | positive | matrixmap.c:2163-2195 | MatrixMap swaps past WinMap to reach merge target | CmpMap(MatrixMap, WinMap, MatrixMap), Series=1 |
| matrixmap-14 | matrix_swap_past_perm.map | cascade | positive | matrixmap.c:2169-2195 | MatrixMap swaps past PermMap to reach merge target | CmpMap(MatrixMap, PermMap, MatrixMap), Series=1 |
| matrixmap-15 | matrix_swap_win_simplifies.map | cascade | positive | matrixmap.c:2201-2264 | Swap produces simpler Mapping (local simplification) | CmpMap(PermMap[dropping axes], MatrixMap), Series=1 |
| matrixmap-16 | neg_matrix_swap_refused.map | cascade | negative | matrixmap.c:2244-2248 | Swap refused: neither swapped Mapping simplifies | CmpMap(WinMap, MatrixMap[full]), Series=1 |
| matrixmap-17 | neg_matrix_parallel_no_merge.map | focused | negative | matrixmap.c:1930+2276 | Parallel mode with non-self-simplifiable MatrixMap: returns -1 | MatrixMap(full) in parallel |

## permmap.c

PermMap's MapMerge: accumulates permutation effects of adjacent PermMaps/
UnitMaps. Series composes sequentially; parallel concatenates side-by-side.
Result may be UnitMap (null permutation), simplified PermMap (cleared invert,
simplified arrays), or unchanged.

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
|---|---|---|---|---|---|---|
| permmap-01 | perm_series_merge.map | cascade | positive | permmap.c:910-1086 | Two or more adjacent PermMaps in series compose into one | CmpMap(PermMap, PermMap), Series=1 |
| permmap-02 | perm_parallel_merge.map | cascade | positive | permmap.c:1091-1280 | Adjacent PermMaps/UnitMaps in parallel compose into one wider PermMap | CmpMap(PermMap, UnitMap), Series=0 |
| permmap-03 | perm_cancel_to_unit.map | focused | positive | permmap.c:1333,1383-1384 | Composed PermMap reduces to UnitMap (both permutations null, nin==nout) | Two inverse PermMaps in series producing identity |
| permmap-04 | perm_invert_normalize.map | focused | positive | permmap.c:1345 | Single PermMap with Invert flag normalized (flag cleared, arrays swapped) | Lone PermMap with invert_list=1 |
| permmap-05 | perm_array_simplify.map | focused | positive | permmap.c:1354-1356 | PermMap simplified: previously-stored array now null after composition | PermMap + UnitMap in series where array becomes identity |
| permmap-06 | perm_inperm_constant_fold.map | focused | positive | permmap.c:1362-1367 | PermMap simplified: inperm array differs after constant folding | PermMap with constants composed with routing PermMap |
| permmap-07 | perm_outperm_constant_fold.map | focused | positive | permmap.c:1371-1376 | PermMap simplified: outperm array differs after re-computation | Similar to permmap-06 affecting outperm |
| permmap-08 | neg_perm_no_merge.map | focused | negative | permmap.c:1379 | No simplification: single canonical PermMap with no mergeable neighbours | Lone forward PermMap, invert=0 |
| permmap-09 | perm_constant_propagation.map | cascade | positive | permmap.c:1050-1078 | Series composition propagates constants through merged PermMap | PermMap with constant outputs + PermMap routing those outputs |
| permmap-10 | perm_bad_propagation.map | cascade | positive | permmap.c:1060-1061 | Series composition propagates AST__BAD through (negative perm index) | PermMap with out-of-range indices in series |

## lutmap.c

LutMap's MapMerge: (1) linear LutMap (detected by GetLinear) replaced by
WinMap; (2) inverse pair cancellation to UnitMap.

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
|---|---|---|---|---|---|---|
| lutmap-01 | lut_linear_to_win.map | focused | positive | lutmap.c:1139-1155 | Linear non-inverted LutMap replaced by WinMap | LutMap with linear table values |
| lutmap-02 | lut_linear_to_win.map | focused | positive | lutmap.c:1139,1156-1157 | Linear inverted LutMap replaced by reversed WinMap | LutMap(linear, Invert=1) |
| lutmap-03 | lut_inverse_cancel.map | focused | positive | lutmap.c:1178-1183,1209-1230 | LutMap cancels with equal upper-neighbour in opposite direction → UnitMap | CmpMap(LutMap, Inverse(LutMap)), Series=1 |
| lutmap-04 | lut_inverse_cancel.map | focused | positive | lutmap.c:1186-1190,1209-1230 | LutMap cancels with equal lower-neighbour in opposite direction → UnitMap | CmpMap(Inverse(LutMap), LutMap), Series=1 |
| lutmap-05 | neg_lut_nonlinear.map | focused | negative | lutmap.c:1139 | LutMap is not linear: WinMap replacement refused | LutMap with non-linear table |
| lutmap-06 | neg_lut_constant.map | focused | negative | lutmap.c:1150 | LutMap is linear but constant (b1==b2): WinMap impossible | LutMap([5,5,5,...]) |
| lutmap-07 | neg_lut_parallel.map | focused | negative | lutmap.c:1174 | Not in series: cancellation skipped | Two LutMaps in parallel |
| lutmap-08 | neg_lut_nonlut_neighbour.map | focused | negative | lutmap.c:1178-1193 | Neither neighbour is a LutMap | CmpMap(LutMap, ZoomMap), Series=1 |
| lutmap-09 | neg_lut_different_tables.map | focused | negative | lutmap.c:1209 | Neighbouring LutMap not inverse-equal (different tables) | Two different LutMaps in opposite directions |

## polymap.c

PolyMap's MapMerge: (1) combine duplicate terms; (2) linearize if only
constant and first-power terms → ShiftMap+MatrixMap; (3) inverse-pair
cancellation.

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
|---|---|---|---|---|---|---|
| polymap-01 | poly_duplicate_terms.map | focused | positive | polymap.c:3607-3717 | Duplicate-power terms combined, reducing coefficient count | PolyMap with coefficients sharing identical power vectors |
| polymap-02 | poly_duplicate_terms.map | focused | positive | polymap.c:3734-3848 | Linear PolyMap (after dedup) replaced by ShiftMap+MatrixMap | PolyMap with nin==nout, only constant and x^1 terms |
| polymap-03 | neg_poly_no_forward.map | focused | negative | polymap.c:3734 | Linearization refused: forward transform undefined (ncoeff_f NULL) | PolyMap defined only with inverse coefficients |
| polymap-04 | neg_poly_nin_ne_nout.map | focused | negative | polymap.c:3734 | Linearization refused: nin != nout | PolyMap(2-in, 3-out) with linear terms |
| polymap-05 | neg_poly_nonlinear.map | focused | negative | polymap.c:3810-3813 | Linearization refused: term has power > 1 or multiple inputs | PolyMap with quadratic term |
| polymap-06 | poly_inverse_cancel.map | focused | positive | polymap.c:3862-3929 | PolyMap + Inverse(PolyMap) cancel to UnitMap | CmpMap(PolyMap, Inverse(PolyMap)), Series=1 |
| polymap-07 | neg_poly_parallel_nonlinear.map | focused | negative | polymap.c:3863 | Inverse-cancel refused: combination is parallel | Two PolyMaps in parallel |
| polymap-08 | neg_poly_nonpoly_neighbour.map | focused | negative | polymap.c:3881 | Inverse-cancel refused: neighbour is not PolyMap | CmpMap(PolyMap, ZoomMap), Series=1 |
| polymap-09 | neg_poly_same_direction.map | focused | negative | polymap.c:3887 | Inverse-cancel refused: neighbour has same invert direction | Two forward PolyMaps in series |
| polymap-10 | neg_poly_different_coeffs.map | focused | negative | polymap.c:3901 | Inverse-cancel refused: astEqual fails (different coefficients) | Two different PolyMaps in opposite directions |

## mathmap.c

MathMap's MapMerge: checks SimpFI/SimpIF permission, then compares function
text of adjacent MathMaps. If forward(first)==inverse(second) and vice versa,
pair cancels to UnitMap.

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
|---|---|---|---|---|---|---|
| mathmap-01 | math_inverse_cancel.map | focused | positive | mathmap.c:4132-4159 | Two MathMaps with matching fwd/inv text cancel to UnitMap | CmpMap(MathMap[SimpFI=1], Inverse(MathMap[SimpIF=1])), Series=1 |
| mathmap-02 | — | focused | negative | mathmap.c:4033 | Parallel mode: refuses | MathMaps in parallel |
| mathmap-03 | neg_math_lone.map | focused | negative | mathmap.c:4041 | No following Mapping (last in list) | Single MathMap |
| mathmap-04 | neg_math_nonmath_neighbour.map | focused | negative | mathmap.c:4047-4051 | Neighbour is not a MathMap | CmpMap(MathMap, ZoomMap), Series=1 |
| mathmap-05 | neg_math_no_simpfi.map | focused | negative | mathmap.c:4064-4067 | SimpFI/SimpIF not set: simplification refused | CmpMap(MathMap[SimpFI=0], Inverse(MathMap)), Series=1 |
| mathmap-06 | — | focused | negative | mathmap.c:4082 | Dimension mismatch: nin(first) != nout(second) | MathMaps of different dimensionality |
| mathmap-07 | — | focused | negative | mathmap.c:4095 | Forward function count mismatch | MathMaps with different output counts |
| mathmap-08 | — | focused | negative | mathmap.c:4107-4110 | Forward function text of first != inverse text of second | Different MathMaps with SimpFI set |
| mathmap-09 | — | focused | negative | mathmap.c:4119 | Inverse function count mismatch | MathMaps whose inv/fwd counts differ |
| mathmap-10 | neg_math_inv_text_mismatch.map | focused | negative | mathmap.c:4124-4127 | Inverse function text of first != forward text of second | MathMaps whose inv(first) != fwd(second) |

## ratemap.c

RateMap's MapMerge: simplifies encapsulated Mapping, then checks both
neighbours for inverse pair cancellation.

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
|---|---|---|---|---|---|---|
| ratemap-01 | ratemap_simplify_interior.map | focused | positive | ratemap.c:723-726 | Encapsulated Mapping simplifies → new RateMap with simplified interior | RateMap(CmpMap(ZoomMap,ZoomMap)) |
| ratemap-02 | ratemap_inverse_cancel.map | focused | positive | ratemap.c:736-761 | RateMap cancels with equal lower-neighbour in opposite direction | CmpMap(Inverse(RateMap), RateMap), Series=1 |
| ratemap-03 | ratemap_inverse_cancel.map | focused | positive | ratemap.c:766-789 | RateMap cancels with equal upper-neighbour in opposite direction | CmpMap(RateMap, Inverse(RateMap)), Series=1 |
| ratemap-04 | neg_ratemap_parallel.map | focused | negative | ratemap.c:730 | Not in series: skipped | RateMaps in parallel |
| ratemap-05 | neg_ratemap_lower_nonrate.map | focused | negative | ratemap.c:736 | Lower neighbour not a RateMap | CmpMap(ZoomMap, RateMap), Series=1 |
| ratemap-06 | neg_ratemap_same_invert.map | focused | negative | ratemap.c:739 | Lower neighbour same invert flag | Two forward RateMaps in series |
| ratemap-07 | neg_ratemap_diff_indices.map | focused | negative | ratemap.c:744-745 | Lower neighbour different iin/iout indices | RateMaps with different axis indices |
| ratemap-08 | neg_ratemap_different_inner.map | focused | negative | ratemap.c:755 | Lower neighbour non-equal encapsulated Mapping | RateMaps wrapping different Mappings |
| ratemap-09 | neg_ratemap_nonratemap_neighbour.map | focused | negative | ratemap.c:767 | Upper neighbour not a RateMap | CmpMap(RateMap, ZoomMap), Series=1 |
| ratemap-10 | neg_ratemap_same_invert.map | focused | negative | ratemap.c:769 | Upper neighbour same invert flag | Two forward RateMaps |
| ratemap-11 | neg_ratemap_diff_indices.map | focused | negative | ratemap.c:772-773 | Upper neighbour different iin/iout | Different axis RateMaps |
| ratemap-12 | neg_ratemap_diff_inner.map | focused | negative | ratemap.c:781 | Upper neighbour non-equal encapsulated Mapping | Different inner Mappings |

## intramap.c

IntraMap's MapMerge: checks following neighbour for same function/IntraFlag
and opposite direction, with AST__SIMPFI/SIMPIF permission.

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
|---|---|---|---|---|---|---|
| intramap-01 | intramap_inverse_cancel.map | focused | positive | intramap.c:1567-1568 | Forward + Inverse IntraMap with SIMPFI cancels to UnitMap | CmpMap(IntraMap, Inverse(IntraMap)), Series=1, SIMPFI set |
| intramap-02 | — | focused | positive | intramap.c:1573-1574 | Inverse + Forward IntraMap with SIMPIF cancels to UnitMap | CmpMap(Inverse(IntraMap), IntraMap), Series=1, SIMPIF set |
| intramap-03 | — | focused | negative | intramap.c:1514 | Not in series or no following Mapping | IntraMaps in parallel |
| intramap-04 | — | focused | negative | intramap.c:1525 | Following Mapping is not IntraMap | CmpMap(IntraMap, ZoomMap), Series=1 |
| intramap-05 | — | focused | negative | intramap.c:1532 | Different transformation functions (ifun differs) | IntraMaps with different registered functions |
| intramap-06 | — | focused | negative | intramap.c:1533-1534 | IntraFlag strings differ | Same function but different flags |
| intramap-07 | — | focused | negative | intramap.c:1562 | Dimension mismatch | Asymmetric IntraMaps |
| intramap-08 | — | focused | negative | intramap.c:1567-1575 | Same direction (both forward or both inverse) | Two forward IntraMaps in series |
| intramap-09 | — | focused | negative | intramap.c:1568 | SIMPFI flag not set on forward-then-inverse pair | IntraMaps without SIMPFI permission |
| intramap-10 | — | focused | negative | intramap.c:1574 | SIMPIF flag not set on inverse-then-forward pair | IntraMaps without SIMPIF permission |

## splinemap.c

SplineMap's MapMerge: checks both neighbours for equal SplineMap in opposite
direction via astEqual.

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
|---|---|---|---|---|---|---|
| splinemap-01 | spline_inverse_cancel.map | focused | positive | splinemap.c:1800-1805 | Lower-neighbour SplineMap equal and opposite → UnitMaps | CmpMap(Inverse(SplineMap), SplineMap), Series=1 |
| splinemap-02 | spline_inverse_cancel.map | focused | positive | splinemap.c:1806-1812 | Upper-neighbour SplineMap equal and opposite → UnitMaps | CmpMap(SplineMap, Inverse(SplineMap)), Series=1 |
| splinemap-03 | neg_spline_parallel.map | focused | negative | splinemap.c:1746 | Not in series: skipped | SplineMaps in parallel |
| splinemap-04 | — | focused | negative | splinemap.c:1763 | No neighbour (boundary) | Single SplineMap |
| splinemap-05 | neg_spline_nonspline_neighbour.map | focused | negative | splinemap.c:1766 | Neighbour is not a SplineMap | CmpMap(SplineMap, ZoomMap), Series=1 |
| splinemap-06 | neg_spline_same_direction.map | focused | negative | splinemap.c:1773 | Same invert flag (same direction) | Two forward SplineMaps in series |
| splinemap-07 | neg_spline_different_coeffs.map | focused | negative | splinemap.c:1794 | astEqual fails (different coefficients) | Different SplineMaps in opposite directions |

## normmap.c

NormMap's MapMerge: (1) simplify encapsulated Frame; (2) basic Frame→UnitMap;
(3) inverse-pair cancellation; (4) duplicate-NormMap elimination.

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
|---|---|---|---|---|---|---|
| normmap-01 | — | focused | positive | normmap.c:626-630 | Encapsulated Frame simplifies → new NormMap with simplified Frame | NormMap wrapping compound Frame that simplifies |
| normmap-02 | normmap_basic_frame_to_unit.map | focused | positive | normmap.c:635-639 | NormMap encapsulating basic Frame replaced by UnitMap (astNorm is no-op) | NormMap wrapping plain Frame |
| normmap-03 | normmap_inverse_cancel.map | focused | positive | normmap.c:650-661 | NormMap cancels with inverse lower-neighbour NormMap | NormMap preceded by Inverse(NormMap) with same Frame |
| normmap-04 | normmap_inverse_cancel_upper.map | focused | positive | normmap.c:666-677 | NormMap cancels with inverse upper-neighbour NormMap | NormMap followed by Inverse(NormMap) with same Frame |
| normmap-05 | normmap_duplicate_elim.map | focused | positive | normmap.c:697-717 | Duplicate adjacent NormMaps (same Frame, same direction) → duplicates become UnitMaps | Two identical NormMaps in series |
| normmap-06 | neg_normmap_diff_frame.map | focused | negative | normmap.c:650-651 | Lower neighbour NormMap: invert flags not opposite | NormMap preceded by same-direction NormMap with different Frame |
| normmap-07 | neg_normmap_different_frames.map | focused | negative | normmap.c:660 | Lower inverse NormMap: Frames not equal | NormMap preceded by Inverse(NormMap) with different Frame |
| normmap-08 | neg_normmap_nonnorm_neighbour.map | focused | negative | normmap.c:666-667 | Upper neighbour is not a NormMap | NormMap followed by non-NormMap |
| normmap-09 | neg_normmap_inv_diff_frame.map | focused | negative | normmap.c:674 | Upper inverse NormMap: Frames differ | NormMap followed by Inverse(NormMap) with different Frame |
| normmap-10 | neg_normmap_diff_frame.map | focused | negative | normmap.c:703,707 | Adjacent same-direction NormMap: Frames differ | Two NormMaps same direction, different Frames |
| normmap-11 | neg_normmap_parallel.map | focused | negative | normmap.c:644 | Parallel mode: no simplification beyond Frame-level | NormMap in parallel |
| normmap-12 | neg_normmap_parallel_sky.map | focused | negative | normmap.c:626+635+644 | Non-basic Frame, doesn't simplify, not in series | NormMap(SkyFrame) in parallel |

## unitnormmap.c

UnitNormMap's MapMerge: merges with adjacent ShiftMap/WinMap (unit-scale) by
adjusting centre, or cancels with inverse UnitNormMap (same centre→UnitMap,
different→ShiftMap).

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
|---|---|---|---|---|---|---|
| unitnormmap-01 | unitnormmap_shift_fwd_merge.map | focused | positive | unitnormmap.c:306-322 | ShiftMap + forward UnitNormMap → new UnitNormMap with adjusted centre | ShiftMap followed by UnitNormMap(fwd) |
| unitnormmap-02 | unitnormmap_winmap_fwd_merge.map | focused | positive | unitnormmap.c:323-343 | WinMap(unit scale) + forward UnitNormMap → new UnitNormMap with adjusted centre | WinMap(scale=1) followed by UnitNormMap(fwd) |
| unitnormmap-03 | neg_unitnormmap_nonunit_scale.map | focused | negative | unitnormmap.c:338,341 | WinMap(non-unit scale) + UnitNormMap: refused | WinMap(scale!=1) followed by UnitNormMap(fwd) |
| unitnormmap-04 | unitnormmap_inv_shift_merge.map | focused | positive | unitnormmap.c:344-360 | Inverse UnitNormMap + ShiftMap → new inverse UnitNormMap with adjusted centre | UnitNormMap(inv) followed by ShiftMap |
| unitnormmap-05 | unitnormmap_inv_winmap_merge.map | focused | positive | unitnormmap.c:361-380 | Inverse UnitNormMap + WinMap(unit scale) → new inverse UnitNormMap | UnitNormMap(inv) followed by WinMap(scale=1) |
| unitnormmap-06 | — | focused | negative | unitnormmap.c:375,379 | UnitNormMap(inv) + WinMap(non-unit scale): refused | UnitNormMap(inv) followed by WinMap(scale!=1) |
| unitnormmap-07 | unitnormmap_inverse_cancel.map | focused | positive | unitnormmap.c:398-401 | Forward + Inverse UnitNormMap with same centre → UnitMap | UnitNormMap(fwd) + Inverse(UnitNormMap), same centre |
| unitnormmap-08 | unitnormmap_inv_fwd_cancel.map | focused | positive | unitnormmap.c:398-401 | Inverse + Forward UnitNormMap with same centre → UnitMap | Inverse(UnitNormMap) + UnitNormMap(fwd), same centre |
| unitnormmap-09 | unitnormmap_diff_centre_to_shift.map | focused | positive | unitnormmap.c:403-415 | Forward + Inverse UnitNormMap with different centres → ShiftMap | UnitNormMap(fwd,c1) + Inverse(UnitNormMap,c2) |
| unitnormmap-10 | — | focused | negative | unitnormmap.c:403 | Inverse + Forward with different centres: no merge (asymmetric) | Inverse(UnitNormMap,c1) + UnitNormMap(fwd,c2) |
| unitnormmap-11 | — | focused | negative | unitnormmap.c:307 | ShiftMap + UnitNormMap(inv): refused | ShiftMap followed by Inverse(UnitNormMap) |
| unitnormmap-12 | — | focused | negative | unitnormmap.c:324 | WinMap + UnitNormMap(inv): refused | WinMap followed by Inverse(UnitNormMap) |
| unitnormmap-13 | — | focused | negative | unitnormmap.c:345 | UnitNormMap(fwd) + ShiftMap: refused | Forward UnitNormMap followed by ShiftMap |
| unitnormmap-14 | — | focused | negative | unitnormmap.c:362 | UnitNormMap(fwd) + WinMap: refused | Forward UnitNormMap followed by WinMap |
| unitnormmap-15 | — | focused | negative | unitnormmap.c:383 | Two UnitNormMaps in same direction: refused | Two forward UnitNormMaps in series |
| unitnormmap-16 | — | focused | negative | unitnormmap.c:302 | Neighbour is not ShiftMap/WinMap/UnitNormMap: refused | UnitNormMap with non-mergeable neighbour |
| unitnormmap-17 | — | focused | negative | unitnormmap.c:791-795 | Parallel mode: never simplifies | UnitNormMap in parallel |
| unitnormmap-18 | — | focused | negative | unitnormmap.c:763 | No neighbour pair produces valid merge | UnitNormMap flanked by non-mergeable classes |

## selectormap.c

SelectorMap's MapMerge: (1) simplify internal regions; (2) inverse-pair
cancellation.

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
|---|---|---|---|---|---|---|
| selectormap-01 | — | focused | positive | selectormap.c:698-714 | Internal regions simplify → new SelectorMap with simplified regions | SelectorMap containing simplifiable Regions |
| selectormap-02 | — | focused | negative | selectormap.c:700-704,727 | No region simplification and no adjacent SelectorMap | SelectorMap with already-simple regions, alone |
| selectormap-03 | — | focused | negative | selectormap.c:731-744 | No adjacent SelectorMap found in series | SelectorMap flanked by non-SelectorMaps |
| selectormap-04 | — | focused | negative | selectormap.c:751-755 | Adjacent SelectorMap not equal-and-opposite | Two SelectorMaps with different regions |
| selectormap-05 | selectormap_inverse_cancel.map | focused | positive | selectormap.c:758-779 | Inverse-pair cancellation → UnitMap | SelectorMap + Inverse(SelectorMap), identical regions |
| selectormap-06 | — | focused | negative | selectormap.c:727 | Parallel mode: Phase 2 skipped | SelectorMap in parallel |

## switchmap.c

SwitchMap's MapMerge: (1) inverse-pair cancellation; (2) invert-flag
normalization (swap selectors, invert routes); (3) simplify internal
selectors/routes.

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
|---|---|---|---|---|---|---|
| switchmap-01 | switchmap_inverse_cancel.map | focused | positive | switchmap.c:1026-1048 | Inverse-pair cancellation → UnitMap | SwitchMap + Inverse(SwitchMap), identical selectors/routes |
| switchmap-02 | — | focused | negative | switchmap.c:999-1013,1021 | Adjacent SwitchMap not equal-and-opposite | Two SwitchMaps with different routes |
| switchmap-03 | switchmap_invert_normalize.map | focused | positive | switchmap.c:1063-1075 | Inverted SwitchMap normalized to non-inverted equivalent | SwitchMap with Invert=1 |
| switchmap-04 | switchmap_internal_simplify.map | focused | positive | switchmap.c:1079-1098 | Internal selectors/routes simplify → new SwitchMap | SwitchMap with simplifiable route maps |
| switchmap-05 | — | focused | negative | switchmap.c:995,1055,1079-1089 | Series: no adjacent match, non-inverted, internals don't simplify | Lone simple SwitchMap in series |
| switchmap-06 | — | focused | negative | switchmap.c:995,1055,1079-1089 | Parallel: same as above | Lone simple SwitchMap in parallel |

## tranmap.c

TranMap's MapMerge: (1) invert normalization (swap+invert components);
(2) simplify each component; (3) equal-component detection (→single Mapping);
(4) adjacent TranMap series merge.

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
|---|---|---|---|---|---|---|
| tranmap-01 | tranmap_invert_normalize.map | focused | positive | tranmap.c:839-849 | Inverted TranMap normalized by swapping and inverting components | TranMap with Invert=1 |
| tranmap-02 | tranmap_component_simplify.map | focused | positive | tranmap.c:852-861 | Component Mappings individually simplified, TranMap rebuilt | TranMap(CmpMap(Z,Z), UnitMap) |
| tranmap-03 | tranmap_equal_components.map | focused | positive | tranmap.c:865-886 | Both components bidirectional and equal → single component Mapping | TranMap(ZoomMap[2], ZoomMap[2]) |
| tranmap-04 | tranmap_adjacent_merge.map | cascade | positive | tranmap.c:902-1001 | Two adjacent TranMaps in series merge by combining fwd/inv legs | CmpMap(TranMap(A,B), TranMap(C,D)), Series=1 |
| tranmap-05 | neg_tranmap_parallel.map | focused | negative | tranmap.c:902 | Parallel mode: adjacent merge skipped | TranMaps in parallel |
| tranmap-06 | neg_tranmap_oneway.map | focused | negative | tranmap.c:865-866 | Equal-component check skipped: components lack bidirectional transforms | TranMap(OneWayFwd, OneWayInv) |
| tranmap-07 | neg_tranmap_unequal_components.map | focused | negative | tranmap.c:881 | CmpMap(fwd,inv(inv)) doesn't simplify to UnitMap: components not equal | TranMap(ShiftMap(1), ShiftMap(2)) |
| tranmap-08 | neg_tranmap_nontranmap_neighbour.map | focused | negative | tranmap.c:905-906 | Higher neighbour is not a TranMap | CmpMap(TranMap, ZoomMap), Series=1 |
| tranmap-09 | neg_tranmap_adjacent_no_merge.map | focused | negative | tranmap.c:959 | Neither fwd nor inv series combination simplified | CmpMap(TranMap(A,B), TranMap(C,D)) with irreducible legs |

## timemap.c

TimeMap's MapMerge: gather adjacent TimeMap steps, eliminate inverse pairs by
argument count (0,1,2,3,5-arg), eliminate no-op steps.

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
|---|---|---|---|---|---|---|
| timemap-01 | time_inverse_cancel.map | focused | positive | timemap.c:2380-2397 | Full cancellation: all steps cancel → UnitMap | Two TimeMaps with inverse steps |
| timemap-02 | time_partial_cancel.map | focused | positive | timemap.c:2401-2406 | Partial cancellation: some steps cancel, result is simplified TimeMap | TimeMap with 3+ steps, one pair cancels |
| timemap-03 | time_merge_no_cancel.map | cascade | positive | timemap.c:2375-2389 | Multi-map merge without step reduction: adjacent TimeMaps merged | Two TimeMaps with non-cancelling steps |
| timemap-04 | time_invert_normalize.map | focused | positive | timemap.c:2388-2389 | Invert-flag clearing: single inverted TimeMap rebuilt with invert=0 | Single TimeMap with Invert=1 |
| timemap-05 | neg_time_parallel.map | focused | negative | timemap.c:2102 | Parallel-mode guard | TimeMap in parallel |
| timemap-06 | neg_time_lone_forward.map | focused | negative | timemap.c:2388-2392 | No simplification: single forward TimeMap, no neighbours | Lone forward TimeMap |
| timemap-07 | time_noop_eliminate.map | focused | positive | timemap.c:2249-2253 | No-op step elimination: MJDTOMJD with zero offset removed | TimeMap with MJDTOMJD(0,0) step |
| timemap-08 | time_inverse_cancel.map | focused | positive | timemap.c:2276-2285 | 1-arg pair cancellation (TAITOTT+TTTOTAI etc.) | Adjacent 1-arg inverse steps with matching arg |
| timemap-09 | time_2arg_swapped_cancel.map | focused | positive | timemap.c:2289-2298 | 2-arg pair cancellation (swapped args: MJDTOJD+JDTOMJD) | Adjacent 2-arg steps with swapped matching args |
| timemap-10 | time_2arg_same_cancel.map | focused | positive | timemap.c:2302-2308 | 2-arg pair cancellation (same order: TAITOUTC+UTCTOTAI) | Adjacent 2-arg steps with same-order args |
| timemap-11 | time_3arg_cancel.map | focused | positive | timemap.c:2311-2321 | 3-arg pair cancellation (GMSTTOLMST+LMSTTOGMST) | Adjacent 3-arg steps with matching args |
| timemap-12 | time_5arg_cancel.map | focused | positive | timemap.c:2324-2336 | 5-arg pair cancellation (TTTOTDB+TDBTOTT) | Adjacent 5-arg steps with matching args |
| timemap-13 | neg_time_arg_mismatch.map | focused | negative | timemap.c:2276-2285 | 1-arg pair with mismatched argument | Steps with different DUT1 values |
| timemap-14 | — | focused | negative | timemap.c:2289-2298 | 2-arg swapped pair with mismatched arguments | MJDTOJD + JDTOMJD with different offsets |
| timemap-15 | neg_time_3arg_mismatch.map | focused | negative | timemap.c:2311-2321 | 3-arg pair with mismatched arguments | GMSTTOLMST + LMSTTOGMST with different lon |
| timemap-16 | neg_time_5arg_mismatch.map | focused | negative | timemap.c:2324-2336 | 5-arg pair with mismatched arguments | TTTOTDB + TDBTOTT with one arg different |
| timemap-17 | time_run_identity.map | cascade | positive | timemap.c:2380-2397 | Run of three TimeMaps whose combined steps cancel to a UnitMap | Three adjacent TimeMaps forming an identity conversion |

## slamap.c

SlaMap's MapMerge: gather adjacent SlaMap steps, eliminate inverse pairs by
argument count (0,1,2,4-arg), eliminate redundant precession, merge adjacent
precession.

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
|---|---|---|---|---|---|---|
| slamap-01 | sla_inverse_cancel.map | focused | positive | slamap.c:3120-3137 | Full cancellation: all steps cancel → UnitMap | Two SlaMaps with inverse steps |
| slamap-02 | sla_partial_cancel.map | focused | positive | slamap.c:3141-3147 | Partial cancellation: some steps cancel, result is simplified SlaMap | SlaMap with 3+ steps, one pair cancels |
| slamap-03 | sla_merge_no_cancel.map | cascade | positive | slamap.c:3115-3128 | Multi-map merge without step reduction: adjacent SlaMaps merged | Two SlaMaps with non-cancelling steps |
| slamap-04 | sla_invert_normalize.map | focused | positive | slamap.c:3128-3129 | Invert-flag clearing: single inverted SlaMap rebuilt with invert=0 | Single SlaMap with Invert=1 |
| slamap-05 | neg_sla_parallel.map | focused | negative | slamap.c:2706 | Parallel-mode guard | SlaMap in parallel |
| slamap-06 | neg_sla_lone_forward.map | focused | negative | slamap.c:3128-3132 | No simplification: single forward SlaMap, no neighbours | Lone forward SlaMap |
| slamap-07 | sla_inverse_cancel.map | focused | positive | slamap.c:3011-3014 | 0-arg pair: galactic (EQGAL+GALEQ) | Adjacent galactic steps |
| slamap-08 | sla_supergalactic_cancel.map | focused | positive | slamap.c:3019-3021 | 0-arg pair: supergalactic (GALSUP+SUPGAL) | Adjacent supergalactic steps |
| slamap-09 | sla_j2000_cancel.map | focused | positive | slamap.c:3065-3068 | 0-arg pair: dynamical J2000 (J2000H+HJ2000) | Adjacent J2000 steps |
| slamap-10 | sla_eterms_cancel.map | focused | positive | slamap.c:2943-2948 | 1-arg pair: E-terms (ADDET+SUBET) | Adjacent E-term steps, same epoch |
| slamap-11 | sla_fk45_cancel.map | focused | positive | slamap.c:2954-2959 | 1-arg pair: FK4/FK5 (FK45Z+FK54Z) | Adjacent FK conversion steps |
| slamap-12 | sla_icrs_cancel.map | focused | positive | slamap.c:2965-2970 | 1-arg pair: ICRS/FK5 (HFK5Z+FK5HZ) | Adjacent ICRS steps |
| slamap-13 | sla_ecliptic_cancel.map | focused | positive | slamap.c:2989-2994 | 1-arg pair: ecliptic (ECLEQ+EQECL) | Adjacent ecliptic steps |
| slamap-14 | sla_helioecl_cancel.map | focused | positive | slamap.c:3056-3061 | 1-arg pair: helio-ecliptic (EQHE+HEEQ) | Adjacent helio-ecliptic steps |
| slamap-15 | sla_ha_cancel.map | focused | positive | slamap.c:3072-3078 | 1-arg pair: HA (R2H+H2R) | Adjacent HA steps |
| slamap-16 | sla_geocentric_cancel.map | focused | positive | slamap.c:2977-2984 | 2-arg pair (cross-matched): geocentric (AMP+MAP) | Adjacent AMP/MAP with crossed args |
| slamap-17 | sla_azel_cancel.map | focused | positive | slamap.c:2998-3005 | 2-arg pair (same order): AzEl (DH2E+DE2H) | Adjacent AzEl steps |
| slamap-18 | sla_hpc_cancel.map | focused | positive | slamap.c:3026-3037 | 4-arg pair: helioprojective-Cartesian (HPCEQ+EQHPC) | Adjacent HPC steps |
| slamap-19 | sla_hpr_cancel.map | focused | positive | slamap.c:3041-3052 | 4-arg pair: helioprojective-Radial (HPREQ+EQHPR) | Adjacent HPR steps |
| slamap-20 | neg_sla_arg_mismatch.map | focused | negative | slamap.c:2943-2948 | 1-arg pair: mismatched argument | Steps with different epochs |
| slamap-21 | neg_sla_2arg_mismatch.map | focused | negative | slamap.c:2977-2984 | 2-arg pair: mismatched arguments | AMP+MAP with non-matching args |
| slamap-22 | neg_sla_4arg_mismatch.map | focused | negative | slamap.c:3026-3037 | 4-arg pair: mismatched arguments | HPC steps with one arg different |
| slamap-23 | sla_prec_redundant.map | focused | positive | slamap.c:2908-2911 | Redundant precession: PREC/PREBN with start==end eliminated | PREC(2000,2000) |
| slamap-24 | sla_prec_merge.map | focused | positive | slamap.c:2929-2936 | Adjacent precession merge: PREC(a,b)+PREC(b,c) → PREC(a,c) | Adjacent PREC steps with common equinox |
| slamap-25 | neg_sla_prec_no_common.map | focused | negative | slamap.c:2929-2931 | Adjacent precession: non-common equinox prevents merge | PREC(1950,1975) + PREC(2000,2025) |
| slamap-26 | sla_run_identity.map | cascade | positive | slamap.c:3120-3137 | Run of three SlaMaps whose combined steps cancel to a UnitMap | Three adjacent SlaMaps forming an identity conversion |
| slamap-27 | sla_run_partial.map | cascade | positive | slamap.c:3141-3147 | Run of three SlaMaps partially cancels to a single SlaMap | Three adjacent SlaMaps, a subset of steps cancel |

## specmap.c

SpecMap's MapMerge: gather adjacent SpecMap steps, eliminate inverse pairs by
argument count (0,1,2,3,6-arg).

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
|---|---|---|---|---|---|---|
| specmap-01 | spec_inverse_cancel.map | focused | positive | specmap.c:2557-2574 | Full cancellation: all steps cancel → UnitMap | Two SpecMaps with inverse steps |
| specmap-02 | spec_partial_cancel.map | focused | positive | specmap.c:2578-2583 | Partial cancellation: some steps cancel, simplified SpecMap | SpecMap with 3+ steps, one pair cancels |
| specmap-03 | spec_merge_no_cancel.map | cascade | positive | specmap.c:2552-2566 | Multi-map merge without step reduction | Two adjacent SpecMaps with non-inverse steps |
| specmap-04 | spec_invert_normalize.map | focused | positive | specmap.c:2565-2566 | Invert-flag clearing | Single SpecMap with Invert=1 |
| specmap-05 | neg_spec_parallel.map | focused | negative | specmap.c:2259 | Parallel-mode guard | SpecMap in parallel |
| specmap-06 | neg_spec_lone_forward.map | focused | negative | specmap.c:2565-2569 | No simplification: single forward SpecMap, no neighbours | Lone forward SpecMap |
| specmap-07 | — | focused | negative | specmap.c:2274-2275 | Nin-mismatch guard: adjacent SpecMap with different nin | SpecMap(nin=1) adjacent to SpecMap(nin=3) |
| specmap-08 | spec_unit_cancel.map | focused | positive | specmap.c:2453-2462 | 0-arg pair: unit conversions (ENTOFR+FRTOEN etc.) | Adjacent 0-arg inverse steps |
| specmap-09 | spec_inverse_cancel.map | focused | positive | specmap.c:2465-2469 | 1-arg pair (FRTOVL+VLTOFR) | Adjacent 1-arg steps, matching arg |
| specmap-10 | spec_lsr_cancel.map | focused | positive | specmap.c:2472-2481 | 2-arg pair: local-standard (LKF2HL+HLF2LK) | Adjacent 2-arg steps, matching args |
| specmap-11 | spec_geocentric_cancel.map | focused | positive | specmap.c:2484-2494 | 3-arg pair: geocentric/barycentric (GEF2HL+HLF2GE) | Adjacent 3-arg steps |
| specmap-12 | spec_topocentric_cancel.map | focused | positive | specmap.c:2498-2512 | 6-arg pair: topocentric (TPF2HL+HLF2TP) | Adjacent 6-arg steps |
| specmap-13 | neg_spec_arg_mismatch.map | focused | negative | specmap.c:2465-2469 | 1-arg pair: mismatched argument | Steps with different rest frequencies |
| specmap-14 | neg_spec_2arg_mismatch.map | focused | negative | specmap.c:2472-2481 | 2-arg pair: mismatched arguments | Steps with different RA/dec |
| specmap-15 | neg_spec_3arg_mismatch.map | focused | negative | specmap.c:2484-2494 | 3-arg pair: mismatched arguments | Steps with different epoch |
| specmap-16 | neg_spec_6arg_mismatch.map | focused | negative | specmap.c:2498-2512 | 6-arg pair: mismatched arguments | Steps with one arg different |
| specmap-17 | spec_run4_identity.map | scenario | positive | specmap.c:2557-2574 | Run of four SpecMaps whose combined steps cancel to a UnitMap | Four adjacent SpecMaps forming an identity conversion |
| specmap-18 | spec_run_inverted_mid.map | cascade | positive | specmap.c:2557-2574 | Run of three SpecMaps with an inverted middle step partially cancels to one SpecMap | Three adjacent SpecMaps, middle one inverted |

## wcsmap.c

WcsMap's MapMerge: (1) AST__WCSBAD→UnitMap; (2) inverse-pair cancellation via
CanMerge; (3) swap past PermMap toward merge target or for local
simplification.

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
|---|---|---|---|---|---|---|
| wcsmap-01 | wcsmap_bad_to_unit.map | focused | positive | wcsmap.c:3088-3099 | WcsMap with AST__WCSBAD type replaced by UnitMap | WcsMap(AST__WCSBAD) |
| wcsmap-02 | wcs_inverse_cancel.map | focused | positive | wcsmap.c:3136-3161 | Adjacent inverse WcsMap pair cancels to UnitMap | CmpMap(WcsMap[TAN], Inverse(WcsMap[TAN])), Series=1 |
| wcsmap-03 | wcsmap_perm_swap_cancel.map | cascade | positive | wcsmap.c:3269-3273 | WcsMap swaps past PermMap to reach inverse merge target | WcsMap + PermMap + Inverse(WcsMap), Series=1 |
| wcsmap-04 | wcsmap_perm_swap_simplify.map | cascade | positive | wcsmap.c:3278-3332 | Swap with PermMap for local simplification (no merge target) | WcsMap + PermMap(adds axes) where swap simplifies |
| wcsmap-05 | — | cascade | negative | wcsmap.c:3307-3311 | Speculative swap refused: neither Mapping simplifies | WcsMap + PermMap(identity) with no merge target |
| wcsmap-06 | neg_wcs_parallel.map | focused | negative | wcsmap.c:3105 | Parallel mode or nmap==1: refused | WcsMap in parallel |
| wcsmap-07 | neg_wcs_different_projection.map | focused | negative | wcsmap.c:3120-3131 | Neighbour not a WcsMap or different projection type | WcsMap[TAN] adjacent to WcsMap[SIN] |
| wcsmap-08 | neg_wcs_same_direction.map | focused | negative | wcsmap.c:CanMerge:880 | Two WcsMaps same invert direction | Two forward WcsMap[TAN] in series |
| wcsmap-09 | — | focused | negative | wcsmap.c:CanMerge:884-885 | Lon/lat axis indices differ | WcsMap pair with different axis assignments |
| wcsmap-10 | neg_wcs_different_params.map | focused | negative | wcsmap.c:CanMerge:896-910 | Projection parameters differ | WcsMap pair with different PV values |
| wcsmap-11 | neg_wcs_nonperm_between.map | cascade | negative | wcsmap.c:3209-3212 | Intervening Mapping not a PermMap: swap blocked | WcsMap + MatrixMap + Inverse(WcsMap) |
| wcsmap-12 | — | cascade | negative | wcsmap.c:CanSwap:1056-1073 | PermMap has non-bidirectional links: swap refused | WcsMap + PermMap(one-way) |
| wcsmap-13 | — | cascade | negative | wcsmap.c:CanSwap:1086-1100 | PermMap doesn't pass through both lon/lat: swap refused | WcsMap + PermMap(disconnects one WCS axis) |

## sphmap.c

SphMap's MapMerge: (1) inverse-pair cancellation (UnitRadius or matching
PolarLong); (2) matrix sandwich — Inv(SphMap)+DiagMatrix/ZoomMap+SphMap →
WinMap (coordinate reflection).

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
|---|---|---|---|---|---|---|
| sphmap-01 | sph_inverse_cancel.map | focused | positive | sphmap.c:751-793 | Inverse(SphMap) + SphMap with matching PolarLong cancels to UnitMap | CmpMap(Inverse(SphMap[UntRd=1]), SphMap), Series=1 |
| sphmap-02 | sph_fwd_inv_unitradius_cancel.map | focused | positive | sphmap.c:759-760 | SphMap(UnitRadius) + Inverse(SphMap) cancels to UnitMap | CmpMap(SphMap(UntRd=1), Inverse(SphMap)), Series=1 |
| sphmap-03 | neg_sph_polarlong_mismatch.map | focused | negative | sphmap.c:751-753 | Inverse+Forward but PolarLong values differ | Inverse(SphMap(PolarLong=0)) + SphMap(PolarLong=pi) |
| sphmap-04 | neg_sph_no_unitradius.map | focused | negative | sphmap.c:759-760 | Forward+Inverse but UnitRadius not set | SphMap(UntRd=0) + Inverse(SphMap) |
| sphmap-05 | — | focused | negative | sphmap.c:751 | Same direction (both forward or both inverse) | Two forward SphMaps in series |
| sphmap-06 | neg_sph_parallel.map | focused | negative | sphmap.c:737 | Parallel mode or last in list | SphMap in parallel |
| sphmap-07 | neg_sph_non_sphmap_neighbour.map | focused | negative | sphmap.c:746 | Following Mapping is not a SphMap | SphMap + ZoomMap in series |
| sphmap-08 | sph_matrix_sandwich.map | cascade | positive | sphmap.c:810-949 | Inv(SphMap) + DiagMatrixMap + SphMap → WinMap | Inv(SphMap) + MatrixMap(diag,equal mag) + SphMap |
| sphmap-09 | sph_zoom_sandwich.map | cascade | positive | sphmap.c:826-833 | Sandwich with ZoomMap instead of MatrixMap | Inv(SphMap) + ZoomMap + SphMap |
| sphmap-10 | — | cascade | negative | sphmap.c:815-816 | Third Mapping not a non-inverted SphMap | Inv(SphMap) + MatrixMap + ZoomMap |
| sphmap-11 | neg_sph_sandwich_wrong_middle.map | cascade | negative | sphmap.c:825-858 | Middle not ZoomMap or diagonal MatrixMap | Inv(SphMap) + ShiftMap + SphMap |
| sphmap-12 | — | cascade | negative | sphmap.c:828-832 | ZoomMap has zero factor | Inv(SphMap) + ZoomMap(0) + SphMap |
| sphmap-13 | neg_sph_sandwich_full_matrix.map | cascade | negative | sphmap.c:837-855 | MatrixMap not diagonal or null | Inv(SphMap) + FullMatrixMap + SphMap |
| sphmap-14 | neg_sph_sandwich_unequal_diag.map | cascade | negative | sphmap.c:843-847 | MatrixMap diagonal: unequal magnitude | Inv(SphMap) + MatrixMap(diag=[1,2,3]) + SphMap |
| sphmap-15 | — | cascade | negative | sphmap.c:849-850 | MatrixMap first diagonal is zero | Inv(SphMap) + MatrixMap(diag=[0,1,1]) + SphMap |
| sphmap-16 | — | cascade | negative | sphmap.c:921 | Adjusted PolarLong doesn't match third SphMap | PolarLong mismatch after sign adjustment |
| sphmap-17 | — | cascade | negative | sphmap.c:810 | Nominated SphMap not inverted | Forward SphMap + MatrixMap + SphMap |
| sphmap-18 | — | cascade | negative | sphmap.c:811 | Fewer than 3 Mappings remain | Inv(SphMap) at end of list |

## pcdmap.c

PcdMap's MapMerge: (1) Disco=0→UnitMap; (2) inverse-pair cancellation;
(3) swap past ZoomMap/PermMap toward merge target or for local simplification.

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
|---|---|---|---|---|---|---|
| pcdmap-01 | pcd_zero_to_unit.map | focused | positive | pcdmap.c:1410-1420 | PcdMap with Disco=0 replaced by UnitMap | PcdMap(Disco=0) |
| pcdmap-02 | pcd_inverse_cancel.map | focused | positive | pcdmap.c:1460-1499 | PcdMap + Inverse(PcdMap) cancel to UnitMap | CmpMap(PcdMap, Inverse(PcdMap)), Series=1 |
| pcdmap-03 | pcd_unit_series_merge.map | focused | positive | pcdmap.c:1460-1499 | PcdMap + UnitMap neighbour: UnitMap eliminated | CmpMap(PcdMap, UnitMap), Series=1 |
| pcdmap-04 | pcd_zoom_swap_cancel.map | cascade | positive | pcdmap.c:1605-1637 | PcdMap swaps with ZoomMap toward merge target | PcdMap + ZoomMap + Inverse(PcdMap) |
| pcdmap-05 | pcd_perm_swap_cancel.map | cascade | positive | pcdmap.c:1605-1637 | PcdMap swaps with axis-swapping PermMap toward merge target | PcdMap + PermMap(swap) + Inverse(PcdMap) |
| pcdmap-06 | pcd_swap_zoom_simplifies.map | cascade | positive | pcdmap.c:1642-1713 | Swap without target if it simplifies one Mapping | PcdMap + ZoomMap(1) where ZoomMap→UnitMap after swap |
| pcdmap-07 | neg_pcd_zoom_no_target.map | cascade | negative | pcdmap.c:1685-1689 | Speculative swap refused: neither simplifies | PcdMap + ZoomMap(2) with no target |
| pcdmap-08 | neg_pcd_parallel.map | focused | negative | pcdmap.c:1432 | Parallel mode: refused | PcdMaps in parallel |
| pcdmap-09 | neg_pcd_nonpcd_neighbour.map | focused | negative | pcdmap.c:1438-1456 | Neighbour not PcdMap/UnitMap/inverse-PcdMap | PcdMap + ShiftMap in series |
| pcdmap-10 | neg_pcd_nonswappable_between.map | cascade | negative | pcdmap.c:1530-1553 | Intervening Mapping not ZoomMap or PermMap: swap blocked | PcdMap + ShiftMap + Inverse(PcdMap) |
| pcdmap-11 | pcd_search_blocked.map | cascade | negative | pcdmap.c:1547-1550 | Non-swappable class blocks the forward swap search | PcdMap + ZoomMap + MatrixMap (Zoom swappable, Matrix blocks) |
| pcdmap-12 | — | focused | negative | pcdmap.c:1514-1522 | CanSwap false: PermMap doesn't simply swap axes | PcdMap + PermMap(identity) |
| pcdmap-13 | neg_pcd_zoom_no_target.map | cascade | negative | pcdmap.c:1558-1586 | Backward (swaplo) swap search: lower neighbour is a swappable ZoomMap but no merge target found below it | CmpMap(ZoomMap, PcdMap), Series=1 |

## dssmap.c

DssMap's MapMerge: absorbs a neighbouring WinMap into modified pixel
parameters (CNPIX, XPIXELSZ, YPIXELSZ) within its FitsChan. Has many guards
(integer CNPIX, non-zero scale, required keywords present).

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
|---|---|---|---|---|---|---|
| dssmap-01 | — | focused | negative | dssmap.c:1007 | Parallel mode: refused | DssMap in parallel |
| dssmap-02 | — | focused | negative | dssmap.c:1015 | No adjacent mapping at expected index | DssMap alone or at boundary |
| dssmap-03 | dssmap_zoom_no_merge.map | focused | negative | dssmap.c:1016 | Adjacent mapping is not a WinMap | DssMap + ZoomMap in series |
| dssmap-04 | — | focused | negative | dssmap.c:1027-1029 | WinMap scale/shift unusable (AST__BAD or zero scale) | WinMap with zero scale preceding DssMap |
| dssmap-05 | — | focused | negative | dssmap.c:1058-1059 | Computed CNPIX values non-integer beyond tolerance | WinMap producing fractional CNPIX |
| dssmap-06 | — | focused | negative | dssmap.c:1068-1096 | Required FITS keywords missing from DssMap FitsChan | DssMap with incomplete FitsChan |
| dssmap-07 | dssmap_winmap_absorb.map | focused | positive | dssmap.c:1042-1046,1100-1123 | Non-inverted DssMap absorbs preceding WinMap | WinMap + DssMap(Invert=0) with valid integer CNPIX |
| dssmap-08 | dssmap_inv_winmap_absorb.map | focused | positive | dssmap.c:1048-1053,1100-1123 | Inverted DssMap absorbs following WinMap | DssMap(Invert=1) + WinMap with valid integer CNPIX |

## grismmap.c

GrismMap's MapMerge: (1) inverse-pair cancellation; (2) ZoomMap absorption
into wavelength parameters (forward GrismMap+Zoom or Zoom+inverse GrismMap).

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
|---|---|---|---|---|---|---|
| grismmap-01 | grism_inverse_cancel.map | focused | positive | grismmap.c:494-505 | Two GrismMaps with matching attributes and opposite invert cancel to UnitMap | CmpMap(GrismMap, Inverse(GrismMap)), Series=1 |
| grismmap-02 | grism_zoom_merge.map | focused | positive | grismmap.c:515-528,563-588 | Forward GrismMap + ZoomMap: zoom absorbed into wavelength params | CmpMap(GrismMap(fwd), ZoomMap), Series=1 |
| grismmap-03 | grism_zoom_inv_merge.map | focused | positive | grismmap.c:540-556,563-588 | ZoomMap + Inverse(GrismMap): zoom absorbed into wavelength params | CmpMap(ZoomMap, GrismMap(inv)), Series=1 |
| grismmap-04 | neg_grism_parallel.map | focused | negative | grismmap.c:1205 | Parallel mode: refused | GrismMaps in parallel |
| grismmap-05 | — | focused | negative | grismmap.c:1212 | No mergeable neighbour found | Single GrismMap with non-mergeable neighbours |
| grismmap-06 | neg_grism_different_attrs.map | focused | negative | grismmap.c:494-501 | Two GrismMaps with differing attributes (NR, NRP, WaveR, etc.) | GrismMaps with different parameters |
| grismmap-07 | — | focused | negative | grismmap.c:499 | GrismM values equal: merge blocked (possible historical bug) | Two GrismMaps with same M |
| grismmap-08 | neg_grism_same_direction.map | focused | negative | grismmap.c:505 | Same invert flag (same direction) | Two forward GrismMaps in series |
| grismmap-09 | neg_grism_inv_then_zoom.map | focused | negative | grismmap.c:515 | First GrismMap is inverted: ZoomMap merge N/A | Inverse(GrismMap) + ZoomMap |
| grismmap-10 | neg_grism_fwd_then_nonzoom.map | focused | negative | grismmap.c:525 | First is forward GrismMap but second not ZoomMap | GrismMap(fwd) + ShiftMap |
| grismmap-11 | neg_grism_zoom_before_fwd.map | focused | negative | grismmap.c:540 | Second GrismMap not inverted: pre-ZoomMap merge N/A | ZoomMap + GrismMap(fwd) |
| grismmap-12 | neg_grism_nonzoom_before_inv.map | focused | negative | grismmap.c:553 | Second is inverted GrismMap but first not ZoomMap | ShiftMap + Inverse(GrismMap) |
| grismmap-13 | — | focused | negative | grismmap.c:563 | ZoomMap zoom factor is zero | GrismMap(fwd) + ZoomMap(0) |

## xphmap.c

XphMap's MapMerge: inverse-pair cancellation only (same order, opposite
direction via astEqual).

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
|---|---|---|---|---|---|---|
| xphmap-01 | — | focused | negative | xphmap.c:524 | Parallel mode: refused | XphMap in parallel |
| xphmap-02 | — | focused | negative | xphmap.c:530-544 | No adjacent XphMap found | XphMap flanked by non-XphMaps |
| xphmap-03 | — | focused | negative | xphmap.c:559 | Adjacent XphMap not equal-and-opposite | XphMaps with different order or same direction |
| xphmap-04 | — | focused | positive | xphmap.c:567-596 | Inverse-pair cancellation → UnitMap | XphMap + Inverse(XphMap), same order |

---

## Region-as-Mapping classes

All four Region classes below share structurally identical MapMerge logic:
(1) self-simplification via astSimplify; (2) parallel merge with adjacent
Region via class-specific MergeXxx helper.

### box.c

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
|---|---|---|---|---|---|---|
| box-01 | — | focused | positive | box.c:1583-1593 | Self-simplification succeeds (astSimplify returns different pointer) | Box with non-trivial base-to-current FrameSet |
| box-02 | — | focused | negative | box.c:1599,1650-1652 | No self-simplification and series mode | Already-simple Box in series |
| box-03 | box_parallel_merge.map | cascade | positive | box.c:1603-1609,1624-1648 | Parallel merge with lower Region via MergeBox | Box in parallel with compatible Region (lower) |
| box-04 | box_parallel_merge.map | cascade | positive | box.c:1614-1621,1624-1648 | Parallel merge with upper Region via MergeBox | Box in parallel with compatible Region (upper) |
| box-05 | — | focused | negative | box.c:1599-1621 | Parallel but no Region neighbour or MergeBox returns NULL | Box in parallel with incompatible Region |
| box-06 | neg_box_asymmetric_2d.map | focused | negative | box.c:1583-1652 | Standalone 2-D Box with asymmetric axis intervals does not self-simplify | Lone Box, differing intervals per axis, no neighbour |

### interval.c

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
|---|---|---|---|---|---|---|
| interval-01 | — | focused | positive | interval.c:1149-1159 | Self-simplification succeeds | Interval with non-trivial FrameSet |
| interval-02 | — | focused | negative | interval.c:1165,1216-1218 | No self-simplification, series mode | Simple Interval in series |
| interval-03 | interval_parallel_merge.map | cascade | positive | interval.c:1169-1175,1190-1213 | Parallel merge with lower Region | Interval + compatible Region (lower) |
| interval-04 | interval_parallel_merge.map | cascade | positive | interval.c:1180-1187,1190-1213 | Parallel merge with upper Region | Interval + compatible Region (upper) |
| interval-05 | — | focused | negative | interval.c:1165-1187 | Parallel but no compatible Region | Interval + non-Region in parallel |

### nullregion.c

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
|---|---|---|---|---|---|---|
| nullregion-01 | — | focused | positive | nullregion.c:496-506 | Self-simplification succeeds | NullRegion with non-trivial FrameSet |
| nullregion-02 | — | focused | negative | nullregion.c:512,563-565 | No self-simplification, series mode | Simple NullRegion in series |
| nullregion-03 | — | cascade | positive | nullregion.c:516-523,537-560 | Parallel merge with lower Region | NullRegion + compatible Region (lower) |
| nullregion-04 | — | cascade | positive | nullregion.c:527-533,537-560 | Parallel merge with upper Region | NullRegion + compatible Region (upper) |
| nullregion-05 | — | focused | negative | nullregion.c:512-534 | Parallel but no compatible Region | NullRegion + non-Region in parallel |

### pointlist.c

| ID | Fixture | Type | Polarity | Lines | Description | Trigger |
|---|---|---|---|---|---|---|
| pointlist-01 | — | focused | positive | pointlist.c:1130-1140 | Self-simplification succeeds | PointList with non-trivial FrameSet |
| pointlist-02 | — | focused | negative | pointlist.c:1146,1197-1199 | No self-simplification, series mode | Simple PointList in series |
| pointlist-03 | — | cascade | positive | pointlist.c:1150-1156,1171-1194 | Parallel merge with lower Region | PointList + compatible Region (lower) |
| pointlist-04 | — | cascade | positive | pointlist.c:1161-1167,1171-1194 | Parallel merge with upper Region | PointList + compatible Region (upper) |
| pointlist-05 | — | focused | negative | pointlist.c:1146-1168 | Parallel but no compatible Region | PointList + non-Region in parallel |
