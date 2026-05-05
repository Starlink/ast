# astSimplify Coverage Matrix

This file is the source-level checklist for simplification coverage. A
simplification is counted as covered only when a test intentionally asserts
the simplified shape or semantics for that rule. Incidental use of
`astSimplify` from FITS, WCS conversion, or region tests is useful regression
coverage, but does not close a rule here unless the test is explicitly tied
to that rule.

## Test Shape

- Prefer focused `.map` / `.simp` fixtures for public `astSimplify` behaviour.
  They are easy to inspect and are already supported by `simplify_tests.txt`.
- Use one primary simplification rule per fixture. Add larger scenario
  fixtures only when the rule being tested is a cascade or reordering rule.
- Use C tests when the expected result is easier to assert by API, requires
  generated data, or needs numerical transformation checks instead of a stable
  serialized form.
- Keep historical end-to-end fixtures such as `brad`, `lsst1`, and `rigby` as
  scenario regressions, not as the only proof for individual rules.

## Rule Inventory

| Owner | Simplification rule or family | Targeted coverage |
| --- | --- | --- |
| Mapping | Single mapping repeatedly delegates to `astMapMerge`; restricted simplify honours `AllowSimplify` | `testresimp.c` |
| CmpMap | Flattens nested CmpMaps, decomposes series/parallel lists, detects repeated simplification loops | `cmpmap_nested_parallel_flatten.map`; scenario fixtures |
| CmpMap | Merges adjacent CmpMaps by reshaping series/parallel combinations | `cmpmap_parallel_series_components.map` |
| CmpMap | Swaps compatible PermMap and parallel CmpMap combinations | `cmpmap_perm_parallel_swap.map` |
| UnitMap | Removes UnitMap from series combinations | `unit_series_elision.map` |
| UnitMap | Merges adjacent UnitMaps in parallel combinations | `unit_parallel_merge.map` |
| ZoomMap | Merges adjacent ZoomMaps in series | `zoom_series_merge.map` |
| ZoomMap | Cancels net unit zoom to UnitMap | `zoom_series_cancel.map` |
| ZoomMap | Merges parallel ZoomMaps/UnitMaps to ZoomMap, MatrixMap, or UnitMap | `zoom_parallel_to_matrix.map` |
| ShiftMap | Normalizes inverted ShiftMap and delegates merge logic through WinMap | `shift_invert_normalize.map`; WinMap merge fixtures |
| WinMap | Converts zero-shift WinMap to diagonal MatrixMap | `win_to_matrix.map` |
| WinMap | Converts unit-scale WinMap to ShiftMap when no better merge exists | `win_to_shift.map` |
| WinMap | Merges with WinMap, ZoomMap, ShiftMap, UnitMap, or diagonal MatrixMap in series and parallel | `win_shift_series_merge.map`, `win_zoom_series_merge.map`, `win_matrix_series_merge.map`, `win_parallel_merge.map`; additional parallel cases planned |
| WinMap | Swaps with compatible MatrixMap, PermMap, or WcsMap to enable later simplification | `win_perm_swap_merge.map`; MatrixMap/WcsMap swap fixtures planned; `rigby` scenario |
| MatrixMap | Converts unit MatrixMap to UnitMap | `matrix_unit_to_unit.map` |
| MatrixMap | Converts equal diagonal MatrixMap to ZoomMap | `matrix_diagonal_to_zoom.map` |
| MatrixMap | Compresses full diagonal matrix storage to diagonal MatrixMap | `matrix_full_to_diagonal.map` |
| MatrixMap | Merges with MatrixMap, ZoomMap, compatible PermMap, UnitMap, or diagonal WinMap | `matrix_matrix_series_merge.map`, `matrix_zoom_series_merge.map`, `matrix_perm_series_merge.map`, `win_matrix_series_merge.map` |
| MatrixMap | Swaps with compatible WinMap or PermMap to enable later simplification | `matrix_perm_series_merge.map`; WinMap swap fixture planned; `rigby` scenario |
| PermMap | Merges adjacent PermMaps/UnitMaps in series and parallel | `perm_series_merge.map`, `perm_parallel_merge.map` |
| LutMap | Converts linear LUT to WinMap | `lut_linear_to_win.map` |
| LutMap | Cancels a LutMap with its inverse | `lut_inverse_cancel.map` |
| PolyMap | Combines duplicate coefficient terms and drops zero terms | `poly_duplicate_terms.map`, `testpolymap.c` |
| PolyMap | Converts linear polynomial to simpler ShiftMap/MatrixMap/WinMap forms | `poly_duplicate_terms.map`, `testpolymap.c` |
| PolyMap | Cancels a PolyMap with its inverse | `poly_inverse_cancel.map` |
| NormMap | Cancels with inverse, collapses repeated NormMaps, simplifies basic Frame to UnitMap | `testnormmap.c` |
| UnitNormMap | Cancels inverse pairs and merges with ShiftMap/WinMap combinations | `testunitnormmap.c` |
| SwitchMap | Normalizes invert flag and cancels inverse pairs | `testswitchmap.c` |
| TimeMap | Collapses inverse time conversions and adjacent compatible conversions | `time_inverse_cancel.map`, `testtime.c`; additional adjacent conversion fixtures planned |
| SlaMap | Collapses adjacent inverse/compatible sky conversions | `sla_inverse_cancel.map`; additional compatible-conversion fixtures planned |
| SpecMap | Collapses adjacent inverse/compatible spectral conversions | `spec_inverse_cancel.map`; additional compatible-conversion fixtures planned |
| WcsMap | Cancels inverse projections and swaps with compatible PermMap/WinMap combinations | `wcs_inverse_cancel.map`; swap fixtures planned; WCS scenario tests |
| SphMap | Cancels spherical/cartesian inverse pairs | `sph_inverse_cancel.map`, `sph_matrix_sandwich.map` |
| PcdMap, DssMap, GrismMap | Cancels inverse pairs and class-specific adjacent simplifications | `pcd_zero_to_unit.map`, `pcd_inverse_cancel.map`, `grism_inverse_cancel.map`, `grism_zoom_merge.map`; DssMap fixture planned |
| MathMap, IntraMap, RateMap, TranMap | Class-specific inverse-pair and no-op simplifications | Planned C/API tests |
| SplineMap | Class-specific inverse-pair and no-op simplifications | Planned C/API tests |
| SelectorMap | Simplifies encapsulated regions and class-specific no-op cases | Planned C/API tests |
| Box, Interval, NullRegion, PointList | Region-as-Mapping `MapMerge` simplifications | Existing region tests; planned targeted fixtures |
| Region, Circle, Ellipse, Polygon, Prism, CmpRegion, Stc | Region `Simplify` overrides | Existing region/STC tests; separate region coverage matrix recommended |

## Adding A Rule Test

1. Add a focused `.map` input and corresponding `.simp` reference.
2. Add one line to `simplify_tests.txt` with a rule-oriented name.
3. If strict serialization varies by platform, set `skip_string_compare` to
   `yes`; the semantic `ast_astequal` check will still run.
4. Update the row above from `Planned` or incidental coverage to the fixture
   or C test that proves the rule.
