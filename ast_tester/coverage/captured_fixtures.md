# Captured differential-coverage fixtures

These `cap_*` fixtures were generated automatically, not hand-authored.
Each is a real top-level Mapping that the wider AST test suite passes to
`astSimplify`, captured by temporary instrumentation, replayed as a
self-contained fixture. Together they close 293 of the 323 differential
coverage gaps (branches the full suite exercised but the simplify fixtures
did not). See `README.md` for the capture method and
`simplify_coverage_gaps.md` for the live gap ledger.

Unlike the hand-authored rule fixtures, each capture exercises several
branches across one or more classes, so they are catalogued here by the
class that dominates their coverage rather than as single-branch inventory
rows.

| Fixture | Dominant class | Target branches covered |
| --- | --- | --- |
| `cap_specmap_01` | specmap.c | 36 |
| `cap_xphmap_01` | xphmap.c | 30 |
| `cap_timemap_02` | timemap.c | 29 |
| `cap_specmap_02` | specmap.c | 24 |
| `cap_wcsmap_01` | wcsmap.c | 21 |
| `cap_interval_01` | interval.c | 20 |
| `cap_pointlist_01` | pointlist.c | 20 |
| `cap_cmpmap_01` | cmpmap.c | 19 |
| `cap_cmpmap_04` | cmpmap.c | 17 |
| `cap_permmap_01` | permmap.c | 17 |
| `cap_pointlist_02` | pointlist.c | 17 |
| `cap_slamap_02` | slamap.c | 17 |
| `cap_matrixmap_01` | matrixmap.c | 16 |
| `cap_specmap_03` | specmap.c | 14 |
| `cap_winmap_01` | winmap.c | 13 |
| `cap_xphmap_02` | xphmap.c | 12 |
| `cap_zoommap_02` | zoommap.c | 11 |
| `cap_slamap_01` | slamap.c | 10 |
| `cap_zoommap_01` | zoommap.c | 10 |
| `cap_specmap_09` | specmap.c | 9 |
| `cap_timemap_01` | timemap.c | 9 |
| `cap_winmap_02` | winmap.c | 9 |
| `cap_matrixmap_02` | matrixmap.c | 8 |
| `cap_slamap_03` | slamap.c | 8 |
| `cap_slamap_05` | slamap.c | 8 |
| `cap_timemap_05` | timemap.c | 8 |
| `cap_cmpmap_02` | cmpmap.c | 7 |
| `cap_cmpmap_03` | cmpmap.c | 7 |
| `cap_matrixmap_03` | matrixmap.c | 7 |
| `cap_specmap_04` | specmap.c | 6 |
| `cap_zoommap_04` | zoommap.c | 6 |
| `cap_slamap_06` | slamap.c | 5 |
| `cap_specmap_06` | specmap.c | 5 |
| `cap_zoommap_03` | zoommap.c | 5 |
| `cap_specmap_05` | specmap.c | 4 |
| `cap_specmap_08` | specmap.c | 4 |
| `cap_specmap_11` | specmap.c | 4 |
| `cap_specmap_12` | specmap.c | 4 |
| `cap_specmap_14` | specmap.c | 4 |
| `cap_timemap_03` | timemap.c | 4 |
| `cap_timemap_04` | timemap.c | 4 |
| `cap_winmap_03` | winmap.c | 4 |
| `cap_slamap_04` | slamap.c | 3 |
| `cap_slamap_08` | slamap.c | 3 |
| `cap_specmap_07` | specmap.c | 3 |
| `cap_timemap_06` | timemap.c | 3 |
| `cap_timemap_07` | timemap.c | 3 |
| `cap_unitmap_01` | unitmap.c | 3 |
| `cap_winmap_04` | winmap.c | 3 |
| `cap_lutmap_01` | lutmap.c | 2 |
| `cap_matrixmap_04` | matrixmap.c | 2 |
| `cap_permmap_02` | permmap.c | 2 |
| `cap_slamap_07` | slamap.c | 2 |
| `cap_specmap_10` | specmap.c | 2 |
| `cap_specmap_13` | specmap.c | 2 |
| `cap_unitnormmap_01` | unitnormmap.c | 2 |
| `cap_mathmap_01` | mathmap.c | 1 |
| `cap_slamap_09` | slamap.c | 1 |
| `cap_sphmap_01` | sphmap.c | 1 |
| `cap_timemap_08` | timemap.c | 1 |

**Total:** 60 fixtures.
