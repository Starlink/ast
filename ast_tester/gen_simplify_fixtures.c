/*
 * gen_simplify_fixtures.c
 *
 * Generates .map and .simp fixture files for the astSimplify branch
 * coverage campaign. Each fixture is created using the public C API,
 * serialized with astToString (.map), simplified with astSimplify,
 * then serialized again (.simp).
 *
 * Build:
 *   cmake --build build-dev --target gen_simplify_fixtures
 * Run:
 *   ./build-dev/ast_tester/gen_simplify_fixtures
 *
 * Output files go to the current working directory (run from ast_tester/).
 */

#include "ast.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void write_fixture(const char *dir, const char *name, AstMapping *map) {
    char path_map[512], path_simp[512];
    AstChannel *chan;
    AstMapping *simp;

    snprintf(path_map, sizeof(path_map), "%s/%s.map", dir, name);
    snprintf(path_simp, sizeof(path_simp), "%s/%s.simp", dir, name);

    chan = astChannel(NULL, NULL, "SinkFile=%s", path_map);
    if (astWrite(chan, map) != 1) {
        fprintf(stderr, "ERROR: astWrite failed for %s.map\n", name);
        chan = astAnnul(chan);
        return;
    }
    chan = astAnnul(chan);

    simp = astSimplify(map);
    chan = astChannel(NULL, NULL, "SinkFile=%s", path_simp);
    if (astWrite(chan, simp) != 1) {
        fprintf(stderr, "ERROR: astWrite failed for %s.simp\n", name);
        chan = astAnnul(chan);
        simp = astAnnul(simp);
        return;
    }
    chan = astAnnul(chan);
    simp = astAnnul(simp);

    printf("  %s\n", name);
}

/* ===== ZoomMap fixtures ===== */

static void gen_zoom_fixtures(const char *dir) {
    printf("ZoomMap fixtures:\n");

    /* zoommap-03: single inverted ZoomMap normalizes */
    {
        AstZoomMap *zm = astZoomMap(2, 2.0, "");
        astInvert(zm);
        write_fixture(dir, "zoom_invert_normalize", (AstMapping*)zm);
        zm = astAnnul(zm);
    }

    /* zoommap-05: parallel all-unit → UnitMap */
    {
        AstUnitMap *um = astUnitMap(2, "");
        AstZoomMap *zm = astZoomMap(3, 1.0, "");
        AstCmpMap *cm = astCmpMap(um, zm, 0, "");
        write_fixture(dir, "zoom_parallel_all_unit", (AstMapping*)cm);
        cm = astAnnul(cm); um = astAnnul(um); zm = astAnnul(zm);
    }

    /* zoommap-06: parallel same factor → single ZoomMap */
    {
        AstZoomMap *z1 = astZoomMap(1, 2.0, "");
        AstZoomMap *z2 = astZoomMap(2, 2.0, "");
        AstCmpMap *cm = astCmpMap(z1, z2, 0, "");
        write_fixture(dir, "zoom_parallel_same_factor", (AstMapping*)cm);
        cm = astAnnul(cm); z1 = astAnnul(z1); z2 = astAnnul(z2);
    }

    /* zoommap-09: absorb into previous MatrixMap */
    {
        double diag[] = {3.0, 5.0};
        AstMatrixMap *mm = astMatrixMap(2, 2, 1, diag, "");
        AstZoomMap *zm = astZoomMap(2, 2.0, "");
        AstCmpMap *cm = astCmpMap(mm, zm, 1, "");
        write_fixture(dir, "zoom_absorb_prev_matrix", (AstMapping*)cm);
        cm = astAnnul(cm); mm = astAnnul(mm); zm = astAnnul(zm);
    }

    /* zoommap-10: absorb into previous WinMap */
    {
        double ina[] = {0, 0}, inb[] = {1, 1};
        double outa[] = {1, 2}, outb[] = {4, 6};
        AstWinMap *wm = astWinMap(2, ina, inb, outa, outb, "");
        AstZoomMap *zm = astZoomMap(2, 2.0, "");
        AstCmpMap *cm = astCmpMap(wm, zm, 1, "");
        write_fixture(dir, "zoom_absorb_prev_win", (AstMapping*)cm);
        cm = astAnnul(cm); wm = astAnnul(wm); zm = astAnnul(zm);
    }

    /* zoommap-11: absorb into previous ShiftMap */
    {
        double shifts[] = {3.0, 5.0};
        AstShiftMap *sm = astShiftMap(2, shifts, "");
        AstZoomMap *zm = astZoomMap(2, 2.0, "");
        AstCmpMap *cm = astCmpMap(sm, zm, 1, "");
        write_fixture(dir, "zoom_absorb_prev_shift", (AstMapping*)cm);
        cm = astAnnul(cm); sm = astAnnul(sm); zm = astAnnul(zm);
    }

    /* zoommap-12: absorb into next MatrixMap */
    {
        double diag[] = {3.0, 5.0};
        AstZoomMap *zm = astZoomMap(2, 2.0, "");
        AstMatrixMap *mm = astMatrixMap(2, 2, 1, diag, "");
        AstCmpMap *cm = astCmpMap(zm, mm, 1, "");
        write_fixture(dir, "zoom_absorb_next_matrix", (AstMapping*)cm);
        cm = astAnnul(cm); zm = astAnnul(zm); mm = astAnnul(mm);
    }

    /* zoommap-13: absorb into next WinMap */
    {
        double ina[] = {0, 0}, inb[] = {1, 1};
        double outa[] = {1, 2}, outb[] = {4, 6};
        AstZoomMap *zm = astZoomMap(2, 2.0, "");
        AstWinMap *wm = astWinMap(2, ina, inb, outa, outb, "");
        AstCmpMap *cm = astCmpMap(zm, wm, 1, "");
        write_fixture(dir, "zoom_absorb_next_win", (AstMapping*)cm);
        cm = astAnnul(cm); zm = astAnnul(zm); wm = astAnnul(wm);
    }

    /* zoommap-14: absorb into next ShiftMap */
    {
        double shifts[] = {3.0, 5.0};
        AstZoomMap *zm = astZoomMap(2, 2.0, "");
        AstShiftMap *sm = astShiftMap(2, shifts, "");
        AstCmpMap *cm = astCmpMap(zm, sm, 1, "");
        write_fixture(dir, "zoom_absorb_next_shift", (AstMapping*)cm);
        cm = astAnnul(cm); zm = astAnnul(zm); sm = astAnnul(sm);
    }
}

/* ===== WinMap fixtures ===== */

static void gen_win_fixtures(const char *dir) {
    printf("WinMap fixtures:\n");

    /* winmap-05: WinMap + WinMap series */
    {
        double ina[] = {0}, inb[] = {1}, outa[] = {3}, outb[] = {5};
        double ina2[] = {0}, inb2[] = {1}, outa2[] = {5}, outb2[] = {9};
        AstWinMap *w1 = astWinMap(1, ina, inb, outa, outb, "");
        AstWinMap *w2 = astWinMap(1, ina2, inb2, outa2, outb2, "");
        AstCmpMap *cm = astCmpMap(w1, w2, 1, "");
        write_fixture(dir, "win_win_series_merge", (AstMapping*)cm);
        cm = astAnnul(cm); w1 = astAnnul(w1); w2 = astAnnul(w2);
    }

    /* winmap-07: ZoomMap + WinMap series (reverse order) */
    {
        double ina[] = {0}, inb[] = {1}, outa[] = {5}, outb[] = {7};
        AstZoomMap *zm = astZoomMap(1, 3.0, "");
        AstWinMap *wm = astWinMap(1, ina, inb, outa, outb, "");
        AstCmpMap *cm = astCmpMap(zm, wm, 1, "");
        write_fixture(dir, "win_zoom_series_merge_rev", (AstMapping*)cm);
        cm = astAnnul(cm); zm = astAnnul(zm); wm = astAnnul(wm);
    }

    /* winmap-09: ShiftMap + WinMap series (reverse order) */
    {
        double shifts[] = {3.0};
        double ina[] = {0}, inb[] = {1}, outa[] = {5}, outb[] = {7};
        AstShiftMap *sm = astShiftMap(1, shifts, "");
        AstWinMap *wm = astWinMap(1, ina, inb, outa, outb, "");
        AstCmpMap *cm = astCmpMap(sm, wm, 1, "");
        write_fixture(dir, "win_shift_series_merge_rev", (AstMapping*)cm);
        cm = astAnnul(cm); sm = astAnnul(sm); wm = astAnnul(wm);
    }

    /* winmap-11: DiagMatrixMap + WinMap series (reverse order) */
    {
        double diag[] = {3.0, 5.0};
        double ina[] = {0, 0}, inb[] = {1, 1}, outa[] = {1, 2}, outb[] = {5, 8};
        AstMatrixMap *mm = astMatrixMap(2, 2, 1, diag, "");
        AstWinMap *wm = astWinMap(2, ina, inb, outa, outb, "");
        AstCmpMap *cm = astCmpMap(mm, wm, 1, "");
        write_fixture(dir, "win_matrix_series_merge_rev", (AstMapping*)cm);
        cm = astAnnul(cm); mm = astAnnul(mm); wm = astAnnul(wm);
    }

    /* winmap-12: WinMap + UnitMap series */
    {
        double ina[] = {0, 0}, inb[] = {1, 1}, outa[] = {3, 5}, outb[] = {5, 9};
        AstWinMap *wm = astWinMap(2, ina, inb, outa, outb, "");
        AstUnitMap *um = astUnitMap(2, "");
        AstCmpMap *cm = astCmpMap(wm, um, 1, "");
        write_fixture(dir, "win_unit_series_merge", (AstMapping*)cm);
        cm = astAnnul(cm); wm = astAnnul(wm); um = astAnnul(um);
    }

    /* winmap-29: WinMap + WinMap parallel */
    {
        double ina[] = {0}, inb[] = {1}, outa[] = {3}, outb[] = {5};
        double ina2[] = {0}, inb2[] = {1}, outa2[] = {7}, outb2[] = {11};
        AstWinMap *w1 = astWinMap(1, ina, inb, outa, outb, "");
        AstWinMap *w2 = astWinMap(1, ina2, inb2, outa2, outb2, "");
        AstCmpMap *cm = astCmpMap(w1, w2, 0, "");
        write_fixture(dir, "win_win_parallel_merge", (AstMapping*)cm);
        cm = astAnnul(cm); w1 = astAnnul(w1); w2 = astAnnul(w2);
    }

    /* winmap-30: WinMap + ZoomMap parallel */
    {
        double ina[] = {0}, inb[] = {1}, outa[] = {3}, outb[] = {5};
        AstWinMap *wm = astWinMap(1, ina, inb, outa, outb, "");
        AstZoomMap *zm = astZoomMap(1, 4.0, "");
        AstCmpMap *cm = astCmpMap(wm, zm, 0, "");
        write_fixture(dir, "win_zoom_parallel_merge", (AstMapping*)cm);
        cm = astAnnul(cm); wm = astAnnul(wm); zm = astAnnul(zm);
    }

    /* winmap-34: WinMap + DiagMatrixMap parallel */
    {
        double ina[] = {0}, inb[] = {1}, outa[] = {3}, outb[] = {5};
        double diag[] = {7.0};
        AstWinMap *wm = astWinMap(1, ina, inb, outa, outb, "");
        AstMatrixMap *mm = astMatrixMap(1, 1, 1, diag, "");
        AstCmpMap *cm = astCmpMap(wm, mm, 0, "");
        write_fixture(dir, "win_diagmatrix_parallel_merge", (AstMapping*)cm);
        cm = astAnnul(cm); wm = astAnnul(wm); mm = astAnnul(mm);
    }

    /* winmap-36: WinMap + UnitMap parallel */
    {
        double ina[] = {0}, inb[] = {1}, outa[] = {3}, outb[] = {5};
        AstWinMap *wm = astWinMap(1, ina, inb, outa, outb, "");
        AstUnitMap *um = astUnitMap(1, "");
        AstCmpMap *cm = astCmpMap(wm, um, 0, "");
        write_fixture(dir, "win_unit_parallel_merge", (AstMapping*)cm);
        cm = astAnnul(cm); wm = astAnnul(wm); um = astAnnul(um);
    }
}

/* ===== UnitMap fixtures ===== */

static void gen_unit_fixtures(const char *dir) {
    printf("UnitMap fixtures:\n");

    /* unitmap-01: single inverted UnitMap → clears invert */
    {
        AstUnitMap *um = astUnitMap(2, "");
        astInvert(um);
        write_fixture(dir, "unit_invert_clear", (AstMapping*)um);
        um = astAnnul(um);
    }
}

/* ===== MatrixMap fixtures ===== */

static void gen_matrix_fixtures(const char *dir) {
    printf("MatrixMap fixtures:\n");

    /* matrixmap-11: DiagMatrixMap + WinMap series → merged via MatWin2 */
    {
        double diag[] = {2.0, 3.0};
        double ina[] = {0, 0}, inb[] = {1, 1}, outa[] = {1, 2}, outb[] = {5, 8};
        AstMatrixMap *mm = astMatrixMap(2, 2, 1, diag, "");
        AstWinMap *wm = astWinMap(2, ina, inb, outa, outb, "");
        AstCmpMap *cm = astCmpMap(mm, wm, 1, "");
        write_fixture(dir, "matrix_diagwin_series_merge", (AstMapping*)cm);
        cm = astAnnul(cm); mm = astAnnul(mm); wm = astAnnul(wm);
    }

    /* matrixmap-12: MatrixMap + UnitMap series → UnitMap eliminated */
    {
        double diag[] = {2.0, 3.0};
        AstMatrixMap *mm = astMatrixMap(2, 2, 1, diag, "");
        AstUnitMap *um = astUnitMap(2, "");
        AstCmpMap *cm = astCmpMap(mm, um, 1, "");
        write_fixture(dir, "matrix_unit_series_merge", (AstMapping*)cm);
        cm = astAnnul(cm); mm = astAnnul(mm); um = astAnnul(um);
    }
}

/* ===== CmpMap fixtures ===== */

static void gen_cmpmap_fixtures(const char *dir) {
    printf("CmpMap fixtures:\n");

    /* cmpmap-01: CmpMap self-simplifies (internal inverse pair cancels) */
    {
        AstZoomMap *z1 = astZoomMap(2, 2.0, "");
        AstZoomMap *z2 = astZoomMap(2, 2.0, "");
        astInvert(z2);
        AstCmpMap *inner = astCmpMap(z1, z2, 1, "");
        double shifts[] = {1.0, 2.0};
        AstShiftMap *sm = astShiftMap(2, shifts, "");
        AstCmpMap *outer = astCmpMap(inner, sm, 1, "");
        write_fixture(dir, "cmpmap_self_simplify", (AstMapping*)outer);
        outer = astAnnul(outer); inner = astAnnul(inner);
        z1 = astAnnul(z1); z2 = astAnnul(z2); sm = astAnnul(sm);
    }
}

/* ===== TranMap fixtures ===== */

static void gen_tranmap_fixtures(const char *dir) {
    printf("TranMap fixtures:\n");

    /* tranmap-01: inverted TranMap → normalizes */
    {
        AstZoomMap *z1 = astZoomMap(1, 2.0, "");
        AstZoomMap *z2 = astZoomMap(1, 0.5, "");
        AstTranMap *tm = astTranMap(z1, z2, "");
        astInvert(tm);
        write_fixture(dir, "tranmap_invert_normalize", (AstMapping*)tm);
        tm = astAnnul(tm); z1 = astAnnul(z1); z2 = astAnnul(z2);
    }

    /* tranmap-02: TranMap with simplifiable components */
    {
        AstZoomMap *za1 = astZoomMap(1, 2.0, "");
        AstZoomMap *za2 = astZoomMap(1, 3.0, "");
        AstCmpMap *fwd = astCmpMap(za1, za2, 1, "");
        AstZoomMap *zb1 = astZoomMap(1, 0.5, "");
        AstZoomMap *zb2 = astZoomMap(1, 2.0, "");
        AstCmpMap *inv = astCmpMap(zb1, zb2, 1, "");
        AstTranMap *tm = astTranMap(fwd, inv, "");
        write_fixture(dir, "tranmap_component_simplify", (AstMapping*)tm);
        tm = astAnnul(tm); fwd = astAnnul(fwd); inv = astAnnul(inv);
        za1 = astAnnul(za1); za2 = astAnnul(za2);
        zb1 = astAnnul(zb1); zb2 = astAnnul(zb2);
    }
}

/* ===== RateMap fixtures ===== */

static void gen_ratemap_fixtures(const char *dir) {
    printf("RateMap fixtures:\n");

    /* ratemap-01: RateMap with simplifiable interior */
    {
        AstZoomMap *z1 = astZoomMap(2, 2.0, "");
        AstZoomMap *z2 = astZoomMap(2, 3.0, "");
        AstCmpMap *inner = astCmpMap(z1, z2, 1, "");
        AstRateMap *rm = astRateMap(inner, 1, 1, "");
        write_fixture(dir, "ratemap_simplify_interior", (AstMapping*)rm);
        rm = astAnnul(rm); inner = astAnnul(inner);
        z1 = astAnnul(z1); z2 = astAnnul(z2);
    }
}

/* ===== SlaMap fixtures ===== */

static void gen_slamap_fixtures(const char *dir) {
    printf("SlaMap fixtures:\n");

    /* slamap-04: single inverted SlaMap normalizes */
    {
        AstSlaMap *sm = astSlaMap(0, "");
        astSlaAdd(sm, "EQGAL", 0, NULL);
        astInvert(sm);
        write_fixture(dir, "sla_invert_normalize", (AstMapping*)sm);
        sm = astAnnul(sm);
    }

    /* slamap-02: partial cancellation (3 steps, middle pair cancels) */
    {
        double beq[] = {1950.0};
        AstSlaMap *sm = astSlaMap(0, "");
        astSlaAdd(sm, "ADDET", 1, beq);
        astSlaAdd(sm, "SUBET", 1, beq);
        astSlaAdd(sm, "EQGAL", 0, NULL);
        write_fixture(dir, "sla_partial_cancel", (AstMapping*)sm);
        sm = astAnnul(sm);
    }

    /* slamap-03: multi-map merge without step reduction */
    {
        double beq[] = {1950.0};
        AstSlaMap *s1 = astSlaMap(0, "");
        astSlaAdd(s1, "ADDET", 1, beq);
        AstSlaMap *s2 = astSlaMap(0, "");
        astSlaAdd(s2, "EQGAL", 0, NULL);
        AstCmpMap *cm = astCmpMap(s1, s2, 1, "");
        write_fixture(dir, "sla_merge_no_cancel", (AstMapping*)cm);
        cm = astAnnul(cm); s1 = astAnnul(s1); s2 = astAnnul(s2);
    }

    /* slamap-08: 0-arg supergalactic */
    {
        AstSlaMap *sm = astSlaMap(0, "");
        astSlaAdd(sm, "GALSUP", 0, NULL);
        astSlaAdd(sm, "SUPGAL", 0, NULL);
        write_fixture(dir, "sla_supergalactic_cancel", (AstMapping*)sm);
        sm = astAnnul(sm);
    }

    /* slamap-09: 0-arg J2000 dynamical */
    {
        AstSlaMap *sm = astSlaMap(0, "");
        astSlaAdd(sm, "J2000H", 0, NULL);
        astSlaAdd(sm, "HJ2000", 0, NULL);
        write_fixture(dir, "sla_j2000_cancel", (AstMapping*)sm);
        sm = astAnnul(sm);
    }

    /* slamap-10: 1-arg E-terms */
    {
        double beq[] = {1950.0};
        AstSlaMap *sm = astSlaMap(0, "");
        astSlaAdd(sm, "ADDET", 1, beq);
        astSlaAdd(sm, "SUBET", 1, beq);
        write_fixture(dir, "sla_eterms_cancel", (AstMapping*)sm);
        sm = astAnnul(sm);
    }

    /* slamap-11: 1-arg FK4/FK5 */
    {
        double beq[] = {1950.0};
        AstSlaMap *sm = astSlaMap(0, "");
        astSlaAdd(sm, "FK45Z", 1, beq);
        astSlaAdd(sm, "FK54Z", 1, beq);
        write_fixture(dir, "sla_fk45_cancel", (AstMapping*)sm);
        sm = astAnnul(sm);
    }

    /* slamap-12: 1-arg ICRS/FK5 */
    {
        double args[] = {2000.0};
        AstSlaMap *sm = astSlaMap(0, "");
        astSlaAdd(sm, "HFK5Z", 1, args);
        astSlaAdd(sm, "FK5HZ", 1, args);
        write_fixture(dir, "sla_icrs_cancel", (AstMapping*)sm);
        sm = astAnnul(sm);
    }

    /* slamap-13: 1-arg ecliptic */
    {
        double args[] = {2000.0};
        AstSlaMap *sm = astSlaMap(0, "");
        astSlaAdd(sm, "EQECL", 1, args);
        astSlaAdd(sm, "ECLEQ", 1, args);
        write_fixture(dir, "sla_ecliptic_cancel", (AstMapping*)sm);
        sm = astAnnul(sm);
    }

    /* slamap-16: 2-arg geocentric (cross-matched AMP+MAP) */
    {
        double amp_args[] = {51544.0, 2000.0};
        double map_args[] = {2000.0, 51544.0};
        AstSlaMap *sm = astSlaMap(0, "");
        astSlaAdd(sm, "AMP", 2, amp_args);
        astSlaAdd(sm, "MAP", 2, map_args);
        write_fixture(dir, "sla_geocentric_cancel", (AstMapping*)sm);
        sm = astAnnul(sm);
    }

    /* slamap-23: redundant precession (start==end) */
    {
        double args[] = {2000.0, 2000.0};
        AstSlaMap *sm = astSlaMap(0, "");
        astSlaAdd(sm, "PREC", 2, args);
        write_fixture(dir, "sla_prec_redundant", (AstMapping*)sm);
        sm = astAnnul(sm);
    }

    /* slamap-24: adjacent precession merge */
    {
        double args1[] = {1950.0, 1975.0};
        double args2[] = {1975.0, 2000.0};
        AstSlaMap *sm = astSlaMap(0, "");
        astSlaAdd(sm, "PREC", 2, args1);
        astSlaAdd(sm, "PREC", 2, args2);
        write_fixture(dir, "sla_prec_merge", (AstMapping*)sm);
        sm = astAnnul(sm);
    }
}

/* ===== SpecMap fixtures ===== */

static void gen_specmap_fixtures(const char *dir) {
    printf("SpecMap fixtures:\n");

    /* specmap-04: single inverted SpecMap normalizes */
    {
        double rf[] = {1.4e9, 0.0};
        AstSpecMap *sm = astSpecMap(1, 0, "");
        astSpecAdd(sm, "FRTOVL", 1, rf);
        astInvert(sm);
        write_fixture(dir, "spec_invert_normalize", (AstMapping*)sm);
        sm = astAnnul(sm);
    }

    /* specmap-02: partial cancellation */
    {
        double rf[] = {1.4e9, 0.0};
        AstSpecMap *sm = astSpecMap(1, 0, "");
        astSpecAdd(sm, "FRTOVL", 1, rf);
        astSpecAdd(sm, "VLTOFR", 1, rf);
        astSpecAdd(sm, "ENTOFR", 0, NULL);
        write_fixture(dir, "spec_partial_cancel", (AstMapping*)sm);
        sm = astAnnul(sm);
    }

    /* specmap-08: 0-arg unit conversion cancel */
    {
        AstSpecMap *sm = astSpecMap(1, 0, "");
        astSpecAdd(sm, "ENTOFR", 0, NULL);
        astSpecAdd(sm, "FRTOEN", 0, NULL);
        write_fixture(dir, "spec_unit_cancel", (AstMapping*)sm);
        sm = astAnnul(sm);
    }

    /* specmap-03: multi-map merge without step reduction */
    {
        double rf[] = {1.4e9, 0.0};
        AstSpecMap *s1 = astSpecMap(1, 0, "");
        astSpecAdd(s1, "FRTOVL", 1, rf);
        AstSpecMap *s2 = astSpecMap(1, 0, "");
        astSpecAdd(s2, "ENTOFR", 0, NULL);
        AstCmpMap *cm = astCmpMap(s1, s2, 1, "");
        write_fixture(dir, "spec_merge_no_cancel", (AstMapping*)cm);
        cm = astAnnul(cm); s1 = astAnnul(s1); s2 = astAnnul(s2);
    }

    /* specmap-10: 2-arg local-standard cancel (szargs=3, pad to be safe) */
    {
        double args[] = {0.5, 1.2, 0.0, 0.0};
        AstSpecMap *sm = astSpecMap(1, 0, "");
        astSpecAdd(sm, "LKF2HL", 2, args);
        astSpecAdd(sm, "HLF2LK", 2, args);
        write_fixture(dir, "spec_lsr_cancel", (AstMapping*)sm);
        sm = astAnnul(sm);
    }

    /* specmap-11: 3-arg geocentric cancel (szargs=4, pad) */
    {
        double args[] = {0.5, 1.2, 51544.0, 0.0, 0.0};
        AstSpecMap *sm = astSpecMap(1, 0, "");
        astSpecAdd(sm, "GEF2HL", 3, args);
        astSpecAdd(sm, "HLF2GE", 3, args);
        write_fixture(dir, "spec_geocentric_cancel", (AstMapping*)sm);
        sm = astAnnul(sm);
    }

    /* specmap-12: 6-arg topocentric cancel (szargs=7, pad) */
    {
        double args[] = {-2.5, 0.9, 1000.0, 51544.0, 0.5, 1.2, 0.0};
        AstSpecMap *sm = astSpecMap(1, 0, "");
        astSpecAdd(sm, "TPF2HL", 6, args);
        astSpecAdd(sm, "HLF2TP", 6, args);
        write_fixture(dir, "spec_topocentric_cancel", (AstMapping*)sm);
        sm = astAnnul(sm);
    }
}

/* ===== TimeMap fixtures ===== */

static void gen_timemap_fixtures(const char *dir) {
    printf("TimeMap fixtures:\n");

    /* timemap-04: single inverted TimeMap normalizes (szargs=2, pad) */
    {
        double dut[] = {0.5, 0.0};
        AstTimeMap *tm = astTimeMap(0, "");
        astTimeAdd(tm, "TAITOTT", 1, dut);
        astInvert(tm);
        write_fixture(dir, "time_invert_normalize", (AstMapping*)tm);
        tm = astAnnul(tm);
    }

    /* timemap-02: partial cancellation */
    {
        double dut[] = {0.5, 0.0};
        double dut2[] = {0.5, 0.0};
        AstTimeMap *tm = astTimeMap(0, "");
        astTimeAdd(tm, "TAITOTT", 1, dut);
        astTimeAdd(tm, "TTTOTAI", 1, dut);
        astTimeAdd(tm, "TTTOTCG", 1, dut2);
        write_fixture(dir, "time_partial_cancel", (AstMapping*)tm);
        tm = astAnnul(tm);
    }

    /* timemap-03: multi-map merge without step reduction */
    {
        double dut[] = {0.5, 0.0};
        double dut2[] = {0.3, 0.0};
        AstTimeMap *t1 = astTimeMap(0, "");
        astTimeAdd(t1, "TAITOTT", 1, dut);
        AstTimeMap *t2 = astTimeMap(0, "");
        astTimeAdd(t2, "TTTOTCG", 1, dut2);
        AstCmpMap *cm = astCmpMap(t1, t2, 1, "");
        write_fixture(dir, "time_merge_no_cancel", (AstMapping*)cm);
        cm = astAnnul(cm); t1 = astAnnul(t1); t2 = astAnnul(t2);
    }

    /* timemap-07: no-op MJDTOMJD step with zero offset (szargs=3, pad) */
    {
        double args[] = {0.0, 0.0, 0.0};
        AstTimeMap *tm = astTimeMap(0, "");
        astTimeAdd(tm, "MJDTOMJD", 2, args);
        write_fixture(dir, "time_noop_eliminate", (AstMapping*)tm);
        tm = astAnnul(tm);
    }

    /* timemap-09: 2-arg swapped pair (MJDTOJD + JDTOMJD) (szargs=3) */
    {
        double args1[] = {0.0, 2400000.5, 0.0};
        double args2[] = {2400000.5, 0.0, 0.0};
        AstTimeMap *tm = astTimeMap(0, "");
        astTimeAdd(tm, "MJDTOJD", 2, args1);
        astTimeAdd(tm, "JDTOMJD", 2, args2);
        write_fixture(dir, "time_2arg_swapped_cancel", (AstMapping*)tm);
        tm = astAnnul(tm);
    }

    /* timemap-10: 2-arg same-order pair (TAITOUTC + UTCTOTAI) (szargs=3) */
    {
        double args[] = {53000.0, 0.5, 0.0};
        AstTimeMap *tm = astTimeMap(0, "");
        astTimeAdd(tm, "TAITOUTC", 2, args);
        astTimeAdd(tm, "UTCTOTAI", 2, args);
        write_fixture(dir, "time_2arg_same_cancel", (AstMapping*)tm);
        tm = astAnnul(tm);
    }

    /* timemap-11: 3-arg pair (GMSTTOLMST + LMSTTOGMST) (szargs=3, exact) */
    {
        double args[] = {-2.5, 1000.0, 1013.0, 0.0};
        AstTimeMap *tm = astTimeMap(0, "");
        astTimeAdd(tm, "GMSTTOLMST", 3, args);
        astTimeAdd(tm, "LMSTTOGMST", 3, args);
        write_fixture(dir, "time_3arg_cancel", (AstMapping*)tm);
        tm = astAnnul(tm);
    }

    /* timemap-12: 5-arg pair (TTTOTDB + TDBTOTT) (szargs=7, pad) */
    {
        double args[] = {53000.0, 0.5, -2.5, 6378.0, 0.0, 0.0, 0.0};
        AstTimeMap *tm = astTimeMap(0, "");
        astTimeAdd(tm, "TTTOTDB", 5, args);
        astTimeAdd(tm, "TDBTOTT", 5, args);
        write_fixture(dir, "time_5arg_cancel", (AstMapping*)tm);
        tm = astAnnul(tm);
    }
}

/* ===== PermMap fixtures ===== */

static void gen_permmap_fixtures(const char *dir) {
    printf("PermMap fixtures:\n");

    /* permmap-03: two PermMaps cancel to UnitMap */
    {
        int inperm[] = {2, 1};
        int outperm[] = {2, 1};
        AstPermMap *p1 = astPermMap(2, inperm, 2, outperm, NULL, "");
        AstPermMap *p2 = astPermMap(2, inperm, 2, outperm, NULL, "");
        AstCmpMap *cm = astCmpMap(p1, p2, 1, "");
        write_fixture(dir, "perm_cancel_to_unit", (AstMapping*)cm);
        cm = astAnnul(cm); p1 = astAnnul(p1); p2 = astAnnul(p2);
    }

    /* permmap-04: single inverted PermMap normalizes */
    {
        int inperm[] = {2, 1};
        int outperm[] = {2, 1};
        AstPermMap *pm = astPermMap(2, inperm, 2, outperm, NULL, "");
        astInvert(pm);
        write_fixture(dir, "perm_invert_normalize", (AstMapping*)pm);
        pm = astAnnul(pm);
    }
}

/* ===== NormMap fixtures ===== */

static void gen_normmap_fixtures(const char *dir) {
    printf("NormMap fixtures:\n");

    /* normmap-02: NormMap wrapping basic Frame → UnitMap */
    {
        AstFrame *f = astFrame(2, "");
        AstNormMap *nm = astNormMap(f, "");
        write_fixture(dir, "normmap_basic_frame_to_unit", (AstMapping*)nm);
        nm = astAnnul(nm); f = astAnnul(f);
    }

    /* normmap-03: inverse NormMap pair cancels (SkyFrame) */
    {
        AstSkyFrame *sf = astSkyFrame("");
        AstNormMap *n1 = astNormMap(sf, "");
        AstNormMap *n2 = astNormMap(sf, "");
        astInvert(n1);
        AstCmpMap *cm = astCmpMap(n1, n2, 1, "");
        write_fixture(dir, "normmap_inverse_cancel", (AstMapping*)cm);
        cm = astAnnul(cm); n1 = astAnnul(n1); n2 = astAnnul(n2);
        sf = astAnnul(sf);
    }

    /* normmap-05: duplicate NormMaps (same direction) → extras become UnitMap */
    {
        AstSkyFrame *sf = astSkyFrame("");
        AstNormMap *n1 = astNormMap(sf, "");
        AstNormMap *n2 = astNormMap(sf, "");
        AstCmpMap *cm = astCmpMap(n1, n2, 1, "");
        write_fixture(dir, "normmap_duplicate_elim", (AstMapping*)cm);
        cm = astAnnul(cm); n1 = astAnnul(n1); n2 = astAnnul(n2);
        sf = astAnnul(sf);
    }
}

/* ===== UnitNormMap fixtures ===== */

static void gen_unitnormmap_fixtures(const char *dir) {
    printf("UnitNormMap fixtures:\n");

    /* unitnormmap-01: ShiftMap + forward UnitNormMap → adjusted centre */
    {
        double centre[] = {1.0, 2.0};
        double shifts[] = {0.5, 0.5};
        AstUnitNormMap *unm = astUnitNormMap(2, centre, "");
        AstShiftMap *sm = astShiftMap(2, shifts, "");
        AstCmpMap *cm = astCmpMap(sm, unm, 1, "");
        write_fixture(dir, "unitnormmap_shift_fwd_merge", (AstMapping*)cm);
        cm = astAnnul(cm); sm = astAnnul(sm); unm = astAnnul(unm);
    }

    /* unitnormmap-04: inverted UnitNormMap + ShiftMap → adjusted centre */
    {
        double centre[] = {1.0, 2.0};
        double shifts[] = {0.5, 0.5};
        AstUnitNormMap *unm = astUnitNormMap(2, centre, "");
        astInvert(unm);
        AstShiftMap *sm = astShiftMap(2, shifts, "");
        AstCmpMap *cm = astCmpMap(unm, sm, 1, "");
        write_fixture(dir, "unitnormmap_inv_shift_merge", (AstMapping*)cm);
        cm = astAnnul(cm); sm = astAnnul(sm); unm = astAnnul(unm);
    }

    /* unitnormmap-07: fwd + inv with same centre → UnitMap */
    {
        double centre[] = {1.0, 2.0};
        AstUnitNormMap *u1 = astUnitNormMap(2, centre, "");
        AstUnitNormMap *u2 = astUnitNormMap(2, centre, "");
        astInvert(u2);
        AstCmpMap *cm = astCmpMap(u1, u2, 1, "");
        write_fixture(dir, "unitnormmap_inverse_cancel", (AstMapping*)cm);
        cm = astAnnul(cm); u1 = astAnnul(u1); u2 = astAnnul(u2);
    }

    /* unitnormmap-09: fwd + inv with different centres → ShiftMap */
    {
        double c1[] = {1.0, 2.0};
        double c2[] = {3.0, 4.0};
        AstUnitNormMap *u1 = astUnitNormMap(2, c1, "");
        AstUnitNormMap *u2 = astUnitNormMap(2, c2, "");
        astInvert(u2);
        AstCmpMap *cm = astCmpMap(u1, u2, 1, "");
        write_fixture(dir, "unitnormmap_diff_centre_to_shift", (AstMapping*)cm);
        cm = astAnnul(cm); u1 = astAnnul(u1); u2 = astAnnul(u2);
    }
}

/* ===== GrismMap fixtures ===== */

static void gen_grismmap_fixtures(const char *dir) {
    printf("GrismMap fixtures:\n");

    /* grismmap-03: ZoomMap + inverted GrismMap → absorbed */
    {
        AstGrismMap *gm = astGrismMap("");
        astInvert(gm);
        AstZoomMap *zm = astZoomMap(1, 2.0, "");
        AstCmpMap *cm = astCmpMap(zm, gm, 1, "");
        write_fixture(dir, "grism_zoom_inv_merge", (AstMapping*)cm);
        cm = astAnnul(cm); gm = astAnnul(gm); zm = astAnnul(zm);
    }
}

/* ===== WcsMap fixtures ===== */

static void gen_wcsmap_fixtures(const char *dir) {
    printf("WcsMap fixtures:\n");

    /* wcsmap-01: AST__WCSBAD → UnitMap */
    {
        AstWcsMap *wm = astWcsMap(2, AST__WCSBAD, 1, 2, "");
        write_fixture(dir, "wcsmap_bad_to_unit", (AstMapping*)wm);
        wm = astAnnul(wm);
    }

    /* wcsmap-03: WcsMap swaps past PermMap to reach inverse merge target */
    {
        int inperm[] = {2, 1};
        int outperm[] = {2, 1};
        AstWcsMap *w1 = astWcsMap(2, AST__TAN, 1, 2, "");
        AstPermMap *pm = astPermMap(2, inperm, 2, outperm, NULL, "");
        AstWcsMap *w2 = astWcsMap(2, AST__TAN, 2, 1, "");
        astInvert(w2);
        AstCmpMap *inner = astCmpMap(pm, w2, 1, "");
        AstCmpMap *outer = astCmpMap(w1, inner, 1, "");
        write_fixture(dir, "wcsmap_perm_swap_cancel", (AstMapping*)outer);
        outer = astAnnul(outer); inner = astAnnul(inner);
        w1 = astAnnul(w1); pm = astAnnul(pm); w2 = astAnnul(w2);
    }
}

/* ===== SphMap fixtures ===== */

static void gen_sphmap_fixtures(const char *dir) {
    printf("SphMap additional fixtures:\n");

    /* sphmap-02: SphMap(UnitRadius=1) + Inverse(SphMap) cancels */
    {
        AstSphMap *s1 = astSphMap("UnitRadius=1");
        AstSphMap *s2 = astSphMap("UnitRadius=1");
        astInvert(s2);
        AstCmpMap *cm = astCmpMap(s1, s2, 1, "");
        write_fixture(dir, "sph_fwd_inv_unitradius_cancel", (AstMapping*)cm);
        cm = astAnnul(cm); s1 = astAnnul(s1); s2 = astAnnul(s2);
    }

    /* sphmap-09: Inv(SphMap) + ZoomMap + SphMap → WinMap (ZoomMap variant) */
    {
        AstSphMap *s1 = astSphMap("UnitRadius=1");
        astInvert(s1);
        AstZoomMap *zm = astZoomMap(3, -1.0, "");
        AstSphMap *s2 = astSphMap("UnitRadius=1");
        AstCmpMap *inner = astCmpMap(zm, s2, 1, "");
        AstCmpMap *outer = astCmpMap(s1, inner, 1, "");
        write_fixture(dir, "sph_zoom_sandwich", (AstMapping*)outer);
        outer = astAnnul(outer); inner = astAnnul(inner);
        s1 = astAnnul(s1); zm = astAnnul(zm); s2 = astAnnul(s2);
    }
}

/* ===== PcdMap fixtures ===== */

static void gen_pcdmap_fixtures(const char *dir) {
    printf("PcdMap fixtures:\n");

    /* pcdmap-03: PcdMap + UnitMap → UnitMap eliminated */
    {
        double pcdcen[] = {100.0, 100.0};
        AstPcdMap *pm = astPcdMap(0.001, pcdcen, "");
        AstUnitMap *um = astUnitMap(2, "");
        AstCmpMap *cm = astCmpMap(pm, um, 1, "");
        write_fixture(dir, "pcd_unit_series_merge", (AstMapping*)cm);
        cm = astAnnul(cm); pm = astAnnul(pm); um = astAnnul(um);
    }

    /* pcdmap-04: PcdMap swaps with ZoomMap to reach inverse PcdMap */
    {
        double pcdcen[] = {0.0, 0.0};
        AstPcdMap *p1 = astPcdMap(0.001, pcdcen, "");
        AstZoomMap *zm = astZoomMap(2, 1.0, "");
        AstPcdMap *p2 = astPcdMap(0.001, pcdcen, "");
        astInvert(p2);
        AstCmpMap *inner = astCmpMap(zm, p2, 1, "");
        AstCmpMap *outer = astCmpMap(p1, inner, 1, "");
        write_fixture(dir, "pcd_zoom_swap_cancel", (AstMapping*)outer);
        outer = astAnnul(outer); inner = astAnnul(inner);
        p1 = astAnnul(p1); zm = astAnnul(zm); p2 = astAnnul(p2);
    }
}

/* ===== SwitchMap fixtures ===== */

static void gen_switchmap_fixtures(const char *dir) {
    printf("SwitchMap fixtures:\n");

    /* switchmap-03: inverted SwitchMap normalizes */
    {
        AstZoomMap *fsel = astZoomMap(1, 1.0, "");
        AstZoomMap *isel = astZoomMap(1, 1.0, "");
        AstZoomMap *r1 = astZoomMap(1, 2.0, "");
        AstMapping *routes[] = {(AstMapping*)r1};
        AstSwitchMap *sw = astSwitchMap(fsel, isel, 1, (void**)routes, "");
        astInvert(sw);
        write_fixture(dir, "switchmap_invert_normalize", (AstMapping*)sw);
        sw = astAnnul(sw); fsel = astAnnul(fsel); isel = astAnnul(isel);
        r1 = astAnnul(r1);
    }

    /* switchmap-04: SwitchMap internal simplification */
    {
        AstZoomMap *fsel = astZoomMap(1, 1.0, "");
        AstZoomMap *isel = astZoomMap(1, 1.0, "");
        AstZoomMap *za = astZoomMap(1, 2.0, "");
        AstZoomMap *zb = astZoomMap(1, 3.0, "");
        AstCmpMap *route = astCmpMap(za, zb, 1, "");
        AstMapping *routes[] = {(AstMapping*)route};
        AstSwitchMap *sw = astSwitchMap(fsel, isel, 1, (void**)routes, "");
        write_fixture(dir, "switchmap_internal_simplify", (AstMapping*)sw);
        sw = astAnnul(sw); fsel = astAnnul(fsel); isel = astAnnul(isel);
        za = astAnnul(za); zb = astAnnul(zb); route = astAnnul(route);
    }
}

/* ===== TranMap additional fixtures ===== */

static void gen_tranmap_extra_fixtures(const char *dir) {
    printf("TranMap extra fixtures:\n");

    /* tranmap-04: adjacent TranMap merge in series */
    {
        AstZoomMap *z2 = astZoomMap(1, 2.0, "");
        AstZoomMap *z05 = astZoomMap(1, 0.5, "");
        AstZoomMap *z3 = astZoomMap(1, 3.0, "");
        AstZoomMap *z033 = astZoomMap(1, 1.0/3.0, "");
        AstTranMap *t1 = astTranMap(z2, z05, "");
        AstTranMap *t2 = astTranMap(z3, z033, "");
        AstCmpMap *cm = astCmpMap(t1, t2, 1, "");
        write_fixture(dir, "tranmap_adjacent_merge", (AstMapping*)cm);
        cm = astAnnul(cm); t1 = astAnnul(t1); t2 = astAnnul(t2);
        z2 = astAnnul(z2); z05 = astAnnul(z05);
        z3 = astAnnul(z3); z033 = astAnnul(z033);
    }
}

/* ===== MatrixMap cascade fixtures ===== */

static void gen_matrix_cascade_fixtures(const char *dir) {
    printf("MatrixMap cascade fixtures:\n");

    /* matrixmap-13: MatrixMap swaps past WinMap to reach merge target */
    {
        double diag1[] = {2.0, 3.0};
        double ina[] = {0, 0}, inb[] = {1, 1}, outa[] = {1, 2}, outb[] = {5, 8};
        double diag2[] = {4.0, 5.0};
        AstMatrixMap *m1 = astMatrixMap(2, 2, 1, diag1, "");
        AstWinMap *wm = astWinMap(2, ina, inb, outa, outb, "");
        AstMatrixMap *m2 = astMatrixMap(2, 2, 1, diag2, "");
        AstCmpMap *inner = astCmpMap(wm, m2, 1, "");
        AstCmpMap *outer = astCmpMap(m1, inner, 1, "");
        write_fixture(dir, "matrix_swap_past_win", (AstMapping*)outer);
        outer = astAnnul(outer); inner = astAnnul(inner);
        m1 = astAnnul(m1); wm = astAnnul(wm); m2 = astAnnul(m2);
    }

    /* matrixmap-15: Swap with WinMap for local simplification.
       A full (non-diagonal) MatrixMap can't directly merge with WinMap
       (MatWin2 requires diagonal). After MatWin swap, the resulting WinMap
       has all scales=1 (absorbed into MatrixMap) and simplifies to ShiftMap
       — a class change that triggers the swap acceptance. */
    {
        double mat[] = {1.0, 0.0, 1.0, 1.0};
        AstMatrixMap *mm = astMatrixMap(2, 2, 0, mat, "");
        double ina[] = {0, 0}, inb[] = {1, 1};
        double outa[] = {1, 2}, outb[] = {4, 6};
        AstWinMap *wm = astWinMap(2, ina, inb, outa, outb, "");
        AstCmpMap *cm = astCmpMap(mm, wm, 1, "");
        write_fixture(dir, "matrix_swap_win_simplifies", (AstMapping*)cm);
        cm = astAnnul(cm); mm = astAnnul(mm); wm = astAnnul(wm);
    }

    /* matrixmap-14: MatrixMap swaps past PermMap to reach merge target */
    {
        double diag1[] = {2.0, 3.0};
        int inperm[] = {2, 1};
        int outperm[] = {2, 1};
        double diag2[] = {4.0, 5.0};
        AstMatrixMap *m1 = astMatrixMap(2, 2, 1, diag1, "");
        AstPermMap *pm = astPermMap(2, inperm, 2, outperm, NULL, "");
        AstMatrixMap *m2 = astMatrixMap(2, 2, 1, diag2, "");
        AstCmpMap *inner = astCmpMap(pm, m2, 1, "");
        AstCmpMap *outer = astCmpMap(m1, inner, 1, "");
        write_fixture(dir, "matrix_swap_past_perm", (AstMapping*)outer);
        outer = astAnnul(outer); inner = astAnnul(inner);
        m1 = astAnnul(m1); pm = astAnnul(pm); m2 = astAnnul(m2);
    }
}

/* ===== SlaMap remaining fixtures ===== */

static void gen_slamap_extra_fixtures(const char *dir) {
    printf("SlaMap extra fixtures:\n");

    /* slamap-14: 1-arg helio-ecliptic (EQHE+HEEQ) */
    {
        double args[] = {51544.0};
        AstSlaMap *sm = astSlaMap(0, "");
        astSlaAdd(sm, "EQHE", 1, args);
        astSlaAdd(sm, "HEEQ", 1, args);
        write_fixture(dir, "sla_helioecl_cancel", (AstMapping*)sm);
        sm = astAnnul(sm);
    }

    /* slamap-15: 1-arg HA (R2H+H2R) */
    {
        double args[] = {3.5};
        AstSlaMap *sm = astSlaMap(0, "");
        astSlaAdd(sm, "R2H", 1, args);
        astSlaAdd(sm, "H2R", 1, args);
        write_fixture(dir, "sla_ha_cancel", (AstMapping*)sm);
        sm = astAnnul(sm);
    }

    /* slamap-17: 2-arg AzEl (H2E+E2H) */
    {
        double args[] = {0.9, 0.001};
        AstSlaMap *sm = astSlaMap(0, "");
        astSlaAdd(sm, "H2E", 2, args);
        astSlaAdd(sm, "E2H", 2, args);
        write_fixture(dir, "sla_azel_cancel", (AstMapping*)sm);
        sm = astAnnul(sm);
    }
}

/* ===== SwitchMap extra fixtures ===== */

static void gen_switchmap_extra_fixtures(const char *dir) {
    printf("SwitchMap extra fixtures:\n");

    /* switchmap-01: SwitchMap + Inverse(SwitchMap) cancel to UnitMap */
    {
        AstZoomMap *fsel = astZoomMap(1, 1.0, "");
        AstZoomMap *isel = astZoomMap(1, 1.0, "");
        AstZoomMap *r1 = astZoomMap(1, 2.0, "");
        AstMapping *routes[] = {(AstMapping*)r1};
        AstSwitchMap *s1 = astSwitchMap(fsel, isel, 1, (void**)routes, "");
        AstSwitchMap *s2 = astSwitchMap(fsel, isel, 1, (void**)routes, "");
        astInvert(s2);
        AstCmpMap *cm = astCmpMap(s1, s2, 1, "");
        write_fixture(dir, "switchmap_inverse_cancel", (AstMapping*)cm);
        cm = astAnnul(cm); s1 = astAnnul(s1); s2 = astAnnul(s2);
        fsel = astAnnul(fsel); isel = astAnnul(isel); r1 = astAnnul(r1);
    }
}

/* ===== PcdMap extra fixtures ===== */

static void gen_pcdmap_extra_fixtures(const char *dir) {
    printf("PcdMap extra fixtures:\n");

    /* pcdmap-06: PcdMap swap without merge target — accepted because
       ZoomMap(1) simplifies to UnitMap after the swap. */
    {
        double pcdcen[] = {0.0, 0.0};
        AstPcdMap *pm = astPcdMap(0.001, pcdcen, "");
        AstZoomMap *zm = astZoomMap(2, 1.0, "");
        AstCmpMap *cm = astCmpMap(pm, zm, 1, "");
        write_fixture(dir, "pcd_swap_zoom_simplifies", (AstMapping*)cm);
        cm = astAnnul(cm); pm = astAnnul(pm); zm = astAnnul(zm);
    }

    /* pcdmap-05: PcdMap swaps with PermMap(axis swap) to reach inverse */
    {
        double pcdcen[] = {0.0, 0.0};
        int inperm[] = {2, 1};
        int outperm[] = {2, 1};
        AstPcdMap *p1 = astPcdMap(0.001, pcdcen, "");
        AstPermMap *pm = astPermMap(2, inperm, 2, outperm, NULL, "");
        AstPcdMap *p2 = astPcdMap(0.001, pcdcen, "");
        astInvert(p2);
        AstCmpMap *inner = astCmpMap(pm, p2, 1, "");
        AstCmpMap *outer = astCmpMap(p1, inner, 1, "");
        write_fixture(dir, "pcd_perm_swap_cancel", (AstMapping*)outer);
        outer = astAnnul(outer); inner = astAnnul(inner);
        p1 = astAnnul(p1); pm = astAnnul(pm); p2 = astAnnul(p2);
    }
}

/* ===== WinMap extra cascade fixtures ===== */

static void gen_win_extra_cascade_fixtures(const char *dir) {
    printf("WinMap extra cascade fixtures:\n");

    /* winmap-26: Swap accepted because swapped Mapping simplifies.
       Setup: PermMap(2→3, passes axes 1,2 and adds constant on axis 3)
       followed by WinMap(3D, a=[0,0,5], b=[1,1,3]).
       After WinPerm swap: new WinMap(2D, a=[0,0], b=[1,1]) → simplifies to UnitMap. */
    {
        int inperm[] = {1, 2};
        int outperm[] = {1, 2, -1};
        double consts[] = {5.0};
        AstPermMap *pm = astPermMap(2, inperm, 3, outperm, consts, "");
        double ina[] = {0, 0, 0}, inb[] = {1, 1, 1};
        double outa[] = {0, 0, 5}, outb[] = {1, 1, 8};
        AstWinMap *wm = astWinMap(3, ina, inb, outa, outb, "");
        AstCmpMap *cm = astCmpMap(pm, wm, 1, "");
        write_fixture(dir, "win_swap_simplifies", (AstMapping*)cm);
        cm = astAnnul(cm); pm = astAnnul(pm); wm = astAnnul(wm);
    }

    /* winmap-27: Swap accepted because outer neighbours can merge.
       Setup: [PermMap(swap)] [WinMap(2D)] [PermMap(swap)]
       WinMap can't merge directly with PermMap. After swap with left PermMap,
       the resulting PermMap is adjacent to the right PermMap and they merge. */
    {
        int inperm[] = {2, 1};
        int outperm[] = {2, 1};
        double ina[] = {0, 0}, inb[] = {1, 1};
        double outa[] = {1, 2}, outb[] = {4, 6};
        AstPermMap *p1 = astPermMap(2, inperm, 2, outperm, NULL, "");
        AstWinMap *wm = astWinMap(2, ina, inb, outa, outb, "");
        AstPermMap *p2 = astPermMap(2, inperm, 2, outperm, NULL, "");
        AstCmpMap *inner = astCmpMap(wm, p2, 1, "");
        AstCmpMap *outer = astCmpMap(p1, inner, 1, "");
        write_fixture(dir, "win_swap_outer_merge", (AstMapping*)outer);
        outer = astAnnul(outer); inner = astAnnul(inner);
        p1 = astAnnul(p1); wm = astAnnul(wm); p2 = astAnnul(p2);
    }

    /* winmap-14: WinMap merges with neighbouring parallel CmpMap (lower) */
    {
        double ina[] = {0, 0}, inb[] = {1, 1};
        double outa[] = {1, 2}, outb[] = {3, 5};
        double shifts[] = {1.0};
        AstShiftMap *sm1 = astShiftMap(1, shifts, "");
        double shifts2[] = {2.0};
        AstShiftMap *sm2 = astShiftMap(1, shifts2, "");
        AstCmpMap *par = astCmpMap(sm1, sm2, 0, "");
        AstWinMap *wm = astWinMap(2, ina, inb, outa, outb, "");
        AstCmpMap *cm = astCmpMap(par, wm, 1, "");
        write_fixture(dir, "win_cmpmap_parallel_merge", (AstMapping*)cm);
        cm = astAnnul(cm); par = astAnnul(par); wm = astAnnul(wm);
        sm1 = astAnnul(sm1); sm2 = astAnnul(sm2);
    }
}

/* ===== SlaMap 4-arg fixtures ===== */

static void gen_slamap_4arg_fixtures(const char *dir) {
    printf("SlaMap 4-arg fixtures:\n");

    /* slamap-18: 4-arg HPC (HPCEQ+EQHPC) */
    {
        double args[] = {51544.0, 0.5, 1.2, 150.0e6};
        AstSlaMap *sm = astSlaMap(0, "");
        astSlaAdd(sm, "HPCEQ", 4, args);
        astSlaAdd(sm, "EQHPC", 4, args);
        write_fixture(dir, "sla_hpc_cancel", (AstMapping*)sm);
        sm = astAnnul(sm);
    }

    /* slamap-19: 4-arg HPR (HPREQ+EQHPR) */
    {
        double args[] = {51544.0, 0.5, 1.2, 150.0e6};
        AstSlaMap *sm = astSlaMap(0, "");
        astSlaAdd(sm, "HPREQ", 4, args);
        astSlaAdd(sm, "EQHPR", 4, args);
        write_fixture(dir, "sla_hpr_cancel", (AstMapping*)sm);
        sm = astAnnul(sm);
    }
}

/* ===== WinMap cascade fixtures ===== */

static void gen_win_cascade_fixtures(const char *dir) {
    printf("WinMap cascade fixtures:\n");

    /* winmap-18: WinMap swaps past MatrixMap to reach merge target */
    {
        double ina[] = {0, 0}, inb[] = {1, 1};
        double outa1[] = {1, 2}, outb1[] = {3, 4};
        double outa2[] = {5, 6}, outb2[] = {7, 8};
        double diag[] = {2.0, 3.0};
        AstWinMap *w1 = astWinMap(2, ina, inb, outa1, outb1, "");
        AstMatrixMap *mm = astMatrixMap(2, 2, 1, diag, "");
        AstWinMap *w2 = astWinMap(2, ina, inb, outa2, outb2, "");
        AstCmpMap *inner = astCmpMap(mm, w2, 1, "");
        AstCmpMap *outer = astCmpMap(w1, inner, 1, "");
        write_fixture(dir, "win_swap_past_matrix", (AstMapping*)outer);
        outer = astAnnul(outer); inner = astAnnul(inner);
        w1 = astAnnul(w1); mm = astAnnul(mm); w2 = astAnnul(w2);
    }

    /* winmap-20: WinMap swaps past WcsMap to reach merge target */
    {
        double ina[] = {0, 0}, inb[] = {1, 1};
        double outa1[] = {0, 0}, outb1[] = {1, 1};
        double outa2[] = {0, 0}, outb2[] = {1, 1};
        AstWinMap *w1 = astWinMap(2, ina, inb, outa1, outb1, "");
        AstWcsMap *wcs = astWcsMap(2, AST__TAN, 1, 2, "");
        AstWinMap *w2 = astWinMap(2, ina, inb, outa2, outb2, "");
        AstCmpMap *inner = astCmpMap(wcs, w2, 1, "");
        AstCmpMap *outer = astCmpMap(w1, inner, 1, "");
        write_fixture(dir, "win_swap_past_wcsmap", (AstMapping*)outer);
        outer = astAnnul(outer); inner = astAnnul(inner);
        w1 = astAnnul(w1); wcs = astAnnul(wcs); w2 = astAnnul(w2);
    }
}

/* ===== SelectorMap fixtures ===== */

static void gen_selectormap_fixtures(const char *dir) {
    printf("SelectorMap fixtures:\n");

    /* selectormap-05: SelectorMap + Inverse(SelectorMap) cancel to UnitMap */
    {
        AstFrame *f = astFrame(2, "");
        double lbnd[] = {0.0, 0.0};
        double ubnd[] = {10.0, 10.0};
        AstBox *box = astBox(f, 1, lbnd, ubnd, NULL, "");
        AstMapping *regs[] = {(AstMapping*)box};
        AstSelectorMap *s1 = astSelectorMap(1, (void**)regs, AST__BAD, "");
        AstSelectorMap *s2 = astSelectorMap(1, (void**)regs, AST__BAD, "");
        astInvert(s2);
        AstCmpMap *cm = astCmpMap(s1, s2, 1, "");
        write_fixture(dir, "selectormap_inverse_cancel", (AstMapping*)cm);
        cm = astAnnul(cm); s1 = astAnnul(s1); s2 = astAnnul(s2);
        box = astAnnul(box); f = astAnnul(f);
    }
}

/* ===== DssMap fixtures ===== */

static void gen_dssmap_fixtures(const char *dir) {
    printf("DssMap fixtures:\n");

    /* dssmap-07: non-inverted DssMap absorbs preceding WinMap.
       Read DSS headers from dss.fits-dss to create a FrameSet, extract
       the Mapping (which contains a DssMap), then compose with a WinMap. */
    {
        FILE *fp = fopen("ast_tester/dss.fits-dss", "r");
        if (fp) {
            char line[256];
            AstFitsChan *fc = astFitsChan(NULL, NULL, "Encoding=DSS");
            while (fgets(line, sizeof(line), fp)) {
                size_t len = strlen(line);
                if (len > 0 && line[len-1] == '\n') line[len-1] = '\0';
                astPutFits(fc, line, 0);
            }
            fclose(fp);
            astClear(fc, "Card");
            AstFrameSet *fs = (AstFrameSet *)astRead(fc);
            if (fs && astIsAFrameSet(fs)) {
                AstMapping *map = astGetMapping(fs, AST__BASE, AST__CURRENT);
                /* Use unsimplified map — it should contain a DssMap. */
                double ina[] = {0, 0}, inb[] = {1, 1};
                double outa[] = {1, 1}, outb[] = {2, 2};
                AstWinMap *wm = astWinMap(2, ina, inb, outa, outb, "");
                AstCmpMap *cm = astCmpMap(wm, map, 1, "");
                write_fixture(dir, "dss_winmap_absorb", (AstMapping*)cm);
                cm = astAnnul(cm); wm = astAnnul(wm);
                map = astAnnul(map);
                fs = astAnnul(fs);
            } else {
                fprintf(stderr, "WARNING: Could not read DSS FrameSet\n");
            }
            fc = astAnnul(fc);
        } else {
            fprintf(stderr, "WARNING: Could not open ast_tester/dss.fits-dss\n");
        }
    }
}

/* ===== Negative (guard-rejection) fixtures ===== */
/* These test that astSimplify does NOT structurally change the input.
   Only the .map file is generated; the test uses the same file as both
   input and reference with skip_string_compare=yes (astequal only). */

static void write_negative_fixture(const char *dir, const char *name, AstMapping *map) {
    char path_map[512];
    AstChannel *chan;
    extern int *astGetStatusPtr_( void );

    if (!astOK) astClearStatus;

    snprintf(path_map, sizeof(path_map), "%s/%s.map", dir, name);

    chan = astChannel(NULL, NULL, "SinkFile=%s", path_map);
    if (astWrite(chan, map) != 1) {
        fprintf(stderr, "ERROR: astWrite failed for %s.map\n", name);
    }
    chan = astAnnul(chan);
    printf("  %s\n", name);
}

static void gen_negative_fixtures(const char *dir) {
    if (!astOK) astClearStatus;
    printf("Negative (guard-rejection) fixtures:\n");

    /* matrixmap-04: diagonal MatrixMap with unequal elements — can't become ZoomMap */
    {
        double diag[] = {2.0, 3.0};
        AstMatrixMap *mm = astMatrixMap(2, 2, 1, diag, "");
        write_negative_fixture(dir, "neg_matrix_diag_unequal", (AstMapping*)mm);
        mm = astAnnul(mm);
    }

    /* matrixmap-06: full MatrixMap with non-zero off-diagonal — stays full */
    {
        double mat[] = {1.0, 2.0, 3.0, 4.0};
        AstMatrixMap *mm = astMatrixMap(2, 2, 0, mat, "");
        write_negative_fixture(dir, "neg_matrix_full_offdiag", (AstMapping*)mm);
        mm = astAnnul(mm);
    }

    /* matrixmap-17: full MatrixMap in parallel — no merge attempted */
    {
        double mat[] = {1.0, 2.0, 3.0, 4.0};
        AstMatrixMap *mm = astMatrixMap(2, 2, 0, mat, "");
        AstZoomMap *zm = astZoomMap(1, 5.0, "");
        AstCmpMap *cm = astCmpMap(mm, zm, 0, "");
        write_negative_fixture(dir, "neg_matrix_parallel_no_merge", (AstMapping*)cm);
        cm = astAnnul(cm); mm = astAnnul(mm); zm = astAnnul(zm);
    }

    /* winmap-13: WinMap with non-mergeable neighbour in series.
       Use a MathMap (not directly mergeable with WinMap) with matching dims. */
    {
        double ina[] = {0, 0}, inb[] = {1, 1};
        double outa[] = {1, 2}, outb[] = {4, 6};
        const char *fwd[] = {"y1 = x1", "y2 = x2"};
        const char *inv[] = {"x1 = y1", "x2 = y2"};
        AstWinMap *wm = astWinMap(2, ina, inb, outa, outb, "");
        AstMathMap *mm = astMathMap(2, 2, 2, fwd, 2, inv, "");
        AstCmpMap *cm = astCmpMap(wm, mm, 1, "");
        write_negative_fixture(dir, "neg_win_nonmergeable_series", (AstMapping*)cm);
        cm = astAnnul(cm); wm = astAnnul(wm); mm = astAnnul(mm);
    }

    /* winmap-38: WinMap with non-mergeable neighbour in parallel */
    {
        double ina[] = {0, 0}, inb[] = {1, 1};
        double outa[] = {1, 2}, outb[] = {4, 6};
        double mat[] = {1.0, 2.0, 3.0, 4.0};
        AstWinMap *wm = astWinMap(2, ina, inb, outa, outb, "");
        AstMatrixMap *mm = astMatrixMap(2, 2, 0, mat, "");
        AstCmpMap *cm = astCmpMap(wm, mm, 0, "");
        write_negative_fixture(dir, "neg_win_nonmergeable_parallel", (AstMapping*)cm);
        cm = astAnnul(cm); wm = astAnnul(wm); mm = astAnnul(mm);
    }

    /* polymap-05: non-linear PolyMap — linearization refused */
    {
        double coeff_f[] = {2.0, 1, 2, 0,
                            3.0, 1, 0, 1};
        int ncoeff_f = 2;
        AstPolyMap *pm = astPolyMap(2, 1, ncoeff_f, coeff_f, 0, NULL, "");
        write_negative_fixture(dir, "neg_poly_nonlinear", (AstMapping*)pm);
        pm = astAnnul(pm);
    }

    /* polymap-09: two forward PolyMaps in series — same direction refuses cancel */
    {
        double coeff_f[] = {2.0, 1, 1};
        AstPolyMap *p1 = astPolyMap(1, 1, 1, coeff_f, 0, NULL, "");
        AstPolyMap *p2 = astPolyMap(1, 1, 1, coeff_f, 0, NULL, "");
        AstCmpMap *cm = astCmpMap(p1, p2, 1, "");
        write_negative_fixture(dir, "neg_poly_same_direction", (AstMapping*)cm);
        cm = astAnnul(cm); p1 = astAnnul(p1); p2 = astAnnul(p2);
    }

    /* mathmap-05: MathMap pair without SimpFI — refuses cancellation */
    {
        const char *fwd[] = {"y = x"};
        const char *inv[] = {"x = y"};
        AstMathMap *m1 = astMathMap(1, 1, 1, fwd, 1, inv, "SimpFI=0,SimpIF=0");
        AstMathMap *m2 = astMathMap(1, 1, 1, fwd, 1, inv, "SimpFI=0,SimpIF=0");
        astInvert(m2);
        AstCmpMap *cm = astCmpMap(m1, m2, 1, "");
        write_negative_fixture(dir, "neg_math_no_simpfi", (AstMapping*)cm);
        cm = astAnnul(cm); m1 = astAnnul(m1); m2 = astAnnul(m2);
    }

    /* slamap-20: 1-arg pair with mismatched argument — prevents cancel */
    {
        double args1[] = {1950.0};
        double args2[] = {2000.0};
        AstSlaMap *sm = astSlaMap(0, "");
        astSlaAdd(sm, "ADDET", 1, args1);
        astSlaAdd(sm, "SUBET", 1, args2);
        write_negative_fixture(dir, "neg_sla_arg_mismatch", (AstMapping*)sm);
        sm = astAnnul(sm);
    }

    /* specmap-13: 1-arg pair with mismatched argument */
    {
        double rf1[] = {1.4e9, 0.0};
        double rf2[] = {2.0e9, 0.0};
        AstSpecMap *sm = astSpecMap(1, 0, "");
        astSpecAdd(sm, "FRTOVL", 1, rf1);
        astSpecAdd(sm, "VLTOFR", 1, rf2);
        write_negative_fixture(dir, "neg_spec_arg_mismatch", (AstMapping*)sm);
        sm = astAnnul(sm);
    }

    /* timemap-13: 1-arg pair with mismatched argument */
    {
        double dut1[] = {0.5, 0.0};
        double dut2[] = {0.3, 0.0};
        AstTimeMap *tm = astTimeMap(0, "");
        astTimeAdd(tm, "TAITOTT", 1, dut1);
        astTimeAdd(tm, "TTTOTAI", 1, dut2);
        write_negative_fixture(dir, "neg_time_arg_mismatch", (AstMapping*)tm);
        tm = astAnnul(tm);
    }

    /* sphmap-07: SphMap followed by non-SphMap — refuses simplification.
       SphMap forward is 3in→2out (Cartesian to spherical). Follow with ZoomMap(2). */
    {
        if (!astOK) astClearStatus;
        AstSphMap *sm = astSphMap("UnitRadius=1");
        AstZoomMap *zm = astZoomMap(2, 2.0, "");
        AstCmpMap *cm = astCmpMap(sm, zm, 1, "");
        write_negative_fixture(dir, "neg_sph_non_sphmap_neighbour", (AstMapping*)cm);
        cm = astAnnul(cm); sm = astAnnul(sm); zm = astAnnul(zm);
    }

    /* grismmap-08: two GrismMaps same direction — refuses cancel */
    {
        if (!astOK) astClearStatus;
        AstGrismMap *g1 = astGrismMap("");
        AstGrismMap *g2 = astGrismMap("");
        AstCmpMap *cm = astCmpMap(g1, g2, 1, "");
        write_negative_fixture(dir, "neg_grism_same_direction", (AstMapping*)cm);
        cm = astAnnul(cm); g1 = astAnnul(g1); g2 = astAnnul(g2);
    }

    /* wcsmap-08: two WcsMaps same direction — refuses cancel */
    {
        if (!astOK) astClearStatus;
        AstWcsMap *w1 = astWcsMap(2, AST__TAN, 1, 2, "");
        AstWcsMap *w2 = astWcsMap(2, AST__TAN, 1, 2, "");
        AstCmpMap *cm = astCmpMap(w1, w2, 1, "");
        write_negative_fixture(dir, "neg_wcs_same_direction", (AstMapping*)cm);
        cm = astAnnul(cm); w1 = astAnnul(w1); w2 = astAnnul(w2);
    }

    /* slamap-25: adjacent precession with non-common equinox */
    {
        double args1[] = {1950.0, 1975.0};
        double args2[] = {2000.0, 2025.0};
        AstSlaMap *sm = astSlaMap(0, "");
        astSlaAdd(sm, "PREC", 2, args1);
        astSlaAdd(sm, "PREC", 2, args2);
        write_negative_fixture(dir, "neg_sla_prec_no_common", (AstMapping*)sm);
        sm = astAnnul(sm);
    }

    /* unitnormmap-03: WinMap(non-unit scale) + UnitNormMap — refuses merge */
    {
        double centre[] = {1.0, 2.0};
        double ina[] = {0, 0}, inb[] = {1, 1};
        double outa[] = {0, 0}, outb[] = {2, 3};
        AstWinMap *wm = astWinMap(2, ina, inb, outa, outb, "");
        AstUnitNormMap *unm = astUnitNormMap(2, centre, "");
        AstCmpMap *cm = astCmpMap(wm, unm, 1, "");
        write_negative_fixture(dir, "neg_unitnormmap_nonunit_scale", (AstMapping*)cm);
        cm = astAnnul(cm); wm = astAnnul(wm); unm = astAnnul(unm);
    }
}

/* ===== IntraMap fixtures ===== */
/* Note: IntraMap requires registered functions, skip for now */

/* ===== Main ===== */

int main(void) {
    int status = 0;
    const char *dir = "ast_tester/simplify_fixtures";

    astWatch(&status);
    astBegin;

    gen_zoom_fixtures(dir);
    gen_win_fixtures(dir);
    gen_unit_fixtures(dir);
    gen_matrix_fixtures(dir);
    gen_cmpmap_fixtures(dir);
    gen_tranmap_fixtures(dir);
    gen_ratemap_fixtures(dir);
    gen_slamap_fixtures(dir);
    gen_specmap_fixtures(dir);
    gen_timemap_fixtures(dir);
    gen_permmap_fixtures(dir);
    gen_normmap_fixtures(dir);
    gen_unitnormmap_fixtures(dir);
    gen_grismmap_fixtures(dir);
    gen_wcsmap_fixtures(dir);
    gen_sphmap_fixtures(dir);
    gen_pcdmap_fixtures(dir);
    gen_switchmap_fixtures(dir);
    gen_tranmap_extra_fixtures(dir);
    gen_matrix_cascade_fixtures(dir);
    gen_slamap_extra_fixtures(dir);
    gen_slamap_4arg_fixtures(dir);
    gen_win_cascade_fixtures(dir);
    gen_selectormap_fixtures(dir);
    gen_switchmap_extra_fixtures(dir);
    gen_pcdmap_extra_fixtures(dir);
    gen_win_extra_cascade_fixtures(dir);
    /* DssMap: skipped — protected constructor, and DSS FitsChan encoding
       no longer creates DssMap objects (decomposes to WcsMap pipeline). */

    gen_negative_fixtures(dir);

    astEnd;

    if (!astOK) {
        fprintf(stderr, "AST error occurred (status=%d)\n", status);
        return 1;
    }

    printf("\nAll fixtures generated successfully.\n");
    return 0;
}
