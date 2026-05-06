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

    /* polymap-09: two forward non-linear PolyMaps in series — same direction refuses */
    {
        double coeff_f[] = {1.0, 1, 2};
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

/* ===== Negative fixtures batch 2 ===== */

static void gen_negative_fixtures_2(const char *dir) {
    if (!astOK) astClearStatus;
    printf("Negative fixtures batch 2:\n");

    /* zoommap-15: ZoomMap can't be absorbed — neighbours not MatrixMap/WinMap/ShiftMap.
       SphMap(forward) is 3in→2out, so pair with ZoomMap(2) for valid dims. */
    {
        if (!astOK) astClearStatus;
        AstSphMap *sm = astSphMap("UnitRadius=1");
        AstZoomMap *zm = astZoomMap(2, 2.0, "");
        AstCmpMap *cm = astCmpMap(sm, zm, 1, "");
        write_negative_fixture(dir, "neg_zoom_no_absorb", (AstMapping*)cm);
        cm = astAnnul(cm); sm = astAnnul(sm); zm = astAnnul(zm);
    }

    /* lutmap-05: non-linear LutMap — WinMap replacement refused */
    {
        if (!astOK) astClearStatus;
        double lut[] = {0.0, 1.0, 4.0, 9.0, 16.0};
        AstLutMap *lm = astLutMap(5, lut, 1.0, 1.0, "");
        write_negative_fixture(dir, "neg_lut_nonlinear", (AstMapping*)lm);
        lm = astAnnul(lm);
    }

    /* lutmap-09: neighbouring LutMap not inverse-equal (different tables) */
    {
        if (!astOK) astClearStatus;
        double lut1[] = {0.0, 1.0, 4.0, 9.0, 16.0};
        double lut2[] = {0.0, 2.0, 6.0, 12.0, 20.0};
        AstLutMap *l1 = astLutMap(5, lut1, 1.0, 1.0, "");
        AstLutMap *l2 = astLutMap(5, lut2, 1.0, 1.0, "");
        astInvert(l2);
        AstCmpMap *cm = astCmpMap(l1, l2, 1, "");
        write_negative_fixture(dir, "neg_lut_different_tables", (AstMapping*)cm);
        cm = astAnnul(cm); l1 = astAnnul(l1); l2 = astAnnul(l2);
    }

    /* pcdmap-08: two PcdMaps in parallel — refuses simplification */
    {
        if (!astOK) astClearStatus;
        double cen[] = {0.0, 0.0};
        AstPcdMap *p1 = astPcdMap(0.001, cen, "");
        AstPcdMap *p2 = astPcdMap(0.002, cen, "");
        AstCmpMap *cm = astCmpMap(p1, p2, 0, "");
        write_negative_fixture(dir, "neg_pcd_parallel", (AstMapping*)cm);
        cm = astAnnul(cm); p1 = astAnnul(p1); p2 = astAnnul(p2);
    }

    /* pcdmap-09: PcdMap + non-PcdMap neighbour in series — refuses merge */
    {
        if (!astOK) astClearStatus;
        double cen[] = {0.0, 0.0};
        double shifts[] = {1.0, 2.0};
        AstPcdMap *pm = astPcdMap(0.001, cen, "");
        AstShiftMap *sm = astShiftMap(2, shifts, "");
        AstCmpMap *cm = astCmpMap(pm, sm, 1, "");
        write_negative_fixture(dir, "neg_pcd_nonpcd_neighbour", (AstMapping*)cm);
        cm = astAnnul(cm); pm = astAnnul(pm); sm = astAnnul(sm);
    }

    /* cmpmap-04: series CmpMap in parallel list — mode mismatch refuses decomposition */
    {
        if (!astOK) astClearStatus;
        AstZoomMap *z1 = astZoomMap(1, 2.0, "");
        AstZoomMap *z2 = astZoomMap(1, 3.0, "");
        AstCmpMap *series_inner = astCmpMap(z1, z2, 1, "");
        AstZoomMap *z3 = astZoomMap(1, 5.0, "");
        AstCmpMap *par = astCmpMap(series_inner, z3, 0, "");
        write_negative_fixture(dir, "neg_cmpmap_mode_mismatch", (AstMapping*)par);
        par = astAnnul(par); series_inner = astAnnul(series_inner);
        z1 = astAnnul(z1); z2 = astAnnul(z2); z3 = astAnnul(z3);
    }

    /* wcsmap-07: WcsMap adjacent to different projection type — refuses merge */
    {
        if (!astOK) astClearStatus;
        AstWcsMap *w1 = astWcsMap(2, AST__TAN, 1, 2, "");
        AstWcsMap *w2 = astWcsMap(2, AST__SIN, 1, 2, "");
        astInvert(w2);
        AstCmpMap *cm = astCmpMap(w1, w2, 1, "");
        write_negative_fixture(dir, "neg_wcs_different_projection", (AstMapping*)cm);
        cm = astAnnul(cm); w1 = astAnnul(w1); w2 = astAnnul(w2);
    }

    /* tranmap-07: TranMap with unequal components — not simplified */
    {
        if (!astOK) astClearStatus;
        AstZoomMap *z1 = astZoomMap(1, 2.0, "");
        AstZoomMap *z2 = astZoomMap(1, 3.0, "");
        AstTranMap *tm = astTranMap(z1, z2, "");
        write_negative_fixture(dir, "neg_tranmap_unequal_components", (AstMapping*)tm);
        tm = astAnnul(tm); z1 = astAnnul(z1); z2 = astAnnul(z2);
    }

    /* splinemap-06: two SplineMaps same direction — refuses cancel */
    {
        if (!astOK) astClearStatus;
        double tx[] = {0, 1}, ty[] = {0, 1};
        double cu[] = {1.5}, cv[] = {2.5};
        AstSplineMap *s1 = astSplineMap(1, 1, 1, 1, tx, ty, cu, cv, "");
        AstSplineMap *s2 = astSplineMap(1, 1, 1, 1, tx, ty, cu, cv, "");
        AstCmpMap *cm = astCmpMap(s1, s2, 1, "");
        write_negative_fixture(dir, "neg_spline_same_direction", (AstMapping*)cm);
        cm = astAnnul(cm); s1 = astAnnul(s1); s2 = astAnnul(s2);
    }

    /* normmap-11: NormMap in parallel — no simplification */
    {
        if (!astOK) astClearStatus;
        AstSkyFrame *sf = astSkyFrame("");
        AstNormMap *nm = astNormMap(sf, "");
        AstZoomMap *zm = astZoomMap(1, 2.0, "");
        AstCmpMap *cm = astCmpMap(nm, zm, 0, "");
        write_negative_fixture(dir, "neg_normmap_parallel", (AstMapping*)cm);
        cm = astAnnul(cm); nm = astAnnul(nm); zm = astAnnul(zm);
        sf = astAnnul(sf);
    }

    /* permmap-08: lone canonical PermMap with no mergeable neighbours */
    {
        if (!astOK) astClearStatus;
        int inperm[] = {2, 1};
        int outperm[] = {2, 1};
        const char *fwd[] = {"y1 = x1", "y2 = x2"};
        const char *inv[] = {"x1 = y1", "x2 = y2"};
        AstPermMap *pm = astPermMap(2, inperm, 2, outperm, NULL, "");
        AstMathMap *mm = astMathMap(2, 2, 2, fwd, 2, inv, "");
        AstCmpMap *cm = astCmpMap(pm, mm, 1, "");
        write_negative_fixture(dir, "neg_perm_no_merge", (AstMapping*)cm);
        cm = astAnnul(cm); pm = astAnnul(pm); mm = astAnnul(mm);
    }
}

/* ===== Negative fixtures batch 3 ===== */

static void gen_negative_fixtures_3(const char *dir) {
    if (!astOK) astClearStatus;
    printf("Negative fixtures batch 3:\n");

    /* winmap-16: WinMap adjacent to series CmpMap (not parallel) — no merge */
    {
        if (!astOK) astClearStatus;
        double ina[] = {0, 0}, inb[] = {1, 1};
        double outa[] = {1, 2}, outb[] = {4, 6};
        AstWinMap *wm = astWinMap(2, ina, inb, outa, outb, "");
        AstZoomMap *z1 = astZoomMap(2, 2.0, "");
        AstZoomMap *z2 = astZoomMap(2, 3.0, "");
        AstCmpMap *series = astCmpMap(z1, z2, 1, "");
        AstCmpMap *cm = astCmpMap(series, wm, 1, "");
        write_negative_fixture(dir, "neg_win_series_cmpmap_neighbour", (AstMapping*)cm);
        cm = astAnnul(cm); series = astAnnul(series);
        wm = astAnnul(wm); z1 = astAnnul(z1); z2 = astAnnul(z2);
    }

    /* winmap-28: WinMap swap refused — neither swapped Mapping simplifies */
    {
        if (!astOK) astClearStatus;
        double ina[] = {0, 0}, inb[] = {1, 1};
        double outa[] = {1, 2}, outb[] = {4, 6};
        double mat[] = {1.0, 2.0, 3.0, 4.0};
        AstWinMap *wm = astWinMap(2, ina, inb, outa, outb, "");
        AstMatrixMap *mm = astMatrixMap(2, 2, 0, mat, "");
        AstCmpMap *cm = astCmpMap(wm, mm, 1, "");
        write_negative_fixture(dir, "neg_win_swap_no_simplify", (AstMapping*)cm);
        cm = astAnnul(cm); wm = astAnnul(wm); mm = astAnnul(mm);
    }

    /* matrixmap-10: PermMap not bidirectional — direct merge blocked */
    {
        if (!astOK) astClearStatus;
        double diag[] = {2.0, 3.0};
        int inperm[] = {1, -1};
        int outperm[] = {1, -1};
        double consts[] = {5.0, 7.0};
        AstMatrixMap *mm = astMatrixMap(2, 2, 1, diag, "");
        AstPermMap *pm = astPermMap(2, inperm, 2, outperm, consts, "");
        AstCmpMap *cm = astCmpMap(mm, pm, 1, "");
        write_negative_fixture(dir, "neg_matrix_perm_not_bidirectional", (AstMapping*)cm);
        cm = astAnnul(cm); mm = astAnnul(mm); pm = astAnnul(pm);
    }

    /* ratemap-09: RateMap where upper neighbour is not a RateMap */
    {
        if (!astOK) astClearStatus;
        AstZoomMap *zm = astZoomMap(1, 2.0, "");
        AstRateMap *rm = astRateMap(zm, 1, 1, "");
        AstZoomMap *z2 = astZoomMap(1, 3.0, "");
        AstCmpMap *cm = astCmpMap(rm, z2, 1, "");
        write_negative_fixture(dir, "neg_ratemap_nonratemap_neighbour", (AstMapping*)cm);
        cm = astAnnul(cm); rm = astAnnul(rm); z2 = astAnnul(z2);
        zm = astAnnul(zm);
    }

    /* ratemap-08: RateMaps with different encapsulated Mappings */
    {
        if (!astOK) astClearStatus;
        AstZoomMap *z1 = astZoomMap(2, 2.0, "");
        AstZoomMap *z2 = astZoomMap(2, 3.0, "");
        AstRateMap *r1 = astRateMap(z1, 1, 1, "");
        AstRateMap *r2 = astRateMap(z2, 1, 1, "");
        astInvert(r2);
        AstCmpMap *cm = astCmpMap(r1, r2, 1, "");
        write_negative_fixture(dir, "neg_ratemap_different_inner", (AstMapping*)cm);
        cm = astAnnul(cm); r1 = astAnnul(r1); r2 = astAnnul(r2);
        z1 = astAnnul(z1); z2 = astAnnul(z2);
    }

    /* cmpmap-06: CmpMap neighbour is not a CmpMap — refuses merge */
    {
        if (!astOK) astClearStatus;
        AstZoomMap *z1 = astZoomMap(1, 2.0, "");
        AstZoomMap *z2 = astZoomMap(1, 3.0, "");
        AstCmpMap *par = astCmpMap(z1, z2, 0, "");
        const char *fwd[] = {"y1 = x1", "y2 = x2"};
        const char *inv[] = {"x1 = y1", "x2 = y2"};
        AstMathMap *mm = astMathMap(2, 2, 2, fwd, 2, inv, "");
        AstCmpMap *cm = astCmpMap(par, mm, 1, "");
        write_negative_fixture(dir, "neg_cmpmap_neighbour_not_cmpmap", (AstMapping*)cm);
        cm = astAnnul(cm); par = astAnnul(par); mm = astAnnul(mm);
        z1 = astAnnul(z1); z2 = astAnnul(z2);
    }

    /* sphmap-03: Inverse(SphMap) + SphMap but PolarLong differs */
    {
        if (!astOK) astClearStatus;
        AstSphMap *s1 = astSphMap("PolarLong=0");
        astInvert(s1);
        AstSphMap *s2 = astSphMap("PolarLong=3.14159");
        AstCmpMap *cm = astCmpMap(s1, s2, 1, "");
        write_negative_fixture(dir, "neg_sph_polarlong_mismatch", (AstMapping*)cm);
        cm = astAnnul(cm); s1 = astAnnul(s1); s2 = astAnnul(s2);
    }

    /* sphmap-04: SphMap(fwd) + Inverse(SphMap) but UnitRadius not set */
    {
        if (!astOK) astClearStatus;
        AstSphMap *s1 = astSphMap("");
        AstSphMap *s2 = astSphMap("");
        astInvert(s2);
        AstCmpMap *cm = astCmpMap(s1, s2, 1, "");
        write_negative_fixture(dir, "neg_sph_no_unitradius", (AstMapping*)cm);
        cm = astAnnul(cm); s1 = astAnnul(s1); s2 = astAnnul(s2);
    }

    /* polymap-10: two different PolyMaps in opposite directions — astEqual fails */
    {
        if (!astOK) astClearStatus;
        double coeff_f1[] = {2.0, 1, 1};
        double coeff_f2[] = {3.0, 1, 1};
        AstPolyMap *p1 = astPolyMap(1, 1, 1, coeff_f1, 0, NULL, "");
        AstPolyMap *p2 = astPolyMap(1, 1, 1, coeff_f2, 0, NULL, "");
        astInvert(p2);
        AstCmpMap *cm = astCmpMap(p1, p2, 1, "");
        write_negative_fixture(dir, "neg_poly_different_coeffs", (AstMapping*)cm);
        cm = astAnnul(cm); p1 = astAnnul(p1); p2 = astAnnul(p2);
    }

    /* lutmap-07: two LutMaps in parallel — cancellation not attempted */
    {
        if (!astOK) astClearStatus;
        double lut[] = {0.0, 1.0, 4.0, 9.0, 16.0};
        AstLutMap *l1 = astLutMap(5, lut, 1.0, 1.0, "");
        AstLutMap *l2 = astLutMap(5, lut, 1.0, 1.0, "");
        AstCmpMap *cm = astCmpMap(l1, l2, 0, "");
        write_negative_fixture(dir, "neg_lut_parallel", (AstMapping*)cm);
        cm = astAnnul(cm); l1 = astAnnul(l1); l2 = astAnnul(l2);
    }
}

/* ===== Negative fixtures batch 4 ===== */

static void gen_negative_fixtures_4(const char *dir) {
    if (!astOK) astClearStatus;
    printf("Negative fixtures batch 4:\n");

    /* slamap-21: 2-arg pair with mismatched arguments (AMP+MAP) */
    {
        if (!astOK) astClearStatus;
        double amp_args[] = {51544.0, 2000.0};
        double map_args[] = {2001.0, 51544.0};
        AstSlaMap *sm = astSlaMap(0, "");
        astSlaAdd(sm, "AMP", 2, amp_args);
        astSlaAdd(sm, "MAP", 2, map_args);
        write_negative_fixture(dir, "neg_sla_2arg_mismatch", (AstMapping*)sm);
        sm = astAnnul(sm);
    }

    /* specmap-14: 2-arg pair with mismatched arguments */
    {
        if (!astOK) astClearStatus;
        double args1[] = {0.5, 1.2, 0.0, 0.0};
        double args2[] = {0.5, 1.3, 0.0, 0.0};
        AstSpecMap *sm = astSpecMap(1, 0, "");
        astSpecAdd(sm, "LKF2HL", 2, args1);
        astSpecAdd(sm, "HLF2LK", 2, args2);
        write_negative_fixture(dir, "neg_spec_2arg_mismatch", (AstMapping*)sm);
        sm = astAnnul(sm);
    }

    /* specmap-15: 3-arg pair with mismatched arguments */
    {
        if (!astOK) astClearStatus;
        double args1[] = {0.5, 1.2, 51544.0, 0.0, 0.0};
        double args2[] = {0.5, 1.2, 51545.0, 0.0, 0.0};
        AstSpecMap *sm = astSpecMap(1, 0, "");
        astSpecAdd(sm, "GEF2HL", 3, args1);
        astSpecAdd(sm, "HLF2GE", 3, args2);
        write_negative_fixture(dir, "neg_spec_3arg_mismatch", (AstMapping*)sm);
        sm = astAnnul(sm);
    }

    /* timemap-14: 2-arg swapped pair with mismatched arguments */
    {
        if (!astOK) astClearStatus;
        double args1[] = {0.0, 2400000.5, 0.0};
        double args2[] = {2400001.0, 0.0, 0.0};
        AstTimeMap *tm = astTimeMap(0, "");
        astTimeAdd(tm, "MJDTOJD", 2, args1);
        astTimeAdd(tm, "JDTOMJD", 2, args2);
        write_negative_fixture(dir, "neg_time_2arg_mismatch", (AstMapping*)tm);
        tm = astAnnul(tm);
    }

    /* timemap-15: 3-arg pair with mismatched arguments */
    {
        if (!astOK) astClearStatus;
        double args1[] = {-2.5, 1000.0, 1013.0, 0.0};
        double args2[] = {-2.6, 1000.0, 1013.0, 0.0};
        AstTimeMap *tm = astTimeMap(0, "");
        astTimeAdd(tm, "GMSTTOLMST", 3, args1);
        astTimeAdd(tm, "LMSTTOGMST", 3, args2);
        write_negative_fixture(dir, "neg_time_3arg_mismatch", (AstMapping*)tm);
        tm = astAnnul(tm);
    }

    /* grismmap-06: two GrismMaps with different attributes */
    {
        if (!astOK) astClearStatus;
        AstGrismMap *g1 = astGrismMap("GrismNR=1.0");
        AstGrismMap *g2 = astGrismMap("GrismNR=1.5");
        astInvert(g2);
        AstCmpMap *cm = astCmpMap(g1, g2, 1, "");
        write_negative_fixture(dir, "neg_grism_different_attrs", (AstMapping*)cm);
        cm = astAnnul(cm); g1 = astAnnul(g1); g2 = astAnnul(g2);
    }

    /* normmap-07: inverse NormMap pair with different Frames */
    {
        if (!astOK) astClearStatus;
        AstSkyFrame *sf1 = astSkyFrame("System=FK5");
        AstSkyFrame *sf2 = astSkyFrame("System=FK4");
        AstNormMap *n1 = astNormMap(sf1, "");
        AstNormMap *n2 = astNormMap(sf2, "");
        astInvert(n1);
        AstCmpMap *cm = astCmpMap(n1, n2, 1, "");
        write_negative_fixture(dir, "neg_normmap_different_frames", (AstMapping*)cm);
        cm = astAnnul(cm); n1 = astAnnul(n1); n2 = astAnnul(n2);
        sf1 = astAnnul(sf1); sf2 = astAnnul(sf2);
    }

    /* unitnormmap-15: two forward UnitNormMaps — same direction refuses */
    {
        if (!astOK) astClearStatus;
        double c1[] = {1.0, 2.0};
        double c2[] = {3.0, 4.0};
        AstUnitNormMap *u1 = astUnitNormMap(2, c1, "");
        AstUnitNormMap *u2 = astUnitNormMap(2, c2, "");
        AstCmpMap *cm = astCmpMap(u1, u2, 1, "");
        write_negative_fixture(dir, "neg_unitnormmap_same_direction", (AstMapping*)cm);
        cm = astAnnul(cm); u1 = astAnnul(u1); u2 = astAnnul(u2);
    }

    /* unitnormmap-13: forward UnitNormMap + ShiftMap — wrong order refuses */
    {
        if (!astOK) astClearStatus;
        double centre[] = {1.0, 2.0};
        double shifts[] = {0.5, 0.5};
        AstUnitNormMap *unm = astUnitNormMap(2, centre, "");
        AstShiftMap *sm = astShiftMap(2, shifts, "");
        AstCmpMap *cm = astCmpMap(unm, sm, 1, "");
        write_negative_fixture(dir, "neg_unitnormmap_fwd_then_shift", (AstMapping*)cm);
        cm = astAnnul(cm); unm = astAnnul(unm); sm = astAnnul(sm);
    }

    /* wcsmap-10: two WcsMaps with different projection parameters */
    {
        if (!astOK) astClearStatus;
        AstWcsMap *w1 = astWcsMap(2, AST__TAN, 1, 2, "PV1_1=0.5");
        AstWcsMap *w2 = astWcsMap(2, AST__TAN, 1, 2, "PV1_1=1.0");
        astInvert(w2);
        AstCmpMap *cm = astCmpMap(w1, w2, 1, "");
        write_negative_fixture(dir, "neg_wcs_different_params", (AstMapping*)cm);
        cm = astAnnul(cm); w1 = astAnnul(w1); w2 = astAnnul(w2);
    }

    /* slamap-06: single forward SlaMap with no neighbours — nothing simplifies */
    {
        if (!astOK) astClearStatus;
        AstSlaMap *sm = astSlaMap(0, "");
        astSlaAdd(sm, "EQGAL", 0, NULL);
        write_negative_fixture(dir, "neg_sla_lone_forward", (AstMapping*)sm);
        sm = astAnnul(sm);
    }
}

/* ===== Negative fixtures batch 5 ===== */

static void gen_negative_fixtures_5(const char *dir) {
    if (!astOK) astClearStatus;
    printf("Negative fixtures batch 5:\n");

    /* slamap-05: SlaMap in parallel — refuses simplification */
    {
        if (!astOK) astClearStatus;
        AstSlaMap *sm = astSlaMap(0, "");
        astSlaAdd(sm, "EQGAL", 0, NULL);
        AstZoomMap *zm = astZoomMap(1, 2.0, "");
        AstCmpMap *cm = astCmpMap(sm, zm, 0, "");
        write_negative_fixture(dir, "neg_sla_parallel", (AstMapping*)cm);
        cm = astAnnul(cm); sm = astAnnul(sm); zm = astAnnul(zm);
    }

    /* specmap-05: SpecMap in parallel — refuses simplification */
    {
        if (!astOK) astClearStatus;
        double rf[] = {1.4e9, 0.0};
        AstSpecMap *sm = astSpecMap(1, 0, "");
        astSpecAdd(sm, "FRTOVL", 1, rf);
        AstZoomMap *zm = astZoomMap(1, 2.0, "");
        AstCmpMap *cm = astCmpMap(sm, zm, 0, "");
        write_negative_fixture(dir, "neg_spec_parallel", (AstMapping*)cm);
        cm = astAnnul(cm); sm = astAnnul(sm); zm = astAnnul(zm);
    }

    /* timemap-05: TimeMap in parallel — refuses simplification */
    {
        if (!astOK) astClearStatus;
        double dut[] = {0.5, 0.0};
        AstTimeMap *tm = astTimeMap(0, "");
        astTimeAdd(tm, "TAITOTT", 1, dut);
        AstZoomMap *zm = astZoomMap(1, 2.0, "");
        AstCmpMap *cm = astCmpMap(tm, zm, 0, "");
        write_negative_fixture(dir, "neg_time_parallel", (AstMapping*)cm);
        cm = astAnnul(cm); tm = astAnnul(tm); zm = astAnnul(zm);
    }

    /* grismmap-04: GrismMap in parallel — refuses */
    {
        if (!astOK) astClearStatus;
        AstGrismMap *gm = astGrismMap("");
        AstZoomMap *zm = astZoomMap(1, 2.0, "");
        AstCmpMap *cm = astCmpMap(gm, zm, 0, "");
        write_negative_fixture(dir, "neg_grism_parallel", (AstMapping*)cm);
        cm = astAnnul(cm); gm = astAnnul(gm); zm = astAnnul(zm);
    }

    /* grismmap-09: inverted GrismMap followed by ZoomMap — wrong order */
    {
        if (!astOK) astClearStatus;
        AstGrismMap *gm = astGrismMap("");
        astInvert(gm);
        AstZoomMap *zm = astZoomMap(1, 2.0, "");
        AstCmpMap *cm = astCmpMap(gm, zm, 1, "");
        write_negative_fixture(dir, "neg_grism_inv_then_zoom", (AstMapping*)cm);
        cm = astAnnul(cm); gm = astAnnul(gm); zm = astAnnul(zm);
    }

    /* grismmap-10: forward GrismMap followed by non-ZoomMap */
    {
        if (!astOK) astClearStatus;
        AstGrismMap *gm = astGrismMap("");
        double shifts[] = {1.0};
        AstShiftMap *sm = astShiftMap(1, shifts, "");
        AstCmpMap *cm = astCmpMap(gm, sm, 1, "");
        write_negative_fixture(dir, "neg_grism_fwd_then_nonzoom", (AstMapping*)cm);
        cm = astAnnul(cm); gm = astAnnul(gm); sm = astAnnul(sm);
    }

    /* mathmap-04: MathMap followed by non-MathMap */
    {
        if (!astOK) astClearStatus;
        const char *fwd[] = {"y = x"};
        const char *inv[] = {"x = y"};
        AstMathMap *mm = astMathMap(1, 1, 1, fwd, 1, inv, "SimpFI=1");
        AstZoomMap *zm = astZoomMap(1, 2.0, "");
        AstCmpMap *cm = astCmpMap(mm, zm, 1, "");
        write_negative_fixture(dir, "neg_math_nonmath_neighbour", (AstMapping*)cm);
        cm = astAnnul(cm); mm = astAnnul(mm); zm = astAnnul(zm);
    }

    /* wcsmap-06: WcsMap in parallel — refuses */
    {
        if (!astOK) astClearStatus;
        AstWcsMap *wm = astWcsMap(2, AST__TAN, 1, 2, "");
        AstZoomMap *zm = astZoomMap(1, 2.0, "");
        AstCmpMap *cm = astCmpMap(wm, zm, 0, "");
        write_negative_fixture(dir, "neg_wcs_parallel", (AstMapping*)cm);
        cm = astAnnul(cm); wm = astAnnul(wm); zm = astAnnul(zm);
    }

    /* polymap-07: two non-linear PolyMaps in parallel — refuses cancel */
    {
        if (!astOK) astClearStatus;
        double coeff_f[] = {1.0, 1, 2};
        AstPolyMap *p1 = astPolyMap(1, 1, 1, coeff_f, 0, NULL, "");
        AstPolyMap *p2 = astPolyMap(1, 1, 1, coeff_f, 0, NULL, "");
        astInvert(p2);
        AstCmpMap *cm = astCmpMap(p1, p2, 0, "");
        write_negative_fixture(dir, "neg_poly_parallel_nonlinear", (AstMapping*)cm);
        cm = astAnnul(cm); p1 = astAnnul(p1); p2 = astAnnul(p2);
    }

    /* splinemap-03: SplineMap in parallel — refuses */
    {
        if (!astOK) astClearStatus;
        double tx[] = {0, 1}, ty[] = {0, 1};
        double cu[] = {1.5}, cv[] = {2.5};
        AstSplineMap *sm = astSplineMap(1, 1, 1, 1, tx, ty, cu, cv, "");
        AstZoomMap *zm = astZoomMap(1, 2.0, "");
        AstCmpMap *cm = astCmpMap(sm, zm, 0, "");
        write_negative_fixture(dir, "neg_spline_parallel", (AstMapping*)cm);
        cm = astAnnul(cm); sm = astAnnul(sm); zm = astAnnul(zm);
    }

    /* lutmap-08: LutMap with non-LutMap neighbour — refuses cancel */
    {
        if (!astOK) astClearStatus;
        double lut[] = {0.0, 1.0, 4.0, 9.0, 16.0};
        AstLutMap *lm = astLutMap(5, lut, 1.0, 1.0, "");
        AstZoomMap *zm = astZoomMap(1, 2.0, "");
        AstCmpMap *cm = astCmpMap(lm, zm, 1, "");
        write_negative_fixture(dir, "neg_lut_nonlut_neighbour", (AstMapping*)cm);
        cm = astAnnul(cm); lm = astAnnul(lm); zm = astAnnul(zm);
    }
}

/* ===== Negative fixtures batch 6 ===== */

static void gen_negative_fixtures_6(const char *dir) {
    if (!astOK) astClearStatus;
    printf("Negative fixtures batch 6:\n");

    /* specmap-06: lone forward SpecMap — nothing simplifies */
    {
        if (!astOK) astClearStatus;
        double rf[] = {1.4e9, 0.0};
        AstSpecMap *sm = astSpecMap(1, 0, "");
        astSpecAdd(sm, "FRTOVL", 1, rf);
        write_negative_fixture(dir, "neg_spec_lone_forward", (AstMapping*)sm);
        sm = astAnnul(sm);
    }

    /* timemap-06: lone forward TimeMap — nothing simplifies */
    {
        if (!astOK) astClearStatus;
        double dut[] = {0.5, 0.0};
        AstTimeMap *tm = astTimeMap(0, "");
        astTimeAdd(tm, "TAITOTT", 1, dut);
        write_negative_fixture(dir, "neg_time_lone_forward", (AstMapping*)tm);
        tm = astAnnul(tm);
    }

    /* grismmap-13: GrismMap + ZoomMap(0) — zero zoom prevents merge */
    {
        if (!astOK) astClearStatus;
        AstGrismMap *gm = astGrismMap("");
        AstZoomMap *zm = astZoomMap(1, 0.0, "");
        AstCmpMap *cm = astCmpMap(gm, zm, 1, "");
        write_negative_fixture(dir, "neg_grism_zoom_zero", (AstMapping*)cm);
        cm = astAnnul(cm); gm = astAnnul(gm); zm = astAnnul(zm);
    }

    /* polymap-08: PolyMap neighbour is not PolyMap — refuses cancel */
    {
        if (!astOK) astClearStatus;
        double coeff_f[] = {1.0, 1, 2};
        AstPolyMap *pm = astPolyMap(1, 1, 1, coeff_f, 0, NULL, "");
        AstZoomMap *zm = astZoomMap(1, 2.0, "");
        AstCmpMap *cm = astCmpMap(pm, zm, 1, "");
        write_negative_fixture(dir, "neg_poly_nonpoly_neighbour", (AstMapping*)cm);
        cm = astAnnul(cm); pm = astAnnul(pm); zm = astAnnul(zm);
    }

    /* splinemap-05: SplineMap with non-SplineMap neighbour */
    {
        if (!astOK) astClearStatus;
        double tx[] = {0, 1}, ty[] = {0, 1};
        double cu[] = {1.5}, cv[] = {2.5};
        AstSplineMap *sm = astSplineMap(1, 1, 1, 1, tx, ty, cu, cv, "");
        AstZoomMap *zm = astZoomMap(2, 2.0, "");
        AstCmpMap *cm = astCmpMap(sm, zm, 1, "");
        write_negative_fixture(dir, "neg_spline_nonspline_neighbour", (AstMapping*)cm);
        cm = astAnnul(cm); sm = astAnnul(sm); zm = astAnnul(zm);
    }

    /* matrixmap-16: MatrixMap swap with WinMap refused (neither simplifies) —
       already covered by neg_win_swap_no_simplify from the WinMap side.
       Test from MatrixMap perspective with different structure. */
    {
        if (!astOK) astClearStatus;
        double mat[] = {1.0, 2.0, 3.0, 4.0};
        double ina[] = {0, 0}, inb[] = {1, 1};
        double outa[] = {1, 2}, outb[] = {3, 5};
        AstMatrixMap *mm = astMatrixMap(2, 2, 0, mat, "");
        AstWinMap *wm = astWinMap(2, ina, inb, outa, outb, "");
        AstCmpMap *cm = astCmpMap(mm, wm, 1, "");
        write_negative_fixture(dir, "neg_matrix_swap_refused", (AstMapping*)cm);
        cm = astAnnul(cm); mm = astAnnul(mm); wm = astAnnul(wm);
    }
}

/* ===== Negative fixtures batch 7 ===== */

static void gen_negative_fixtures_7(const char *dir) {
    if (!astOK) astClearStatus;
    printf("Negative fixtures batch 7:\n");

    /* polymap-04: PolyMap with nin != nout — linearization refused.
       2 inputs, 1 output, linear terms — but nin(2)!=nout(1) blocks linearize. */
    {
        if (!astOK) astClearStatus;
        double coeff_f[] = {3.0, 1, 1, 0,
                            5.0, 1, 0, 1};
        AstPolyMap *pm = astPolyMap(2, 1, 2, coeff_f, 0, NULL, "");
        write_negative_fixture(dir, "neg_poly_nin_ne_nout", (AstMapping*)pm);
        pm = astAnnul(pm);
    }

    /* tranmap-05: TranMaps in parallel — refuses adjacent merge */
    {
        if (!astOK) astClearStatus;
        AstZoomMap *z1 = astZoomMap(1, 2.0, "");
        AstZoomMap *z2 = astZoomMap(1, 0.5, "");
        AstZoomMap *z3 = astZoomMap(1, 3.0, "");
        AstZoomMap *z4 = astZoomMap(1, 1.0/3.0, "");
        AstTranMap *t1 = astTranMap(z1, z2, "");
        AstTranMap *t2 = astTranMap(z3, z4, "");
        AstCmpMap *cm = astCmpMap(t1, t2, 0, "");
        write_negative_fixture(dir, "neg_tranmap_parallel", (AstMapping*)cm);
        cm = astAnnul(cm); t1 = astAnnul(t1); t2 = astAnnul(t2);
        z1 = astAnnul(z1); z2 = astAnnul(z2);
        z3 = astAnnul(z3); z4 = astAnnul(z4);
    }

    /* tranmap-08: TranMap where higher neighbour is not a TranMap */
    {
        if (!astOK) astClearStatus;
        AstZoomMap *z1 = astZoomMap(1, 2.0, "");
        AstZoomMap *z2 = astZoomMap(1, 0.5, "");
        AstTranMap *tm = astTranMap(z1, z2, "");
        AstZoomMap *z3 = astZoomMap(1, 5.0, "");
        AstCmpMap *cm = astCmpMap(tm, z3, 1, "");
        write_negative_fixture(dir, "neg_tranmap_nontranmap_neighbour", (AstMapping*)cm);
        cm = astAnnul(cm); tm = astAnnul(tm);
        z1 = astAnnul(z1); z2 = astAnnul(z2); z3 = astAnnul(z3);
    }

    /* sphmap-11: SphMap sandwich middle is not ZoomMap/diagonal MatrixMap */
    {
        if (!astOK) astClearStatus;
        double shifts[] = {1.0, 2.0, 3.0};
        AstSphMap *s1 = astSphMap("UnitRadius=1");
        astInvert(s1);
        AstShiftMap *sm = astShiftMap(3, shifts, "");
        AstSphMap *s2 = astSphMap("UnitRadius=1");
        AstCmpMap *inner = astCmpMap(sm, s2, 1, "");
        AstCmpMap *outer = astCmpMap(s1, inner, 1, "");
        write_negative_fixture(dir, "neg_sph_sandwich_wrong_middle", (AstMapping*)outer);
        outer = astAnnul(outer); inner = astAnnul(inner);
        s1 = astAnnul(s1); sm = astAnnul(sm); s2 = astAnnul(s2);
    }

    /* sphmap-14: SphMap sandwich MatrixMap with unequal magnitude diagonals */
    {
        if (!astOK) astClearStatus;
        double diag[] = {1.0, 2.0, 3.0};
        AstSphMap *s1 = astSphMap("UnitRadius=1");
        astInvert(s1);
        AstMatrixMap *mm = astMatrixMap(3, 3, 1, diag, "");
        AstSphMap *s2 = astSphMap("UnitRadius=1");
        AstCmpMap *inner = astCmpMap(mm, s2, 1, "");
        AstCmpMap *outer = astCmpMap(s1, inner, 1, "");
        write_negative_fixture(dir, "neg_sph_sandwich_unequal_diag", (AstMapping*)outer);
        outer = astAnnul(outer); inner = astAnnul(inner);
        s1 = astAnnul(s1); mm = astAnnul(mm); s2 = astAnnul(s2);
    }

    /* pcdmap-10: PcdMap intervening mapping not ZoomMap/PermMap — blocks swap */
    {
        if (!astOK) astClearStatus;
        double cen[] = {0.0, 0.0};
        double shifts[] = {1.0, 2.0};
        AstPcdMap *p1 = astPcdMap(0.001, cen, "");
        AstShiftMap *sm = astShiftMap(2, shifts, "");
        AstPcdMap *p2 = astPcdMap(0.001, cen, "");
        astInvert(p2);
        AstCmpMap *inner = astCmpMap(sm, p2, 1, "");
        AstCmpMap *outer = astCmpMap(p1, inner, 1, "");
        write_negative_fixture(dir, "neg_pcd_nonswappable_between", (AstMapping*)outer);
        outer = astAnnul(outer); inner = astAnnul(inner);
        p1 = astAnnul(p1); sm = astAnnul(sm); p2 = astAnnul(p2);
    }

    /* wcsmap-11: WcsMap with non-PermMap intervening — blocks swap.
       Use a MathMap as the non-swappable blocker (non-zero shifts would
       self-simplify away, and ShiftMap(0,0) becomes UnitMap). */
    {
        if (!astOK) astClearStatus;
        const char *fwd[] = {"y1 = x1", "y2 = x2"};
        const char *inv[] = {"x1 = y1", "x2 = y2"};
        AstWcsMap *w1 = astWcsMap(2, AST__TAN, 1, 2, "");
        AstMathMap *mm = astMathMap(2, 2, 2, fwd, 2, inv, "");
        AstWcsMap *w2 = astWcsMap(2, AST__TAN, 1, 2, "");
        astInvert(w2);
        AstCmpMap *inner = astCmpMap(mm, w2, 1, "");
        AstCmpMap *outer = astCmpMap(w1, inner, 1, "");
        write_negative_fixture(dir, "neg_wcs_nonperm_between", (AstMapping*)outer);
        outer = astAnnul(outer); inner = astAnnul(inner);
        w1 = astAnnul(w1); mm = astAnnul(mm); w2 = astAnnul(w2);
    }
}

/* ===== Negative fixtures batch 8 ===== */

static void gen_negative_fixtures_8(const char *dir) {
    if (!astOK) astClearStatus;
    printf("Negative fixtures batch 8:\n");

    /* sphmap-06: SphMap in parallel — refuses */
    {
        if (!astOK) astClearStatus;
        AstSphMap *sm = astSphMap("UnitRadius=1");
        AstZoomMap *zm = astZoomMap(1, 2.0, "");
        AstCmpMap *cm = astCmpMap(sm, zm, 0, "");
        write_negative_fixture(dir, "neg_sph_parallel", (AstMapping*)cm);
        cm = astAnnul(cm); sm = astAnnul(sm); zm = astAnnul(zm);
    }

    /* sphmap-13: sandwich MatrixMap not diagonal (full matrix) */
    {
        if (!astOK) astClearStatus;
        double mat[] = {1.0, 0.5, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
        AstSphMap *s1 = astSphMap("UnitRadius=1");
        astInvert(s1);
        AstMatrixMap *mm = astMatrixMap(3, 3, 0, mat, "");
        AstSphMap *s2 = astSphMap("UnitRadius=1");
        AstCmpMap *inner = astCmpMap(mm, s2, 1, "");
        AstCmpMap *outer = astCmpMap(s1, inner, 1, "");
        write_negative_fixture(dir, "neg_sph_sandwich_full_matrix", (AstMapping*)outer);
        outer = astAnnul(outer); inner = astAnnul(inner);
        s1 = astAnnul(s1); mm = astAnnul(mm); s2 = astAnnul(s2);
    }

    /* normmap-08: NormMap followed by non-NormMap — no cancel */
    {
        if (!astOK) astClearStatus;
        AstSkyFrame *sf = astSkyFrame("");
        AstNormMap *nm = astNormMap(sf, "");
        AstZoomMap *zm = astZoomMap(2, 2.0, "");
        AstCmpMap *cm = astCmpMap(nm, zm, 1, "");
        write_negative_fixture(dir, "neg_normmap_nonnorm_neighbour", (AstMapping*)cm);
        cm = astAnnul(cm); nm = astAnnul(nm); zm = astAnnul(zm);
        sf = astAnnul(sf);
    }

    /* grismmap-11: ZoomMap before forward (non-inverted) GrismMap — wrong config */
    {
        if (!astOK) astClearStatus;
        AstZoomMap *zm = astZoomMap(1, 2.0, "");
        AstGrismMap *gm = astGrismMap("");
        AstCmpMap *cm = astCmpMap(zm, gm, 1, "");
        write_negative_fixture(dir, "neg_grism_zoom_before_fwd", (AstMapping*)cm);
        cm = astAnnul(cm); zm = astAnnul(zm); gm = astAnnul(gm);
    }

    /* grismmap-12: non-ZoomMap before inverted GrismMap */
    {
        if (!astOK) astClearStatus;
        double shifts[] = {1.0};
        AstShiftMap *sm = astShiftMap(1, shifts, "");
        AstGrismMap *gm = astGrismMap("");
        astInvert(gm);
        AstCmpMap *cm = astCmpMap(sm, gm, 1, "");
        write_negative_fixture(dir, "neg_grism_nonzoom_before_inv", (AstMapping*)cm);
        cm = astAnnul(cm); sm = astAnnul(sm); gm = astAnnul(gm);
    }

    /* unitnormmap-16: UnitNormMap with non-mergeable neighbour.
       UnitNormMap(2) forward: Nin=2, Nout=3. Follow with ZoomMap(3). */
    {
        if (!astOK) astClearStatus;
        double centre[] = {1.0, 2.0};
        AstUnitNormMap *unm = astUnitNormMap(2, centre, "");
        AstZoomMap *zm = astZoomMap(3, 2.0, "");
        AstCmpMap *cm = astCmpMap(unm, zm, 1, "");
        write_negative_fixture(dir, "neg_unitnormmap_nonmergeable", (AstMapping*)cm);
        cm = astAnnul(cm); unm = astAnnul(unm); zm = astAnnul(zm);
    }
}

/* ===== Positive cascade fixtures batch 2 ===== */

static void gen_cascade_positives_2(const char *dir) {
    if (!astOK) astClearStatus;
    printf("Cascade positives batch 2:\n");

    /* cmpmap-09: two parallel CmpMaps in series, components pair and simplify.
       Use MathMaps (non-simplifiable as parallel) so the parallel CmpMaps
       don't self-simplify before the pairing code runs.
       CmpMap(MathMap(2x)||MathMap(3x)) then CmpMap(MathMap(x/2)||MathMap(x/3))
       with SimpFI=1 so the pairs cancel. */
    {
        if (!astOK) astClearStatus;
        const char *fwd1[] = {"y = 2*x"};
        const char *inv1[] = {"x = 0.5*y"};
        const char *fwd2[] = {"y = 3*x"};
        const char *inv2[] = {"x = y/3"};
        AstMathMap *ma1 = astMathMap(1, 1, 1, fwd1, 1, inv1, "SimpFI=1,SimpIF=1");
        AstMathMap *mb1 = astMathMap(1, 1, 1, fwd2, 1, inv2, "SimpFI=1,SimpIF=1");
        AstCmpMap *par1 = astCmpMap(ma1, mb1, 0, "");

        AstMathMap *ma2 = astMathMap(1, 1, 1, fwd1, 1, inv1, "SimpFI=1,SimpIF=1");
        AstMathMap *mb2 = astMathMap(1, 1, 1, fwd2, 1, inv2, "SimpFI=1,SimpIF=1");
        astInvert(ma2);
        astInvert(mb2);
        AstCmpMap *par2 = astCmpMap(ma2, mb2, 0, "");

        AstCmpMap *cm = astCmpMap(par1, par2, 1, "");
        write_fixture(dir, "cmpmap_parallel_in_series_merge", (AstMapping*)cm);
        cm = astAnnul(cm); par1 = astAnnul(par1); par2 = astAnnul(par2);
        ma1 = astAnnul(ma1); mb1 = astAnnul(mb1);
        ma2 = astAnnul(ma2); mb2 = astAnnul(mb2);
    }

    /* cmpmap-18: PermMap swap with aconstants — first component gets constants.
       Use MathMaps (non-simplifiable) as CmpMap components so the parallel
       CmpMap doesn't self-simplify before the PermMap swap code runs.
       Square PermMap(2→2): output 1 = constant, output 2 = input 1. */
    {
        if (!astOK) astClearStatus;
        int inperm[] = {2, -1};
        int outperm[] = {-1, 1};
        double consts[] = {99.0};
        const char *fwda[] = {"y = 2*x"};
        const char *inva[] = {"x = 0.5*y"};
        const char *fwdb[] = {"y = 3*x"};
        const char *invb[] = {"x = y/3"};
        AstPermMap *pm = astPermMap(2, inperm, 2, outperm, consts, "");
        AstMathMap *ma = astMathMap(1, 1, 1, fwda, 1, inva, "");
        AstMathMap *mb = astMathMap(1, 1, 1, fwdb, 1, invb, "");
        AstCmpMap *par = astCmpMap(ma, mb, 0, "");
        AstCmpMap *cm = astCmpMap(pm, par, 1, "");
        write_fixture(dir, "cmpmap_perm_swap_aconstants", (AstMapping*)cm);
        cm = astAnnul(cm); pm = astAnnul(pm); par = astAnnul(par);
        ma = astAnnul(ma); mb = astAnnul(mb);
    }

    /* cmpmap-19: PermMap swap with bconstants — second component gets constants.
       For canswap: outperm[0]=from input at index npin-nin2a=1 (public: 2),
       outperm[1]<0 (constant), inperm[1]=from output 0 (public: 1),
       inperm[0]<0 (constant, avoids bidirectional check failure). */
    {
        if (!astOK) astClearStatus;
        int inperm[] = {-1, 1};
        int outperm[] = {2, -1};
        double consts[] = {99.0};
        const char *fwda[] = {"y = 2*x"};
        const char *inva[] = {"x = 0.5*y"};
        const char *fwdb[] = {"y = 3*x"};
        const char *invb[] = {"x = y/3"};
        AstPermMap *pm = astPermMap(2, inperm, 2, outperm, consts, "");
        AstMathMap *ma = astMathMap(1, 1, 1, fwda, 1, inva, "");
        AstMathMap *mb = astMathMap(1, 1, 1, fwdb, 1, invb, "");
        AstCmpMap *par = astCmpMap(ma, mb, 0, "");
        AstCmpMap *cm = astCmpMap(pm, par, 1, "");
        write_fixture(dir, "cmpmap_perm_swap_bconstants", (AstMapping*)cm);
        cm = astAnnul(cm); pm = astAnnul(pm); par = astAnnul(par);
        ma = astAnnul(ma); mb = astAnnul(mb);
    }

    /* cmpmap-03 (where>0): CmpMap decomposition when not first in list.
       ZoomMap || CmpMap(Zoom||Zoom) || ZoomMap in parallel — the inner
       CmpMap decomposes at where=1 (not position 0). */
    {
        if (!astOK) astClearStatus;
        AstZoomMap *z1 = astZoomMap(1, 2.0, "");
        AstZoomMap *z2 = astZoomMap(1, 3.0, "");
        AstZoomMap *z3 = astZoomMap(1, 4.0, "");
        AstCmpMap *inner = astCmpMap(z2, z3, 0, "");
        AstCmpMap *left = astCmpMap(z1, inner, 0, "");
        AstZoomMap *z4 = astZoomMap(1, 5.0, "");
        AstCmpMap *outer = astCmpMap(left, z4, 0, "");
        write_fixture(dir, "cmpmap_decompose_middle", (AstMapping*)outer);
        outer = astAnnul(outer); left = astAnnul(left);
        inner = astAnnul(inner);
        z1 = astAnnul(z1); z2 = astAnnul(z2);
        z3 = astAnnul(z3); z4 = astAnnul(z4);
    }

    /* cmpmap-07 with invert: series CmpMaps in parallel, one inverted.
       Use MathMaps so they don't self-simplify before the merge fires.
       Tests line 1567-1568 (invert1 flag flip). */
    {
        if (!astOK) astClearStatus;
        const char *fwd[] = {"y = 2*x"};
        const char *inv[] = {"x = 0.5*y"};
        AstMathMap *m1a = astMathMap(1, 1, 1, fwd, 1, inv, "SimpFI=1,SimpIF=1");
        AstMathMap *m1b = astMathMap(1, 1, 1, fwd, 1, inv, "SimpFI=1,SimpIF=1");
        AstCmpMap *ser1 = astCmpMap(m1a, m1b, 1, "");
        astInvert(ser1);

        AstMathMap *m2a = astMathMap(1, 1, 1, fwd, 1, inv, "SimpFI=1,SimpIF=1");
        AstMathMap *m2b = astMathMap(1, 1, 1, fwd, 1, inv, "SimpFI=1,SimpIF=1");
        AstCmpMap *ser2 = astCmpMap(m2a, m2b, 1, "");

        AstCmpMap *par = astCmpMap(ser1, ser2, 0, "");
        write_fixture(dir, "cmpmap_series_in_parallel_inverted", (AstMapping*)par);
        par = astAnnul(par); ser1 = astAnnul(ser1); ser2 = astAnnul(ser2);
        m1a = astAnnul(m1a); m1b = astAnnul(m1b);
        m2a = astAnnul(m2a); m2b = astAnnul(m2b);
    }

    /* winmap-15: WinMap merges with UPPER parallel CmpMap neighbour.
       Parallel CmpMap follows the WinMap in the series. */
    {
        if (!astOK) astClearStatus;
        double ina[] = {0, 0}, inb[] = {1, 1};
        double outa[] = {1, 2}, outb[] = {3, 5};
        AstWinMap *wm = astWinMap(2, ina, inb, outa, outb, "");
        double shifts1[] = {1.0};
        double shifts2[] = {2.0};
        AstShiftMap *sm1 = astShiftMap(1, shifts1, "");
        AstShiftMap *sm2 = astShiftMap(1, shifts2, "");
        AstCmpMap *par = astCmpMap(sm1, sm2, 0, "");
        AstCmpMap *cm = astCmpMap(wm, par, 1, "");
        write_fixture(dir, "win_upper_cmpmap_parallel_merge", (AstMapping*)cm);
        cm = astAnnul(cm); wm = astAnnul(wm); par = astAnnul(par);
        sm1 = astAnnul(sm1); sm2 = astAnnul(sm2);
    }

    /* permmap-09: constant propagation through series composition.
       First PermMap has a constant output, second PermMap routes that
       constant through — the combined PermMap propagates the constant. */
    {
        if (!astOK) astClearStatus;
        int in1[] = {1, -1};
        int out1[] = {1, 2};
        double c1[] = {42.0};
        AstPermMap *p1 = astPermMap(2, in1, 2, out1, c1, "");
        int in2[] = {2, 1};
        int out2[] = {2, 1};
        AstPermMap *p2 = astPermMap(2, in2, 2, out2, NULL, "");
        AstCmpMap *cm = astCmpMap(p1, p2, 1, "");
        write_fixture(dir, "perm_constant_propagation", (AstMapping*)cm);
        cm = astAnnul(cm); p1 = astAnnul(p1); p2 = astAnnul(p2);
    }

    /* unitnormmap-02: WinMap(unit scale) + forward UnitNormMap → adjusted centre */
    {
        if (!astOK) astClearStatus;
        double centre[] = {1.0, 2.0};
        double ina[] = {0, 0}, inb[] = {1, 1};
        double outa[] = {0.5, 0.5}, outb[] = {1.5, 1.5};
        AstWinMap *wm = astWinMap(2, ina, inb, outa, outb, "");
        AstUnitNormMap *unm = astUnitNormMap(2, centre, "");
        AstCmpMap *cm = astCmpMap(wm, unm, 1, "");
        write_fixture(dir, "unitnormmap_winmap_fwd_merge", (AstMapping*)cm);
        cm = astAnnul(cm); wm = astAnnul(wm); unm = astAnnul(unm);
    }

    /* unitnormmap-05: inverted UnitNormMap + WinMap(unit scale) → adjusted centre.
       Inverted UnitNormMap(2): Nin=3, Nout=2. WinMap must have Nin=2. */
    {
        if (!astOK) astClearStatus;
        double centre[] = {1.0, 2.0};
        double ina[] = {0, 0}, inb[] = {1, 1};
        double outa[] = {0.5, 0.5}, outb[] = {1.5, 1.5};
        AstUnitNormMap *unm = astUnitNormMap(2, centre, "");
        astInvert(unm);
        AstWinMap *wm = astWinMap(2, ina, inb, outa, outb, "");
        AstCmpMap *cm = astCmpMap(unm, wm, 1, "");
        write_fixture(dir, "unitnormmap_inv_winmap_merge", (AstMapping*)cm);
        cm = astAnnul(cm); unm = astAnnul(unm); wm = astAnnul(wm);
    }
}

/* ===== Positive cascade fixtures batch 3 ===== */

static void gen_cascade_positives_3(const char *dir) {
    if (!astOK) astClearStatus;
    printf("Cascade positives batch 3:\n");

    /* wcsmap-04: WcsMap + PermMap swap for local simplification.
       WcsMap(3,TAN,lonax=1,latax=2) has internal lonax=0, latax=1.
       PermMap(3→1) where inverse assigns constants to both celestial
       axes (inperm[0]<0,inperm[1]<0). After WcsPerm: the WcsMap
       becomes a UnitMap (class change triggers swap acceptance). */
    {
        if (!astOK) astClearStatus;
        int inperm[] = {-1, -2, 1};
        int outperm[] = {3};
        double consts[] = {0.5, 0.5};
        AstWcsMap *wm = astWcsMap(3, AST__TAN, 1, 2, "");
        AstPermMap *pm = astPermMap(3, inperm, 1, outperm, consts, "");
        AstCmpMap *cm = astCmpMap(wm, pm, 1, "");
        write_fixture(dir, "wcsmap_perm_swap_simplify", (AstMapping*)cm);
        cm = astAnnul(cm); wm = astAnnul(wm); pm = astAnnul(pm);
    }

    /* permmap-05: PermMap with redundant outperm simplifies after
       composition with UnitMap. A 3-axis PermMap that passes axes 1,2
       through and assigns axis 3 a constant — composed with a UnitMap(3)
       should merge and potentially simplify the stored arrays. */
    {
        if (!astOK) astClearStatus;
        int inperm[] = {1, 2, -1};
        int outperm[] = {1, 2, 0};
        double consts[] = {42.0};
        AstPermMap *pm = astPermMap(3, inperm, 3, outperm, consts, "");
        AstUnitMap *um = astUnitMap(3, "");
        AstCmpMap *cm = astCmpMap(pm, um, 1, "");
        write_fixture(dir, "perm_array_simplify", (AstMapping*)cm);
        cm = astAnnul(cm); pm = astAnnul(pm); um = astAnnul(um);
    }

    /* permmap-05 variant: PermMap with explicit identity outperm=[1,2]
       composed with UnitMap. The composition should detect identity and
       null out the array. */
    {
        if (!astOK) astClearStatus;
        int inperm[] = {2, 1};
        int outperm[] = {1, 2};
        AstPermMap *pm = astPermMap(2, inperm, 2, outperm, NULL, "");
        AstUnitMap *um = astUnitMap(2, "");
        AstCmpMap *cm = astCmpMap(pm, um, 1, "");
        write_fixture(dir, "perm_identity_array_null", (AstMapping*)cm);
        cm = astAnnul(cm); pm = astAnnul(pm); um = astAnnul(um);
    }
}

/* ===== Negative fixtures batch 9 ===== */

static void gen_negative_fixtures_9(const char *dir) {
    if (!astOK) astClearStatus;
    printf("Negative fixtures batch 9:\n");

    /* winmap-04: standalone WinMap with mixed shifts/scales — no self-simplify */
    {
        if (!astOK) astClearStatus;
        double ina[] = {0, 0}, inb[] = {1, 1};
        double outa[] = {1, 2}, outb[] = {4, 6};
        AstWinMap *wm = astWinMap(2, ina, inb, outa, outb, "");
        write_negative_fixture(dir, "neg_win_mixed_no_selfsimplify", (AstMapping*)wm);
        wm = astAnnul(wm);
    }

    /* ratemap-04: RateMaps in parallel — refuses */
    {
        if (!astOK) astClearStatus;
        AstZoomMap *zm = astZoomMap(1, 2.0, "");
        AstRateMap *r1 = astRateMap(zm, 1, 1, "");
        AstRateMap *r2 = astRateMap(zm, 1, 1, "");
        AstCmpMap *cm = astCmpMap(r1, r2, 0, "");
        write_negative_fixture(dir, "neg_ratemap_parallel", (AstMapping*)cm);
        cm = astAnnul(cm); r1 = astAnnul(r1); r2 = astAnnul(r2);
        zm = astAnnul(zm);
    }

    /* timemap-16: 5-arg pair with mismatched arguments */
    {
        if (!astOK) astClearStatus;
        double args1[] = {53000.0, 0.5, -2.5, 6378.0, 0.0, 0.0, 0.0};
        double args2[] = {53000.0, 0.5, -2.6, 6378.0, 0.0, 0.0, 0.0};
        AstTimeMap *tm = astTimeMap(0, "");
        astTimeAdd(tm, "TTTOTDB", 5, args1);
        astTimeAdd(tm, "TDBTOTT", 5, args2);
        write_negative_fixture(dir, "neg_time_5arg_mismatch", (AstMapping*)tm);
        tm = astAnnul(tm);
    }

    /* specmap-16: 6-arg pair with mismatched arguments */
    {
        if (!astOK) astClearStatus;
        double args1[] = {-2.5, 0.9, 1000.0, 51544.0, 0.5, 1.2, 0.0};
        double args2[] = {-2.5, 0.9, 1000.0, 51544.0, 0.5, 1.3, 0.0};
        AstSpecMap *sm = astSpecMap(1, 0, "");
        astSpecAdd(sm, "TPF2HL", 6, args1);
        astSpecAdd(sm, "HLF2TP", 6, args2);
        write_negative_fixture(dir, "neg_spec_6arg_mismatch", (AstMapping*)sm);
        sm = astAnnul(sm);
    }

    /* slamap-22: 4-arg pair with mismatched arguments (HPCEQ+EQHPC) */
    {
        if (!astOK) astClearStatus;
        double args1[] = {51544.0, 0.5, 1.2, 150.0e6};
        double args2[] = {51544.0, 0.5, 1.2, 151.0e6};
        AstSlaMap *sm = astSlaMap(0, "");
        astSlaAdd(sm, "HPCEQ", 4, args1);
        astSlaAdd(sm, "EQHPC", 4, args2);
        write_negative_fixture(dir, "neg_sla_4arg_mismatch", (AstMapping*)sm);
        sm = astAnnul(sm);
    }

    /* polymap-03: PolyMap with only inverse coefficients — no forward → linearize refused */
    {
        if (!astOK) astClearStatus;
        double coeff_i[] = {1.0, 1, 1};
        AstPolyMap *pm = astPolyMap(1, 1, 0, NULL, 1, coeff_i, "");
        write_negative_fixture(dir, "neg_poly_no_forward", (AstMapping*)pm);
        pm = astAnnul(pm);
    }

    /* splinemap-07: two different SplineMaps in opposite directions — astEqual fails */
    {
        if (!astOK) astClearStatus;
        double tx[] = {0, 1}, ty[] = {0, 1};
        double cu1[] = {1.5}, cv1[] = {2.5};
        double cu2[] = {3.0}, cv2[] = {4.0};
        AstSplineMap *s1 = astSplineMap(1, 1, 1, 1, tx, ty, cu1, cv1, "");
        AstSplineMap *s2 = astSplineMap(1, 1, 1, 1, tx, ty, cu2, cv2, "");
        astInvert(s2);
        AstCmpMap *cm = astCmpMap(s1, s2, 1, "");
        write_negative_fixture(dir, "neg_spline_different_coeffs", (AstMapping*)cm);
        cm = astAnnul(cm); s1 = astAnnul(s1); s2 = astAnnul(s2);
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
    gen_negative_fixtures_2(dir);
    gen_negative_fixtures_3(dir);
    gen_negative_fixtures_4(dir);
    gen_negative_fixtures_5(dir);
    gen_negative_fixtures_6(dir);
    gen_negative_fixtures_7(dir);
    gen_negative_fixtures_8(dir);
    gen_cascade_positives_2(dir);
    gen_cascade_positives_3(dir);
    gen_negative_fixtures_9(dir);

    astEnd;

    if (!astOK) {
        fprintf(stderr, "AST error occurred (status=%d)\n", status);
        return 1;
    }

    printf("\nAll fixtures generated successfully.\n");
    return 0;
}
