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

    astEnd;

    if (!astOK) {
        fprintf(stderr, "AST error occurred (status=%d)\n", status);
        return 1;
    }

    printf("\nAll fixtures generated successfully.\n");
    return 0;
}
