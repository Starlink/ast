/*
 * gen_native_fixtures.c
 *
 * Reads AST objects from input files and writes them back using a default
 * Channel (Comment=1, Full=0, Indent=3) to produce byte-exact reference
 * fixtures for the Rust port's native serialization round-trip tests.
 *
 * For each input file, produces a .native output in the specified directory
 * that represents the canonical serialization of that object.
 *
 * Also generates simple standalone objects (UnitMap, ZoomMap, ShiftMap,
 * MatrixMap, PermMap, CmpMap, Frame, CmpFrame, FrameSet with SkyFrame)
 * as byte-exact reference fixtures.
 *
 * Build:
 *   cmake --build build --target gen_native_fixtures
 * Run:
 *   ./build/ast_tester/gen_native_fixtures <output_dir>
 *   (reads from ast_tester/ input files, writes to <output_dir>)
 */

#include "ast.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Write an AST object to a file using default Channel settings. */
static int write_object(const char *path, AstObject *obj) {
    AstChannel *chan;

    chan = astChannel(NULL, NULL, "SinkFile=%s", path);
    if (!astOK) {
        fprintf(stderr, "ERROR: cannot create Channel for %s\n", path);
        return 0;
    }
    if (astWrite(chan, obj) != 1) {
        fprintf(stderr, "ERROR: astWrite failed for %s\n", path);
        chan = astAnnul(chan);
        return 0;
    }
    chan = astAnnul(chan);
    return 1;
}

/* Read an AST object from a native-format file. */
static AstObject *read_object(const char *path) {
    AstChannel *chan;
    AstObject *obj;

    chan = astChannel(NULL, NULL, "SourceFile=%s", path);
    if (!astOK) {
        fprintf(stderr, "ERROR: cannot create Channel for %s\n", path);
        return NULL;
    }
    obj = astRead(chan);
    chan = astAnnul(chan);
    if (!obj) {
        fprintf(stderr, "ERROR: astRead returned NULL for %s\n", path);
    }
    return obj;
}

/* Normalize an existing .ast file: read it and write it back with default
 * Channel settings. This produces the canonical form for byte-exact
 * comparison. */
static int normalize_file(const char *input_path, const char *output_path) {
    AstObject *obj = read_object(input_path);
    if (!obj) return 0;
    int ok = write_object(output_path, obj);
    obj = astAnnul(obj);
    return ok;
}

/* ===== Generate standalone fixtures ===== */

static void gen_unitmap(const char *dir) {
    char path[512];
    AstUnitMap *m = astUnitMap(3, "");
    snprintf(path, sizeof(path), "%s/unitmap.ast", dir);
    write_object(path, (AstObject *)m);
    m = astAnnul(m);
    printf("  unitmap.ast\n");
}

static void gen_zoommap(const char *dir) {
    char path[512];
    AstZoomMap *m = astZoomMap(2, 3.5, "");
    snprintf(path, sizeof(path), "%s/zoommap.ast", dir);
    write_object(path, (AstObject *)m);
    m = astAnnul(m);
    printf("  zoommap.ast\n");
}

static void gen_shiftmap(const char *dir) {
    char path[512];
    double shift[] = {1.5, -2.7};
    AstShiftMap *m = astShiftMap(2, shift, "");
    snprintf(path, sizeof(path), "%s/shiftmap.ast", dir);
    write_object(path, (AstObject *)m);
    m = astAnnul(m);
    printf("  shiftmap.ast\n");
}

static void gen_matrixmap(const char *dir) {
    char path[512];
    /* 2x2 rotation by 30 degrees */
    double cos30 = 0.86602540378443864676;
    double sin30 = 0.5;
    double matrix[] = {cos30, -sin30, sin30, cos30};
    AstMatrixMap *m = astMatrixMap(2, 2, 0, matrix, "");
    snprintf(path, sizeof(path), "%s/matrixmap.ast", dir);
    write_object(path, (AstObject *)m);
    m = astAnnul(m);
    printf("  matrixmap.ast\n");
}

static void gen_matrixmap_diagonal(const char *dir) {
    char path[512];
    double diag[] = {2.0, 3.0};
    AstMatrixMap *m = astMatrixMap(2, 2, 1, diag, "");
    snprintf(path, sizeof(path), "%s/matrixmap_diagonal.ast", dir);
    write_object(path, (AstObject *)m);
    m = astAnnul(m);
    printf("  matrixmap_diagonal.ast\n");
}

static void gen_permmap(const char *dir) {
    char path[512];
    int outperm[] = {2, 1};
    int inperm[] = {2, 1};
    double constants[] = {42.0};
    AstPermMap *m = astPermMap(2, inperm, 2, outperm, constants, "");
    snprintf(path, sizeof(path), "%s/permmap.ast", dir);
    write_object(path, (AstObject *)m);
    m = astAnnul(m);
    printf("  permmap.ast\n");
}

static void gen_permmap_with_constant(const char *dir) {
    char path[512];
    /* 2 in → 3 out; third output is constant 99.0 */
    int outperm[] = {1, 2, -1};
    int inperm[] = {1, 2};
    double constants[] = {99.0};
    AstPermMap *m = astPermMap(2, inperm, 3, outperm, constants, "");
    snprintf(path, sizeof(path), "%s/permmap_constant.ast", dir);
    write_object(path, (AstObject *)m);
    m = astAnnul(m);
    printf("  permmap_constant.ast\n");
}

static void gen_winmap(const char *dir) {
    char path[512];
    double ina[] = {0, 0}, inb[] = {1, 1};
    double outa[] = {10, 20}, outb[] = {30, 60};
    AstWinMap *m = astWinMap(2, ina, inb, outa, outb, "");
    snprintf(path, sizeof(path), "%s/winmap.ast", dir);
    write_object(path, (AstObject *)m);
    m = astAnnul(m);
    printf("  winmap.ast\n");
}

static void gen_cmpmap_series(const char *dir) {
    char path[512];
    AstZoomMap *z = astZoomMap(2, 2.0, "");
    double shift[] = {5.0, -3.0};
    AstShiftMap *s = astShiftMap(2, shift, "");
    AstCmpMap *m = astCmpMap(z, s, 1, "");
    snprintf(path, sizeof(path), "%s/cmpmap_series.ast", dir);
    write_object(path, (AstObject *)m);
    m = astAnnul(m);
    z = astAnnul(z);
    s = astAnnul(s);
    printf("  cmpmap_series.ast\n");
}

static void gen_cmpmap_parallel(const char *dir) {
    char path[512];
    AstZoomMap *z = astZoomMap(1, 2.0, "");
    double shift[] = {5.0};
    AstShiftMap *s = astShiftMap(1, shift, "");
    AstCmpMap *m = astCmpMap(z, s, 0, "");
    snprintf(path, sizeof(path), "%s/cmpmap_parallel.ast", dir);
    write_object(path, (AstObject *)m);
    m = astAnnul(m);
    z = astAnnul(z);
    s = astAnnul(s);
    printf("  cmpmap_parallel.ast\n");
}

static void gen_basicframe(const char *dir) {
    char path[512];
    AstFrame *f = astFrame(2, "Title=Pixel Coordinates,Domain=PIXEL");
    astSetC(f, "Label(1)", "Column");
    astSetC(f, "Label(2)", "Row");
    astSetC(f, "Unit(1)", "pixel");
    astSetC(f, "Unit(2)", "pixel");
    snprintf(path, sizeof(path), "%s/basicframe.ast", dir);
    write_object(path, (AstObject *)f);
    f = astAnnul(f);
    printf("  basicframe.ast\n");
}

static void gen_skyframe(const char *dir) {
    char path[512];
    AstSkyFrame *f = astSkyFrame("System=FK5,Equinox=J2000.0,Epoch=2000.0");
    snprintf(path, sizeof(path), "%s/skyframe.ast", dir);
    write_object(path, (AstObject *)f);
    f = astAnnul(f);
    printf("  skyframe.ast\n");
}

static void gen_skyframe_galactic(const char *dir) {
    char path[512];
    AstSkyFrame *f = astSkyFrame("System=Galactic");
    snprintf(path, sizeof(path), "%s/skyframe_galactic.ast", dir);
    write_object(path, (AstObject *)f);
    f = astAnnul(f);
    printf("  skyframe_galactic.ast\n");
}

static void gen_specframe(const char *dir) {
    char path[512];
    AstSpecFrame *f = astSpecFrame("System=FREQ,Unit=GHz,StdOfRest=Heliocentric");
    snprintf(path, sizeof(path), "%s/specframe.ast", dir);
    write_object(path, (AstObject *)f);
    f = astAnnul(f);
    printf("  specframe.ast\n");
}

static void gen_timeframe(const char *dir) {
    char path[512];
    AstTimeFrame *f = astTimeFrame("TimeScale=UTC,System=MJD");
    snprintf(path, sizeof(path), "%s/timeframe.ast", dir);
    write_object(path, (AstObject *)f);
    f = astAnnul(f);
    printf("  timeframe.ast\n");
}

static void gen_cmpframe(const char *dir) {
    char path[512];
    AstSkyFrame *sky = astSkyFrame("System=FK5");
    AstSpecFrame *spec = astSpecFrame("System=FREQ,Unit=GHz");
    AstCmpFrame *f = astCmpFrame(sky, spec, "");
    snprintf(path, sizeof(path), "%s/cmpframe.ast", dir);
    write_object(path, (AstObject *)f);
    f = astAnnul(f);
    sky = astAnnul(sky);
    spec = astAnnul(spec);
    printf("  cmpframe.ast\n");
}

static void gen_frameset_simple(const char *dir) {
    char path[512];
    /* FrameSet: base=pixel Frame(2), current=sky SkyFrame, mapping=ZoomMap */
    AstFrame *base = astFrame(2, "Title=Pixel,Domain=GRID");
    AstSkyFrame *sky = astSkyFrame("System=FK5,Equinox=J2000.0");
    AstZoomMap *zm = astZoomMap(2, 0.0002777777777777778, "");
    AstFrameSet *fs = astFrameSet(base, "");
    astAddFrame(fs, AST__BASE, zm, sky);
    snprintf(path, sizeof(path), "%s/frameset_simple.ast", dir);
    write_object(path, (AstObject *)fs);
    fs = astAnnul(fs);
    base = astAnnul(base);
    sky = astAnnul(sky);
    zm = astAnnul(zm);
    printf("  frameset_simple.ast\n");
}

static void gen_frameset_3frame(const char *dir) {
    char path[512];
    /* FrameSet with 3 frames: GRID → PIXEL → SKY */
    AstFrame *grid = astFrame(2, "Title=Grid,Domain=GRID");
    AstFrame *pixel = astFrame(2, "Title=Pixel,Domain=PIXEL");
    AstSkyFrame *sky = astSkyFrame("System=FK5,Equinox=J2000.0");
    double shift[] = {-0.5, -0.5};
    AstShiftMap *sm = astShiftMap(2, shift, "");
    AstZoomMap *zm = astZoomMap(2, 0.0002777777777777778, "");
    AstFrameSet *fs = astFrameSet(grid, "");
    astAddFrame(fs, AST__BASE, sm, pixel);
    astAddFrame(fs, 2, zm, sky);
    snprintf(path, sizeof(path), "%s/frameset_3frame.ast", dir);
    write_object(path, (AstObject *)fs);
    fs = astAnnul(fs);
    grid = astAnnul(grid);
    pixel = astAnnul(pixel);
    sky = astAnnul(sky);
    sm = astAnnul(sm);
    zm = astAnnul(zm);
    printf("  frameset_3frame.ast\n");
}

static void gen_frameset_skyspec(const char *dir) {
    char path[512];
    /* FrameSet: base=3D Frame, current=SkyFrame+SpecFrame via CmpFrame */
    AstFrame *base = astFrame(3, "Domain=GRID");
    AstSkyFrame *sky = astSkyFrame("System=FK5");
    AstSpecFrame *spec = astSpecFrame("System=FREQ,Unit=GHz");
    AstCmpFrame *current = astCmpFrame(sky, spec, "");
    /* Use a simple 3x3 diagonal mapping */
    double diag[] = {0.001, 0.001, 1e9};
    AstMatrixMap *mm = astMatrixMap(3, 3, 1, diag, "");
    AstFrameSet *fs = astFrameSet(base, "");
    astAddFrame(fs, AST__BASE, mm, current);
    snprintf(path, sizeof(path), "%s/frameset_skyspec.ast", dir);
    write_object(path, (AstObject *)fs);
    fs = astAnnul(fs);
    base = astAnnul(base);
    sky = astAnnul(sky);
    spec = astAnnul(spec);
    current = astAnnul(current);
    mm = astAnnul(mm);
    printf("  frameset_skyspec.ast\n");
}

/* ===== Normalize existing .ast files ===== */

static void normalize_existing(const char *ast_tester_dir, const char *outdir) {
    /* List of existing .ast files to normalize */
    static const char *files[] = {
        "timj.ast",
        "lsst.ast",
        "degen1.ast",
        "splittest1.ast",
        "specflux.ast",
        NULL
    };

    char inpath[512], outpath[512];
    printf("\nNormalizing existing .ast files:\n");
    for (int i = 0; files[i]; i++) {
        snprintf(inpath, sizeof(inpath), "%s/%s", ast_tester_dir, files[i]);
        snprintf(outpath, sizeof(outpath), "%s/%s", outdir, files[i]);
        if (normalize_file(inpath, outpath)) {
            printf("  %s\n", files[i]);
        }
    }
}

int main(int argc, char *argv[]) {
    int status_value = 0;
    int *status = &status_value;

    if (argc < 2) {
        fprintf(stderr,
            "Usage: gen_native_fixtures <output_dir> [ast_tester_dir]\n"
            "\n"
            "  output_dir     : where to write generated .ast fixture files\n"
            "  ast_tester_dir : path to ast_tester/ for normalizing existing files\n"
            "                   (defaults to the directory containing this program)\n");
        return 1;
    }

    const char *outdir = argv[1];
    const char *ast_tester_dir = argc > 2 ? argv[2] : ".";

    astWatch(status);

    printf("Generating standalone fixtures in %s:\n", outdir);

    /* Mappings */
    gen_unitmap(outdir);
    gen_zoommap(outdir);
    gen_shiftmap(outdir);
    gen_matrixmap(outdir);
    gen_matrixmap_diagonal(outdir);
    gen_permmap(outdir);
    gen_permmap_with_constant(outdir);
    gen_winmap(outdir);
    gen_cmpmap_series(outdir);
    gen_cmpmap_parallel(outdir);

    /* Frames */
    gen_basicframe(outdir);
    gen_skyframe(outdir);
    gen_skyframe_galactic(outdir);
    gen_specframe(outdir);
    gen_timeframe(outdir);
    gen_cmpframe(outdir);

    /* FrameSets */
    gen_frameset_simple(outdir);
    gen_frameset_3frame(outdir);
    gen_frameset_skyspec(outdir);

    /* Normalize existing complex .ast files */
    normalize_existing(ast_tester_dir, outdir);

    if (!astOK) {
        fprintf(stderr, "ERROR: AST status is bad at end\n");
        return 1;
    }

    printf("\nDone. All fixtures written to %s\n", outdir);
    return 0;
}
