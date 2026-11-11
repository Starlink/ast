/*
 * check_transform_oracle: read a transform oracle file (written by
 * gen_transform_oracle), reload each named fixture, re-transform the
 * recorded inputs, and compare outputs within tolerance.  The oracle
 * file is the sole manifest; this program does no corpus scanning.
 *
 * Usage:
 *   check_transform_oracle <root> <oracle_file> [<overrides>]
 *   check_transform_oracle --selftest
 */
#include "ast.h"
#include "transform_oracle.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

typedef struct {
    char relpath[1024];
    int nin, nout, forward;
    int npoint;
    double **in;
    double **out;
} Section;

static double **alloc_cols( int ncol, int nrow ) {
    double **c = malloc( sizeof(double *) * (size_t) ncol );
    for ( int a = 0; a < ncol; a++ ) c[a] = malloc( sizeof(double) * (size_t) nrow );
    return c;
}

static void free_sections( Section *secs, int nsec );

static int oracle_read_file( const char *path, Section **out_secs, int *out_n ) {
    FILE *fp = fopen( path, "r" );
    if ( !fp ) { fprintf( stderr, "cannot open oracle %s\n", path ); return 1; }

    size_t cap = 64, n = 0;
    Section *secs = malloc( sizeof(Section) * cap );
    char line[8192];
    Section *cur = NULL;
    int row = 0, rowcap = 0;

    while ( fgets( line, (int) sizeof line, fp ) ) {
        if ( line[0] == '#' || line[0] == '\n' ) continue;
        if ( line[0] == '[' ) {
            if ( n == cap ) { cap *= 2; secs = realloc( secs, sizeof(Section)*cap ); }
            cur = &secs[n++];
            char dir[16] = "forward";
            /* header: [<relpath>  nin=N nout=M dir=word] */
            if ( sscanf( line, "[%1023s nin=%d nout=%d dir=%15[^]]",
                         cur->relpath, &cur->nin, &cur->nout, dir ) < 3 ) {
                fprintf( stderr, "bad section header: %s", line );
                fclose( fp ); return 1;
            }
            cur->forward = ( strncmp( dir, "inverse", 7 ) != 0 );
            cur->npoint = 0;
            rowcap = oracle_sample_axis_count() + 8;
            cur->in  = alloc_cols( cur->nin,  rowcap );
            cur->out = alloc_cols( cur->nout, rowcap );
            row = 0;
            continue;
        }
        if ( !cur ) continue;
        /* a data row: nin inputs then nout outputs */
        char *p = line, *tok;
        int col = 0, ok = 1;
        if ( row >= rowcap ) {
            rowcap *= 2;
            for ( int a = 0; a < cur->nin;  a++ )
                cur->in[a]  = realloc( cur->in[a],  sizeof(double)*(size_t)rowcap );
            for ( int a = 0; a < cur->nout; a++ )
                cur->out[a] = realloc( cur->out[a], sizeof(double)*(size_t)rowcap );
        }
        while ( ( tok = strtok( p, " \t\r\n" ) ) ) {
            p = NULL;
            int good = 0;
            double v = oracle_parse_double( tok, &good );
            if ( !good ) { ok = 0; break; }
            if ( col < cur->nin ) cur->in[col][row] = v;
            else if ( col < cur->nin + cur->nout ) cur->out[col - cur->nin][row] = v;
            col++;
        }
        /* A malformed data row (unparseable token or wrong column count) is a
           corrupt/truncated oracle, not something to skip silently: fail so
           the gate cannot pass on reduced coverage. */
        if ( !ok || col != cur->nin + cur->nout ) {
            fprintf( stderr,
                "%s: malformed data row in section [%s] (row %d: got %d of "
                "%d columns)\n", path, cur->relpath, row, col,
                cur->nin + cur->nout );
            fclose( fp );
            free_sections( secs, (int) n );
            return 1;
        }
        row++; cur->npoint = row;
    }
    fclose( fp );
    *out_secs = secs; *out_n = (int) n;
    return 0;
}

/* Re-transform a section's recorded inputs with the live mapping into
   freshly allocated columns the caller frees.  Returns NULL on arity
   mismatch (a structural change to the fixture). */
static double **transform_section( AstMapping *map, const Section *s ) {
    if ( astGetI( map, "Nin" ) != ( s->forward ? s->nin : s->nout ) ) return NULL;
    if ( astGetI( map, "Nout" ) != ( s->forward ? s->nout : s->nin ) ) return NULL;
    double **live = alloc_cols( s->nout, s->npoint );
    astTranP( map, s->npoint, s->nin, (const double **) s->in,
              s->forward, s->nout, live );
    return live;
}

static void free_cols( double **c, int ncol ) {
    if ( !c ) return;
    for ( int a = 0; a < ncol; a++ ) free( c[a] );
    free( c );
}

static void free_sections( Section *secs, int nsec ) {
    if ( !secs ) return;
    for ( int i = 0; i < nsec; i++ ) {
        free_cols( secs[i].in,  secs[i].nin );
        free_cols( secs[i].out, secs[i].nout );
    }
    free( secs );
}

/* Is output axis `axis` of `frame` cyclic with period 2*pi?  Tested
   generically via astNorm: a point and the same point shifted by 2*pi on
   that axis normalize to the same value iff the axis wraps (e.g. a SkyFrame
   longitude).  Latitude reflects rather than wraps, and basic frames
   (GRID pixels, SpecFrame) are not cyclic, so their pixel/linear axes are
   compared without wrap -- which is what catches a genuine 2*pi drift in a
   non-angular output. */
static int axis_is_cyclic( AstFrame *frame, int axis ) {
    if ( !frame ) return 0;
    int naxes = astGetI( frame, "Naxes" );
    if ( !astOK || axis >= naxes ) { astClearStatus; return 0; }
    double twopi = 2.0 * acos( -1.0 );
    double *p = malloc( sizeof(double) * (size_t) naxes );
    double *q = malloc( sizeof(double) * (size_t) naxes );
    for ( int i = 0; i < naxes; i++ ) { p[i] = 0.1; q[i] = 0.1; }
    q[axis] += twopi;
    astNorm( frame, p );
    astNorm( frame, q );
    int cyclic = astOK && fabs( p[axis] - q[axis] ) < 1e-6;
    if ( !astOK ) astClearStatus;
    free( p ); free( q );
    return cyclic;
}

/* Match one output value: wrap-aware (2*pi periodic) only on axes flagged
   angular, plain otherwise. */
static int axis_match( double got, double ref, double rtol, double atol,
                       int angular ) {
    return angular ? oracle_within_tol_wrap( got, ref, rtol, atol )
                   : oracle_within_tol( got, ref, rtol, atol );
}

/* Compare two equal-shape output sets; report and count mismatches.
   The absolute tolerance on each axis is scaled to that axis's recorded
   values (see oracle_axis_atol).
   is_ang[a] != 0 marks output axis a as angular (wrap-aware comparison).
   skip_rows is a bitmap of row indices to leave uncompared, for rows whose
   recorded value sits on a mathematical singularity (see load_overrides). */
static int compare_outputs( const char *label, const char *relpath,
                            double **got, double **ref, int ncol, int npoint,
                            double rtol, double atol, const int *is_ang,
                            uint64_t skip_rows ) {
    int fails = 0;
    double *axis_atol = malloc( sizeof(double) * (size_t) ncol );
    for ( int a = 0; a < ncol; a++ )
        axis_atol[a] = oracle_axis_atol( atol, oracle_axis_scale( ref[a], npoint ) );
    for ( int p = 0; p < npoint; p++ ) {
        if ( p < 64 && ( skip_rows >> p ) & 1u ) {
            fprintf( stderr, "SKIP [%s] %s row=%d (overrides file)\n",
                     label, relpath, p );
            continue;
        }
        for ( int a = 0; a < ncol; a++ ) {
            if ( !axis_match( got[a][p], ref[a][p], rtol, axis_atol[a],
                              is_ang[a] ) ) {
                double d = fabs( got[a][p] - ref[a][p] );
                double rel = ref[a][p] != 0.0 ? d / fabs( ref[a][p] ) : d;
                fprintf( stderr,
                    "MISMATCH [%s] %s row=%d axis=%d ref=%.17g got=%.17g rel=%.3g\n",
                    label, relpath, p, a, ref[a][p], got[a][p], rel );
                fails++;
            }
        }
    }
    free( axis_atol );
    return fails;
}

/* A Mapping and its simplification nearly always agree at one point, the
   origin, whatever the simplification did: a Mapping that sends every input
   to zero agrees there with a UnitMap.  So a fixture counts as checked only
   if it gives a value at this many distinct sampled inputs. */
#define MIN_CHECKED_INPUTS 2

/* Number of distinct input points among the rows p of in (nin x npoint)
   for which keep[p] is non-zero. */
static int distinct_inputs( double **in, int nin, int npoint, const int *keep ) {
    int n = 0;
    for ( int p = 0; p < npoint; p++ ) {
        if ( !keep[p] ) continue;
        int seen = 0;
        for ( int q = 0; q < p && !seen; q++ ) {
            if ( !keep[q] ) continue;
            int same = 1;
            for ( int a = 0; a < nin && same; a++ )
                if ( in[a][q] != in[a][p] ) same = 0;
            seen = same;
        }
        if ( !seen ) n++;
    }
    return n;
}

/* Number of distinct sampled inputs at which a section's recorded outputs
   include at least one value that is not AST__BAD. */
static int valued_inputs( const Section *s ) {
    int *keep = calloc( (size_t) s->npoint, sizeof(int) );
    for ( int p = 0; p < s->npoint; p++ )
        for ( int a = 0; a < s->nout; a++ )
            if ( s->out[a][p] != AST__BAD ) { keep[p] = 1; break; }
    int n = distinct_inputs( s->in, s->nin, s->npoint, keep );
    free( keep );
    return n;
}

/* Equivalence check C.  Simplification must not change a Mapping's values,
   so the .map and .simp forward outputs must agree at every sampled point,
   both as recorded and as computed now.  Values that are both AST__BAD
   agree, since a Region must stay BAD outside itself.

   A row whose recorded outputs disagree is reported as a simplification
   that changes the transform, unless it is one of accept_rows (a bitmap of
   0-based rows from the overrides file), in which case the row is not
   compared at all.  A listed row whose recorded outputs agree is stale;
   *nstale is set to the number of them.

   in (nin x npoint) holds the sampled inputs, which the .map and .simp
   sections share.  *nagree is set to the number of distinct inputs at which
   the recorded outputs agree on every axis and include a value that is not
   AST__BAD: zero when the recorded values agree nowhere, or only where both
   forms are BAD.  The absolute tolerance on each axis is scaled to the
   larger of the two forms' recorded values. */
static int compare_equiv( const char *relpath, double **in, int nin,
                          double **map_live, double **simp_live,
                          double **map_ref,  double **simp_ref,
                          int ncol, int npoint, double rtol, double atol,
                          const int *is_ang, uint64_t accept_rows,
                          int *nagree, int *nstale ) {
    int fails = 0;
    int *agree = calloc( (size_t) npoint, sizeof(int) );
    double *axis_atol = malloc( sizeof(double) * (size_t) ncol );
    *nstale = 0;
    for ( int a = 0; a < ncol; a++ )
        axis_atol[a] = oracle_axis_atol( atol,
                                         fmax( oracle_axis_scale( map_ref[a], npoint ),
                                               oracle_axis_scale( simp_ref[a], npoint ) ) );
    for ( int p = 0; p < npoint; p++ ) {
        int bad_axis = -1, row_valued = 0;
        for ( int a = 0; a < ncol && bad_axis < 0; a++ ) {
            if ( !axis_match( map_ref[a][p], simp_ref[a][p], rtol, axis_atol[a],
                              is_ang[a] ) ) bad_axis = a;
            else if ( map_ref[a][p] != AST__BAD ) row_valued = 1;
        }
        int accepted = p < 64 && ( ( accept_rows >> p ) & 1u );
        if ( accepted ) {
            if ( bad_axis < 0 ) {
                fprintf( stderr,
                    "MISMATCH [equiv-rows] %s row=%d: the recorded outputs "
                    "agree, so the row's \"equiv-rows\" override is stale\n",
                    relpath, p );
                ( *nstale )++;
            }
            continue;
        }
        if ( bad_axis >= 0 ) {
            char mbuf[64], sbuf[64];
            oracle_format_double( mbuf, sizeof mbuf, map_ref[bad_axis][p] );
            oracle_format_double( sbuf, sizeof sbuf, simp_ref[bad_axis][p] );
            fprintf( stderr,
                "MISMATCH [equiv-recorded] %s row=%d axis=%d map=%s simp=%s: "
                "the simplification changes the transform\n",
                relpath, p, bad_axis, mbuf, sbuf );
            fails++;
            continue;
        }
        agree[p] = row_valued;
        for ( int a = 0; a < ncol; a++ ) {
            if ( !axis_match( map_live[a][p], simp_live[a][p], rtol, axis_atol[a],
                              is_ang[a] ) ) {
                double d = fabs( map_live[a][p] - simp_live[a][p] );
                double rel = simp_live[a][p] != 0.0 ? d / fabs( simp_live[a][p] ) : d;
                fprintf( stderr,
                    "MISMATCH [equiv] %s row=%d axis=%d map=%.17g simp=%.17g rel=%.3g\n",
                    relpath, p, a, map_live[a][p], simp_live[a][p], rel );
                fails++;
            }
        }
    }
    *nagree = distinct_inputs( in, nin, npoint, agree );
    free( agree );
    free( axis_atol );
    return fails;
}

/* Per-fixture check overrides, keyed by relpath.  Round-trip lines:
   "relpath off" disables round-trip, "relpath <rtol> <atol>" loosens it.
   Golden lines: "relpath golden off" / "relpath golden-fwd off" /
   "relpath golden-inv off" skip the golden comparison for the named
   direction(s), and "relpath golden-fwd-rows <list>" /
   "relpath golden-inv-rows <list>" skip only the listed rows, where <list>
   is comma-separated 0-based row indices below 64.  A row list is the
   preferred form: it leaves the rest of the section pinned.  Equivalence
   lines: "relpath equiv off", keyed on a .map relpath, accepts a .map/.simp
   pair whose recorded outputs agree, with a value that is not AST__BAD, at
   fewer than MIN_CHECKED_INPUTS sampled inputs, so that its equivalence is
   barely checked, "relpath equiv-rows <list>" accepts a recorded
   disagreement between the .map and .simp at only the listed rows, which
   are then not compared, and "relpath equiv-tol <rtol> <atol>" compares
   the pair to a looser tolerance, for a Mapping whose inverse is itself an
   approximation.  Value lines: "relpath values off" accepts a .map with no
   .simp partner whose recorded outputs have a value that is not AST__BAD at
   fewer than MIN_CHECKED_INPUTS sampled inputs.  A fixture may have lines of
   each kind; they are looked up independently. */
typedef struct {
    char relpath[1024];
    int rtrip;                  /* round-trip entry?  else golden entry */
    double rtol, atol;          /* round-trip: loosened tolerances */
    int off;                    /* round-trip: check disabled */
    int golden_fwd_off;         /* golden: skip forward sections */
    int golden_inv_off;         /* golden: skip inverse sections */
    uint64_t golden_fwd_rows;   /* golden: forward rows to leave uncompared */
    uint64_t golden_inv_rows;   /* golden: inverse rows to leave uncompared */
    int equiv_off;              /* equivalence: barely-compared pair accepted */
    uint64_t equiv_rows;        /* equivalence: rows accepted to disagree */
    double equiv_rtol, equiv_atol;  /* equivalence: tolerances, 0 = default */
    int values_off;             /* values: mostly-BAD lone .map accepted */
} Override;
static Override *g_over = NULL;
static int g_nover = 0;

/* Parse a comma-separated list of 0-based row indices into a bitmap.  Rows
   at or beyond 64 are ignored, as is any token that is not a number. */
static uint64_t parse_row_list( const char *list ) {
    uint64_t rows = 0;
    const char *p = list;
    while ( *p ) {
        char *end;
        long v = strtol( p, &end, 10 );
        if ( end == p ) break;
        if ( v >= 0 && v < 64 ) rows |= (uint64_t) 1u << v;
        p = end;
        while ( *p == ',' || *p == ' ' ) p++;
    }
    return rows;
}

static void load_overrides( const char *path ) {
    if ( !path ) return;
    FILE *fp = fopen( path, "r" );
    if ( !fp ) return;
    size_t cap = 16; g_over = malloc( sizeof(Override) * cap );
    char line[1024];
    while ( fgets( line, (int) sizeof line, fp ) ) {
        if ( line[0] == '#' || line[0] == '\n' ) continue;
        char rel[1024], what[64];
        if ( sscanf( line, "%1023s %63s", rel, what ) < 2 ) continue;
        if ( g_nover == (int) cap ) { cap *= 2; g_over = realloc( g_over, sizeof(Override)*cap ); }
        Override *o = &g_over[g_nover++];
        memset( o, 0, sizeof *o );
        snprintf( o->relpath, sizeof o->relpath, "%s", rel );
        if ( strcmp( what, "golden-fwd-rows" ) == 0 ||
             strcmp( what, "golden-inv-rows" ) == 0 ) {
            char list[1024] = "";
            sscanf( line, "%*s %*s %1023s", list );
            uint64_t rows = parse_row_list( list );
            if ( strcmp( what, "golden-fwd-rows" ) == 0 ) o->golden_fwd_rows = rows;
            else                                         o->golden_inv_rows = rows;
        } else if ( strcmp( what, "equiv-tol" ) == 0 ) {
            sscanf( line, "%*s %*s %lf %lf", &o->equiv_rtol, &o->equiv_atol );
        } else if ( strcmp( what, "equiv-rows" ) == 0 ) {
            char list[1024] = "";
            sscanf( line, "%*s %*s %1023s", list );
            o->equiv_rows = parse_row_list( list );
        } else if ( strcmp( what, "equiv" ) == 0 ) {
            char arg[64] = "";
            sscanf( line, "%*s %*s %63s", arg );
            o->equiv_off = ( strcmp( arg, "off" ) == 0 );
        } else if ( strcmp( what, "values" ) == 0 ) {
            char arg[64] = "";
            sscanf( line, "%*s %*s %63s", arg );
            o->values_off = ( strcmp( arg, "off" ) == 0 );
        } else if ( strncmp( what, "golden", 6 ) == 0 ) {
            o->golden_fwd_off = ( strcmp( what, "golden" ) == 0 ||
                                  strcmp( what, "golden-fwd" ) == 0 );
            o->golden_inv_off = ( strcmp( what, "golden" ) == 0 ||
                                  strcmp( what, "golden-inv" ) == 0 );
        } else {
            o->rtrip = 1;
            if ( strcmp( what, "off" ) == 0 ) { o->off = 1; }
            else { o->rtol = ORACLE_DEF_RTRIP_RTOL; o->atol = ORACLE_DEF_RTRIP_ATOL;
                   sscanf( line, "%*s %lf %lf", &o->rtol, &o->atol ); }
        }
    }
    fclose( fp );
}

/* Fill rtol/atol for a fixture's round-trip; return 0 if it is disabled. */
static int rtrip_tol( const char *relpath, double *rtol, double *atol ) {
    *rtol = ORACLE_DEF_RTRIP_RTOL; *atol = ORACLE_DEF_RTRIP_ATOL;
    for ( int i = 0; i < g_nover; i++ ) {
        if ( g_over[i].rtrip && strcmp( g_over[i].relpath, relpath ) == 0 ) {
            if ( g_over[i].off ) return 0;
            *rtol = g_over[i].rtol; *atol = g_over[i].atol; return 1;
        }
    }
    return 1;
}

/* Is a .map fixture accepted although its recorded .map and .simp outputs
   agree, with a value, at fewer than MIN_CHECKED_INPUTS sampled inputs? */
static int equiv_off( const char *relpath ) {
    for ( int i = 0; i < g_nover; i++ ) {
        if ( g_over[i].equiv_off && strcmp( g_over[i].relpath, relpath ) == 0 )
            return 1;
    }
    return 0;
}

/* Bitmap of rows at which a .map fixture's recorded .map and .simp outputs
   are accepted to disagree. */
static uint64_t equiv_accept_rows( const char *relpath ) {
    uint64_t rows = 0;
    for ( int i = 0; i < g_nover; i++ )
        if ( strcmp( g_over[i].relpath, relpath ) == 0 ) rows |= g_over[i].equiv_rows;
    return rows;
}

/* Fill rtol/atol for a .map fixture's equivalence check. */
static void equiv_tol( const char *relpath, double *rtol, double *atol ) {
    *rtol = ORACLE_DEF_EQUIV_RTOL;
    *atol = ORACLE_DEF_EQUIV_ATOL;
    for ( int i = 0; i < g_nover; i++ ) {
        if ( strcmp( g_over[i].relpath, relpath ) == 0 &&
             ( g_over[i].equiv_rtol > 0.0 || g_over[i].equiv_atol > 0.0 ) ) {
            *rtol = g_over[i].equiv_rtol;
            *atol = g_over[i].equiv_atol;
        }
    }
}

/* Is a .map fixture with no .simp partner accepted although its recorded
   outputs have a value at fewer than MIN_CHECKED_INPUTS sampled inputs? */
static int values_off( const char *relpath ) {
    for ( int i = 0; i < g_nover; i++ ) {
        if ( g_over[i].values_off && strcmp( g_over[i].relpath, relpath ) == 0 )
            return 1;
    }
    return 0;
}

/* Is the golden comparison disabled for this fixture and direction? */
static int golden_off( const char *relpath, int forward ) {
    for ( int i = 0; i < g_nover; i++ ) {
        if ( !g_over[i].rtrip && strcmp( g_over[i].relpath, relpath ) == 0 ) {
            if ( forward ? g_over[i].golden_fwd_off
                         : g_over[i].golden_inv_off ) return 1;
        }
    }
    return 0;
}

/* Bitmap of rows the golden comparison should leave uncompared for this
   fixture and direction. */
static uint64_t golden_skip_rows( const char *relpath, int forward ) {
    uint64_t rows = 0;
    for ( int i = 0; i < g_nover; i++ ) {
        if ( !g_over[i].rtrip && strcmp( g_over[i].relpath, relpath ) == 0 ) {
            rows |= forward ? g_over[i].golden_fwd_rows
                            : g_over[i].golden_inv_rows;
        }
    }
    return rows;
}

/* Round-trip accuracy check (GRID-domain corpora only).  fwd is a fixture's
   forward section; inv_live are the live outputs of its inverse section
   (i.e. inverse(forward(P))).  Assert those recover the forward inputs P,
   skipping points where the forward output was BAD (e.g. projection corners
   that have no inverse).

   Alias-aware: a recovered pixel that differs from P is still accepted if it
   forward-maps to the same sky as P.  This happens when an image is wider
   than the projection period (e.g. a quad-cube image spanning >360 deg),
   where the inverse legitimately returns a different but equivalent
   preimage.  Only a recovered pixel that maps to a *different* sky is a real
   inverse failure. */
static int compare_rtrip( AstMapping *map, const Section *fwd,
                          double **inv_live, double rtol, double atol ) {
    int fails = 0;

    /* Forward-transform the recovered pixels so aliased preimages can be
       detected by comparing their sky against the original forward output. */
    double **sky2 = alloc_cols( fwd->nout, fwd->npoint );
    astTranP( map, fwd->npoint, fwd->nin, (const double **) inv_live,
              1, fwd->nout, sky2 );
    int have_sky2 = astOK;
    if ( !astOK ) astClearStatus;
    double *sky_atol = malloc( sizeof(double) * (size_t) fwd->nout );
    for ( int b = 0; b < fwd->nout; b++ )
        sky_atol[b] = oracle_axis_atol( ORACLE_DEF_ATOL,
                                        oracle_axis_scale( fwd->out[b], fwd->npoint ) );

    for ( int p = 0; p < fwd->npoint; p++ ) {
        int fbad = 0;
        for ( int b = 0; b < fwd->nout; b++ )
            if ( fwd->out[b][p] == AST__BAD ) { fbad = 1; break; }
        if ( fbad ) continue;                 /* nothing to invert */

        int pixel_ok = 1;
        for ( int a = 0; a < fwd->nin; a++ )
            if ( !oracle_within_tol( inv_live[a][p], fwd->in[a][p], rtol, atol ) )
                { pixel_ok = 0; break; }
        if ( pixel_ok ) continue;

        /* Different pixel: accept if it forward-maps to the same sky. */
        if ( have_sky2 ) {
            int alias = 1;
            for ( int b = 0; b < fwd->nout; b++ )
                if ( !oracle_within_tol_wrap( sky2[b][p], fwd->out[b][p],
                                              ORACLE_DEF_RTOL, sky_atol[b] ) )
                    { alias = 0; break; }
            if ( alias ) continue;
        }

        /* Genuine round-trip failure: report the offending axes. */
        for ( int a = 0; a < fwd->nin; a++ ) {
            if ( !oracle_within_tol( inv_live[a][p], fwd->in[a][p], rtol, atol ) ) {
                double d = fabs( inv_live[a][p] - fwd->in[a][p] );
                double rel = fwd->in[a][p] != 0.0 ? d / fabs( fwd->in[a][p] ) : d;
                fprintf( stderr,
                    "MISMATCH [round-trip] %s row=%d axis=%d pixel=%.17g "
                    "recovered=%.17g rel=%.3g\n",
                    fwd->relpath, p, a, fwd->in[a][p], inv_live[a][p], rel );
                fails++;
            }
        }
    }
    free_cols( sky2, fwd->nout );
    free( sky_atol );
    return fails;
}

static int selftest( void ) {
    const char *tmp = "/tmp/oracle_selftest.txt";
    FILE *fp = fopen( tmp, "w" );
    fputs( "# header\n\n"
           "[foo.map  nin=2 nout=1 dir=forward]\n"
           "  1 2  3\n"
           "  4 5  BAD\n\n", fp );
    fclose( fp );
    Section *s = NULL; int n = 0;
    if ( oracle_read_file( tmp, &s, &n ) || n != 1 ) {
        fprintf( stderr, "selftest: parse failed\n" ); return 1;
    }
    int ok = ( strcmp( s[0].relpath, "foo.map" ) == 0 && s[0].nin == 2 &&
               s[0].nout == 1 && s[0].forward == 1 && s[0].npoint == 2 &&
               s[0].in[0][0] == 1 && s[0].in[1][1] == 5 &&
               s[0].out[0][0] == 3 && s[0].out[0][1] == AST__BAD );
    free_sections( s, n );

    /* Overrides parsing and lookup: round-trip lines and golden lines are
       independent even for the same fixture. */
    const char *tmpov = "/tmp/oracle_selftest_overrides.txt";
    fp = fopen( tmpov, "w" );
    fputs( "# comment\n"
           "rtoff.head off\n"
           "loose.head 2e-3 4e-2\n"
           "noisy.map golden-inv off\n"
           "both.map golden off\n"
           "fwdonly.map golden-fwd off\n"
           "rows.map golden-inv-rows 2,3\n"
           "rows.map golden-fwd-rows 0\n"
           "waived.map equiv off\n"
           "notoff.map equiv on\n"
           "allbad.map values off\n"
           "differs.map equiv-rows 1,2\n"
           "loose.map equiv-tol 1e-7 2e-7\n", fp );
    fclose( fp );
    load_overrides( tmpov );
    double rtol, atol;
    ok = ok && !rtrip_tol( "rtoff.head", &rtol, &atol );
    ok = ok && rtrip_tol( "loose.head", &rtol, &atol ) &&
         rtol == 2e-3 && atol == 4e-2;
    ok = ok && rtrip_tol( "noisy.map", &rtol, &atol );   /* golden line only */
    ok = ok && !golden_off( "rtoff.head", 1 ) && !golden_off( "rtoff.head", 0 );
    ok = ok && !golden_off( "noisy.map", 1 ) && golden_off( "noisy.map", 0 );
    ok = ok && golden_off( "both.map", 1 ) && golden_off( "both.map", 0 );
    ok = ok && golden_off( "fwdonly.map", 1 ) && !golden_off( "fwdonly.map", 0 );

    /* A row list narrows the skip instead of disabling the direction. */
    ok = ok && !golden_off( "rows.map", 0 ) && !golden_off( "rows.map", 1 );
    ok = ok && golden_skip_rows( "rows.map", 0 ) == 0xcu;   /* rows 2 and 3 */
    ok = ok && golden_skip_rows( "rows.map", 1 ) == 0x1u;   /* row 0        */
    ok = ok && golden_skip_rows( "noisy.map", 0 ) == 0u;
    ok = ok && parse_row_list( "" ) == 0u;
    ok = ok && parse_row_list( "0,63" ) == ( 1ull | ( 1ull << 63 ) );
    ok = ok && parse_row_list( "64,-1" ) == 0u;

    /* An equivalence line affects neither round-trip nor golden lookups,
       and only "off" waives the check. */
    ok = ok && equiv_off( "waived.map" ) && !equiv_off( "noisy.map" );
    ok = ok && !equiv_off( "notoff.map" );
    ok = ok && rtrip_tol( "waived.map", &rtol, &atol ) &&
         rtol == ORACLE_DEF_RTRIP_RTOL && atol == ORACLE_DEF_RTRIP_ATOL;
    ok = ok && !golden_off( "waived.map", 1 ) && !golden_off( "waived.map", 0 );

    /* A value line is looked up independently of an equivalence line, and
       a row list for disagreements independently of both. */
    ok = ok && values_off( "allbad.map" ) && !equiv_off( "allbad.map" );
    ok = ok && !values_off( "waived.map" );
    ok = ok && equiv_accept_rows( "differs.map" ) == 0x6u &&
         !equiv_off( "differs.map" ) && equiv_accept_rows( "waived.map" ) == 0u;
    ok = ok && rtrip_tol( "differs.map", &rtol, &atol ) &&
         rtol == ORACLE_DEF_RTRIP_RTOL;
    equiv_tol( "loose.map", &rtol, &atol );
    ok = ok && rtol == 1e-7 && atol == 2e-7;
    equiv_tol( "differs.map", &rtol, &atol );
    ok = ok && rtol == ORACLE_DEF_EQUIV_RTOL && atol == ORACLE_DEF_EQUIV_ATOL;
    ok = ok && rtrip_tol( "loose.map", &rtol, &atol ) &&
         rtol == ORACLE_DEF_RTRIP_RTOL;

    /* compare_equiv fails a row whose recorded .map and .simp outputs
       disagree, unless the row is accepted, and counts the distinct inputs
       at which they agree with a value: not rows that disagree, not rows
       that are BAD in both, and not a repeated input twice.  An accepted
       row that agrees is stale. */
    {
        double x[] = { 0.0, 1.0, 2.0 }, *in[] = { x };
        double r0[] = { 1.0, 2.0, 3.0 }, r1[] = { 1.0, 9.0, AST__BAD };
        double *mref[] = { r0 }, *sref[] = { r1 };
        double *mlive[] = { r0 }, *slive[] = { r0 };
        int ang[] = { 0 }, nagree = -1, nstale = -1;
        fprintf( stderr, "selftest: expect two MISMATCH [equiv-recorded], "
                         "rows 1 and 2\n" );
        ok = ok && compare_equiv( "selftest", in, 1, mlive, slive, mref, sref,
                                  1, 3, 1e-9, 1e-9, ang, 0u, &nagree,
                                  &nstale ) == 2 && nagree == 1 && nstale == 0;
        ok = ok && compare_equiv( "selftest", in, 1, mlive, slive, mref, sref,
                                  1, 3, 1e-9, 1e-9, ang, 0x6u, &nagree,
                                  &nstale ) == 0 && nagree == 1 && nstale == 0;
        double *sref_none[] = { ( double[] ) { 5.0, 9.0, AST__BAD } };
        ok = ok && compare_equiv( "selftest", in, 1, mlive, slive, mref,
                                  sref_none, 1, 3, 1e-9, 1e-9, ang, 0x7u,
                                  &nagree, &nstale ) == 0 && nagree == 0;
        double *bad2[] = { ( double[] ) { AST__BAD, AST__BAD, 3.0 } };
        ok = ok && compare_equiv( "selftest", in, 1, bad2, bad2, bad2, bad2,
                                  1, 3, 1e-9, 1e-9, ang, 0u, &nagree,
                                  &nstale ) == 0 && nagree == 1;
        fprintf( stderr, "selftest: expect one MISMATCH [equiv-rows], row 0\n" );
        ok = ok && compare_equiv( "selftest", in, 1, bad2, bad2, bad2, bad2,
                                  1, 3, 1e-9, 1e-9, ang, 0x1u, &nagree,
                                  &nstale ) == 0 && nstale == 1;
        double *rep_in[] = { ( double[] ) { 0.0, 0.0, 1.0 } };
        ok = ok && compare_equiv( "selftest", rep_in, 1, mlive, slive, mref,
                                  mref, 1, 3, 1e-9, 1e-9, ang, 0u, &nagree,
                                  &nstale ) == 0 && nagree == 2;

        /* valued_inputs applies the same rules to one section: a row
           counts if any of its outputs is not BAD. */
        double o0[] = { AST__BAD, 1.0, AST__BAD, 4.0 };
        double o1[] = { AST__BAD, AST__BAD, 2.0, 5.0 };
        double xi[] = { 0.0, 1.0, 2.0, 1.0 };
        Section sec = { .nin = 1, .nout = 2, .forward = 1, .npoint = 4,
                        .in = ( double *[] ) { xi },
                        .out = ( double *[] ) { o0, o1 } };
        ok = ok && valued_inputs( &sec ) == 2;
    }

    /* compare_outputs scales the absolute tolerance per axis.  Each axis
       below is out by 1e-9: that is a mismatch on the axis in metres, where
       it is 0.2% of the value, but not on the axis in radians. */
    {
        double rad_ref[] = { 0.0, 1.5 }, rad_got[] = { 1e-9, 1.5 };
        double m_ref[] = { 5.0e-7, 5.05e-7 }, m_got[] = { 5.01e-7, 5.05e-7 };
        double *ref[] = { rad_ref, m_ref }, *got[] = { rad_got, m_got };
        int ang[] = { 1, 0 };
        fprintf( stderr, "selftest: expect one MISMATCH, on axis 1\n" );
        ok = ok && compare_outputs( "selftest", "selftest", got, ref, 2, 2,
                                    ORACLE_DEF_RTOL, ORACLE_DEF_ATOL, ang,
                                    0 ) == 1;
    }

    printf( "selftest: %s\n", ok ? "ok" : "FAIL" );
    return ok ? 0 : 1;
}

/* Stem of a relpath (path minus final extension). */
static void path_stem( const char *relpath, char *stem, size_t n ) {
    const char *ext = strrchr( relpath, '.' );
    size_t len = ext ? (size_t)( ext - relpath ) : strlen( relpath );
    snprintf( stem, n, "%.*s", (int) len, relpath );
}

static int has_ext( const char *relpath, const char *ext ) {
    const char *dot = strrchr( relpath, '.' );
    return dot && strcmp( dot, ext ) == 0;
}

/* Does the oracle hold a forward section for the .simp partner of the .map
   named by relpath? */
static int has_simp_partner( const Section *secs, int nsec, const char *relpath ) {
    char simp[1024];
    path_stem( relpath, simp, sizeof simp );
    size_t len = strlen( simp );
    if ( len + 6 > sizeof simp ) return 0;
    strcpy( simp + len, ".simp" );
    for ( int i = 0; i < nsec; i++ )
        if ( secs[i].forward && strcmp( secs[i].relpath, simp ) == 0 ) return 1;
    return 0;
}

/* Were two sections sampled at exactly the same inputs? */
static int same_inputs( const Section *a, const Section *b ) {
    if ( a->nin != b->nin || a->npoint != b->npoint ) return 0;
    for ( int i = 0; i < a->nin; i++ )
        for ( int p = 0; p < a->npoint; p++ )
            if ( a->in[i][p] != b->in[i][p] ) return 0;
    return 1;
}

int main( int argc, char *argv[] ) {
    int status_value = 0; int *status = &status_value;
    astWatch( status );

    oracle_register_intramaps( status );
    if ( argc == 2 && strcmp( argv[1], "--selftest" ) == 0 ) return selftest();
    if ( argc < 3 ) {
        fprintf( stderr,
            "Usage: check_transform_oracle <root> <oracle_file> [<overrides>]\n" );
        return 2;
    }
    const char *root   = argv[1];
    const char *oracle = argv[2];
    if ( argc >= 4 ) load_overrides( argv[3] );

    Section *secs = NULL; int nsec = 0;
    if ( oracle_read_file( oracle, &secs, &nsec ) ) return 2;

    int total_fail = 0, checked = 0;

    /* One-slot cache of the most recent forward section, used to pair a
       fixture's inverse section (emitted immediately after) for the
       round-trip accuracy check. */
    Section *prev_fwd = NULL;

    /* One-slot cache of the most recent .map forward section and its live
       outputs, so check C can compare the matching .simp forward section
       against it.  The generator emits a stem's .map sections immediately
       before its .simp sections; only .map forward sections update the
       cache and only .simp forward sections consume it, so intervening
       inverse sections do not disturb it. */
    Section *prev_map = NULL;
    double **prev_map_live = NULL;

    for ( int i = 0; i < nsec; i++ ) {
        Section *s = &secs[i];
        astBegin;
        AstFrame *base = NULL, *cur = NULL;
        AstMapping *map = oracle_load_mapping( root, s->relpath, &base, &cur );
        if ( !map ) {
            fprintf( stderr, "LOADFAIL %s (fixture missing or unreadable)\n",
                     s->relpath );
            total_fail++;
            astEnd;
            continue;
        }
        double **live = transform_section( map, s );
        if ( !live ) {
            fprintf( stderr, "ARITY %s dir=%s (Nin/Nout changed)\n",
                     s->relpath, s->forward ? "forward" : "inverse" );
            total_fail++;
            map = astAnnul( map ); astEnd; continue;
        }

        /* Which output axes are angular?  The output frame for this direction
           is the current frame (forward) or base frame (inverse) of a
           FrameSet fixture.  A bare Mapping has no frame, so fall back to
           treating every axis as angular (these are coordinate mappings with
           no pixel axes to protect). */
        AstFrame *outfr = s->forward ? cur : base;
        int *is_ang = malloc( sizeof(int) * (size_t) s->nout );
        for ( int a = 0; a < s->nout; a++ )
            is_ang[a] = outfr ? axis_is_cyclic( outfr, a ) : 1;

        /* Golden check: applies uniformly to forward and inverse sections,
           unless disabled for this fixture/direction by the overrides file
           (used for sections whose recorded values sit on a mathematical
           singularity and are not architecture-stable).  The section's live
           outputs are still computed and still feed the round-trip and
           equivalence checks below. */
        char label[32];
        snprintf( label, sizeof label, "golden-%s",
                  s->forward ? "fwd" : "inv" );
        if ( golden_off( s->relpath, s->forward ) ) {
            fprintf( stderr, "SKIP [%s] %s (overrides file)\n",
                     label, s->relpath );
        } else {
            total_fail += compare_outputs( label, s->relpath, live, s->out,
                                           s->nout, s->npoint,
                                           ORACLE_DEF_RTOL, ORACLE_DEF_ATOL,
                                           is_ang,
                                           golden_skip_rows( s->relpath,
                                                             s->forward ) );
            checked++;
        }

        /* Round-trip accuracy (GRID-domain corpora): an inverse section
           paired with the immediately preceding forward section of the same
           fixture verifies inverse(forward(P)) ~= P. */
        int grid_domain = has_ext( s->relpath, ".head" ) ||
                          has_ext( s->relpath, ".ast" );
        if ( !s->forward && grid_domain && prev_fwd &&
             strcmp( prev_fwd->relpath, s->relpath ) == 0 &&
             prev_fwd->nin == s->nout && prev_fwd->npoint == s->npoint ) {
            double rtol, atol;
            if ( rtrip_tol( s->relpath, &rtol, &atol ) )
                total_fail += compare_rtrip( map, prev_fwd, live, rtol, atol );
        }
        if ( s->forward ) prev_fwd = s;

        /* Equivalence check C: this section is a .simp forward whose stem
           matches the cached preceding .map forward. */
        int is_fwd = s->forward;
        if ( is_fwd && has_ext( s->relpath, ".simp" ) && prev_map &&
             prev_map->nout == s->nout && prev_map->npoint == s->npoint ) {
            char stem_s[1024], stem_m[1024];
            path_stem( s->relpath, stem_s, sizeof stem_s );
            path_stem( prev_map->relpath, stem_m, sizeof stem_m );
            if ( strcmp( stem_s, stem_m ) == 0 &&
                 !same_inputs( prev_map, s ) ) {
                /* The comparison pairs rows by index, so it means nothing
                   unless the two forms were sampled at the same points. */
                fprintf( stderr,
                    "MISMATCH [equiv-inputs] %s: the .map and .simp sections "
                    "were sampled at different inputs\n", prev_map->relpath );
                total_fail++;
            } else if ( strcmp( stem_s, stem_m ) == 0 ) {
                int nagree, nstale;
                double ertol, eatol;
                equiv_tol( prev_map->relpath, &ertol, &eatol );
                total_fail += compare_equiv( prev_map->relpath, s->in, s->nin,
                                             prev_map_live, live,
                                             prev_map->out, s->out,
                                             s->nout, s->npoint,
                                             ertol, eatol, is_ang,
                                             equiv_accept_rows( prev_map->relpath ),
                                             &nagree, &nstale );
                total_fail += nstale;

                /* A pair whose recorded outputs agree, with a value, at too
                   few inputs is barely compared, so it fails unless the
                   overrides file accepts it, and an override for a pair
                   that is compared is stale.  Both depend only on the
                   recorded values. */
                int waived = equiv_off( prev_map->relpath );
                if ( nagree < MIN_CHECKED_INPUTS && !waived ) {
                    fprintf( stderr,
                        "MISMATCH [equiv-never] %s: the recorded .map and "
                        ".simp outputs agree with a value at %d sampled "
                        "input(s), fewer than %d, so their equivalence is "
                        "barely checked\n", prev_map->relpath, nagree,
                        MIN_CHECKED_INPUTS );
                    total_fail++;
                } else if ( nagree < MIN_CHECKED_INPUTS ) {
                    fprintf( stderr, "SKIP [equiv] %s (overrides file)\n",
                             prev_map->relpath );
                } else if ( waived ) {
                    fprintf( stderr,
                        "MISMATCH [equiv-override] %s: agrees with a value at "
                        "%d sampled inputs, so its \"equiv off\" override is "
                        "stale\n", prev_map->relpath, nagree );
                    total_fail++;
                }
            }
        }

        /* A .map with no .simp partner is checked only against its recorded
           outputs, so those must give a value at enough inputs to pin
           something; the same rules apply as to a pair. */
        if ( is_fwd && has_ext( s->relpath, ".map" ) &&
             !has_simp_partner( secs, nsec, s->relpath ) ) {
            int nvalued = valued_inputs( s );
            int waived = values_off( s->relpath );
            if ( nvalued < MIN_CHECKED_INPUTS && !waived ) {
                fprintf( stderr,
                    "MISMATCH [values-never] %s: the recorded outputs have a "
                    "value at %d sampled input(s), fewer than %d, so the "
                    "golden check pins almost nothing but BAD\n",
                    s->relpath, nvalued, MIN_CHECKED_INPUTS );
                total_fail++;
            } else if ( nvalued < MIN_CHECKED_INPUTS ) {
                fprintf( stderr, "SKIP [values] %s (overrides file)\n",
                         s->relpath );
            } else if ( waived ) {
                fprintf( stderr,
                    "MISMATCH [values-override] %s: has a value at %d "
                    "sampled inputs, so its \"values off\" override is "
                    "stale\n", s->relpath, nvalued );
                total_fail++;
            }
        }

        free( is_ang );
        if ( base ) base = astAnnul( base );
        if ( cur )  cur  = astAnnul( cur );

        /* Refresh the cache on a .map forward section; otherwise free live. */
        if ( is_fwd && has_ext( s->relpath, ".map" ) ) {
            free_cols( prev_map_live, prev_map ? prev_map->nout : 0 );
            prev_map = s;
            prev_map_live = live;           /* hand ownership to the cache */
        } else {
            free_cols( live, s->nout );
        }

        map = astAnnul( map );
        astEnd;
    }
    free_cols( prev_map_live, prev_map ? prev_map->nout : 0 );
    free_sections( secs, nsec );

    printf( "check_transform_oracle: %d sections checked, %d mismatch(es)\n",
            checked, total_fail );
    return ( total_fail == 0 && astOK ) ? 0 : 1;
}
