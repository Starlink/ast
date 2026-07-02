/* Enable POSIX 2008 (clock_gettime) and glibc extensions (drand48). */
#define _DEFAULT_SOURCE

/*
 * bench_simd.c
 *
 * Benchmark AST transforms at varying N to evaluate SIMD vectorization gains.
 *
 * Benchmarks (each run scalar, and SIMD when the library supports it):
 *   sphmap_fwd: SphMap (x,y,z) -> (lon,lat): atan2/sqrt via libmvec
 *   sphmap_inv: SphMap (lon,lat) -> (x,y,z): sin/cos via libmvec
 *   poly5_fwd: degree-5 2-D PolyMap (Roman WCS SIP-like distortion)
 *   poly1_fwd: degree-1 2-D PolyMap (Roman WCS rotation/scale)
 *   matdiag2_fwd: 2x2 diagonal MatrixMap (scale)
 *   matfull2_fwd: 2x2 full MatrixMap (rotation)
 *
 * CSV output (transform column uses _simd / _scalar suffix):
 *   transform,n_points,rep,time_s
 * where time_s is the amortised time for one transform of n_points points
 * (each rep times a calibrated batch of repeated calls; see run_sweep).
 *
 * Usage:
 *   bench_simd [-o output.csv] [-r nreps] [--markdown]
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "ast.h"

/* N sweep: 1 through 4088x4088 (full Roman WFI science area). */
static const size_t N_SWEEP[] = {
    1, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 4096, 16384, 65536,
    262144, 1048576, 4194304, 16711744
};
static const size_t N_N_SWEEP = sizeof(N_SWEEP) / sizeof(N_SWEEP[0]);

#define RAND_SEED 42UL
#define DEFAULT_REPS 5
#define IMAGE_HALF 2048.0   /* pixel coord half-range for poly/matrix */

/*
 * At small N a single astTranP call is dominated by call overhead and clock
 * resolution rather than transform throughput. Each measurement therefore
 * times a batch of repeated calls, growing the batch until it runs for at
 * least MIN_BATCH_S, and reports the amortised per-call time. MAX_BATCH_CALLS
 * bounds the batch so the tiniest N values do not run unboundedly.
 */
#define MIN_BATCH_S 0.02
#define MAX_BATCH_CALLS ( 1 << 22 )

/*
 * Degree-5 2-D PolyMap (Roman WCS SIP-like distortion, 10 terms).
 * Row format: { coeff, out_coord (1-based), pow_x, pow_y }.
 */
#define POLY5_NCOEFF 10
static const double POLY5_COEFF_F[POLY5_NCOEFF * 4] = {
     0.11034133, 1, 1, 0,
    -3.0e-8,     1, 2, 0,
     3.4168e-4,  1, 0, 1,
    -1.0e-8,     1, 0, 2,
     1.4e-7,     1, 1, 1,
     3.1436e-4,  2, 1, 0,
     7.0e-8,     2, 2, 0,
     0.10828278, 2, 0, 1,
     2.1e-7,     2, 0, 2,
    -2.0e-8,     2, 1, 1,
};

/* Degree-1 2-D PolyMap (Roman WCS linear transform, 4 terms). */
#define POLY1_NCOEFF 4
static const double POLY1_COEFF_F[POLY1_NCOEFF * 4] = {
    -0.5,       1, 1, 0,
    -0.8660254, 1, 0, 1,
    -0.8660254, 2, 1, 0,
     0.5,       2, 0, 1,
};

/* 2x2 diagonal scale */
static const double DIAG2X2[2] = { 2.0, 0.5 };

/*
 * 2x2 full rotation matrix (same linear transform as poly1_fwd, row-major):
 */
static const double ROT2X2[4] = {
    -0.5,       -0.8660254,
    -0.8660254,  0.5
};

#define MAX_SUM_ENTRIES 512

typedef struct {
    char base[64]; /* transform name without _simd/_scalar suffix */
    size_t N;
    int is_simd;   /* 0 = scalar, 1 = SIMD */
    double mpxs;   /* median throughput in Mpx/s */
} SumEntry;

static SumEntry g_sum[MAX_SUM_ENTRIES];
static int g_nsum = 0;

/*
 * Whether the loaded libast supports runtime SIMD. Determined at startup by
 * detect_simd() rather than a compile-time macro, since UseSIMD is a runtime
 * tuning knob: the benchmark need not be rebuilt to match the library.
 */
static int g_have_simd = 0;

static void record_sum( const char *base, size_t N, int is_simd,
                        double median_s ) {
    if( g_nsum >= MAX_SUM_ENTRIES )
        return;
    strncpy( g_sum[g_nsum].base, base, sizeof(g_sum[0].base) - 1 );
    g_sum[g_nsum].base[ sizeof(g_sum[0].base) - 1 ] = '\0';
    g_sum[g_nsum].N = N;
    g_sum[g_nsum].is_simd = is_simd;
    g_sum[g_nsum].mpxs = ( median_s > 0.0 ) ? (double)N / median_s / 1e6 : 0.0;
    g_nsum++;
}

static double lookup_mpxs( const char *base, size_t N, int is_simd ) {
    for( int i = 0; i < g_nsum; i++ ) {
        if( strcmp( g_sum[i].base, base ) == 0
            && g_sum[i].N == N
            && g_sum[i].is_simd == is_simd )
            return g_sum[i].mpxs;
    }
    return -1.0;
}

/* Return the set of unique base labels in insertion order. */
static void unique_bases( char bases[][64], int *n ) {
    *n = 0;
    for( int i = 0; i < g_nsum; i++ ) {
        int found = 0;
        for( int j = 0; j < *n; j++ ) {
            if( strcmp( bases[j], g_sum[i].base ) == 0 ) {
                found = 1;
                break;
            }
        }
        if( !found && *n < MAX_SUM_ENTRIES )
            strncpy( bases[(*n)++], g_sum[i].base, 63 );
    }
}

static double now_s( void ) {
    struct timespec ts;
    clock_gettime( CLOCK_MONOTONIC, &ts );
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

/* Run `ncalls` transforms of N points and return the elapsed seconds, or
   -1.0 if a transform failed. */
static double time_batch( AstMapping *map, int ncalls, int forward, size_t N,
                          int nin, const double **ptr_in,
                          int nout, double **ptr_out ) {
    double t0 = now_s();
    for( int k = 0; k < ncalls; k++ ) {
        astTranP( map, (int)N, nin, ptr_in, forward, nout, ptr_out );
        if( !astOK )
            return -1.0;
    }
    return now_s() - t0;
}

/*
 * Run the N-sweep for one direction.
 * Records medians in g_sum for the summary table.
 */
static int run_sweep( FILE *fout, AstMapping *map,
                      const char *base_label, int forward,
                      int nin, const double **ptr_in,
                      int nout, double **ptr_out,
                      int nreps, double *rep_times,
                      int use_simd ) {
    const char *suffix = use_simd ? "_simd" : "_scalar";
    char full_label[80];
    snprintf( full_label, sizeof(full_label), "%s%s", base_label, suffix );

    /* Apply UseSIMD setting (use_simd=1 is only ever requested when the
       library reports SIMD support, so the astSet always succeeds). */
    astSet( (AstObject *) map, use_simd ? "UseSIMD=1" : "UseSIMD=0" );
    if( !astOK ) {
        fprintf( stderr, "error: astSet(UseSIMD) failed\n" );
        return 1;
    }

    fprintf( stderr, "=== %s ===\n", full_label );

    for( size_t ni = 0; ni < N_N_SWEEP; ni++ ) {
        size_t N = N_SWEEP[ni];

        /* Calibrate the batch size so a batch runs for at least MIN_BATCH_S.
           This also serves as a warm-up for the rep measurements below. */
        int ncalls = 1;
        for( ;; ) {
            double dt = time_batch( map, ncalls, forward, N,
                                    nin, ptr_in, nout, ptr_out );
            if( dt < 0.0 ) {
                fprintf( stderr, "error: astTranP failed\n" );
                return 1;
            }
            if( dt >= MIN_BATCH_S || ncalls >= MAX_BATCH_CALLS )
                break;
            ncalls *= 2;
        }

        for( int rep = 0; rep < nreps; rep++ ) {
            double dt = time_batch( map, ncalls, forward, N,
                                    nin, ptr_in, nout, ptr_out );
            if( dt < 0.0 ) {
                fprintf( stderr, "error: astTranP failed\n" );
                return 1;
            }

            /* Amortised per-call time (one call transforms N points). */
            rep_times[rep] = dt / ncalls;
            fprintf( fout, "%s,%zu,%d,%.9e\n", full_label, N, rep, rep_times[rep] );
        }
        fflush( fout );

        /* Sort reps */
        for( int a = 0; a < nreps - 1; a++ ) {
            for( int b = a + 1; b < nreps; b++ ) {
                if( rep_times[b] < rep_times[a] ) {
                    double tmp = rep_times[a];
                    rep_times[a] = rep_times[b];
                    rep_times[b] = tmp;
                }
            }
        }
        double median = rep_times[ nreps / 2 ];

        record_sum( base_label, N, use_simd, median );

        fprintf( stderr,
                 "  %-30s  N=%-8zu  x%-8d  median=%10.6f ms  (%6.1f Mpx/s)\n",
                 full_label, N, ncalls, median * 1e3, (double)N / median / 1e6 );
    }
    return 0;
}

static int run_both( FILE *fout, AstMapping *map,
                     const char *base_label, int forward,
                     int nin, const double **ptr_in,
                     int nout, double **ptr_out,
                     int nreps, double *rep_times ) {
    int rc;

    /* Scalar run always. */
    rc = run_sweep( fout, map, base_label, forward,
                    nin, ptr_in, nout, ptr_out, nreps, rep_times, 0 );
    if( rc )
        return rc;

    /* SIMD run only when the library reports SIMD support. */
    if( g_have_simd )
        rc = run_sweep( fout, map, base_label, forward,
                        nin, ptr_in, nout, ptr_out, nreps, rep_times, 1 );
    return rc;
}

static int bench_sphmap( FILE *fout, AstSphMap *sphmap,
                         const double *x, const double *y, const double *z,
                         const double *lon, const double *lat,
                         double *out0, double *out1, double *out2,
                         int nreps, double *rep_times ) {
    const double *fwd_in[3] = { x, y, z };
    double *fwd_out[2] = { out0, out1 };
    const double *inv_in[2] = { lon, lat };
    double *inv_out[3] = { out0, out1, out2 };
    int rc;

    rc = run_both( fout, (AstMapping *)sphmap, "sphmap_fwd", 1,
                   3, fwd_in, 2, fwd_out, nreps, rep_times );
    if( rc )
        return rc;

    return run_both( fout, (AstMapping *)sphmap, "sphmap_inv", 0,
                     2, inv_in, 3, inv_out, nreps, rep_times );
}

static int bench_polymap( FILE *fout, const char *base_label,
                          int ncoeff, const double *coeff_f,
                          const double *in0, const double *in1,
                          double *out0, double *out1,
                          int nreps, double *rep_times ) {
    AstPolyMap *polymap = astPolyMap( 2, 2, ncoeff, coeff_f, 0, NULL, " " );
    if( !astOK ) {
        fprintf( stderr, "error: failed to create PolyMap '%s'\n", base_label );
        return 1;
    }

    const double *ptr_in[2] = { in0, in1 };
    double *ptr_out[2] = { out0, out1 };

    int rc = run_both( fout, (AstMapping *)polymap, base_label, 1,
                       2, ptr_in, 2, ptr_out, nreps, rep_times );
    astAnnul( polymap );
    return rc;
}

static int bench_matrixmap( FILE *fout,
                            const double *in0, const double *in1,
                            double *out0, double *out1,
                            int nreps, double *rep_times ) {
    const double *ptr_in[2] = { in0, in1 };
    double *ptr_out[2] = { out0, out1 };
    int rc;

    /* Diagonal 2x2. */
    AstMatrixMap *diag = astMatrixMap( 2, 2, 1, DIAG2X2, " " );
    if( !astOK ) {
        fprintf( stderr, "error: failed to create diagonal MatrixMap\n" );
        return 1;
    }
    rc = run_both( fout, (AstMapping *)diag, "matdiag2_fwd", 1,
                   2, ptr_in, 2, ptr_out, nreps, rep_times );
    astAnnul( diag );
    if( rc )
        return rc;

    /* Full 2x2 rotation. */
    AstMatrixMap *rot = astMatrixMap( 2, 2, 0, ROT2X2, " " );
    if( !astOK ) {
        fprintf( stderr, "error: failed to create rotation MatrixMap\n" );
        return 1;
    }
    rc = run_both( fout, (AstMapping *)rot, "matfull2_fwd", 1,
                   2, ptr_in, 2, ptr_out, nreps, rep_times );
    astAnnul( rot );
    return rc;
}

/*
 * Probe whether the loaded libast supports runtime SIMD. The default value of
 * the UseSIMD attribute is 1 when SIMD is compiled in and 0 otherwise, so
 * reading it from a fresh mapping reports support without triggering an error.
 */
static int detect_simd( void ) {
    int have = 0;
    AstMatrixMap *probe = astMatrixMap( 2, 2, 1, DIAG2X2, " " );
    if( astOK )
        have = astGetI( (AstObject *) probe, "UseSIMD" );
    if( probe )
        probe = astAnnul( probe );
    return have;
}

static void print_summary( void ) {
    char bases[MAX_SUM_ENTRIES][64];
    int nbases;

    unique_bases( bases, &nbases );

    printf( "\n=== Summary ===\n" );
    if( g_have_simd ) {
        printf( "%-22s  %8s  %14s  %12s  %7s\n",
                "Transform", "N", "Scalar Mpx/s", "SIMD Mpx/s", "Speedup" );
        printf( "%-22s  %8s  %14s  %12s  %7s\n",
                "----------------------", "--------", "--------------",
                "------------", "-------" );
    } else {
        printf( "%-22s  %8s  %14s\n",
                "Transform", "N", "Scalar Mpx/s" );
        printf( "%-22s  %8s  %14s\n",
                "----------------------", "--------", "--------------" );
    }

    for( int b = 0; b < nbases; b++ ) {
        for( size_t ni = 0; ni < N_N_SWEEP; ni++ ) {
            size_t N = N_SWEEP[ni];
            double sc = lookup_mpxs( bases[b], N, 0 );
            if( sc < 0.0 )
                continue;

            printf( "%-22s  %8zu  %14.1f", bases[b], N, sc );
            if( g_have_simd ) {
                double sm = lookup_mpxs( bases[b], N, 1 );
                double sp = ( sm > 0.0 && sc > 0.0 ) ? sm / sc : 0.0;
                printf( "  %12.1f  %6.2fx", sm, sp );
            }
            printf( "\n" );
        }
    }
}

static void print_summary_markdown( void ) {
    char bases[MAX_SUM_ENTRIES][64];
    int nbases;

    unique_bases( bases, &nbases );

    printf( "\n### SIMD benchmark summary\n\n" );
    if( g_have_simd )
        printf( "| Transform | N | Scalar Mpx/s | SIMD Mpx/s | Speedup |\n"
                "|-----------|--:|-------------:|-----------:|--------:|\n" );
    else
        printf( "| Transform | N | Scalar Mpx/s |\n"
                "|-----------|--:|-------------:|\n" );

    for( int b = 0; b < nbases; b++ ) {
        for( size_t ni = 0; ni < N_N_SWEEP; ni++ ) {
            size_t N = N_SWEEP[ni];
            double sc = lookup_mpxs( bases[b], N, 0 );
            if( sc < 0.0 )
                continue;

            printf( "| %s | %zu | %.1f |", bases[b], N, sc );
            if( g_have_simd ) {
                double sm = lookup_mpxs( bases[b], N, 1 );
                double sp = ( sm > 0.0 && sc > 0.0 ) ? sm / sc : 0.0;
                printf( " %.1f | %.2fx |", sm, sp );
            }
            printf( "\n" );
        }
    }
}

int main( int argc, char *argv[] ) {
    const char *outpath = NULL;
    int nreps = DEFAULT_REPS;
    int markdown = 0;
    FILE *fout = stdout;
    int status = 0;
    int rc = 0;
    size_t max_n = N_SWEEP[N_N_SWEEP - 1];
    double *x = NULL;
    double *y = NULL;
    double *z = NULL;
    double *lon = NULL;
    double *lat = NULL;
    double *out0 = NULL;
    double *out1 = NULL;
    double *out2 = NULL;
    double *rep_times = NULL;
    AstSphMap *sphmap = NULL;

    for( int i = 1; i < argc; i++ ) {
        if( strcmp( argv[i], "-o" ) == 0 && i + 1 < argc ) {
            outpath = argv[++i];
        } else if( strcmp( argv[i], "-r" ) == 0 && i + 1 < argc ) {
            nreps = atoi( argv[++i] );
            if( nreps < 1 )
                nreps = 1;
        } else if( strcmp( argv[i], "--markdown" ) == 0 ) {
            markdown = 1;
        } else {
            fprintf( stderr, "Usage: %s [-o output.csv] [-r nreps] [--markdown]\n",
                     argv[0] );
            return 1;
        }
    }

    if( outpath ) {
        fout = fopen( outpath, "w" );
        if( !fout ) {
            fprintf( stderr, "error: cannot open %s for writing\n", outpath );
            return 1;
        }
    }

    astWatch( &status );
    astBegin;

    g_have_simd = detect_simd();
    fprintf( stderr, "SIMD support: %s\n", g_have_simd ? "yes" : "no" );

    sphmap = astSphMap( "UnitRadius=1" );
    if( !astOK ) {
        fprintf( stderr, "error: failed to create SphMap (status=%d)\n", status );
        rc = 1;
        goto done;
    }

    x = malloc( max_n * sizeof(double) );
    y = malloc( max_n * sizeof(double) );
    z = malloc( max_n * sizeof(double) );
    lon = malloc( max_n * sizeof(double) );
    lat = malloc( max_n * sizeof(double) );
    out0 = malloc( max_n * sizeof(double) );
    out1 = malloc( max_n * sizeof(double) );
    out2 = malloc( max_n * sizeof(double) );
    rep_times = malloc( (size_t)nreps * sizeof(double) );

    if( !x || !y || !z || !lon || !lat || !out0 || !out1 || !out2 || !rep_times ) {
        fprintf( stderr, "error: out of memory\n" );
        rc = 1;
        goto done;
    }

    /* Random unit-sphere Cartesian points and matching spherical coords. */
    srand48( (long)RAND_SEED );
    for( size_t i = 0; i < max_n; i++ ) {
        double lo = drand48() * 2.0 * M_PI;
        double la = (drand48() - 0.5) * M_PI;
        double cp = cos( la );
        x[i] = cos( lo ) * cp;
        y[i] = sin( lo ) * cp;
        z[i] = sin( la );
        lon[i] = lo;
        lat[i] = la;
    }

    fprintf( fout, "transform,n_points,rep,time_s\n" );

    rc = bench_sphmap( fout, sphmap, x, y, z, lon, lat,
                       out0, out1, out2, nreps, rep_times );
    if( rc )
        goto done;

    /* Reinitialise x/y as centred pixel coordinates. */
    srand48( (long)RAND_SEED );
    for( size_t i = 0; i < max_n; i++ ) {
        x[i] = (drand48() - 0.5) * 2.0 * IMAGE_HALF;
        y[i] = (drand48() - 0.5) * 2.0 * IMAGE_HALF;
    }

    rc = bench_polymap( fout, "poly5_fwd", POLY5_NCOEFF, POLY5_COEFF_F,
                        x, y, out0, out1, nreps, rep_times );
    if( rc )
        goto done;

    rc = bench_polymap( fout, "poly1_fwd", POLY1_NCOEFF, POLY1_COEFF_F,
                        x, y, out0, out1, nreps, rep_times );
    if( rc )
        goto done;

    rc = bench_matrixmap( fout, x, y, out0, out1, nreps, rep_times );

done:
    free( rep_times );
    free( x );
    free( y );
    free( z );
    free( lon );
    free( lat );
    free( out0 );
    free( out1 );
    free( out2 );

    astEnd;

    if( outpath )
        fclose( fout );

    if( !rc && !status ) {
        if( markdown )
            print_summary_markdown();
        else
            print_summary();
    }

    return (rc || status) ? 1 : 0;
}
