/*
 * transform_oracle_util.c
 *
 * Shared helpers for gen_transform_oracle and check_transform_oracle:
 * tolerance comparison, low-discrepancy sampling, fixture loading, and
 * full-precision double formatting/parsing.  Public AST API only.
 */
#include "transform_oracle.h"
#include <math.h>
#include <stddef.h>

int oracle_within_tol( double got, double ref, double rtol, double atol ) {
    int got_bad = ( got == AST__BAD );
    int ref_bad = ( ref == AST__BAD );
    if ( got_bad || ref_bad ) return ( got_bad && ref_bad );
    if ( isnan( got ) || isnan( ref ) ) return 0;
    return fabs( got - ref ) <= atol + rtol * fabs( ref );
}

/* Number of sampled points per fixture (tunable). */
#define ORACLE_NPOINT 24

/* Distinct prime bases per axis; extended cyclically beyond the table. */
static const unsigned oracle_primes[] = {
    2u, 3u, 5u, 7u, 11u, 13u, 17u, 19u, 23u, 29u, 31u, 37u
};

double oracle_halton( unsigned index, unsigned base ) {
    double f = 1.0, r = 0.0;
    while ( index > 0u ) {
        f /= (double) base;
        r += f * (double) ( index % base );
        index /= base;
    }
    return r;
}

int oracle_sample_axis_count( void ) { return ORACLE_NPOINT; }

void oracle_sample_points( int naxis, const double *lo, const double *hi,
                           int npoint, double **out ) {
    /* First four rows are deterministic edge/structure points. */
    for ( int a = 0; a < naxis; a++ ) {
        double zero = 0.0 < lo[a] ? lo[a] : ( 0.0 > hi[a] ? hi[a] : 0.0 );
        if ( npoint > 0 ) out[a][0] = lo[a];                 /* all-lo  */
        if ( npoint > 1 ) out[a][1] = hi[a];                 /* all-hi  */
        if ( npoint > 2 ) out[a][2] = zero;                  /* zero    */
        if ( npoint > 3 ) out[a][3] = 0.5 * ( lo[a] + hi[a] ); /* mid   */
    }
    /* Remaining rows: Halton spread, distinct prime base per axis. */
    int nedge = npoint < 4 ? npoint : 4;
    for ( int a = 0; a < naxis; a++ ) {
        unsigned base = oracle_primes[ (size_t) a % ( sizeof(oracle_primes)
                                       / sizeof(oracle_primes[0]) ) ];
        for ( int i = nedge; i < npoint; i++ ) {
            double u = oracle_halton( (unsigned)( i - nedge + 1 ), base );
            out[a][i] = lo[a] + u * ( hi[a] - lo[a] );
        }
    }
}
