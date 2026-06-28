/*
 * transform_oracle_util.c
 *
 * Shared helpers for gen_transform_oracle and check_transform_oracle:
 * tolerance comparison, low-discrepancy sampling, fixture loading, and
 * full-precision double formatting/parsing.  Public AST API only.
 */
#include "transform_oracle.h"
#include <math.h>

int oracle_within_tol( double got, double ref, double rtol, double atol ) {
    int got_bad = ( got == AST__BAD );
    int ref_bad = ( ref == AST__BAD );
    if ( got_bad || ref_bad ) return ( got_bad && ref_bad );
    if ( isnan( got ) || isnan( ref ) ) return 0;
    return fabs( got - ref ) <= atol + rtol * fabs( ref );
}
