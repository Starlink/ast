/* Shared declarations for the transform-output oracle programs. */
#ifndef TRANSFORM_ORACLE_H
#define TRANSFORM_ORACLE_H

#include "ast.h"
#include <stddef.h>

#define ORACLE_BAD_TOKEN "BAD"

/* Golden comparison tolerance.  Set to absorb cross-architecture
   floating-point differences while still flagging real algorithmic change
   (orders of magnitude larger, >= ~1e-6 relative).  ARM vs x86 `libm`
   differences reach ~5e-9 for transcendental inverses sampled near a
   singularity (e.g. a relativistic spectral inverse at beta ~ 1), well
   above the ~2e-12 seen merely across compilers on one machine. */
#define ORACLE_DEF_RTOL        1e-7
#define ORACLE_DEF_ATOL        1e-7
#define ORACLE_DEF_EQUIV_RTOL  1e-9
#define ORACLE_DEF_EQUIV_ATOL  1e-9

/* Round-trip accuracy (inverse(forward(P)) vs P) for GRID-domain corpora.
   The inputs are pixel coordinates, so the criterion is an *absolute*
   pixel tolerance (a relative bound is meaningless near pixel 1): the
   inverse must recover the pixel to better than this many pixels.  Loose
   enough to admit iterative distortion inverses (SIP ~0.01 px), tight
   enough to flag inverses that do not recover the pixel at all. */
#define ORACLE_DEF_RTRIP_RTOL  1e-6
#define ORACLE_DEF_RTRIP_ATOL  5e-2

/* Return 1 if got and ref agree within |got-ref| <= atol + rtol*|ref|.
   AST__BAD matches only AST__BAD; any NaN never matches. */
int oracle_within_tol( double got, double ref, double rtol, double atol );

/* As oracle_within_tol, but also matches angles (radians) that are equal
   modulo 2*pi -- a longitude may emerge as 0 vs 2*pi or +pi vs -pi across
   architectures or normalization conventions.  Use for angle-bearing
   outputs (golden, equivalence); not for pixel comparisons (round-trip). */
int oracle_within_tol_wrap( double got, double ref, double rtol, double atol );

double oracle_halton( unsigned index, unsigned base );
int    oracle_sample_axis_count( void );
void   oracle_sample_points( int naxis, const double *lo, const double *hi,
                             int npoint, double **out );

AstMapping *oracle_load_mapping( const char *root, const char *relpath );
void   oracle_format_double( char *buf, size_t buflen, double v );
double oracle_parse_double( const char *tok, int *ok );

#endif
