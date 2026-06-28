/* Shared declarations for the transform-output oracle programs. */
#ifndef TRANSFORM_ORACLE_H
#define TRANSFORM_ORACLE_H

#include "ast.h"

#define ORACLE_BAD_TOKEN "BAD"

/* Default comparison tolerances (tuned in the tolerance-tuning task). */
#define ORACLE_DEF_RTOL        1e-12
#define ORACLE_DEF_ATOL        1e-12
#define ORACLE_DEF_EQUIV_RTOL  1e-9
#define ORACLE_DEF_EQUIV_ATOL  1e-9

/* Return 1 if got and ref agree within |got-ref| <= atol + rtol*|ref|.
   AST__BAD matches only AST__BAD; any NaN never matches. */
int oracle_within_tol( double got, double ref, double rtol, double atol );

#endif
