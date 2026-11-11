/*
 * testslamap: SlaMap regression tests (C only; no Fortran original).
 *
 * An SlaMap conversion whose stored argument is AST__BAD (e.g. an
 * ecliptic conversion with an undefined epoch) cannot be evaluated: its
 * transform must return AST__BAD for every point rather than propagating
 * NaN through the underlying PAL routines, and must not raise an AST
 * error.
 *
 * The FK4/FK5 conversions are checked against ERFA 2.0.1's eraFk54z and
 * eraFk45z, whose outputs for the same arguments are the wanted values.
 *
 * PREC and EQECL are checked against values formed from ERFA 2.0.1's
 * IAU 2006 precession-only matrices (the rp output of eraBp06, as
 * returned by pyerfa's bp06): PREC from EP0 to EP1 is rp(EP1)*rp(EP0)^T,
 * and EQECL at a date is Rx(eraObl06(date)) * rp(epoch of date), the
 * mean ecliptic of date with PAL's IAU 2006 obliquity.  These include no frame bias, so PREC is the identity
 * when its epochs are equal, and an epoch of exactly 2000.0 is not a
 * special case.
 */
#include "ast.h"
#include <math.h>
#include <stdio.h>

static int failed = 0;

/* Transform (ra, dec) radians forward through a one-step SlaMap with a
   single argument and compare the result with the wanted values. */
static void check_fk( const char *cvt, double bepoch, double ra, double dec,
                      double ra_want, double dec_want, double tol ) {
   AstSlaMap *sm;
   double args[ 1 ];
   double raout, decout;

   astBegin;
   sm = astSlaMap( 0, " " );
   args[ 0 ] = bepoch;
   astSlaAdd( sm, cvt, 1, args );
   astTran2( sm, 1, &ra, &dec, 1, &raout, &decout );
   if( !astOK ) {
      astClearStatus;
      printf( "Error: %s at B%g raised an AST error\n", cvt, bepoch );
      failed = 1;
   } else if( fabs( raout - ra_want ) > tol ||
              fabs( decout - dec_want ) > tol ) {
      printf( "Error: %s at B%g gives (%.17g,%.17g), expected "
              "(%.17g,%.17g)\n", cvt, bepoch, raout, decout, ra_want,
              dec_want );
      failed = 1;
   }
   astEnd;
}

/* Transform (ra, dec) radians forward through a one-step SlaMap with the
   supplied arguments and compare the result with the wanted values. */
static void check_cvt( const char *cvt, int nargs, const double *args,
                       double ra, double dec, double ra_want,
                       double dec_want, double tol ) {
   AstSlaMap *sm;
   double raout, decout;

   astBegin;
   sm = astSlaMap( 0, " " );
   astSlaAdd( sm, cvt, nargs, args );
   astTran2( sm, 1, &ra, &dec, 1, &raout, &decout );
   if( !astOK ) {
      astClearStatus;
      printf( "Error: %s (%g, ...) raised an AST error\n", cvt, args[ 0 ] );
      failed = 1;
   } else if( fabs( raout - ra_want ) > tol ||
              fabs( decout - dec_want ) > tol ) {
      printf( "Error: %s (%g, %g) gives (%.17g,%.17g), expected "
              "(%.17g,%.17g)\n", cvt, args[ 0 ], nargs > 1 ? args[ 1 ] : 0.0,
              raout, decout, ra_want, dec_want );
      failed = 1;
   }
   astEnd;
}

static void check_prec( double ep0, double ep1, double ra_want,
                        double dec_want ) {
   double args[ 2 ];
   args[ 0 ] = ep0;
   args[ 1 ] = ep1;
   check_cvt( "PREC", 2, args, 1.0, 0.5, ra_want, dec_want, 1.0E-14 );
}

int main( void ) {
   int status_value = 0;
   int *status = &status_value;
   double args[ 1 ];
   double xin[ 1 ] = { 1.0 };
   double yin[ 1 ] = { 0.5 };
   double xout[ 1 ], yout[ 1 ];

   astWatch( status );
   astBegin;

   AstSlaMap *sm = astSlaMap( 0, " " );
   args[ 0 ] = AST__BAD;
   astSlaAdd( sm, "EQECL", 1, args );

   astTran2( sm, 1, xin, yin, 1, xout, yout );

   if( !astOK ) {
      printf( "Error 1: transform with BAD argument raised an AST error\n" );
   } else if( xout[ 0 ] != AST__BAD || yout[ 0 ] != AST__BAD ) {
      printf( "Error 2: got (%.17g,%.17g), expected AST__BAD\n",
              xout[ 0 ], yout[ 0 ] );
      *status = 1;
   }

   astEnd;

/* FK5 to FK4 away from B1950.  ERFA, like the Fortran SLALIB original,
   applies the fictitious proper motion to the Cartesian unit vector, where
   PAL's C version applied it linearly in (RA, Dec).  The two differ most
   for points carried close to a pole: at (0.6, -1.55) and B1980 by 7.1e-10
   radians in RA, and at (1.0, 0.5) by only 1e-14.  At B1950 they agree to
   rounding, and FK4 to FK5 differs only by rounding. */
   if( *status == 0 ) {
      check_fk( "FK54Z", 1980.0, 0.6, -1.55, 0.74967952758676781,
                -1.5538073132385521, 1.0E-12 );
      check_fk( "FK54Z", 1950.0, 1.0, 0.5, 0.98660250110964187,
                0.49734772100081265, 1.0E-14 );
      check_fk( "FK45Z", 1980.0, 1.0, 0.5, 1.0134307040454991,
                0.50259711081913927, 1.0E-14 );
      check_fk( "FK45Z", 1950.0, 1.0, 0.5, 1.0134306526650936,
                0.50259743655625932, 1.0E-14 );
      if( failed ) *status = 1;
   }

/* IAU 2006 precession.  A frame bias of up to 23 mas (1.1e-7 radians) or
   a reversed matrix product (6 mas between 1950 and 2050) would fail
   these by far more than the tolerance. */
   if( *status == 0 ) {
      check_prec( 2000.0, 2000.0, 1.0000000000000002, 0.5000000000000001 );
      check_prec( 2000.0, 2025.0, 1.0067112795780082, 0.5013055113074502 );
      check_prec( 2010.0, 2000.0, 0.9973178941015486, 0.4994739358685401 );
      check_prec( 1950.0, 2050.0, 1.0268929459157925, 0.5051393318180102 );

/* J2000 equatorial to the mean ecliptic of MJD 58849, which precesses
   from J2000 to the epoch of the date. */
      args[ 0 ] = 58849.0;
      check_cvt( "EQECL", 1, args, 1.0, 0.5, 1.0758075590444474,
                 0.14668872937582086, 1.0E-14 );
      if( failed ) *status = 1;
   }

   if( *status == 0 ) {
      printf( " All SlaMap tests passed\n" );
   } else {
      printf( "SlaMap tests failed\n" );
   }
   return *status;
}
