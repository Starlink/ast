/*
*  Purpose:
*     Exercise WcsMap projections at the points where the wcslib code AST
*     carries has differed from upstream wcslib: face centres, poles,
*     longitude edges and non-default projection parameters.
*
*  Notes:
*     - Public C API only, with astWatch/astOK for error checking.
*     - A WcsMap works in radians on both sides.  The cases here are written
*       in degrees, as wcslib's own tests are, and converted at the call.
*     - Unless a comment says otherwise, every wanted value is the output of
*       upstream wcslib 8.10's prjs2x or prjx2s for the same projection and
*       parameters, printed with %.17g.
*/
#include "ast.h"

#include <math.h>
#include <stdio.h>

/* "status" is the AST status watched by astWatch; "failed" records a test
   failure, so one failed case does not stop AST working for the next. */
static int status = 0;
static int failed = 0;

static void fail( const char *text ) {
   failed = 1;
   printf( "%s\n", text );
}

/* Compare one output (radians, or AST__BAD) with a wanted value (degrees,
   or AST__BAD to require a bad output). */
static void check_value( const char *label, const char *which, double got,
                         double want, double tol ) {
   char text[ 200 ];

   if( want == AST__BAD ) {
      if( got != AST__BAD ) {
         snprintf( text, sizeof text, "testwcsmap: %s: %s is %.17g, wanted "
                   "bad", label, which, got*AST__DR2D );
         fail( text );
      }
   } else if( got == AST__BAD ) {
      snprintf( text, sizeof text, "testwcsmap: %s: %s is bad, wanted %.17g",
                label, which, want );
      fail( text );
   } else if( fabs( got*AST__DR2D - want ) > tol ) {
      snprintf( text, sizeof text, "testwcsmap: %s: %s is %.17g, wanted "
                "%.17g (tolerance %g)", label, which, got*AST__DR2D, want,
                tol );
      fail( text );
   }
}

/* Transform one point through a 2-axis WcsMap, forward when "forward" is
   non-zero, and check both outputs. */
static void check_tran( const char *label, int type, const char *pv,
                        int forward, double a_deg, double b_deg,
                        double a_want, double b_want, double tol,
                        const char *a_name, const char *b_name ) {
   AstWcsMap *map;
   double ain, bin, aout, bout;

   astBegin;
   map = astWcsMap( 2, type, 1, 2, "%s", pv );
   if( !astOK ) {
      astClearStatus;
      fail( "testwcsmap: could not create the WcsMap" );
      printf( "   (case: %s)\n", label );
      astEnd;
      return;
   }

   ain = a_deg*AST__DD2R;
   bin = b_deg*AST__DD2R;
   astTran2( map, 1, &ain, &bin, forward, &aout, &bout );
   if( !astOK ) {
      astClearStatus;
      fail( "testwcsmap: astTran2 failed" );
      printf( "   (case: %s)\n", label );
   } else {
      check_value( label, a_name, aout, a_want, tol );
      check_value( label, b_name, bout, b_want, tol );
   }
   astEnd;
}

/* Native spherical (phi, theta) to projection plane (x, y). */
static void check_fwd( const char *label, int type, const char *pv,
                       double phi_deg, double theta_deg, double x_want,
                       double y_want, double tol ) {
   check_tran( label, type, pv, 1, phi_deg, theta_deg, x_want, y_want, tol,
               "x", "y" );
}

/* Projection plane (x, y) to native spherical (phi, theta). */
static void check_rev( const char *label, int type, const char *pv,
                       double x_deg, double y_deg, double phi_want,
                       double theta_want, double tol ) {
   check_tran( label, type, pv, 0, x_deg, y_deg, phi_want, theta_want, tol,
               "phi", "theta" );
}

int main( void ) {
   astWatch( &status );

/* QSC, face 4.  Within about 0.008 degrees of the face centre the forward
   projection takes a small-angle branch, which squared the longitude
   offset where it should have taken it once. */
   check_fwd( "qsc face 4 centre", AST__QSC, "", -90.005, 0.0,
              269.99486913923511, 0.0, 1.0E-10 );
   check_rev( "qsc face 4 centre reverse", AST__QSC, "",
              269.99486913923511, 0.0, -90.005, 0.0, 1.0E-9 );
   check_fwd( "qsc face 4 near", AST__QSC, "", -89.995, 0.004,
              270.00573448754045, 0.0047579971900078627, 1.0E-10 );
   check_rev( "qsc face 4 near reverse", AST__QSC, "",
              270.00573448754045, 0.0047579971900078627, -89.995, 0.004,
              1.0E-9 );

/* HPX polar zone with an odd H: the facet centre was computed with
   floor(... + 1/2), an integer division that added nothing. */
   check_fwd( "hpx H=3 K=3 south pole", AST__HPX, "PV2_1=3,PV2_2=3",
              -180.0, -90.0, -120.0, -120.0, 1.0E-10 );
   check_fwd( "hpx H=4 K=4 south pole", AST__HPX, "PV2_1=4,PV2_2=4",
              -131.9, -90.0, -90.0, -112.5, 1.0E-10 );

/* HPX reverse: the polar zone extends to (K+1)/2 facets of 180/H degrees,
   not to a fixed 90 degrees, which is only right for H = 4. */
   check_rev( "hpx H=3 K=3 polar zone", AST__HPX, "PV2_1=3,PV2_2=3",
              -0.4, -119.5, -48.000000000000171, -89.610150847356522,
              1.0E-9 );

/* HPX reverse: a point on a facet edge, where (x - x_c)/sigma reaches
   180/H only up to rounding. */
   check_fwd( "hpx facet edge", AST__HPX, "", 90.0, -89.26,
              134.28818979161747, -89.288189791617469, 1.0E-10 );
   check_rev( "hpx facet edge reverse", AST__HPX, "",
              134.28818979161747, -89.288189791617469, 90.0, -89.26,
              1.0E-9 );

/* HPX at the pole: the forward projection puts it at the apex of the
   facet, and the reverse must accept it rather than reject sigma = 0. */
   check_fwd( "hpx pole phi=-180", AST__HPX, "", -180.0, -90.0,
              -135.0, -90.0, 1.0E-10 );
   check_rev( "hpx pole phi=-180 reverse", AST__HPX, "", -135.0, -90.0,
              -135.0, -90.0, 1.0E-10 );
   check_fwd( "hpx pole phi=-45", AST__HPX, "", -45.0, -90.0,
              -45.0, -90.0, 1.0E-10 );
   check_rev( "hpx pole phi=-45 reverse", AST__HPX, "", -45.0, -90.0,
              -45.0, -90.0, 1.0E-10 );
   check_fwd( "hpx pole phi=0", AST__HPX, "", 0.0, -90.0,
              45.0, -90.0, 1.0E-10 );
   check_rev( "hpx pole phi=0 reverse", AST__HPX, "", 45.0, -90.0,
              45.0, -90.0, 1.0E-10 );
   check_fwd( "hpx pole phi=90", AST__HPX, "", 90.0, -90.0,
              135.0, -90.0, 1.0E-10 );
   check_rev( "hpx pole phi=90 reverse", AST__HPX, "", 135.0, -90.0,
              135.0, -90.0, 1.0E-10 );

/* ZPN with a degree-2 polynomial: beyond the turning point of r(theta)
   the projection folds back on itself, and the point is on the far side.
   wcslib 8.10 reports it as bad at every degree. */
   check_fwd( "zpn degree 2 far side", AST__ZPN,
              "PV2_1=1.27118,PV2_2=-0.210643", 5.4085552, -89.2063216,
              AST__BAD, AST__BAD, 0.0 );

/* PCO within 1e-4 degrees of the equator: cot(theta) is huge and
   1 - cos(phi*sin(theta)) cancels, so y must come from a series instead
   (upstream 4.22). */
   check_fwd( "pco near equator", AST__PCO, "", 90.0, 1.0E-9,
              90.0, 2.2337005501361701e-09, 1.0E-15 );
   check_fwd( "pco near equator, far from the meridian", AST__PCO, "",
              -170.0, -5.0E-5, -169.99999999993526, -0.0002700860857957293,
              1.0E-12 );

   if( !failed && astOK ) {
      printf( " All WcsMap tests passed\n" );
      return 0;
   }
   printf( "WcsMap tests failed\n" );
   return 1;
}
