/*
*  Name:
*     testsplinemapedge

*  Purpose:
*     Check the values of a SplineMap on the edges of the area it covers.

*  Description:
*     A SplineMap covers a rectangle, and positions on the edges of that
*     rectangle are inside it.  The spline value on the upper edge of each
*     axis is the limit approached from below, which is found by stepping
*     back from the upper limit to the last knot interval of non-zero
*     width.  This test checks that the step back works when there is only
*     one knot interval, so that the last interval is also the first, and
*     that it is abandoned, giving bad values, only if no interval has a
*     non-zero width.  It checks:
*
*       - a bilinear SplineMap with a single knot interval on each axis, at
*         its corners, on each of its edges and inside it;
*       - a SplineMap of order one, which is constant over its single knot
*         interval, at the same positions;
*       - a SplineMap with two knot intervals on one axis and one on the
*         other, in both arrangements;
*       - positions just outside the area covered, which are bad;
*       - a SplineMap whose knots are all equal on one axis, which covers
*         no area and so is bad everywhere.

*  Notes:
*     - This is a native C test with no Fortran original.
*     - Public C API only.
*/

#include <math.h>
#include <stdio.h>
#include "sae_par.h"
#include "ast.h"

/* The first and second outputs of the bilinear SplineMaps used here. */
static double uval( double x, double y ) {
   return 1.0 + 2.0*x + 3.0*y + 4.0*x*y;
}
static double vval( double x, double y ) {
   return 5.0 - x + 2.0*y - 3.0*x*y;
}

/* Transform a grid of positions covering the edges, corners and interior
   of the rectangle (0:xhi,0:yhi), plus positions just outside it, and
   compare with the wanted values.  If "bilinear" is non-zero the wanted
   values are those of "uval" and "vval", otherwise they are the constants
   "ucon" and "vcon".  If "allbad" is non-zero, every output must be bad.
   Returns the number of failures. */
static int check_map( const char *what, AstSplineMap *map, double xhi,
                      double yhi, int bilinear, double ucon, double vcon,
                      int allbad ) {
   double frac[ 7 ] = { -1.0e-6, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0 + 1.0e-6 };
   double xin[ 49 ];
   double yin[ 49 ];
   double uout[ 49 ];
   double vout[ 49 ];
   double uexp;
   double vexp;
   int i;
   int ix;
   int iy;
   int inside;
   int nfail = 0;

   if ( !astOK || !map ) {
      fprintf( stderr, "testsplinemapedge: could not build %s\n", what );
      if ( !astOK ) astClearStatus;
      return 1;
   }

   for ( iy = 0; iy < 7; iy++ ) {
      for ( ix = 0; ix < 7; ix++ ) {
         xin[ ix + 7*iy ] = frac[ ix ]*xhi;
         yin[ ix + 7*iy ] = frac[ iy ]*yhi;
      }
   }

   astTran2( map, 49, xin, yin, 1, uout, vout );
   if ( !astOK ) {
      fprintf( stderr, "testsplinemapedge: could not transform %s\n", what );
      astClearStatus;
      return 1;
   }

   for ( iy = 0; iy < 7; iy++ ) {
      for ( ix = 0; ix < 7; ix++ ) {
         i = ix + 7*iy;
         inside = ( ix > 0 && ix < 6 && iy > 0 && iy < 6 ) && !allbad;
         if ( !inside ) {
            uexp = AST__BAD;
            vexp = AST__BAD;
         } else if ( bilinear ) {
            uexp = uval( xin[ i ], yin[ i ] );
            vexp = vval( xin[ i ], yin[ i ] );
         } else {
            uexp = ucon;
            vexp = vcon;
         }
         if ( inside ? ( uout[ i ] == AST__BAD || vout[ i ] == AST__BAD ||
                         fabs( uout[ i ] - uexp ) > 1.0e-12 ||
                         fabs( vout[ i ] - vexp ) > 1.0e-12 )
                     : ( uout[ i ] != AST__BAD || vout[ i ] != AST__BAD ) ) {
            fprintf( stderr, "testsplinemapedge: %s at (%g,%g) gave (%g,%g), "
                     "expected (%g,%g)\n", what, xin[ i ], yin[ i ],
                     uout[ i ], vout[ i ], uexp, vexp );
            nfail++;
         }
      }
   }
   return nfail;
}

int main( void ) {
   int status = SAI__OK;
   int nfail = 0;

/* Knots for a single knot interval covering (0:1), for orders two and one,
   and for two knot intervals covering (0:2) for order two. */
   double t2[ 4 ] = { 0.0, 0.0, 1.0, 1.0 };
   double t1[ 2 ] = { 0.0, 1.0 };
   double t22[ 5 ] = { 0.0, 0.0, 1.0, 2.0, 2.0 };

/* Knots that are all equal, so that they cover no area. */
   double tzero[ 4 ] = { 1.0, 1.0, 1.0, 1.0 };

/* Coefficients of an order two spline are its values at the knots.  The
   coefficients are stored with the Y index varying fastest. */
   double cu[ 4 ];
   double cv[ 4 ];
   double cu32[ 6 ];
   double cv32[ 6 ];
   double cu23[ 6 ];
   double cv23[ 6 ];
   double c1u[ 1 ] = { 7.0 };
   double c1v[ 1 ] = { -3.0 };
   int ix;
   int iy;

   astWatch( &status );
   astBegin;

   for ( ix = 0; ix < 2; ix++ ) {
      for ( iy = 0; iy < 2; iy++ ) {
         cu[ iy + 2*ix ] = uval( ix, iy );
         cv[ iy + 2*ix ] = vval( ix, iy );
      }
   }
   for ( ix = 0; ix < 3; ix++ ) {
      for ( iy = 0; iy < 2; iy++ ) {
         cu32[ iy + 2*ix ] = uval( ix, iy );
         cv32[ iy + 2*ix ] = vval( ix, iy );
      }
   }
   for ( ix = 0; ix < 2; ix++ ) {
      for ( iy = 0; iy < 3; iy++ ) {
         cu23[ iy + 3*ix ] = uval( ix, iy );
         cv23[ iy + 3*ix ] = vval( ix, iy );
      }
   }

/* A single knot interval on each axis. */
   nfail += check_map( "a bilinear SplineMap with one interval per axis",
                       astSplineMap( 2, 2, 2, 2, t2, t2, cu, cv, " " ),
                       1.0, 1.0, 1, 0.0, 0.0, 0 );
   nfail += check_map( "a constant SplineMap with one interval per axis",
                       astSplineMap( 1, 1, 1, 1, t1, t1, c1u, c1v, " " ),
                       1.0, 1.0, 0, 7.0, -3.0, 0 );

/* Two knot intervals on one axis and one on the other.  A bilinear
   function is reproduced exactly by a piecewise bilinear spline. */
   nfail += check_map( "a SplineMap with two X intervals and one Y interval",
                       astSplineMap( 2, 2, 3, 2, t22, t2, cu32, cv32, " " ),
                       2.0, 1.0, 1, 0.0, 0.0, 0 );
   nfail += check_map( "a SplineMap with one X interval and two Y intervals",
                       astSplineMap( 2, 2, 2, 3, t2, t22, cu23, cv23, " " ),
                       1.0, 2.0, 1, 0.0, 0.0, 0 );

/* Knots that are all equal on one axis cover no area, so there is no
   interval from which to find a value, even at the knots themselves. */
   nfail += check_map( "a SplineMap with equal X knots",
                       astSplineMap( 2, 2, 2, 2, tzero, t2, cu, cv, " " ),
                       1.0, 1.0, 1, 0.0, 0.0, 1 );
   nfail += check_map( "a SplineMap with equal Y knots",
                       astSplineMap( 2, 2, 2, 2, t2, tzero, cu, cv, " " ),
                       1.0, 1.0, 1, 0.0, 0.0, 1 );

   astEnd;
   if ( !astOK ) {
      fprintf( stderr, "testsplinemapedge: error on exit\n" );
      astClearStatus;
      nfail++;
   }

   if ( nfail == 0 ) {
      printf( " All SplineMap edge tests passed\n" );
      return 0;
   }
   printf( "SplineMap edge tests failed (%d)\n", nfail );
   return 1;
}
