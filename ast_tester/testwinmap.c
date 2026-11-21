/*
 * testwinmap: WinMap simplification tests (C only; no Fortran original).
 *
 * A WinMap's terms carry the units of its axes, so a legitimate shift may be
 * far smaller than 1.0E-12: a spectral axis in joules has shifts of order
 * 1e-26.  These tests check that astSimplify:
 *
 *  - keeps a small non-zero shift, rather than replacing the WinMap with a
 *    diagonal MatrixMap that has none;
 *  - maps CRPIX to CRVAL for the pixel-to-energy chain of a FITS axis in
 *    joules, ShiftMap(-CRPIX) then ZoomMap(CDELT) then ShiftMap(CRVAL);
 *  - does not cancel a WinMap against the inverse of one whose shift
 *    differs, however small both shifts are; and that astEqual tells such
 *    WinMaps apart;
 *  - does not replace a WinMap or a diagonal MatrixMap whose small scale
 *    terms differ with a ZoomMap; and that astEqual tells apart ShiftMaps,
 *    ZoomMaps and MatrixMaps whose small terms differ, while still finding
 *    two matrices equal if they differ by rounding noise;
 *  - still sets to zero a shift that cancels to within its rounding error
 *    when two WinMaps are merged, including when the terms are large
 *    (x + 1.1 followed by 7e6*x - 7.7e6 is 7e6*x), so that the result has
 *    no shift;
 *  - does the same for a series of ShiftMaps in which the larger shifts
 *    cancel first, so that what the smaller ones leave is the rounding
 *    error of a larger value, and replaces the series with a UnitMap;
 *  - does the same for a shift that is moved past a full MatrixMap, so
 *    that the dump of the result does not depend on whether the compiler
 *    contracts a multiply and an add;
 *  - keeps the difference between two large shifts that are nearly equal
 *    and opposite, if it is larger than their rounding error.
 *
 * They also check that a WinMap made from the corners of two windows that
 * differ only in scale has no shift, and so simplifies to a ZoomMap.
 *
 * Every simplified Mapping must also transform as the original does.
 */
#include "ast.h"
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>

static int status = 0;
static int failed = 0;

static void fail( const char *fmt, ... )
   __attribute__((format(printf, 1, 2)));

static void fail( const char *fmt, ... ) {
   va_list ap;
   printf( "Error: " );
   va_start( ap, fmt );
   vprintf( fmt, ap );
   va_end( ap );
   printf( "\n" );
   failed = 1;
}

/* Check that "simp" transforms "n" 1-D inputs as "map" does, to a relative
   tolerance. */
static void check_same( AstMapping *map, AstMapping *simp, int n,
                        const double *in, const char *what ) {
   double out1[ 8 ], out2[ 8 ];
   int i;
   astTran1( map, n, in, 1, out1 );
   astTran1( simp, n, in, 1, out2 );
   for( i = 0; i < n; i++ ) {
      if( fabs( out2[ i ] - out1[ i ] ) > 1.0E-12*fabs( out1[ i ] ) ) {
         fail( "%s: simplified Mapping gives %.17g at %g, original %.17g",
               what, out2[ i ], in[ i ], out1[ i ] );
         return;
      }
   }
}

/* Check that "simp" transforms "n" 1-D inputs as "map" does, to an absolute
   tolerance. */
static void check_same_within( AstMapping *map, AstMapping *simp, int n,
                               const double *in, double tol,
                               const char *what ) {
   double out1[ 8 ], out2[ 8 ];
   int i;
   astTran1( map, n, in, 1, out1 );
   astTran1( simp, n, in, 1, out2 );
   for( i = 0; i < n; i++ ) {
      if( !( fabs( out2[ i ] - out1[ i ] ) <= tol ) ) {
         fail( "%s: simplified Mapping gives %.17g at %g, original %.17g",
               what, out2[ i ], in[ i ], out1[ i ] );
         return;
      }
   }
}

/* A WinMap for x -> shift + scale*x, given as windows (0,1) -> (shift,
   shift + scale), so that the shift is stored exactly. */
static AstWinMap *winmap( double shift, double scale ) {
   double ina = 0.0, inb = 1.0, outa = shift, outb = shift + scale;
   return astWinMap( 1, &ina, &inb, &outa, &outb, " " );
}

static void check_small_shift( void ) {
   AstMapping *map, *simp;
   double in[ 3 ] = { 0.0, 50.0, 100.0 };
   double out;

   astBegin;

/* x -> 1.5e-14 + 3e-16*x, a grism parameter for a ruling density of
   3e-6 m^-1 (the shift that places the reference pixel). */
   map = (AstMapping *) astCmpMap( astZoomMap( 1, 3.0E-16, " " ),
                                   astShiftMap( 1, (double[]){ 1.5E-14 },
                                                " " ), 1, " " );
   simp = astSimplify( map );
   astTran1( simp, 1, in, 1, &out );
   if( out != 1.5E-14 ) {
      fail( "a shift of 1.5e-14 simplifies to %.17g", out );
   }
   check_same( map, simp, 3, in, "small shift" );

   astEnd;
}

static void check_joules( void ) {
   AstMapping *map, *simp;
   double crpix = 512.5, cdelt = 9.425E-29, crval = 9.425E-25;
   double in[ 3 ] = { 1.0, 512.5, 1024.0 };
   double out;

   astBegin;
   map = (AstMapping *) astCmpMap(
            astCmpMap( astShiftMap( 1, (double[]){ -crpix }, " " ),
                       astZoomMap( 1, cdelt, " " ), 1, " " ),
            astShiftMap( 1, &crval, " " ), 1, " " );
   simp = astSimplify( map );
   astTran1( simp, 1, &crpix, 1, &out );
   if( fabs( out - crval ) > 1.0E-14*crval ) {
      fail( "energy axis in joules maps CRPIX to %.17g J, not CRVAL %.17g J",
            out, crval );
   }
   check_same( map, simp, 3, in, "energy axis in joules" );
   astEnd;
}

static void check_distinct_shifts( void ) {
   AstWinMap *w1, *w2;
   AstMapping *map, *simp;
   double in[ 3 ] = { 0.0, 1.0, -2.0 };
   double out;

   astBegin;
   w1 = winmap( 1.0E-14, 1.0 );
   w2 = winmap( 5.0E-14, 1.0 );

   if( astEqual( w1, w2 ) ) {
      fail( "WinMaps with shifts 1e-14 and 5e-14 are equal" );
   }

/* x + 1e-14 followed by the inverse of x + 5e-14 is x - 4e-14, not the
   identity. */
   astInvert( w2 );
   map = (AstMapping *) astCmpMap( w1, w2, 1, " " );
   simp = astSimplify( map );
   astTran1( simp, 1, in, 1, &out );
   if( fabs( out + 4.0E-14 ) > 1.0E-12*4.0E-14 ) {
      fail( "x + 1e-14 then the inverse of x + 5e-14 maps 0 to %.17g, not"
            " -4e-14", out );
   }
   check_same( map, simp, 3, in, "distinct small shifts" );
   astEnd;
}

static void check_residue( void ) {
   AstMapping *map, *simp;
   double in[ 3 ] = { 1.0, 2.0, -3.0 };
   double zero = 0.0, out;

   astBegin;

/* x + 1.1 followed by 7e6*x - 7.7e6 is 7e6*x.  Merging leaves a shift of
   about 1e-9, within the rounding error of the 7.7e6 terms that cancel. */
   map = (AstMapping *) astCmpMap( astShiftMap( 1, (double[]){ 1.1 }, " " ),
                                   winmap( -7.7E6, 7.0E6 ), 1, " " );
   simp = astSimplify( map );
   astTran1( simp, 1, &zero, 1, &out );
   if( out != 0.0 ) {
      fail( "x + 1.1 then 7e6*x - 7.7e6 maps 0 to %.17g, not 0", out );
   }
   if( astIsAWinMap( simp ) ) {
      fail( "x + 1.1 then 7e6*x - 7.7e6 simplifies to a WinMap with a"
            " shift" );
   }
   check_same( map, simp, 3, in, "cancelling shifts" );
   astEnd;
}

/* Shifts that cancel in a series of ShiftMaps.  The first two are added
   before either is canceled, and their sum is rounded.  Once the larger
   has been canceled, that rounding error is much larger than the rounding
   error of the smaller shift it is left in, so it can only be recognized
   by looking at every value that was added. */
static void check_identity_chain( void ) {
   static const double shifts[ 3 ][ 2 ] = { { 1023.5, 0.3 }, { -10.0, -0.1 },
                                            { 512.5, 1.0E-3 } };
   AstMapping *map, *part, *simp;
   AstShiftMap *sm[ 4 ];
   double in[ 3 ] = { 0.0, 1.0, -37.5 };
   double sh;
   int i, j, k;

   astBegin;
   for( i = 0; i < 3; i++ ) {

/* The larger shift, the smaller, and then either the inverse of each or
   a ShiftMap with the opposite shift. */
      for( j = 0; j < 2; j++ ) {
         for( k = 0; k < 4; k++ ) {
            sh = shifts[ i ][ k % 2 ];
            if( k > 1 && j ) sh = -sh;
            sm[ k ] = astShiftMap( 1, &sh, " " );
            if( k > 1 && !j ) astInvert( sm[ k ] );
         }
         map = (AstMapping *) sm[ 0 ];
         for( k = 1; k < 4; k++ ) {
            part = (AstMapping *) astCmpMap( map, sm[ k ], 1, " " );
            map = part;
         }
         simp = astSimplify( map );
         if( !astIsAUnitMap( simp ) ) {
            fail( "shifts of %g and %g followed by %s simplify to a %s, not"
                  " a UnitMap", shifts[ i ][ 0 ], shifts[ i ][ 1 ],
                  j ? "the opposite shifts" : "their inverses",
                  astGetC( simp, "Class" ) );
         }
         check_same_within( map, simp, 3, in, 4.0E-13,
                            "shifts that cancel" );
      }
   }
   astEnd;
}

/* A shift moved past a full MatrixMap is the sum of the old shifts, each
   multiplied by a matrix element.  Here the two products that make up the
   first new shift cancel, but only if both are rounded in the same way,
   which they are not if the compiler contracts a multiply and an add.  The
   ShiftMap in the result must have no shift on that axis either way. */
static void check_matrix_swap( void ) {
   AstMapping *map, *simp;
   AstMatrixMap *m1, *m2;
   AstWinMap *wm;
   char *text;
   double mat1[ 4 ] = { 0.3, -0.3, 1.0, 1.0 };
   double mat2[ 4 ] = { 1.5, 0.25, -0.5, 2.0 };
   double ina[ 2 ] = { 0.0, 0.0 }, inb[ 2 ] = { 1.0, 1.0 };
   double outa[ 2 ] = { 0.1, 0.1 }, outb[ 2 ] = { 1.1, 1.1 };
   double xin[ 3 ] = { 0.0, 1.0, -3.0 }, yin[ 3 ] = { 0.0, 2.0, 5.0 };
   double x1[ 3 ], y1[ 3 ], x2[ 3 ], y2[ 3 ];
   int i;

   astBegin;
   m1 = astMatrixMap( 2, 2, 0, mat1, "Invert=1" );
   wm = astWinMap( 2, ina, inb, outa, outb, " " );
   m2 = astMatrixMap( 2, 2, 0, mat2, " " );
   map = (AstMapping *) astCmpMap( astCmpMap( m1, wm, 1, " " ), m2, 1, " " );
   simp = astSimplify( map );

   text = astToString( simp );
   if( text && strstr( text, "Sft1" ) ) {
      fail( "a shift moved past a MatrixMap leaves a shift on an axis"
            " where its terms cancel" );
   }
   text = astFree( text );

   astTran2( map, 3, xin, yin, 1, x1, y1 );
   astTran2( simp, 3, xin, yin, 1, x2, y2 );
   for( i = 0; i < 3; i++ ) {
      if( fabs( x2[ i ] - x1[ i ] ) > 1.0E-13 ||
          fabs( y2[ i ] - y1[ i ] ) > 1.0E-13 ) {
         fail( "shift moved past a MatrixMap: simplified Mapping gives"
               " (%.17g,%.17g) at (%g,%g), original (%.17g,%.17g)",
               x2[ i ], y2[ i ], xin[ i ], yin[ i ], x1[ i ], y1[ i ] );
         break;
      }
   }
   astEnd;
}

/* Two shifts that are nearly equal and opposite leave a difference that
   is exact.  It is kept if it is four or more units in the last place of
   the shifts, as each of these is: two Julian Dates 161 microseconds
   apart, two Modified Julian Dates 2.5 microseconds apart, and two
   frequencies of 230 GHz that are 0.12 mHz apart. */
static void check_exact_difference( void ) {
   static const double pairs[ 3 ][ 2 ] = {
      { 2460000.5, 2460000.499999998 },
      { 60000.0, 59999.99999999997 },
      { 2.3E11, 2.3E11 - 1.2E-4 } };
   AstMapping *map, *simp;
   double zero = 0.0, out, sh, want;
   int i;

   astBegin;
   for( i = 0; i < 3; i++ ) {
      want = pairs[ i ][ 0 ] - pairs[ i ][ 1 ];
      sh = -pairs[ i ][ 1 ];
      map = (AstMapping *) astCmpMap( astShiftMap( 1, &pairs[ i ][ 0 ], " " ),
                                      astShiftMap( 1, &sh, " " ), 1, " " );
      simp = astSimplify( map );
      astTran1( simp, 1, &zero, 1, &out );
      if( out != want ) {
         fail( "shifts of %.17g and %.17g simplify to a %s that maps 0 to"
               " %.17g, not %.17g", pairs[ i ][ 0 ], sh,
               astGetC( simp, "Class" ), out, want );
      }
   }
   astEnd;
}

/* The corners of two windows that differ only in scale give a shift that
   is zero but for rounding, so the WinMap is a ZoomMap if the scale is the
   same on every axis.  A window that is shifted as well keeps its shift. */
static void check_corners( void ) {
   AstMapping *simp;
   AstWinMap *wm;
   double ina[ 2 ] = { 0.1, 0.7 }, inb[ 2 ] = { 0.3, 1.9 };
   double outa[ 2 ] = { 0.3, 2.1 }, outb[ 2 ] = { 0.9, 5.7 };
   double a, b, oa, ob, scale, zero = 0.0, out;
   int i, j, nbad = 0;

   astBegin;
   wm = astWinMap( 2, ina, inb, outa, outb, " " );
   simp = astSimplify( wm );
   if( !astIsAZoomMap( simp ) ) {
      fail( "a WinMap that scales both axes by 3 simplifies to a %s, not a"
            " ZoomMap", astGetC( simp, "Class" ) );
   }

   for( i = 1; i < 40; i++ ) {
      for( j = 1; j < 25; j++ ) {
         scale = 1.0 + 0.1*j;
         a = 0.01*i;
         b = a + 0.37*j;
         oa = scale*a;
         ob = scale*b;
         wm = astWinMap( 1, &a, &b, &oa, &ob, " " );
         simp = astSimplify( wm );
         if( !astIsAZoomMap( simp ) ) nbad++;
         simp = astAnnul( simp );
         wm = astAnnul( wm );
      }
   }
   if( nbad ) {
      fail( "%d of 936 WinMaps that only scale an axis do not simplify to a"
            " ZoomMap", nbad );
   }

/* x -> 3*x + 1e-13, from the windows (0.1,0.3) and (0.3 + 1e-13, 0.9 +
   1e-13), has a shift that is far larger than the rounding error of the
   corners, so it is kept. */
   a = 0.1;
   b = 0.3;
   oa = 0.3 + 1.0E-13;
   ob = 0.9 + 1.0E-13;
   wm = astWinMap( 1, &a, &b, &oa, &ob, " " );
   simp = astSimplify( wm );
   astTran1( simp, 1, &zero, 1, &out );
   if( astIsAZoomMap( simp ) || fabs( out - 1.0E-13 ) > 1.0E-15 ) {
      fail( "a window shifted by 1e-13 gives a %s that maps 0 to %.17g",
            astGetC( simp, "Class" ), out );
   }
   astEnd;
}

/* Terms that carry the units of an axis are compared using a relative
   tolerance alone by each of the classes that a WinMap is simplified to. */
static void check_small_terms( void ) {
   AstMapping *simp;
   AstWinMap *wm;
   double scl[ 2 ] = { 1.0E-26, 3.0E-26 };
   double ina[ 2 ] = { 0.0, 0.0 }, inb[ 2 ] = { 1.0, 1.0 };
   double xin = 50.0, yin = 50.0, xout, yout;
   double sh1 = 1.5E-14, sh2 = 2.5E-14;
   double full1[ 4 ] = { 1.0, 1.0E-7, 0.0, 1.0 };
   double full2[ 4 ] = { 1.0, 1.000005E-7, 0.0, 1.0 };
   double tiny1[ 4 ] = { 1.0E-26, 1.0E-33, 0.0, 1.0E-26 };
   double tiny2[ 4 ] = { 1.0E-26, 1.000005E-33, 0.0, 1.0E-26 };
   double tiny3[ 4 ] = { 1.0E-26, 3.0E-33, 0.0, 1.0E-26 };
   int pass;

   astBegin;

/* x -> 1e-26*x, y -> 3e-26*y is not a ZoomMap, whether it is given as a
   WinMap or as a diagonal MatrixMap. */
   for( pass = 0; pass < 2; pass++ ) {
      if( pass ) {
         simp = astSimplify( astMatrixMap( 2, 2, 1, scl, " " ) );
      } else {
         wm = astWinMap( 2, ina, inb, ina, scl, " " );
         simp = astSimplify( wm );
      }
      astTran2( simp, 1, &xin, &yin, 1, &xout, &yout );
      if( astIsAZoomMap( simp ) || fabs( xout - 5.0E-25 ) > 1.0E-37 ||
          fabs( yout - 1.5E-24 ) > 1.0E-36 ) {
         fail( "scales of 1e-26 and 3e-26 in a %s simplify to a %s that"
               " maps (50,50) to (%.17g,%.17g)", pass ? "MatrixMap" :
               "WinMap", astGetC( simp, "Class" ), xout, yout );
      }
   }

   if( astEqual( astZoomMap( 1, scl[ 0 ], " " ),
                 astZoomMap( 1, scl[ 1 ], " " ) ) ) {
      fail( "ZoomMaps with factors 1e-26 and 3e-26 are equal" );
   }
   if( astEqual( astShiftMap( 1, &sh1, " " ), astShiftMap( 1, &sh2, " " ) ) ) {
      fail( "ShiftMaps with shifts 1.5e-14 and 2.5e-14 are equal" );
   }
   if( astEqual( astMatrixMap( 1, 1, 1, &scl[ 0 ], " " ),
                 astMatrixMap( 1, 1, 1, &scl[ 1 ], " " ) ) ) {
      fail( "MatrixMaps with elements 1e-26 and 3e-26 are equal" );
   }

/* Two matrices are still equal if they differ by an amount that is
   negligible beside the other elements in the same row and column,
   whatever the size of the elements. */
   if( !astEqual( astMatrixMap( 2, 2, 0, full1, " " ),
                  astMatrixMap( 2, 2, 0, full2, " " ) ) ) {
      fail( "matrices of order one that differ by 5e-13 are not equal" );
   }
   if( !astEqual( astMatrixMap( 2, 2, 0, tiny1, " " ),
                  astMatrixMap( 2, 2, 0, tiny2, " " ) ) ) {
      fail( "matrices of order 1e-26 that differ by 5e-39 are not equal" );
   }
   if( astEqual( astMatrixMap( 2, 2, 0, tiny1, " " ),
                 astMatrixMap( 2, 2, 0, tiny3, " " ) ) ) {
      fail( "matrices of order 1e-26 that differ by 2e-33 are equal" );
   }
   astEnd;
}

int main( void ) {
   astWatch( &status );
   astBegin;

   check_small_shift();
   check_joules();
   check_distinct_shifts();
   check_residue();
   check_identity_chain();
   check_matrix_swap();
   check_exact_difference();
   check_corners();
   check_small_terms();

   if( !astOK ) {
      fail( "AST error" );
      astClearStatus;
   }
   astEnd;

   if( !failed ) {
      printf( " All WinMap tests passed\n" );
      return 0;
   }
   printf( "WinMap tests failed\n" );
   return 1;
}
