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
 *  - still sets to zero a shift that cancels to within its rounding error
 *    when two WinMaps are merged, including when the terms are large
 *    (x + 1.1 followed by 7e6*x - 7.7e6 is 7e6*x), so that the result has
 *    no shift.
 *
 * Every simplified Mapping must also transform as the original does.
 */
#include "ast.h"
#include <math.h>
#include <stdarg.h>
#include <stdio.h>

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

int main( void ) {
   astWatch( &status );
   astBegin;

   check_small_shift();
   check_joules();
   check_distinct_shifts();
   check_residue();

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
