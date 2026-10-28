/*
 * testspecmap: SpecMap regression tests (C only; no Fortran original).
 *
 * astSpecAdd must read only the arguments the caller supplies.  The
 * conversions to and from heliocentric frequency store one argument more
 * than the caller passes, so an implementation that copies the stored
 * length reads past the end of the caller's array; AddressSanitizer reports
 * that as a stack-buffer-overflow in the call below.
 */
#include "ast.h"
#include <math.h>
#include <stdio.h>

static int status = 0;
static int failed = 0;

/* Topocentric to heliocentric frequency, with the six documented TPF2HL
   arguments in an array of exactly that size. */
static void check_topocentric_frequency( void ) {
   AstSpecMap *map;
   double args[ 6 ] = { -2.71, 0.34, 4200.0, 60000.0, 1.570796326794896,
                        0.1 };
   double in = 1.0E9, out;
   const double want = 1000085270.7870382;

   astBegin;
   map = astSpecMap( 1, 0, " " );
   astSpecAdd( map, "TPF2HL", 6, args );
   astTran1( map, 1, &in, 1, &out );
   if( !astOK ) {
      astClearStatus;
      printf( "Error: TPF2HL transformation failed\n" );
      failed = 1;
   } else if( fabs( out - want ) > 2.0E-6 ) {
      printf( "Error: TPF2HL gives %.17g Hz, expected %.17g\n", out, want );
      failed = 1;
   }
   astEnd;
}

int main( void ) {
   astWatch( &status );

   check_topocentric_frequency();

   if( !failed && astOK ) {
      printf( " All SpecMap tests passed\n" );
      return 0;
   }
   printf( "SpecMap tests failed\n" );
   return 1;
}
