/*
 * compare_dumps_selftest.c
 *
 * Tests the tolerance compare_dumps applies, rather than assuming it.  A
 * comparator nobody checks is how a gate quietly stops being one.
 *
 * compare_dumps.c is included rather than linked so that its helpers can stay
 * static.  COMPARE_DUMPS_NO_MAIN leaves out its main, so that each program's
 * main is in its own source file.
 */
#define COMPARE_DUMPS_NO_MAIN
#include "compare_dumps.c"

static int check( int expect, int got, const char *what ) {
   if( expect != got ) {
      printf( "compare_dumps selftest: %s -- expected %d, got %d\n",
              what, expect, got );
      return 1;
   }
   return 0;
}

int main( void ) {
   double up = 1.0;
   int fails = 0;
   int i;

   fails += check( 1, within_ulps( 1.0, 1.0 ), "exact equality matches" );
   fails += check( 1, within_ulps( 1.0, nextafter( 1.0, 2.0 ) ),
                   "one ULP matches" );

   for( i = 0; i < 5; i++ ) up = nextafter( up, 2.0 );
   fails += check( 0, within_ulps( 1.0, up ), "five ULP does not match" );

   fails += check( 0, within_ulps( nextafter( 0.0, 1.0 ), nextafter( 0.0, -1.0 ) ),
                   "opposite signs never match" );
   fails += check( 0, within_ulps( 1.0, 1.0000001 ),
                   "a real change does not match" );

   /* Token comparison: a non-numeric token is compared exactly, and a differing
      token count is a difference. */
   fails += check( 1, tokens_match( "\"FK5\"", "\"FK5\"" ), "string token matches" );
   fails += check( 0, tokens_match( "\"FK5\"", "\"ICRS\"" ),
                   "differing string tokens do not match" );
   fails += check( 0, tokens_match( "1.0x", "1.0" ),
                   "a partly numeric token is not a number" );
   fails += check( 0, lines_match( "Zoom = 2", "Zoom = 2 2" ),
                   "a dropped value is a difference" );
   fails += check( 1, lines_match( "   Zoom = 1", "   Zoom = 1" ),
                   "identical lines match" );

   /* astequal widens the tolerance, but only when asked. */
   use_astequal = 1;
   /* astEQUAL is |a-b| <= max( 1e5*(|a|+|b|)*DBL_EPSILON, 1e-12 ), so at 1e6
      it reaches about 4.4e-5 and at 1.0 about 4.4e-11. */
   fails += check( 1, within_ulps( 1.0e6, 1.0e6 + 1.0e-5 ),
                   "astequal absorbs a cancellation-scale difference" );
   fails += check( 0, within_ulps( 1.0e6, 1.0e6 + 1.0 ),
                   "astequal still rejects a real change" );
   fails += check( 1, within_ulps( 1.0, 1.0 + 1.0e-11 ),
                   "astequal reaches 1e-11 at unit magnitude" );
   fails += check( 0, within_ulps( 1.0, 1.0 + 1.0e-6 ),
                   "astequal stops well below 1e-6" );
   use_astequal = 0;
   fails += check( 0, within_ulps( 1.0e6, 1.0e6 + 1.0e-5 ),
                   "without astequal that difference is rejected" );

   if( fails ) {
      printf( "compare_dumps selftest: %d failures\n", fails );
      return 1;
   }
   printf( "compare_dumps selftest: all checks passed\n" );
   return 0;
}
