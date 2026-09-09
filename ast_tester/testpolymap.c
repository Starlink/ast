/*
 *  Test the PolyMap class (PolyTran, PolyCoeffs, IterInverse).
 *  Converted from the Fortran test testpolymap.f.
 *  Direct conversion; no material differences from the Fortran original.
 */
#include "ast.h"
#include <stdio.h>
#include <math.h>

static void stopit( int i, int *status ) {
   if( *status != 0 ) return;
   printf( "Error %d\n", i );
   *status = 1;
}

int main( void ) {
   int status_value = 0;
   int *status = &status_value;
   int i, nco;
   double lbnd[] = { -1000.0, -1000.0 };
   double ubnd[] = { 1000.0, 1000.0 };
   double acc = 1.0e-7, errlim, maxacc = 1.0e-3;
   int maxord = 10;
   double xin[3], yin[3], xout[3], yout[3], xin2[3], yin2[3];
   double cofs[20];
   AstPolyMap *pm, *pm2;
   AstMapping *simplified;

   double coeff[] = { 1.0, 1, 0, 0,
                      2.0, 1, 1, 0,
                      1.0, 2, 0, 0,
                      3.0, 2, 0, 1 };

   double coeffb[] = { 1.0, 1, 0, 0,
                       0.5, 1, 0, 0,
                       2.0, 1, 1, 0,
                       1.0, 2, 0, 0,
                       3.0, 2, 0, 1,
                      -0.5, 2, 0, 1 };

   double coeff2[] = { 1.0, 1, 0, 0,
                       2.0, 1, 1, 0,
                       1.0, 1, 0, 1,
                       1.0, 2, 0, 0,
                       1.0, 2, 1, 0,
                       2.0, 2, 0, 1 };

   double coeff3[] = { -0.1,   1, 0, 0,
                        0.99,  1, 1, 0,
                        1.0e-4,1, 1, 1,
                       -0.1,   2, 0, 0,
                        0.99,  2, 0, 1,
                        1.0e-4,2, 1, 1 };

   double coeff_1d[] = { 1.0, 1, 0,
                         2.0, 1, 1 };

   double coeff_2in1out[] = { 2.0, 1, 2, 0,
                              3.0, 1, 0, 1 };

   errlim = 1000 * acc;

   astWatch( status );
   astBegin;

   /* Basic 2D PolyMap and PolyCoeffs. */
   pm = astPolyMap( 2, 2, 4, coeff, 0, coeff, " " );

   astPolyCoeffs( pm, 1, 0, NULL, &nco );
   if( nco != 4 ) stopit( -1, status );

   astPolyCoeffs( pm, 0, 0, NULL, &nco );
   if( nco != 0 ) stopit( -2, status );

   astPolyCoeffs( pm, 1, 20, cofs, &nco );
   if( nco != 4 ) stopit( -3, status );
   for( i = 0; i < 16; i++ ) {
      if( cofs[i] != coeff[i] ) stopit( -4, status );
   }

   /* PolyTran: fit inverse and test round-trip. */
   pm2 = astPolyTran( pm, 0, acc, maxacc, maxord, lbnd, ubnd );
   xin[0] = 1.0;  xin[1] = 100.0;  xin[2] = -50.0;
   yin[0] = 1.0;  yin[1] = 100.0;  yin[2] = -50.0;

   astTran2( pm2, 3, xin, yin, 1, xout, yout );
   astTran2( pm2, 3, xout, yout, 0, xin2, yin2 );
   for( i = 0; i < 3; i++ ) {
      if( fabs( xin[i] - xin2[i] ) > errlim ) stopit( 1, status );
      if( fabs( yin[i] - yin2[i] ) > errlim ) stopit( 2, status );
   }

   /* IterInverse round-trip. */
   astSetL( pm2, "IterInverse", 1 );
   astTran2( pm2, 3, xout, yout, 0, xin2, yin2 );
   for( i = 0; i < 3; i++ ) {
      if( fabs( xin[i] - xin2[i] ) > errlim ) stopit( 1001, status );
      if( fabs( yin[i] - yin2[i] ) > errlim ) stopit( 1002, status );
   }

   /* Linear PolyMap should simplify to WinMap. */
   simplified = astSimplify( pm );
   if( !astIsAWinMap( simplified ) ) stopit( 1003, status );
   xin[0] = 1.0; yin[0] = 2.0;
   astTran2( simplified, 1, xin, yin, 1, xout, yout );
   if( xout[0] != 3.0 || yout[0] != 7.0 ) stopit( 1004, status );

   /* PolyMap with extra terms should also simplify to WinMap. */
   pm = astPolyMap( 2, 2, 6, coeffb, 0, coeffb, " " );
   simplified = astSimplify( pm );
   if( !astIsAWinMap( simplified ) ) stopit( 1005, status );
   xin[0] = 1.0; yin[0] = 2.0;
   astTran2( simplified, 1, xin, yin, 1, xout, yout );
   if( xout[0] != 3.5 || yout[0] != 6.0 ) stopit( 1006, status );

   /* 1D PolyMap. */
   pm = astPolyMap( 1, 1, 2, coeff_1d, 0, coeff_1d, " " );
   pm2 = astPolyTran( pm, 0, acc, maxacc, maxord, lbnd, ubnd );
   xin[0] = 1.0; xin[1] = 100.0; xin[2] = -50.0;
   astTran1( pm2, 3, xin, 1, xout );
   astTran1( pm2, 3, xout, 0, xin2 );
   for( i = 0; i < 3; i++ ) {
      if( fabs( xin[i] - xin2[i] ) > errlim ) stopit( 3, status );
   }

   astSetL( pm2, "IterInverse", 1 );
   astTran1( pm2, 3, xout, 0, xin2 );
   for( i = 0; i < 3; i++ ) {
      if( fabs( xin[i] - xin2[i] ) > errlim ) stopit( 3001, status );
   }

   /* Non-linear 2D: simplify should remain PolyMap. */
   pm = astPolyMap( 2, 2, 6, coeff2, 0, coeff2, " " );
   simplified = astSimplify( pm );
   if( !astIsAPolyMap( simplified ) ) stopit( 3002, status );
   if( !astEqual( simplified, pm ) ) stopit( 3003, status );

   pm2 = astPolyTran( pm, 0, acc, maxacc, maxord, lbnd, ubnd );
   xin[0] = 1.0;  xin[1] = 100.0;  xin[2] = -50.0;
   yin[0] = 1.0;  yin[1] = 100.0;  yin[2] = -50.0;
   astTran2( pm2, 3, xin, yin, 1, xout, yout );
   astTran2( pm2, 3, xout, yout, 0, xin2, yin2 );
   for( i = 0; i < 3; i++ ) {
      if( fabs( xin[i] - xin2[i] ) > errlim ) stopit( 4, status );
      if( fabs( yin[i] - yin2[i] ) > errlim ) stopit( 5, status );
   }

   astSetL( pm2, "IterInverse", 1 );
   astTran2( pm2, 3, xout, yout, 0, xin2, yin2 );
   for( i = 0; i < 3; i++ ) {
      if( fabs( xin[i] - xin2[i] ) > errlim ) stopit( 4001, status );
      if( fabs( yin[i] - yin2[i] ) > errlim ) stopit( 5001, status );
   }

   /* Another 2D non-linear PolyMap. */
   pm = astPolyMap( 2, 2, 6, coeff3, 0, coeff3, " " );
   pm2 = astPolyTran( pm, 0, acc, maxacc, maxord, lbnd, ubnd );
   xin[0] = 1.0;  xin[1] = 100.0;  xin[2] = -50.0;
   yin[0] = 1.0;  yin[1] = 100.0;  yin[2] = -50.0;
   astTran2( pm2, 3, xin, yin, 1, xout, yout );
   astTran2( pm2, 3, xout, yout, 0, xin2, yin2 );
   for( i = 0; i < 3; i++ ) {
      if( fabs( xin[i] - xin2[i] ) > errlim ) stopit( 6, status );
      if( fabs( yin[i] - yin2[i] ) > errlim ) stopit( 7, status );
   }

   astSetL( pm2, "IterInverse", 1 );
   astTran2( pm2, 3, xout, yout, 0, xin2, yin2 );
   for( i = 0; i < 3; i++ ) {
      if( fabs( xin[i] - xin2[i] ) > errlim ) stopit( 6001, status );
      if( fabs( yin[i] - yin2[i] ) > errlim ) stopit( 7001, status );
   }

   /* IterInverse attribute behaviour. */
   if( !astGetL( pm, "TranForward" ) ) stopit( 8001, status );
   if( !astGetL( pm, "IterInverse" ) ) stopit( 8002, status );
   if( !astGetL( pm, "TranInverse" ) ) stopit( 8003, status );

   astSetL( pm, "IterInverse", 0 );
   if( !astGetL( pm, "TranForward" ) ) stopit( 8004, status );
   if( astGetL( pm, "IterInverse" ) ) stopit( 8005, status );
   if( astGetL( pm, "TranInverse" ) ) stopit( 8006, status );

   astInvert( pm );
   if( astGetL( pm, "TranForward" ) ) stopit( 8007, status );
   if( astGetL( pm, "IterInverse" ) ) stopit( 8008, status );
   if( !astGetL( pm, "TranInverse" ) ) stopit( 8009, status );

   astSetL( pm, "IterInverse", 1 );
   if( !astGetL( pm, "TranForward" ) ) stopit( 8010, status );
   if( !astGetL( pm, "IterInverse" ) ) stopit( 8011, status );
   if( !astGetL( pm, "TranInverse" ) ) stopit( 8012, status );

   astInvert( pm );
   if( !astGetL( pm, "TranForward" ) ) stopit( 8013, status );
   if( !astGetL( pm, "IterInverse" ) ) stopit( 8014, status );
   if( !astGetL( pm, "TranInverse" ) ) stopit( 8015, status );

   /* A PolyMap with unequal Nin and Nout cannot use an iterative inverse,
      so IterInverse must default to zero for a forward-only 2-in 1-out
      PolyMap, and the map must not claim to define an inverse
      transformation. */
   pm2 = astPolyMap( 2, 1, 2, coeff_2in1out, 0, NULL, " " );
   if( astGetL( pm2, "IterInverse" ) ) stopit( 8016, status );
   if( astGetL( pm2, "TranInverse" ) ) stopit( 8017, status );
   if( !astGetL( pm2, "TranForward" ) ) stopit( 8018, status );

   /* astEqual must compare the inverse transformation's own coefficients, not
      the forward ones a second time, and must index every coefficient array by
      the axis count it was allocated with.

      The forward transformation here is y = x^2, chosen because it is not
      linear: a linear forward transformation is rebuilt as a MatrixMap and a
      ShiftMap by MapMerge, which discards the explicit inverse and would
      cancel the pair whatever astEqual said. */
   {
      double cf_sq[]  = { 1.0, 1, 2 };        /* y = x^2               */
      double ci_half[] = { 0.5, 1, 1 };       /* x = 0.5 y             */
      double ci_seven[] = { 0.7, 1, 1 };      /* x = 0.7 y             */
      AstPolyMap *pa, *pb, *pbi;
      AstCmpMap *series;
      AstMapping *simp;

      pa = astPolyMap( 1, 1, 1, cf_sq, 1, ci_half, " " );
      pb = astPolyMap( 1, 1, 1, cf_sq, 1, ci_half, " " );
      if( !astEqual( pa, pb ) ) stopit( 8019, status );

      pb = astPolyMap( 1, 1, 1, cf_sq, 1, ci_seven, " " );
      if( astEqual( pa, pb ) ) stopit( 8020, status );

      /* The consequence in MapMerge: a PolyMap and a neighbour used in the
         opposite direction are replaced by a UnitMap when astEqual says they
         match, so a pair whose inverses differ must not cancel. */
      pbi = astCopy( pb );
      astInvert( pbi );
      series = astCmpMap( pa, pbi, 1, " " );
      simp = astSimplify( series );
      if( astIsAUnitMap( simp ) ) stopit( 8021, status );

      /* The same pair with equal inverses must still cancel. */
      pb = astPolyMap( 1, 1, 1, cf_sq, 1, ci_half, " " );
      pbi = astCopy( pb );
      astInvert( pbi );
      series = astCmpMap( pa, pbi, 1, " " );
      simp = astSimplify( series );
      if( !astIsAUnitMap( simp ) ) stopit( 8022, status );
   }

   /* A PolyMap with more outputs than inputs exercises the array bounds: the
      inverse arrays have one element per input and "mxpow_f" one per input
      too, while the forward arrays have one per output. Comparing such a
      PolyMap with a copy of itself used to read past the end of both, which
      is a heap overflow and could report two identical PolyMaps as unequal. */
   {
      double cf_wide[ 24 ];
      double ci_wide[ 10 ];
      AstPolyMap *pw, *pwcopy;
      int io, nwide = 8;

      for( io = 0; io < nwide; io++ ) {
         cf_wide[ 3*io + 0 ] = io + 1.0;    /* coefficient           */
         cf_wide[ 3*io + 1 ] = io + 1.0;    /* one-based output index */
         cf_wide[ 3*io + 2 ] = 1.0;         /* power of the input    */
      }
      for( io = 0; io < 10; io++ ) ci_wide[ io ] = 0.0;
      ci_wide[ 0 ] = 1.0;                   /* coefficient           */
      ci_wide[ 1 ] = 1.0;                   /* one-based input index */
      ci_wide[ 2 ] = 1.0;                   /* power of output 1     */

      pw = astPolyMap( 1, nwide, nwide, cf_wide, 1, ci_wide, " " );
      pwcopy = astCopy( pw );
      if( !astEqual( pw, pwcopy ) ) stopit( 8023, status );
      if( !astEqual( pwcopy, pw ) ) stopit( 8024, status );
   }

   astEnd;
   astFlushMemory( 1 );

   if( *status == 0 ) {
      printf( " All PolyMap tests passed\n" );
   } else {
      printf( "PolyMap tests failed\n" );
   }
   return *status;
}
