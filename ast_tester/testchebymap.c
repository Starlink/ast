/*
 *  Test the ChebyMap class.
 *  Converted from the Fortran test testchebymap.f.
 *
 *  The checkdump round-trip uses Channel SinkFile/SourceFile instead of
 *  Fortran channel source/sink callbacks (functionally equivalent).
 */
#include "ast.h"
#include <stdio.h>
#include <math.h>
#include "ast_err.h"

static void inverse_tests( int *status );

static void stopit( int i, int *status ) {
   if( *status != 0 ) return;
   printf( "Error %d\n", i );
   *status = 1;
}

static void checkdump( AstObject *obj, int *status ) {
   AstChannel *ch;
   AstObject *result;
   if( *status != 0 ) return;
   ch = astChannel( NULL, NULL, " " );
   astSet( ch, "SinkFile=fred.tmp" );
   if( astWrite( ch, obj ) != 1 ) { stopit( -1, status ); return; }
   astClear( ch, "SinkFile" );
   astSet( ch, "SourceFile=fred.tmp" );
   result = astRead( ch );
   if( !result ) { stopit( -2, status ); return; }
   astClear( ch, "SourceFile" );
   if( !astEqual( result, obj ) ) stopit( -3, status );
   remove( "fred.tmp" );
}

int main( void ) {
   int status_value = 0;
   int *status = &status_value;
   int i, j, nco;
   double lbnd[2], ubnd[2], dval;
   double dlbnd[2], dubnd[2], tlbnd[2], tubnd[2];
   double xin[5], xout[5], xrec[5], yrec[5], work[5];
   double yin[5], yout[5], xi, yi, xv, yv;
   double cofs[100];
   AstChebyMap *cm, *cm2, *cm3;

   /*  f(x) = 1.5*T0(x') - 1.0*T2(x') + 2.0*T3(x') + 1.3*T4(x') */
   double coeffs_1[] = { 1.5, 1, 0,
                         -1.0, 1, 2,
                          2.0, 1, 3,
                          1.3, 1, 4 };

   /* f(x) = 1.0*T0(x') - 2.0*T1(x') */
   double coeffs_2[] = { 1.0, 1, 0,
                         -2.0, 1, 1 };

   /* fx(x,y) = 1.0*T0(x')T0(y') - 2.0*T1(x')T2(y') + T1(y')
      fy(x,y) = 1.5*T0(x')T0(y') - 2.5*T1(x')T2(y') */
   double coeffs_3[] = { 1.0, 1, 0, 0,
                         -2.0, 1, 1, 2,
                          1.0, 1, 0, 1,
                          1.5, 2, 0, 0,
                         -2.5, 2, 1, 2 };

   /* fx(x,y) = T1(x') + T1(y')
      fy(x,y) = T1(x') - T1(y') */
   double coeffs_4[] = { 1.0, 1, 1, 0,
                         1.0, 1, 0, 1,
                         1.0, 2, 1, 0,
                        -1.0, 2, 0, 1 };

   astWatch( status );
   astBegin;

   lbnd[0] = -1.0; lbnd[1] = -1.0;
   ubnd[0] = 1.0;  ubnd[1] = 1.0;

   xin[0] = -1.0; xin[1] = -0.5; xin[2] = 0.0; xin[3] = 0.5; xin[4] = 1.0;

   /* 1D order 1: constant = 1.5 */
   cm = astChebyMap( 1, 1, 1, coeffs_1, 0, NULL, lbnd, ubnd, lbnd, ubnd, " " );
   astTran1( cm, 5, xin, 1, xout );
   for( i = 0; i < 5; i++ )
      if( xout[i] != 1.5 ) stopit( 0, status );

   /* 1D order 3: 2.5 - 2*x^2 */
   cm = astChebyMap( 1, 1, 2, coeffs_1, 0, NULL, lbnd, ubnd, lbnd, ubnd, " " );
   astTran1( cm, 5, xin, 1, xout );
   for( i = 0; i < 5; i++ ) {
      dval = 2.5 - 2.0*xin[i]*xin[i];
      if( xout[i] != dval ) stopit( 1, status );
   }

   /* 1D order 4: 2.5 - 6x - 2x^2 + 8x^3 */
   cm = astChebyMap( 1, 1, 3, coeffs_1, 0, NULL, lbnd, ubnd, lbnd, ubnd, " " );
   astTran1( cm, 5, xin, 1, xout );
   for( i = 0; i < 5; i++ ) {
      dval = 2.5 - 6.0*xin[i] - 2.0*xin[i]*xin[i] + 8.0*xin[i]*xin[i]*xin[i];
      if( xout[i] != dval ) stopit( 2, status );
   }

   /* 1D order 5: full Chebyshev evaluation */
   cm = astChebyMap( 1, 1, 4, coeffs_1, 0, NULL, lbnd, ubnd, lbnd, ubnd, " " );
   astTran1( cm, 5, xin, 1, xout );
   for( i = 0; i < 5; i++ ) {
      work[0] = 1.0;
      work[1] = xin[i];
      for( j = 2; j < 5; j++ )
         work[j] = 2.0*xin[i]*work[j-1] - work[j-2];
      dval = 1.5*work[0] - 1.0*work[2] + 2.0*work[3] + 1.3*work[4];
      if( xout[i] != dval ) stopit( 3, status );
   }

   if( !astGetL( cm, "IterInverse" ) ) stopit( 4, status );

   /* PolyTran on 1D ChebyMap, order 2: 1 - 2*x */
   lbnd[0] = -1.0; ubnd[0] = 1.0;
   cm = astChebyMap( 1, 1, 2, coeffs_2, 0, NULL, lbnd, ubnd, lbnd, ubnd, " " );
   cm2 = astPolyTran( cm, 0, 0.01, 0.01, 5, lbnd, ubnd );
   if( !cm2 ) {
      stopit( 5, status );
   } else {
      xin[0] = -1.0; xin[1] = -0.5; xin[2] = 0.0; xin[3] = 0.5; xin[4] = 1.0;
      astTran1( cm2, 5, xin, 1, xout );
      astTran1( cm2, 5, xout, 0, xrec );
      for( i = 0; i < 5; i++ )
         if( fabs( xrec[i] - xin[i] ) > 1.0e-3*fabs( xin[i] ) )
            stopit( 6, status );

      astChebyDomain( cm2, 0, dlbnd, dubnd );
      if( dlbnd[0] != -1.0 ) stopit( 501, status );
      if( dubnd[0] != 3.0 ) stopit( 502, status );

      astPolyCoeffs( cm2, 0, 100, cofs, &nco );
      if( nco != 1 ) stopit( 503, status );
      if( fabs( cofs[0] + 1.0 ) > 1.0e-10 ) stopit( 504, status );
      if( fabs( cofs[1] - 1.0 ) > 1.0e-10 ) stopit( 505, status );
      if( fabs( cofs[2] - 1.0 ) > 1.0e-10 ) stopit( 506, status );
   }

   /* PolyTran on 1D ChebyMap, order 5 with wider domain */
   lbnd[0] = -100.0; ubnd[0] = 100.0;
   cm = astChebyMap( 1, 1, 4, coeffs_1, 0, NULL, lbnd, ubnd, lbnd, ubnd, " " );
   { double fit_lbnd = -5.0, fit_ubnd = 50.0;
   cm2 = astPolyTran( cm, 0, 0.01, 0.01, 10, &fit_lbnd, &fit_ubnd ); }
   if( !cm2 ) {
      stopit( 7, status );
   } else {
      xin[0] = 0.0; xin[1] = 10.0; xin[2] = 20.0; xin[3] = 30.0; xin[4] = 40.0;
      astTran1( cm2, 5, xin, 1, xout );
      astTran1( cm2, 5, xout, 0, xrec );
      for( i = 0; i < 5; i++ )
         if( fabs( xrec[i] - xin[i] ) > 0.01 ) stopit( 8, status );
   }

   /* astEqual and astCopy */
   cm3 = astCopy( cm2 );
   if( !astEqual( cm2, cm3 ) ) stopit( 9, status );

   checkdump( (AstObject *)cm2, status );

   /* Simple 2D ChebyMap: fx = T1(x')+T1(y'), fy = T1(x')-T1(y') */
   lbnd[0] = -1.0; lbnd[1] = -1.0;
   ubnd[0] = 1.0;  ubnd[1] = 1.0;
   cm = astChebyMap( 2, 2, 4, coeffs_4, 0, NULL, lbnd, ubnd, lbnd, ubnd, " " );

   cm2 = astCopy( cm );
   astInvert( cm2 );
   cm3 = (AstChebyMap *)astSimplify( astCmpMap( cm, cm2, 1, " " ) );
   if( !astIsAUnitMap( cm3 ) ) stopit( 1000, status );

   xin[0] = 0.5; xin[1] = 0.0; xin[2] = -0.5; xin[3] = 0.0;
   yin[0] = 0.0; yin[1] = 0.5; yin[2] = 0.0;  yin[3] = -0.5;
   astTran2( cm, 4, xin, yin, 1, xout, yout );
   for( i = 0; i < 4; i++ ) {
      xv = xin[i] + yin[i];
      yv = xin[i] - yin[i];
      if( fabs( xout[i] - xv ) > 1.0e-6*fabs(xv) ||
          fabs( yout[i] - yv ) > 1.0e-6*fabs(yv) )
         stopit( 101, status );
   }

   cm2 = astPolyTran( cm, 0, 0.01, 0.01, 10, lbnd, ubnd );
   if( !cm2 ) {
      stopit( 102, status );
   } else {
      astTran2( cm2, 4, xout, yout, 0, xrec, yrec );
      for( i = 0; i < 4; i++ )
         if( fabs( xrec[i] - xin[i] ) > 0.01 ||
             fabs( yrec[i] - yin[i] ) > 0.01 )
            stopit( 103, status );
   }

   astPolyCoeffs( cm2, 0, 100, cofs, &nco );
   if( nco != 4 ) {
      stopit( 104, status );
   } else {
      for( i = 0; i < 16; i++ )
         if( fabs( cofs[i] - coeffs_4[i] ) > 0.01 ) stopit( 105, status );
   }

   astChebyDomain( cm2, 0, dlbnd, dubnd );
   if( dlbnd[0] != -2.0 ) stopit( 106, status );
   if( dlbnd[1] != -2.0 ) stopit( 107, status );
   if( dubnd[0] != 2.0 ) stopit( 108, status );
   if( dubnd[1] != 2.0 ) stopit( 109, status );

   /* 2D ChebyMap with non-unit domain */
   lbnd[0] = 0.0; lbnd[1] = 0.0;
   ubnd[0] = 10.0; ubnd[1] = 10.0;
   cm = astChebyMap( 2, 2, 5, coeffs_3, 0, NULL, lbnd, ubnd, lbnd, ubnd, " " );

   xin[0] = 0.0; xin[1] = 2.0; xin[2] = 6.0; xin[3] = 10.0;
   yin[0] = 2.0; yin[1] = 5.0; yin[2] = 8.0; yin[3] = 0.0;
   astTran2( cm, 4, xin, yin, 1, xout, yout );
   for( i = 0; i < 4; i++ ) {
      xi = 2.0*(xin[i] - lbnd[0])/(ubnd[0] - lbnd[0]) - 1.0;
      yi = 2.0*(yin[i] - lbnd[1])/(ubnd[1] - lbnd[1]) - 1.0;
      xv = 1 - 2*xi*(2*yi*yi - 1) + yi;
      yv = 1.5 - 2.5*xi*(2*yi*yi - 1);
      if( fabs( xout[i] - xv ) > 1.0e-6*fabs(xv) ||
          fabs( yout[i] - yv ) > 1.0e-6*fabs(yv) )
         stopit( 10, status );
   }

   /* 2D PolyTran inverse */
   tlbnd[0] = 4.0; tlbnd[1] = 4.0;
   tubnd[0] = 6.0; tubnd[1] = 6.0;
   cm2 = astPolyTran( cm, 0, 0.01, 0.01, 10, tlbnd, tubnd );
   if( !cm2 ) {
      stopit( 11, status );
   } else {
      xin[0] = 4.0; xin[1] = 4.5; xin[2] = 5.0; xin[3] = 5.5;
      yin[0] = 6.0; yin[1] = 5.5; yin[2] = 5.0; yin[3] = 4.5;
      astTran2( cm2, 4, xin, yin, 1, xout, yout );
      astTran2( cm2, 4, xout, yout, 0, xrec, yrec );
      for( i = 0; i < 4; i++ )
         if( fabs( xrec[i] - xin[i] ) > 0.01 ||
             fabs( yrec[i] - yin[i] ) > 0.01 )
            stopit( 12, status );
   }

   /* Test recovery of forward coefficients */
   astPolyCoeffs( cm2, 1, 0, NULL, &nco );
   if( nco != 5 ) stopit( 13, status );

   astPolyCoeffs( cm2, 1, 100, cofs, &nco );
   if( nco != 5 ) {
      stopit( 14, status );
   } else {
      for( i = 0; i < 20; i++ )
         if( cofs[i] != coeffs_3[i] ) stopit( 15, status );
   }

   /* Test recovery of inverse coefficients */
   astPolyCoeffs( cm2, 0, 0, NULL, &nco );
   if( nco != 9 ) stopit( 16, status );

   astPolyCoeffs( cm2, 0, 100, cofs, &nco );
   if( nco != 9 ) {
      stopit( 17, status );
   } else {
      if( fabs( cofs[0] - 5.0000000000000018 ) > 1.0e-6 ) stopit( 18, status );
      if( fabs( cofs[12] - 0.35096188953505458 ) > 1.0e-6 ) stopit( 19, status );
      if( cofs[14] != 2.0 ) stopit( 20, status );
   }

   /* Domain bounding box */
   astChebyDomain( cm, 1, dlbnd, dubnd );
   if( dlbnd[0] != lbnd[0] ) stopit( 21, status );
   if( dlbnd[1] != lbnd[1] ) stopit( 22, status );
   if( dubnd[0] != ubnd[0] ) stopit( 23, status );
   if( dubnd[1] != ubnd[1] ) stopit( 24, status );

   astChebyDomain( cm, 0, dlbnd, dubnd );
   if( dlbnd[0] != -2.0 ) stopit( 25, status );
   if( dlbnd[1] != -1.0 ) stopit( 26, status );
   if( dubnd[0] != 4.0 ) stopit( 27, status );
   if( dubnd[1] != 4.0 ) stopit( 28, status );

   astChebyDomain( cm2, 1, dlbnd, dubnd );
   if( dlbnd[0] != lbnd[0] ) stopit( 29, status );
   if( dlbnd[1] != lbnd[1] ) stopit( 30, status );
   if( dubnd[0] != ubnd[0] ) stopit( 31, status );
   if( dubnd[1] != ubnd[1] ) stopit( 32, status );

   astChebyDomain( cm2, 0, dlbnd, dubnd );
   if( fabs( dlbnd[0] - 0.432 ) > 1.0e-6 ) stopit( 33, status );
   if( fabs( dlbnd[1] - 1.000816 ) > 1.0e-6 ) stopit( 34, status );
   if( fabs( dubnd[0] - 1.568 ) > 1.0e-6 ) stopit( 35, status );
   if( fabs( dubnd[1] - 1.9991836 ) > 1.0e-6 ) stopit( 36, status );

   /* astRate at x0=0 for a ChebyMap defined only over [-1,1]. The
    * derivative is well defined there, but a too-large initial search
    * interval used to push the sample points outside the domain, making
    * astRate return AST__BAD at x0=0 (while x0=0.5 worked). f(x)=3*T1(x)
    * is the linear map f(x)=3x, so the gradient is 3 everywhere. */
   {
      double rate_coeffs[] = { 3.0, 1, 1 };
      double rlbnd = -1.0, rubnd = 1.0;
      double zero_at[1] = { 0.0 }, half_at[1] = { 0.5 }, rr;
      AstChebyMap *rcm = astChebyMap( 1, 1, 1, rate_coeffs, 0, NULL,
                                      &rlbnd, &rubnd, &rlbnd, &rubnd, " " );
      rr = astRate( (AstMapping *)rcm, half_at, 1, 1 );
      if( fabs( rr - 3.0 ) > 1.0e-6 ) stopit( 600, status );
      rr = astRate( (AstMapping *)rcm, zero_at, 1, 1 );
      if( rr == AST__BAD || fabs( rr - 3.0 ) > 1.0e-6 ) stopit( 601, status );
   }

/* A bounding box is required whenever coefficients are supplied for that
   direction. Omitting it must be reported, not read. */
   if( *status == 0 ) {
      double box_coeffs[] = { 1.0, 1, 1 };
      double blo = -1.0, bhi = 1.0;

      cm = astChebyMap( 1, 1, 1, box_coeffs, 0, NULL, NULL, &bhi, NULL, NULL,
                        " " );
      if( *status != AST__NOBOX || cm ) stopit( 610, status );
      astClearStatus;

      cm = astChebyMap( 1, 1, 1, box_coeffs, 0, NULL, &blo, NULL, NULL, NULL,
                        " " );
      if( *status != AST__NOBOX || cm ) stopit( 611, status );
      astClearStatus;

      cm = astChebyMap( 1, 1, 0, NULL, 1, box_coeffs, NULL, NULL, NULL, &bhi,
                        " " );
      if( *status != AST__NOBOX || cm ) stopit( 612, status );
      astClearStatus;

/* Omitting the box for a direction with no coefficients remains valid. */
      cm = astChebyMap( 1, 1, 1, box_coeffs, 0, NULL, &blo, &bhi, NULL, NULL,
                        " " );
      if( !cm || *status != 0 ) stopit( 613, status );
      cm = astAnnul( cm );
   }

   inverse_tests( status );

   astEnd;
   astFlushMemory( 1 );

   if( *status == 0 ) {
      printf( " All ChebyMap tests passed\n" );
   } else {
      printf( "ChebyMap tests failed\n" );
   }
   return *status;
}

static void inverse_near( double got, double want, double tol, int code,
                          int *status ) {
   if( *status != 0 ) return;
   if( got == AST__BAD || !isfinite(got) || fabs(got-want) > tol ) {
      fprintf( stderr, "Inverse test %d: got %.17g, expected %.17g\n",
               code, got, want );
      stopit( code, status );
   }
}

static void inverse_tests( int *status ) {
   AstChebyMap *cm;
   double lo[] = { 0, -3, 10 }, hi[] = { 10, 5, 12 };
   double linear[] = { 2, 1, 1 };
   double quad[] = { 1, 1, 1, .125, 1, 2 };
   double x[441], y[441], u[441], v[441], xr[441], yr[441];
   int i, j, n;
   if( *status != 0 ) return;

   cm = astChebyMap( 1, 1, 1, linear, 0, NULL, lo, hi, NULL, NULL, "" );
   if( !astGetI(cm,"TranInverse") || !astGetI(cm,"IterInverse") ||
       astTest(cm,"IterInverse") ) stopit( 700, status );
   u[0] = -2; u[1] = -1.6; u[2] = 2;
   astTran1( cm, 3, u, 0, xr );
   inverse_near( xr[0], 0, 1e-12, 701, status );
   inverse_near( xr[1], 1, 1e-12, 702, status );
   inverse_near( xr[2], 10, 1e-12, 703, status );
   astSet( cm, "IterInverse=0" );
   if( astGetI(cm,"TranInverse") || astGetI(cm,"IterInverse") ) stopit(704,status);
   astClear( cm, "IterInverse" );
   if( !astGetI(cm,"TranInverse") || astTest(cm,"IterInverse") ) stopit(705,status);
   astInvert( cm );
   astTran1( cm, 3, u, 1, xr );
   inverse_near( xr[1], 1, 1e-12, 706, status );
   cm = astAnnul( cm );

/* Sample targets independently of AST's forward transformation. */
   cm = astChebyMap( 1, 1, 2, quad, 0, NULL, lo, hi, NULL, NULL,
                     "IterInverse=1,NiterInverse=20,TolInverse=1e-12" );
   for( i = 0; i <= 100; i++ ) {
      u[i] = -.875 + 2.0*i/100;
      x[i] = 5*(1 + 2*(u[i]+.125)/(1+sqrt(1+u[i]+.125)));
   }
   astTran1( cm, 101, u, 0, xr );
   astTran1( cm, 101, x, 1, v );
   for( i = 0; i <= 100; i++ ) {
      inverse_near( xr[i], x[i], 1e-10, 710, status );
      inverse_near( v[i], u[i], 1e-12, 711, status );
   }
   u[0] = -1; u[1] = 2; u[2] = AST__BAD; u[3] = NAN; u[4] = INFINITY;
   u[5] = -.125;
   astTran1( cm, 6, u, 0, xr );
   for( i = 0; i < 5; i++ ) if( xr[i] != AST__BAD ) stopit(712,status);
   inverse_near( xr[5], 5, 1e-10, 713, status );
   astSet( cm, "NiterInverse=1" );
   u[0] = .3;
   astTran1( cm, 1, u, 0, xr );
   if( xr[0] != AST__BAD ) stopit(714,status);
   cm = astAnnul( cm );

/* Two polynomial shears have a unique analytic inverse and determinant
   one in normalized coordinates. This expanded single-map representation
   has both a sixth-order term and a mixed T1(z)*T3(w) term. */
   {
      const double a = .125, b = .0625;
      double coeffs[] = { 1, 1, 1, 0, a, 1, 0, 3,
                         1, 2, 0, 1, b, 2, 2, 0, 4*a*b, 2, 1, 3,
                         a*a*b, 2, 0, 6, a*a*b, 2, 0, 0 };
      cm = astChebyMap( 2, 2, 7, coeffs, 0, NULL, lo, hi, NULL, NULL,
                        "IterInverse=1,NiterInverse=30,TolInverse=1e-12" );
      n = 0;
      for( i = 0; i <= 20; i++ ) {
         for( j = 0; j <= 20; j++, n++ ) {
            double z = -1 + i/10.0, w = -1 + j/10.0;
            x[n] = 5*(z+1); y[n] = 4*(w+1)-3;
            u[n] = z + a*(4*w*w*w-3*w);
            v[n] = w + b*(2*u[n]*u[n]-1);
         }
      }
      astTran2( cm, n, x, y, 1, xr, yr );
      for( i = 0; i < n; i++ ) {
         inverse_near( xr[i], u[i], 1e-12, 720, status );
         inverse_near( yr[i], v[i], 1e-12, 721, status );
      }
      astTran2( cm, n, u, v, 0, xr, yr );
      for( i = 0; i < n; i++ ) {
         double w = v[i] - b*(2*u[i]*u[i]-1);
         double z = u[i] - a*(4*w*w*w-3*w);
         inverse_near( xr[i], 5*(z+1), 1e-10, 722, status );
         inverse_near( yr[i], 4*(w+1)-3, 1e-10, 723, status );
      }
      u[0] = 4; v[0] = 2;
      astTran2( cm, 1, u, v, 0, xr, yr );
      if( xr[0] != AST__BAD || yr[0] != AST__BAD ) stopit(724,status);
      cm = astAnnul( cm );
   }

/* A monotonic cubic whose centre linearization overshoots the box. */
   {
      double coeffs[] = { 1, 1, 1, .25, 1, 3 };
      cm = astChebyMap( 1, 1, 2, coeffs, 0, NULL, lo, hi, NULL, NULL,
                        "NiterInverse=30,TolInverse=1e-12" );
      for( i = 0; i <= 100; i++ ) {
         double z = -1 + i/50.0;
         u[i] = z*z*z + .25*z;
         x[i] = 5*(z+1);
      }
      astTran1( cm, 101, u, 0, xr );
      for( i = 0; i <= 100; i++ ) inverse_near(xr[i],x[i],1e-10,730,status);
      cm = astAnnul( cm );
   }

/* A constant forward map: accept an exact solution, reject other targets. */
   {
      double coeffs[] = { 7, 1, 0 };
      cm = astChebyMap( 1, 1, 1, coeffs, 0, NULL, lo, hi, NULL, NULL, "" );
      u[0] = 7; u[1] = 8;
      astTran1( cm, 2, u, 0, xr );
      inverse_near( xr[0], 5, 1e-12, 740, status );
      if( xr[1] != AST__BAD ) stopit(741,status);
      cm = astAnnul( cm );
   }

/* Explicit inverse coefficients remain preferred until iteration is set. */
   {
      double inverse[] = { 42, 1, 0 };
      cm = astChebyMap( 1, 1, 1, linear, 1, inverse, lo, hi, lo, hi, "" );
      if( astGetI(cm,"IterInverse") ) stopit(750,status);
      u[0] = 0;
      astTran1( cm, 1, u, 0, xr );
      inverse_near( xr[0], 42, 0, 751, status );
      astSet( cm, "IterInverse=1" );
      astTran1( cm, 1, u, 0, xr );
      inverse_near( xr[0], 5, 1e-12, 752, status );
      astClear( cm, "IterInverse" );
      astTran1( cm, 1, u, 0, xr );
      inverse_near( xr[0], 42, 0, 753, status );
      cm = astAnnul( cm );
   }

/* Three dimensions, with one domain far from zero. Dyadic samples make
   the reference coordinates representable even on the translated axis. */
   {
      double blo[] = { 0, -3, 1e9 }, bhi[] = { 10, 5, 1e9+16 };
      double coeffs[] = { 1, 1, 1, 0, 0, .125, 1, 0, 2, 0,
                         1, 2, 0, 1, 0, .125, 2, 0, 0, 3,
                         1, 3, 0, 0, 1 };
      double input[3*125], output[3*125], recovered[3*125];
      int k;
      cm = astChebyMap( 3, 3, 5, coeffs, 0, NULL, blo, bhi, NULL, NULL,
                        "NiterInverse=20,TolInverse=1e-12" );
      n = 0;
      for( i = 0; i < 5; i++ ) {
         for( j = 0; j < 5; j++ ) {
            for( k = 0; k < 5; k++, n++ ) {
               double z = -1 + .5*i, w = -1 + .5*j, t = -1 + .5*k;
               input[n] = 5*(z+1);
               input[125+n] = 4*(w+1)-3;
               input[250+n] = 1e9+8*(t+1);
               output[n] = z + .125*(2*w*w-1);
               output[125+n] = w + .125*(4*t*t*t-3*t);
               output[250+n] = t;
            }
         }
      }
      astTranN( cm, n, 3, n, input, 1, 3, n, recovered );
      for( i = 0; i < 3*n; i++ ) inverse_near(recovered[i],output[i],1e-12,760,status);
      astTranN( cm, n, 3, n, output, 0, 3, n, recovered );
      for( i = 0; i < 3*n; i++ ) inverse_near(recovered[i],input[i],1e-10,761,status);
      cm = astAnnul( cm );
   }

/* Branch selection is local; any returned root must reproduce the target. */
   {
      double coeffs[] = { 1, 1, 1, 1, 1, 2 };
      cm = astChebyMap( 1, 1, 2, coeffs, 0, NULL, lo, hi, NULL, NULL,
                        "NiterInverse=30,TolInverse=1e-12" );
      u[0] = -1; u[1] = -.875; u[2] = -1.25;
      astTran1( cm, 3, u, 0, xr );
      astTran1( cm, 2, xr, 1, v );
      inverse_near( v[0], u[0], 1e-11, 770, status );
      inverse_near( v[1], u[1], 1e-11, 771, status );
      if( xr[2] != AST__BAD ) stopit(772,status);
      cm = astAnnul( cm );
   }

/* Reconstructing an endpoint from scale/offset can round just outside
   the forward domain. The inverse should find the adjacent valid value. */
   {
      double blo = -45.701777348010886, bhi = -40.27813046921121;
      cm = astChebyMap( 1, 1, 1, linear, 0, NULL, &blo, &bhi, NULL, NULL,
                        "NiterInverse=20,TolInverse=1e-12" );
      u[0] = -2; u[1] = 2;
      astTran1( cm, 2, u, 0, xr );
      inverse_near( xr[0], blo, 1e-11, 775, status );
      inverse_near( xr[1], bhi, 1e-11, 776, status );
      cm = astAnnul( cm );
   }

/* Unsupported configurations must never advertise an iterative inverse. */
   if( *status == 0 ) {
      double coeffs[] = { 1, 1, 1, 0 };
      cm = astChebyMap( 2, 1, 1, coeffs, 0, NULL, lo, hi, NULL, NULL, "" );
      if( astGetI(cm,"IterInverse") || astGetI(cm,"TranInverse") ) stopit(780,status);
      if( *status == 0 ) {
         astSet( cm, "IterInverse=1" );
         int expected = *status == AST__ATTIN;
         astClearStatus;
         if( !expected ) stopit(781,status);
      }
      cm = astAnnul( cm );
      cm = astChebyMap( 1, 1, 0, NULL, 1, linear, NULL, NULL, lo, hi, "" );
      if( astGetI(cm,"IterInverse") || astGetI(cm,"TranForward") ) stopit(782,status);
      if( *status == 0 ) {
         astSet( cm, "IterInverse=1" );
         int expected = *status == AST__ATTIN;
         astClearStatus;
         if( !expected ) stopit(783,status);
      }
      cm = astAnnul( cm );
   }

/* A loaded normalization can describe a box whose width is below the
   resolution of its own centre. Inverting such a ChebyMap yields bad
   values, but must not report an error from within astTransform. */
   if( *status == 0 ) {
      const char *degenerate = " Begin ChebyMap\n"
                               " Nin = 1\n"
                               " IsA Mapping\n"
                               " MPF1 = 1\n"
                               " NCF1 = 1\n"
                               " CF1 = 1\n"
                               " PF1 = 1\n"
                               " IsA PolyMap\n"
                               " FSCL1 = 1\n"
                               " FOFF1 = 1e17\n"
                               " End ChebyMap\n";
      cm = (AstChebyMap *) astFromString( degenerate );
      if( !cm || !astOK ) {
         astClearStatus;
         stopit( 796, status );
      } else {
         u[ 0 ] = 0.5;
         xr[ 0 ] = 0.0;
         astTran1( cm, 1, u, 0, xr );
         if( *status != 0 ) {
            astClearStatus;
            stopit( 797, status );
         } else if( xr[ 0 ] != AST__BAD ) {
            stopit( 798, status );
         }
         cm = astAnnul( cm );
      }
   }

/* Refitting the forward transformation samples the iterative inverse,
   which evaluates the forward transformation. The fit must therefore see
   the supplied transformation unchanged at every polynomial order, even
   though the fit normalizes its own coefficients differently. */
   if( *status == 0 ) {
      /* f(x) = 500*T1(x/1000) + 100*T2(x/1000) = 0.5x + 2e-4 x^2 - 100,
         monotonic on [-1000,1000] with output range [-400,600]. */
      double fwd[] = { 500.0, 1, 1, 100.0, 1, 2 };
      double inv[] = { 1.0, 1, 1 };
      double ilo[] = { -1000.0 }, ihi[] = { 1000.0 };
      double olo[] = { -400.0 }, ohi[] = { 600.0 };
      double slo[] = { -350.0 }, shi[] = { 0.0 };
      AstChebyMap *refit;

      cm = astChebyMap( 1, 1, 2, fwd, 1, inv, ilo, ihi, olo, ohi,
                        "IterInverse=1" );

/* The sampled sub-range of the output box covers inputs [-691,186]. */
      refit = (AstChebyMap *) astPolyTran( cm, 1, 1e-7, 1e-3, 8, slo, shi );
      if( !refit ) {
         stopit( 790, status );
      } else {
         for( i = 0; i < 5; i++ ) {
            double xi = -600.0 + 150.0*i;
            double want = 0.5*xi + 2.0e-4*xi*xi - 100.0;
            astTran1( refit, 1, &xi, 1, xr );
            if( fabs( xr[ 0 ] - want ) > 1e-3 ) stopit( 791 + i, status );
         }
         refit = astAnnul( refit );
      }
      cm = astAnnul( cm );
   }

/* Earlier versions of AST allowed IterInverse to be set on an
   inverse-only ChebyMap, so dumps recording that combination exist. They
   must still load, with the recorded value having no effect. */
   if( *status == 0 ) {
      const char *legacy = " Begin ChebyMap\n"
                           " Nin = 1\n"
                           " IsA Mapping\n"
                           " MPI1 = 1\n"
                           " NCI1 = 1\n"
                           " CI1 = 1\n"
                           " PI1 = 1\n"
                           " IterInv = 1\n"
                           " IsA PolyMap\n"
                           " ISCL1 = 0.1\n"
                           " IOFF1 = 0\n"
                           " End ChebyMap\n";
      AstObject *obj = astFromString( legacy );
      if( !obj || !astOK ) {
         astClearStatus;
         stopit( 784, status );
      } else {
         cm = (AstChebyMap *) obj;
         if( astGetI( cm, "TranForward" ) ) stopit( 785, status );
         if( astGetI( cm, "IterInverse" ) ) stopit( 786, status );
         if( !astGetI( cm, "TranInverse" ) ) stopit( 787, status );

/* The inverse coefficients survive the load and still evaluate. */
         u[ 0 ] = 5.0;
         xr[ 0 ] = AST__BAD;
         astTran1( cm, 1, u, 0, xr );
         if( fabs( xr[ 0 ] - 0.5 ) > 1e-12 ) stopit( 788, status );
         cm = astAnnul( cm );
      }
   }
}
