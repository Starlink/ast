/* Test the protected building blocks of ChebyMap's iterative inverse. */
#define astCLASS testchebyinverse
#define THREAD_SAFE 1
#include "error.h"
#include "object.h"
#include "mapping.h"
#include "polymap.h"
#include "chebymap.h"
#include "memory.h"
#include <math.h>
#include <stdio.h>

void astBegin_( void );
void astEnd_( int * );

static int failures = 0;

static void check( int ok, const char *message ) {
   if( !ok ) {
      fprintf( stderr, "%s\n", message );
      failures++;
   }
}

static void near( double got, double want, const char *message ) {
   if( got == AST__BAD || !isfinite( got ) ||
       fabs( got - want ) > 1.0e-11*(1.0 + fabs( want )) ) {
      fprintf( stderr, "%s: got %.17g, expected %.17g\n", message, got, want );
      failures++;
   }
}

/* Independent derivative reference: T_n' = n U_(n-1). */
static double tprime( int n, double z ) {
   double prev = 1.0, value = 2*z, next;
   int i;
   if( n == 0 ) return 0.0;
   if( n == 1 ) return 1.0;
   for( i = 2; i < n; i++ ) {
      next = 2*z*value - prev;
      prev = value;
      value = next;
   }
   return n*value;
}

static void derivatives( int *status ) {
   double lbnd = -3.0, ubnd = 13.0;
   double x[] = { -3.0, -1.0, 5.0, 9.0, 13.0 };
   double got[5];
   int n, i;

   for( n = 0; n <= 12 && astOK; n++ ) {
      double coeffs[] = { 2.5, 1, n };
      AstChebyMap *cm = astChebyMap( 1, 1, 1, coeffs, 0, NULL,
                                     &lbnd, &ubnd, NULL, NULL, "", status );
      AstPolyMap **jac = astGetJacobian( cm );
      if( !astOK ) break;
      check( astGetTranForward( jac[0] ), "Zero derivative must be defined" );
      check( !astGetIterInverse( jac[0] ), "Derivative inverse must be disabled" );
      astTran1( jac[0], 5, x, 1, got );
      for( i = 0; i < 5; i++ ) {
         near( got[i], 2.5/8*tprime( n, (x[i]-5)/8 ), "Scaled derivative" );
      }
      astInvert( cm );
      check( astGetJacobian( cm ) == jac, "Jacobian describes original forward" );
      cm = astAnnul( cm );
   }

/* Cross terms, duplicate orders, cancellation and a constant output.
   F = (3*T3(z)*T2(w), 2*T4(t), 7), with different axis normalizations. */
   {
      double lo[] = { 0, -3, 10 }, hi[] = { 10, 5, 12 };
      double coeffs[] = { 4, 1, 3, 2, 0, -1, 1, 3, 2, 0,
                         2, 2, 0, 0, 4,  7, 3, 0, 0, 0,
                         5, 1, 1, 0, 0, -5, 1, 1, 0, 0 };
      double at[] = { 6, 2, 10.5 }, value[3];
      double z = .2, w = .25, t = -.5;
      AstChebyMap *cm = astChebyMap( 3, 3, 6, coeffs, 0, NULL,
                                     lo, hi, NULL, NULL, "", status );
      AstPolyMap **jac = astGetJacobian( cm );
      for( i = 0; i < 3 && astOK; i++ ) {
         astTranN( jac[i], 1, 3, 1, at, 1, 3, 1, value );
         near( value[0], i == 0 ? 3*.2*tprime(3,z)*(2*w*w-1) :
                        i == 1 ? 3*.25*(4*z*z*z-3*z)*tprime(2,w) : 0,
               "Mixed derivative" );
         near( value[1], i == 2 ? 2*tprime(4,t) : 0, "Independent axis" );
         near( value[2], 0, "Constant output" );
      }
      cm = astAnnul( cm );
   }
}

int main( void ) {
   int status_value = 0;
   int *status = &status_value;
   astWatch( status );
   astBegin_();
   derivatives( status );
   astEnd_( status );
   astFlushMemory( 1 );
   check( astOK, "Unexpected AST error" );
   if( !failures ) puts( "All ChebyMap inverse tests passed" );
   return failures ? 1 : 0;
}
