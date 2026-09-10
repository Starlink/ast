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

static void seeds( int *status ) {
   double lo = 0, hi = 10;
   double coeffs[] = { 2, 1, 1 };
   double target[] = { -2, -1.6, 2 }, got[3], dlo, dhi;
   AstChebyMap *cm = astChebyMap( 1, 1, 1, coeffs, 0, NULL,
                                  &lo, &hi, NULL, NULL, "", status );
   AstMapping *guess = astLinearGuess( cm );
   astTran1( guess, 3, target, 0, got );
   near( got[0], 0, "Affine seed lower endpoint" );
   near( got[1], 1, "Affine seed physical normalization" );
   near( got[2], 10, "Affine seed upper endpoint" );
   guess = astAnnul( guess );
   astInvert( cm );
   astChebyDomain( cm, 0, &dlo, &dhi );
   near( dlo, lo, "Original forward lower bound" );
   near( dhi, hi, "Original forward upper bound" );
   guess = astLinearGuess( cm );
   astTran1( guess, 3, target, 0, got );
   near( got[1], 1, "Seed ignores Invert" );
   guess = astAnnul( guess );
   cm = astAnnul( cm );

/* T3 has slope -3 at the centre although it has no T1 coefficient. */
   coeffs[2] = 3;
   cm = astChebyMap( 1, 1, 1, coeffs, 0, NULL,
                     &lo, &hi, NULL, NULL, "", status );
   guess = astLinearGuess( cm );
   target[0] = 1.2;
   astTran1( guess, 1, target, 0, got );
   near( got[0], 4, "Higher-order contribution to centre Jacobian" );
   guess = astAnnul( guess );
   cm = astAnnul( cm );

/* A singular constant map has a usable midpoint seed, not an AST error. */
   coeffs[2] = 0;
   cm = astChebyMap( 1, 1, 1, coeffs, 0, NULL,
                     &lo, &hi, NULL, NULL, "", status );
   guess = astLinearGuess( cm );
   astTran1( guess, 1, target, 0, got );
   near( got[0], 5, "Midpoint fallback" );
   guess = astAnnul( guess );
   cm = astAnnul( cm );
}

static void caches( int *status ) {
   double lo = -1, hi = 1, ilo = -3, ihi = 3;
   double forward[] = { 2, 1, 1 }, inverse[] = { 1, 1, 1 };
   double target = .6, got;
   size_t before;
   char *dump;
   AstChebyMap *cm = astChebyMap( 1, 1, 1, forward, 1, inverse,
                                  &lo, &hi, &ilo, &ihi, "", status );
   AstChebyMap *copy, *fitted, *loaded;
   AstMapping *guess;
   (void) astGetJacobian( cm );
   before = astGetObjSize( cm );
   guess = astLinearGuess( cm );
   check( astGetObjSize(cm) > before, "Object size includes the affine cache" );
   guess = astAnnul( guess );

   copy = astCopy( cm );
   astSet( copy, "IterInverse=1,NiterInverse=0", status );
   astTran1( copy, 1, &target, 0, &got );
   near( got, .3, "Copy retains original forward transformation" );
   dump = astToString( copy );
   loaded = astFromString( dump );
   dump = astFree( dump );
   copy = astAnnul( copy );
   astTran1( loaded, 1, &target, 0, &got );
   near( got, .3, "Reload reconstructs inverse caches" );
   loaded = astAnnul( loaded );

/* Fitting the forward transformation to the explicit inverse changes
   2*x to 3*x. Zero iterations expose any stale affine initial guess. */
   fitted = astPolyTran( cm, 1, 1e-12, 1e-12, 2, &ilo, &ihi );
   check( fitted != NULL, "Forward fit succeeded" );
   if( fitted ) {
      astSet( fitted, "IterInverse=1,NiterInverse=0", status );
      astTran1( fitted, 1, &target, 0, &got );
      near( got, .2, "Forward fit invalidates the copied seed" );
      fitted = astAnnul( fitted );
   }
   cm = astAnnul( cm );

/* Legacy ChebyMap dumps can define ordinary polynomials by omitting
   normalization coefficients. Their inverse uses the parent hooks. */
   loaded = astFromString( " Begin ChebyMap\n Nin = 1\n IsA Mapping\n"
                          " MPF1 = 1\n NCF1 = 1\n CF1 = 2\n PF1 = 1\n"
                          " IsA PolyMap\n End ChebyMap\n" );
   check( loaded != NULL, "Load ordinary-polynomial ChebyMap" );
   if( loaded ) {
      target = 4;
      astTran1( loaded, 1, &target, 0, &got );
      near( got, 2, "Ordinary-polynomial fallback outside Chebyshev interval" );
      loaded = astAnnul( loaded );
   }

/* A normalization whose evaluable interval is narrower than the spacing of
   the representable values around it has no usable domain: the bounds
   recovered from it cannot be moved to positions the forward evaluator
   accepts. */
   loaded = astFromString( " Begin ChebyMap\n Nin = 1\n IsA Mapping\n"
                          " MPF1 = 1\n NCF1 = 1\n CF1 = 2\n PF1 = 1\n"
                          " IsA PolyMap\n"
                          " FSCL1 = 368727204320.32996\n"
                          " FOFF1 = -15386044096693470\n"
                          " End ChebyMap\n" );
   check( loaded != NULL, "Load ChebyMap with an unresolvable box" );
   if( loaded ) {
      double dlo = -99, dhi = 99;
      astChebyDomain( loaded, 1, &dlo, &dhi );
      check( dlo == AST__BAD && dhi == AST__BAD,
             "Unresolvable box has no evaluable domain" );
      loaded = astAnnul( loaded );
   }

/* A zero scale describes no bounding box, so there is no finite domain to
   restrict iteration to. */
   loaded = astFromString( " Begin ChebyMap\n Nin = 1\n IsA Mapping\n"
                          " MPF1 = 1\n NCF1 = 1\n CF1 = 2\n PF1 = 1\n"
                          " IsA PolyMap\n FSCL1 = 0\n FOFF1 = 0\n"
                          " End ChebyMap\n" );
   check( loaded != NULL, "Load ChebyMap with a zero scale" );
   if( loaded ) {
      double dlo = -99, dhi = 99;
      astChebyDomain( loaded, 1, &dlo, &dhi );
      check( dlo == AST__BAD && dhi == AST__BAD,
             "Zero scale has no evaluable domain" );
      loaded = astAnnul( loaded );
   }
}

/* The bounded solver, driven through a ChebyMap on [-1,1]. The monomial
   -0.125 + x + 0.25 x^2 is 0 T0 + 1 T1 + 0.125 T2 because x^2 = (T0+T2)/2. */
static void bounded_solver( int *status ) {
   double lo = -1, hi = 1;
   double coeffs[] = { 1, 1, 1, .125, 1, 2 };
   double target[] = { -.875, -.125, 0, .3, 1.125, -1, 2, AST__BAD, NAN, INFINITY };
   double mono[] = { -.125, 1, 0, 1, 1, 1, .25, 1, 2 };
   double got[10];
   AstChebyMap *cm = astChebyMap( 1, 1, 2, coeffs, 0, NULL, &lo, &hi,
                                  NULL, NULL, "", status );
   AstPolyMap *pm;
   int i;
   astSet( cm, "NiterInverse=20,TolInverse=1e-12", status );
   astTran1( cm, 10, target, 0, got );
   for( i = 0; i < 5; i++ ) {
      near( got[i], 2*(target[i]+.125)/(1+sqrt(1+target[i]+.125)),
            "Bounded quadratic inverse" );
   }
   for( i = 5; i < 10; i++ ) {
      check( got[i] == AST__BAD, "Unsolved bounded input must return BAD" );
   }
   astSet( cm, "NiterInverse=1", status );
   astTran1( cm, 1, target+3, 0, got );
   check( got[0] == AST__BAD, "Exhaustion must not return the last iterate" );
   cm = astAnnul( cm );

/* The same polynomial as an unbounded PolyMap retains its historical
   last-iterate behavior. */
   pm = astPolyMap( 1, 1, 3, mono, 0, NULL, "NiterInverse=1", status );
   astTran1( pm, 1, target+3, 0, got );
   check( got[0] != AST__BAD && isfinite(got[0]), "Legacy PolyMap exhaustion" );
   pm = astAnnul( pm );
}

/* A forward series with no linear term has a singular Jacobian at the
   domain midpoint, which is where the fallback seed lands. The solver must
   step off that point rather than report every target as unsolved. */
static void singular_seed( int *status ) {
   double lo = 0, hi = 10;
   double coeffs[] = { 1, 1, 2 };
   double target[] = { 0, .5, -.5, .9, -.99 };
   double got[5], back[5];
   double lo2[] = { 0, 0 }, hi2[] = { 10, 10 };
   double coeffs2[] = { 1, 1, 2, 0,  1, 2, 0, 1 };
   double tx[] = { 0, .5, -.5 }, ty[] = { .2, -.4, .9 };
   double gx[3], gy[3], bx[3], by[3];
   AstChebyMap *cm;
   int i;

   cm = astChebyMap( 1, 1, 1, coeffs, 0, NULL, &lo, &hi, NULL, NULL, "", status );
   astSet( cm, "TolInverse=1e-12", status );
   check( astGetI( cm, "TranInverse" ) == 1, "Pure T2 offers an inverse" );
   astTran1( cm, 5, target, 0, got );
   astTran1( cm, 5, got, 1, back );
   for( i = 0; i < 5; i++ ) {
      check( got[i] != AST__BAD, "Pure T2 target solved" );
      near( back[i], target[i], "Pure T2 round trip" );
   }
   cm = astAnnul( cm );

   cm = astChebyMap( 2, 2, 2, coeffs2, 0, NULL, lo2, hi2, NULL, NULL, "", status );
   astSet( cm, "TolInverse=1e-12", status );
   astTran2( cm, 3, tx, ty, 0, gx, gy );
   astTran2( cm, 3, gx, gy, 1, bx, by );
   for( i = 0; i < 3; i++ ) {
      check( gx[i] != AST__BAD && gy[i] != AST__BAD, "T2,T1 target solved" );
      near( bx[i], tx[i], "T2,T1 round trip x" );
      near( by[i], ty[i], "T2,T1 round trip y" );
   }
   cm = astAnnul( cm );
}

int main( void ) {
   int status_value = 0;
   int *status = &status_value;
   astWatch( status );
   astBegin_();
   derivatives( status );
   seeds( status );
   bounded_solver( status );
   singular_seed( status );
   caches( status );
   astEnd_( status );
   astFlushMemory( 1 );
   check( astOK, "Unexpected AST error" );
   if( !failures ) puts( "All ChebyMap inverse tests passed" );
   return failures ? 1 : 0;
}
