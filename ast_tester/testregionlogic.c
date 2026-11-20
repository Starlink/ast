/*
*  Name:
*     testregionlogic

*  Purpose:
*     Check results that depend on how Regions keep and combine their state.

*  Description:
*     Each case here gave a wrong answer without any error, because a
*     Region kept information about itself that it did not bring up to date
*     when it was moved, or combined the results for its parts wrongly.
*     The cases cover:
*
*       - Moving: a Region centres its uncertainty Region on itself.  An
*         Interval so moved must report its new bounds.
*       - Prisms: the overlap of a Prism with a Box and with another Prism,
*         where the second Region has been negated.  The negation was
*         allowed for twice.
*       - CmpRegions: the intersection of two Regions that are not bounded,
*         one inside the other, is the inner Region and so is not bounded.
*         It was reported to be bounded, and asking for its bounds was an
*         error.

*  Notes:
*     - This is a native C test with no Fortran original.
*     - Public C API only.  Failures are counted and reported at the end.
*/

#include <math.h>
#include <stdio.h>
#include "sae_par.h"
#include "ast.h"

/* Check that the status is good after the operation described by "what".
   Returns 0 if so, and 1 otherwise, after clearing the status. */
static int check_ok( const char *what ) {
   int code = astStatus;
   if ( astOK ) return 0;
   astClearStatus;
   fprintf( stderr, "testregionlogic: %s gave status %d\n", what, code );
   return 1;
}

/* Check that "value" is within "tol" of "expected".  Returns 0 if so, and
   1 otherwise. */
static int check_double( const char *what, double value, double expected,
                         double tol ) {
   if ( !( fabs( value - expected ) <= tol ) ) {
      fprintf( stderr, "testregionlogic: %s is %.17g but %.17g was "
               "expected\n", what, value, expected );
      return 1;
   }
   return 0;
}

/* Check that the overlap of "reg1" and "reg2" is "forward", and that the
   overlap of "reg2" and "reg1" is "reverse".  Returns the number of
   failures. */
static int check_overlap( const char *what, void *reg1, void *reg2,
                          int forward, int reverse ) {
   int f = astOverlap( reg1, reg2 );
   int r = astOverlap( reg2, reg1 );
   if ( check_ok( what ) ) return 1;
   if ( f != forward || r != reverse ) {
      fprintf( stderr, "testregionlogic: %s gave %d and %d but %d and %d "
               "were expected\n", what, f, r, forward, reverse );
      return 1;
   }
   return 0;
}

/* The uncertainty Region given to astBox is copied and centred on the Box.
   The Interval class did not change the bounds it keeps when it was
   moved. */
static int test_moves( void ) {
   AstFrame *frm;
   AstRegion *box;
   AstRegion *unc;
   double c1[ 2 ] = { 0.0, 0.0 };
   double c2[ 2 ] = { 10.0, 10.0 };
   double hi[ 2 ] = { 0.1, 0.1 };
   double lbnd[ 2 ];
   double lo[ 2 ] = { -0.1, -0.1 };
   double ubnd[ 2 ];
   int i;
   int nfail = 0;

   astBegin;
   frm = astFrame( 2, " " );

   box = (AstRegion *) astBox( frm, 1, c1, c2,
                               (AstRegion *) astInterval( frm, lo, hi, NULL,
                                                          " " ), " " );
   unc = astGetUnc( box, 0 );
   astGetRegionBounds( unc, lbnd, ubnd );
   nfail += check_ok( "The bounds of an Interval used as an uncertainty" );
   for ( i = 0; i < 2; i++ ) {
      nfail += check_double( "The lower bound of an Interval moved to the "
                             "centre of a Box", lbnd[ i ], 4.9, 1.0E-6 );
      nfail += check_double( "The upper bound of an Interval moved to the "
                             "centre of a Box", ubnd[ i ], 5.1, 1.0E-6 );
   }

   astEnd;
   return nfail;
}

/* A Prism made from a Circle of radius 5 and the Interval 2 to 8. */
static AstRegion *cylinder( void ) {
   double centre[ 2 ] = { 0.0, 0.0 };
   double hi[ 1 ] = { 8.0 };
   double lo[ 1 ] = { 2.0 };
   double radius[ 1 ] = { 5.0 };
   return (AstRegion *) astPrism( astCircle( astFrame( 2, " " ), 1, centre,
                                             radius, NULL, " " ),
                                  astInterval( astFrame( 1, " " ), lo, hi,
                                               NULL, " " ), " " );
}

/* The overlap of a Prism with another Region is found from the overlaps
   of its two components with the matching axes of the other Region, and
   then changed to allow for either Region being negated.  The axes picked
   from the other Region were themselves negated if it was, so that its
   negation was allowed for twice. */
static int test_prisms( void ) {
   AstRegion *box;
   AstRegion *negbox;
   AstRegion *negprism;
   double c1[ 3 ] = { -10.0, -10.0, 0.0 };
   double c2[ 3 ] = { 10.0, 10.0, 10.0 };
   int nfail = 0;

   astBegin;
   box = (AstRegion *) astBox( astFrame( 3, " " ), 1, c1, c2, NULL, " " );
   negbox = astCopy( box );
   astNegate( negbox );
   negprism = cylinder();
   astNegate( negprism );

   nfail += check_overlap( "The overlap of a Prism and a Box that holds it",
                           cylinder(), box, 2, 3 );
   nfail += check_overlap( "The overlap of a Prism and the negation of a "
                           "Box that holds it", cylinder(), negbox, 1, 1 );
   nfail += check_overlap( "The overlap of two equal Prisms", cylinder(),
                           cylinder(), 5, 5 );
   nfail += check_overlap( "The overlap of a Prism and its negation",
                           cylinder(), negprism, 6, 6 );

   astEnd;
   return nfail;
}

/* The intersection of the half-planes x >= 0 and x >= 1 is the half-plane
   x >= 1. */
static int test_cmpregions( void ) {
   AstFrame *frm;
   AstRegion *and;
   double bad[ 2 ] = { AST__BAD, AST__BAD };
   double lbnd[ 2 ];
   double one[ 2 ] = { 1.0, AST__BAD };
   double ubnd[ 2 ];
   double zero[ 2 ] = { 0.0, AST__BAD };
   int bounded;
   int nfail = 0;

   astBegin;
   frm = astFrame( 2, " " );
   and = (AstRegion *) astCmpRegion( astInterval( frm, zero, bad, NULL,
                                                  " " ),
                                     astInterval( frm, one, bad, NULL, " " ),
                                     AST__AND, " " );
   bounded = astGetI( and, "Bounded" );
   nfail += check_ok( "Asking if the intersection of two half-planes is "
                      "bounded" );
   if ( bounded ) {
      fprintf( stderr, "testregionlogic: The intersection of two "
               "half-planes, one inside the other, is reported to be "
               "bounded\n" );
      nfail++;
   }
   astGetRegionBounds( and, lbnd, ubnd );
   nfail += check_ok( "The bounds of the intersection of two half-planes" );

   astEnd;
   return nfail;
}

int main( void ) {
   int status = SAI__OK;
   int nfail = 0;

   astWatch( &status );

   nfail += test_moves();
   nfail += test_prisms();
   nfail += test_cmpregions();

   if ( !astOK ) astClearStatus;

   if ( nfail == 0 ) {
      printf( " All Region logic tests passed\n" );
      return 0;
   }
   printf( "Region logic tests failed (%d)\n", nfail );
   return 1;
}
