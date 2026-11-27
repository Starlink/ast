/*
*  Name:
*     testpolymapinv

*  Purpose:
*     Check that a PolyMap only offers an iterative inverse that it can
*     evaluate.

*  Description:
*     The iterative inverse of a PolyMap finds an input position by
*     repeatedly applying the forward transformation, so it needs the
*     forward transformation to be defined as well as equal numbers of
*     inputs and outputs.  This test checks that:
*
*       - a PolyMap with no forward transformation does not offer an
*         iterative inverse by default, so that one with no coefficients at
*         all has no inverse transformation, and one with only inverse
*         coefficients uses them;
*       - asking such a PolyMap for an iterative inverse through the
*         IterInverse attribute is an error, while a dump written by an
*         earlier version of AST that records one still loads, with the
*         value treated as unset;
*       - a PolyMap with only forward coefficients still has an iterative
*         inverse by default, and it inverts the forward transformation.

*  Notes:
*     - This is a native C test with no Fortran original.
*     - Public C API only.
*/

#include <math.h>
#include <stdio.h>
#include "sae_par.h"
#include "ast.h"
#include "ast_err.h"

/* Check that an integer attribute of a Mapping has the wanted value.
   Returns 0 if it does, 1 otherwise. */
static int check_attr( const char *what, AstMapping *map, const char *attr,
                       int wanted ) {
   int value = astGetI( map, attr );
   if ( !astOK || value != wanted ) {
      fprintf( stderr, "testpolymapinv: %s has %s=%d, expected %d\n", what,
               attr, value, wanted );
      if ( !astOK ) astClearStatus;
      return 1;
   }
   return 0;
}

/* Check that the last operation reported the error "expected", and clear
   it.  Returns 0 if it did, 1 otherwise. */
static int check_status( const char *what, int expected ) {
   int code = astStatus;
   if ( !astOK ) astClearStatus;
   if ( code != expected ) {
      fprintf( stderr, "testpolymapinv: %s gave status %d, expected %d\n",
               what, code, expected );
      return 1;
   }
   return 0;
}

/* Check that transforming "in" with the Mapping, in the direction given by
   "forward", gives "wanted".  Returns 0 if it does, 1 otherwise. */
static int check_tran( const char *what, AstMapping *map, int forward,
                       double in, double wanted ) {
   double out = AST__BAD;
   astTran1( map, 1, &in, forward, &out );
   if ( !astOK || fabs( out - wanted ) > 1.0e-9 ) {
      fprintf( stderr, "testpolymapinv: %s gave %g, expected %g\n", what,
               out, wanted );
      if ( !astOK ) astClearStatus;
      return 1;
   }
   return 0;
}

int main( void ) {
   int status = SAI__OK;
   int nfail = 0;
   double x = 0.5;
   double y;

/* The polynomial 1 + 2*x, as a forward or as an inverse transformation. */
   double coeff[] = { 1.0, 1, 0,
                      2.0, 1, 1 };

/* A polynomial with two inputs and one output. */
   double coeff_2in[] = { 2.0, 1, 1, 0,
                          3.0, 1, 0, 1 };

   AstPolyMap *pm;
   AstObject *obj;

   astWatch( &status );
   astBegin;

/* A PolyMap with no coefficients has no transformation in either
   direction, so it does not offer an iterative inverse, and an attempt to
   use the inverse transformation is refused. */
   pm = astPolyMap( 1, 1, 0, NULL, 0, NULL, " " );
   nfail += check_attr( "an empty PolyMap", (AstMapping *) pm,
                        "TranForward", 0 );
   nfail += check_attr( "an empty PolyMap", (AstMapping *) pm,
                        "IterInverse", 0 );
   nfail += check_attr( "an empty PolyMap", (AstMapping *) pm,
                        "TranInverse", 0 );
   astTran1( pm, 1, &x, 0, &y );
   nfail += check_status( "the inverse of an empty PolyMap", AST__TRNND );

/* The same is true once it has been inverted. */
   astInvert( pm );
   nfail += check_attr( "an inverted empty PolyMap", (AstMapping *) pm,
                        "TranForward", 0 );
   astTran1( pm, 1, &x, 1, &y );
   nfail += check_status( "the forward transformation of an inverted empty "
                          "PolyMap", AST__TRNND );

/* It cannot be given an iterative inverse. */
   astSetI( pm, "IterInverse", 1 );
   nfail += check_status( "setting IterInverse for an empty PolyMap",
                          AST__ATTIN );
   nfail += check_attr( "an empty PolyMap", (AstMapping *) pm,
                        "IterInverse", 0 );

/* A PolyMap with only inverse coefficients uses them for its inverse
   transformation. */
   pm = astPolyMap( 1, 1, 0, NULL, 2, coeff, " " );
   nfail += check_attr( "an inverse-only PolyMap", (AstMapping *) pm,
                        "TranForward", 0 );
   nfail += check_attr( "an inverse-only PolyMap", (AstMapping *) pm,
                        "IterInverse", 0 );
   nfail += check_attr( "an inverse-only PolyMap", (AstMapping *) pm,
                        "TranInverse", 1 );
   nfail += check_tran( "the inverse of an inverse-only PolyMap",
                        (AstMapping *) pm, 0, 0.5, 2.0 );

/* It cannot be given an iterative inverse either, and still uses its
   inverse coefficients after the attempt. */
   astSetI( pm, "IterInverse", 1 );
   nfail += check_status( "setting IterInverse for an inverse-only PolyMap",
                          AST__ATTIN );
   nfail += check_attr( "an inverse-only PolyMap", (AstMapping *) pm,
                        "IterInverse", 0 );
   nfail += check_tran( "the inverse of an inverse-only PolyMap after "
                        "setting IterInverse", (AstMapping *) pm, 0, 0.5,
                        2.0 );

/* The same request made through the options of the constructor is also an
   error. */
   pm = astPolyMap( 1, 1, 0, NULL, 2, coeff, "IterInverse=1" );
   nfail += check_status( "constructing an inverse-only PolyMap with "
                          "IterInverse=1", AST__ATTIN );

/* A PolyMap with only forward coefficients has an iterative inverse by
   default, which inverts the forward transformation. */
   pm = astPolyMap( 1, 1, 2, coeff, 0, NULL, " " );
   nfail += check_attr( "a forward-only PolyMap", (AstMapping *) pm,
                        "IterInverse", 1 );
   nfail += check_attr( "a forward-only PolyMap", (AstMapping *) pm,
                        "TranInverse", 1 );
   nfail += check_tran( "the forward transformation of a forward-only "
                        "PolyMap", (AstMapping *) pm, 1, 0.5, 2.0 );
   nfail += check_tran( "the iterative inverse of a forward-only PolyMap",
                        (AstMapping *) pm, 0, 2.0, 0.5 );

/* With unequal numbers of inputs and outputs it has none. */
   pm = astPolyMap( 2, 1, 2, coeff_2in, 0, NULL, " " );
   nfail += check_attr( "a forward-only PolyMap with two inputs and one "
                        "output", (AstMapping *) pm, "IterInverse", 0 );
   nfail += check_attr( "a forward-only PolyMap with two inputs and one "
                        "output", (AstMapping *) pm, "TranInverse", 0 );

/* A dump of a PolyMap with no coefficients loads, and has no inverse
   transformation. */
   obj = astFromString( " Begin PolyMap\n Nin = 1\n IsA Mapping\n"
                        " End PolyMap\n" );
   if ( astOK && obj ) {
      nfail += check_attr( "a loaded empty PolyMap", (AstMapping *) obj,
                           "TranInverse", 0 );
   } else {
      fprintf( stderr, "testpolymapinv: an empty PolyMap failed to load\n" );
      if ( !astOK ) astClearStatus;
      nfail++;
   }

/* Earlier versions of AST let IterInverse be set without a forward
   transformation, and wrote the value to dumps. Such a dump still loads,
   with the recorded value treated as unset, so the PolyMap uses its
   inverse coefficients and reports no iterative inverse. */
   obj = astFromString( " Begin PolyMap\n Nin = 1\n IsA Mapping\n"
                        " MPI1 = 1\n NCI1 = 2\n CI1 = 1\n PI1 = 0\n"
                        " CI2 = 2\n PI2 = 1\n IterInv = 1\n End PolyMap\n" );
   if ( astOK && obj ) {
      if ( astTest( obj, "IterInverse" ) ) {
         fprintf( stderr, "testpolymapinv: a loaded unusable IterInverse "
                  "value is set\n" );
         nfail++;
      }
      nfail += check_attr( "a loaded inverse-only PolyMap with an unusable "
                           "IterInverse value", (AstMapping *) obj,
                           "IterInverse", 0 );
      nfail += check_tran( "the inverse of a loaded inverse-only PolyMap "
                           "with an unusable IterInverse value",
                           (AstMapping *) obj, 0, 0.5, 2.0 );
   } else {
      fprintf( stderr, "testpolymapinv: a PolyMap with an unusable "
               "IterInverse value failed to load\n" );
      if ( !astOK ) astClearStatus;
      nfail++;
   }

/* Earlier versions of AST accepted any TolInverse and NiterInverse value,
   so a dump may record values the setters now refuse. Such a dump still
   loads, with the unusable values treated as unset, so the defaults apply
   and the iterative inverse still works. */
   obj = astFromString( " Begin PolyMap\n Nin = 1\n IsA Mapping\n"
                        " MPF1 = 1\n NCF1 = 2\n CF1 = 1\n PF1 = 0\n"
                        " CF2 = 2\n PF2 = 1\n NiterInv = -1\n TolInv = 0\n"
                        " End PolyMap\n" );
   if ( astOK && obj ) {
      if ( astTest( obj, "NiterInverse" ) || astTest( obj, "TolInverse" ) ) {
         fprintf( stderr, "testpolymapinv: a loaded unusable NiterInverse "
                  "or TolInverse value is set\n" );
         nfail++;
      }
      nfail += check_attr( "a loaded PolyMap with unusable iteration "
                           "controls", (AstMapping *) obj, "NiterInverse", 4 );
      if ( astGetD( obj, "TolInverse" ) != 1.0E-6 ) {
         fprintf( stderr, "testpolymapinv: a loaded PolyMap with an unusable "
                  "TolInverse value has TolInverse=%g, expected 1e-6\n",
                  astGetD( obj, "TolInverse" ) );
         if ( !astOK ) astClearStatus;
         nfail++;
      }
      nfail += check_tran( "the iterative inverse of a loaded PolyMap with "
                           "unusable iteration controls", (AstMapping *) obj,
                           0, 2.0, 0.5 );
   } else {
      fprintf( stderr, "testpolymapinv: a PolyMap with unusable NiterInverse "
               "and TolInverse values failed to load\n" );
      if ( !astOK ) astClearStatus;
      nfail++;
   }

   astEnd;
   if ( !astOK ) {
      fprintf( stderr, "testpolymapinv: error on exit\n" );
      astClearStatus;
      nfail++;
   }

   if ( nfail == 0 ) {
      printf( " All PolyMap iterative inverse tests passed\n" );
      return 0;
   }
   printf( "PolyMap iterative inverse tests failed (%d)\n", nfail );
   return 1;
}
