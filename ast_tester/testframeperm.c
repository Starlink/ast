/*
*  Name:
*     testframeperm

*  Purpose:
*     Check that astLoadFrame validates the axis permutation array.

*  Description:
*     A Frame dump carries one "Prm" value per axis, each an index into
*     the axes.  astLoadFrame stores these without the library applying any
*     range check, and Dump later inverts the array into a buffer sized to
*     the number of axes.  An out-of-range or repeated index therefore has
*     to be rejected at load time, or it becomes an out-of-bounds write
*     when the loaded object is written back out.

*     This test confirms that an invalid permutation fails the load cleanly
*     and that a valid one - including a non-identity permutation - still
*     loads and round-trips unchanged.  It also loads a CmpFrame and a
*     FrameSet, whose own axis counts are not available while their base
*     Frame is being loaded, to confirm the validation does not disturb the
*     loading of a class that over-rides astGetNaxes.

*     A CmpFrame carries its own axis permutation array (keyword "Axp"),
*     read by the same unchecked path and used as a write index by methods
*     such as astNorm.  The test checks that an invalid CmpFrame array is
*     rejected and that a valid permuted CmpFrame loads and stays usable.

*  Notes:
*     - This is a native C test with no Fortran original.
*     - Public C API only.  Errors expected from a deliberately invalid
*       dump are caught through astOK and cleared with astClearStatus so the
*       remaining cases can run; failures of the test itself are counted and
*       reported at the end.
*/

#include <stdio.h>
#include <string.h>
#include "sae_par.h"
#include "ast.h"
#include "ast_err.h"

/* Attempt to load a dump that should be rejected.  Returns 0 if it was
   rejected cleanly (NULL object, error status raised), 1 otherwise.  Any
   raised status is cleared before returning so the caller can continue. */
static int check_rejected( const char *what, const char *dump ) {
   AstObject *obj = astFromString( dump );
   int bad_status = !astOK;
   if ( obj ) obj = astAnnul( obj );
   if ( !astOK ) astClearStatus;
   if ( obj != NULL || !bad_status ) {
      fprintf( stderr, "testframeperm: %s was accepted but should have "
               "been rejected\n", what );
      return 1;
   }
   return 0;
}

/* Load a dump that should succeed, check its class, and return the object
   (NULL on failure, with a message).  The caller annuls it. */
static AstObject *check_accepted( const char *what, const char *dump,
                                  const char *class ) {
   AstObject *obj = astFromString( dump );
   if ( !astOK || !obj ) {
      fprintf( stderr, "testframeperm: %s failed to load\n", what );
      if ( !astOK ) astClearStatus;
      if ( obj ) obj = astAnnul( obj );
      return NULL;
   }
   if ( strcmp( astGetC( obj, "Class" ), class ) != 0 ) {
      fprintf( stderr, "testframeperm: %s loaded as %s, expected %s\n",
               what, astGetC( obj, "Class" ), class );
      obj = astAnnul( obj );
      return NULL;
   }
   return obj;
}

int main( void ) {
   int status = SAI__OK;
   int nfail = 0;
   astWatch( &status );

/* Invalid permutations must be rejected rather than reaching Dump. */
   nfail += check_rejected( "out-of-range Prm",
      " Begin Frame\n Naxes = 2\n Prm1 = 6\n Prm2 = 1\n End Frame\n" );
   nfail += check_rejected( "large-offset Prm",
      " Begin Frame\n Naxes = 2\n Prm1 = 100000\n Prm2 = 1\n End Frame\n" );
   nfail += check_rejected( "zero (below-range) Prm",
      " Begin Frame\n Naxes = 2\n Prm1 = 0\n Prm2 = 2\n End Frame\n" );
   nfail += check_rejected( "repeated Prm",
      " Begin Frame\n Naxes = 2\n Prm1 = 1\n Prm2 = 1\n End Frame\n" );

/* A valid, non-identity permutation must load and round-trip unchanged. */
   {
      AstFrame *f = astFrame( 3, " " );
      int perm[ 3 ] = { 3, 1, 2 };
      char *dump1, *dump2;
      AstObject *obj;

      astPermAxes( f, perm );
      dump1 = astToString( f );
      f = astAnnul( f );

      obj = check_accepted( "permuted Frame", dump1, "Frame" );
      if ( obj ) {
         dump2 = astToString( obj );
         if ( astOK && strcmp( dump1, dump2 ) != 0 ) {
            fprintf( stderr, "testframeperm: permuted Frame did not "
                     "round-trip unchanged\n" );
            nfail++;
         }
         dump2 = astFree( dump2 );
         obj = astAnnul( obj );
      } else {
         nfail++;
      }
      dump1 = astFree( dump1 );
   }

/* Loading a class whose axis count is not yet available during the base
   Frame load must still work. */
   {
      AstCmpFrame *cmp = astCmpFrame( astFrame( 1, " " ), astFrame( 1, " " ),
                                      " " );
      char *dump = astToString( cmp );
      AstObject *obj;
      cmp = astAnnul( cmp );
      obj = check_accepted( "CmpFrame", dump, "CmpFrame" );
      if ( obj ) obj = astAnnul( obj ); else nfail++;
      dump = astFree( dump );
   }

/* A CmpFrame carries its own axis permutation array (keyword "Axp"), read
   by the same unchecked path as the Frame class.  A valid non-identity
   permutation must load, round-trip unchanged and remain usable - astNorm
   scatters coordinates through the array as a write index - while an
   out-of-range entry must be rejected rather than producing an object that
   writes out of bounds when used. */
   {
      AstCmpFrame *cmp = astCmpFrame( astFrame( 1, "Domain=A" ),
                                      astFrame( 1, "Domain=B" ), " " );
      int perm[ 2 ] = { 2, 1 };
      char *dump1;
      AstObject *obj;

      astPermAxes( cmp, perm );
      dump1 = astToString( cmp );
      cmp = astAnnul( cmp );

      obj = check_accepted( "permuted CmpFrame", dump1, "CmpFrame" );
      if ( obj ) {
         char *dump2 = astToString( obj );
         double point[ 2 ] = { 0.1, 0.2 };
         if ( astOK && strcmp( dump1, dump2 ) != 0 ) {
            fprintf( stderr, "testframeperm: permuted CmpFrame did not "
                     "round-trip unchanged\n" );
            nfail++;
         }
         astNorm( (AstFrame *) obj, point );
         if ( !astOK ) {
            fprintf( stderr, "testframeperm: astNorm failed on a valid "
                     "permuted CmpFrame\n" );
            nfail++;
         }
         dump2 = astFree( dump2 );
         obj = astAnnul( obj );
      } else {
         nfail++;
      }

/* Corrupt the first permutation value to an out-of-range index. */
      {
         char *bad = astStore( NULL, dump1, strlen( dump1 ) + 1 );
         if ( astOK ) {
            char *axp = strstr( bad, "Axp1 = " );
            if ( axp ) {
               axp[ strlen( "Axp1 = " ) ] = '9';
               nfail += check_rejected( "out-of-range CmpFrame Axp", bad );
            } else {
               fprintf( stderr, "testframeperm: could not find Axp1 in the "
                        "CmpFrame dump\n" );
               nfail++;
            }
         }
         bad = astFree( bad );
      }
      dump1 = astFree( dump1 );
   }
   {
      AstFrameSet *fs = astFrameSet( astFrame( 2, " " ), " " );
      char *dump = astToString( fs );
      AstObject *obj;
      fs = astAnnul( fs );
      obj = check_accepted( "FrameSet", dump, "FrameSet" );
      if ( obj ) obj = astAnnul( obj ); else nfail++;
      dump = astFree( dump );
   }

   if ( !astOK ) astClearStatus;

   if ( nfail == 0 ) {
      printf( " All Frame permutation tests passed\n" );
      return 0;
   }
   printf( "Frame permutation tests failed (%d)\n", nfail );
   return 1;
}
