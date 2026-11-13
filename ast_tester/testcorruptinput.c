/*
*  Name:
*     testcorruptinput

*  Purpose:
*     Check that corrupt or malicious input is rejected without harm.

*  Description:
*     An object dump or a FITS header can come from an untrusted source, so
*     reading one must never crash, corrupt memory, or run for a time or
*     use an amount of memory out of proportion to its size.  Each case
*     here is a small input of the kind produced by a fuzzer, which changes
*     one value in a valid input, and each made the library do one of those
*     things.  The test reads every input and requires either the stated
*     error or, where the input is merely unusual, a usable object.  It is
*     intended to be run with the address and undefined-behavior sanitizers
*     enabled, which is how most of the faults are detected.  The cases
*     cover:
*
*       - Integers: a value too large for an int, which used to wrap round
*         to another value.
*       - Axis limit: a Mapping or Frame cannot have more than AST__MXAXES
*         inputs, outputs or axes, whether it comes from a constructor or
*         from a dump.  The limit exists so that a file that declares an
*         axis count of billions cannot make the library allocate memory
*         and loop in proportion to that count.

*  Notes:
*     - This is a native C test with no Fortran original.
*     - Public C API only.  Errors expected from a deliberately invalid
*       input are caught through astStatus and cleared with astClearStatus
*       so the remaining cases can run; failures of the test itself are
*       counted and reported at the end.
*/

#include <stdio.h>
#include <string.h>
#include "sae_par.h"
#include "ast.h"
#include "ast_err.h"

/* Check that the error "expected" has been reported, for the operation
   described by "what" that was meant to fail and return "obj" (which is
   annulled if it is not NULL).  Returns 0 if so, and 1 otherwise.  Any
   raised status is cleared before returning. */
static int check_failed( const char *what, void *obj, int expected ) {
   int code = astStatus;
   if ( !astOK ) astClearStatus;
   if ( obj ) {
      obj = astAnnul( obj );
      fprintf( stderr, "testcorruptinput: %s gave an object but should "
               "have failed\n", what );
      return 1;
   }
   if ( code != expected ) {
      fprintf( stderr, "testcorruptinput: %s gave status %d but status %d "
               "was expected\n", what, code, expected );
      return 1;
   }
   return 0;
}

/* Attempt to load a dump that should be rejected with the error code
   "expected".  Returns 0 if no object was returned and that error was
   reported, 1 otherwise. */
static int check_rejected( const char *what, const char *dump, int expected ) {
   return check_failed( what, astFromString( dump ), expected );
}

/* Check that the operation described by "what" gave an object, "obj", and
   no error.  The object is annulled.  Returns 0 if so, and 1 otherwise. */
static int check_made( const char *what, void *obj ) {
   int code = astStatus;
   if ( !astOK ) astClearStatus;
   if ( obj ) obj = astAnnul( obj );
   if ( code != SAI__OK ) {
      fprintf( stderr, "testcorruptinput: %s gave status %d but should "
               "have worked\n", what, code );
      return 1;
   }
   return 0;
}

/* Attempt to load a dump that should be accepted.  Returns 0 if an object
   was returned without error, 1 otherwise. */
static int check_accepted( const char *what, const char *dump ) {
   AstObject *obj = astFromString( dump );
   if ( !obj && astOK ) {
      fprintf( stderr, "testcorruptinput: %s gave no object\n", what );
      return 1;
   }
   return check_made( what, obj );
}

/* Integers that cannot be held in an int. */
static int test_integers( void ) {
   int nfail = 0;

   nfail += check_accepted( "Frame with the largest int as MaxAxes",
      " Begin Frame\n Naxes = 1\n MaxAx = 2147483647\n End Frame\n" );
   nfail += check_accepted( "Frame with the smallest int as MatchEnd",
      " Begin Frame\n Naxes = 1\n MchEnd = -2147483648\n End Frame\n" );
   nfail += check_rejected( "UnitMap with Nin one more than the largest int",
      " Begin UnitMap\n Nin = 2147483648\n IsA Mapping\n End UnitMap\n",
      AST__BADIN );
   nfail += check_rejected( "UnitMap with Nin that wraps to a small int",
      " Begin UnitMap\n Nin = 4294967298\n IsA Mapping\n End UnitMap\n",
      AST__BADIN );
   nfail += check_rejected( "UnitMap with Nin beyond a 64 bit int",
      " Begin UnitMap\n Nin = 99999999999999999999\n IsA Mapping\n"
      " End UnitMap\n", AST__BADIN );
   nfail += check_rejected( "UnitMap with Nin one less than the smallest int",
      " Begin UnitMap\n Nin = -2147483649\n IsA Mapping\n End UnitMap\n",
      AST__BADIN );
   nfail += check_rejected( "UnitMap with a floating point Nin",
      " Begin UnitMap\n Nin = 1e308\n IsA Mapping\n End UnitMap\n",
      AST__BADIN );
   nfail += check_rejected( "UnitMap with text after Nin",
      " Begin UnitMap\n Nin = 2 3\n IsA Mapping\n End UnitMap\n",
      AST__BADIN );
   return nfail;
}

/* The limit on the number of axes. */
static int test_axis_limit( void ) {
   AstFrame *big;
   AstFrame *half;
   char dump[ 200 ];
   int nfail = 0;

/* Constructors accept the limit and reject one more. */
   nfail += check_made( "astUnitMap with AST__MXAXES inputs",
                        astUnitMap( AST__MXAXES, " " ) );
   nfail += check_failed( "astUnitMap with too many inputs",
                          astUnitMap( AST__MXAXES + 1, " " ), AST__BADNI );
   nfail += check_failed( "astUnitMap with INT_MAX inputs",
                          astUnitMap( 2147483647, " " ), AST__BADNI );
   nfail += check_failed( "astFrame with too many axes",
                          astFrame( AST__MXAXES + 1, " " ), AST__NAXIN );
   nfail += check_failed( "astFrame with INT_MAX axes",
                          astFrame( 2147483647, " " ), AST__NAXIN );

/* Two Frames that are each within the limit cannot be combined into a
   CmpFrame that is beyond it. */
   big = astFrame( AST__MXAXES, " " );
   half = astFrame( AST__MXAXES/2, " " );
   nfail += check_made( "astCmpFrame with AST__MXAXES axes",
                        astCmpFrame( half, half, " " ) );
   nfail += check_failed( "astCmpFrame with too many axes",
                          astCmpFrame( big, half, " " ), AST__NAXIN );
   if ( big ) big = astAnnul( big );
   if ( half ) half = astAnnul( half );
   if ( !astOK ) astClearStatus;

/* Loaders do the same. */
   sprintf( dump, " Begin UnitMap\n Nin = %d\n IsA Mapping\n End UnitMap\n",
            AST__MXAXES );
   nfail += check_accepted( "UnitMap dump with AST__MXAXES inputs", dump );
   sprintf( dump, " Begin UnitMap\n Nin = %d\n IsA Mapping\n End UnitMap\n",
            AST__MXAXES + 1 );
   nfail += check_rejected( "UnitMap dump with too many inputs", dump,
                            AST__BADNI );
   sprintf( dump, " Begin PermMap\n Nin = 1\n Nout = %d\n IsA Mapping\n"
            " End PermMap\n", AST__MXAXES + 1 );
   nfail += check_rejected( "PermMap dump with too many outputs", dump,
                            AST__BADNO );
   nfail += check_rejected( "ShiftMap dump with INT_MAX inputs",
      " Begin ShiftMap\n Nin = 2147483647\n IsA Mapping\n End ShiftMap\n",
      AST__BADNI );
   sprintf( dump, " Begin Frame\n Naxes = %d\n End Frame\n",
            AST__MXAXES + 1 );
   nfail += check_rejected( "Frame dump with too many axes", dump,
                            AST__NAXIN );
   nfail += check_rejected( "Frame dump with INT_MAX axes",
      " Begin Frame\n Naxes = 2147483647\n End Frame\n", AST__NAXIN );
   nfail += check_rejected( "TimeFrame dump with INT_MAX axes",
      " Begin TimeFrame\n Naxes = 2147483647\n IsA Frame\n"
      " End TimeFrame\n", AST__NAXIN );
   return nfail;
}

int main( void ) {
   int status = SAI__OK;
   int nfail = 0;
   astWatch( &status );
   astBegin;

   nfail += test_integers();
   nfail += test_axis_limit();

   astEnd;
   if ( !astOK ) astClearStatus;

   if ( nfail == 0 ) {
      printf( " All corrupt input tests passed\n" );
      return 0;
   }
   printf( "Corrupt input tests failed (%d)\n", nfail );
   return 1;
}
