/*
*  Name:
*     testcast

*  Purpose:
*     Check that a Frame can be compared with a Frame of an ancestral
*     class before any Frame of an intermediate class has been used.

*  Description:
*     When astConvert or astFindFrame finds that a template Frame does not
*     match a target, it casts the template into the class of the target
*     and tries again.  The CmpFrame, Region and FrameSet classes each do
*     the cast by comparing the class of the target with their own class,
*     which needs the virtual function table of that class.  That table
*     used to be set up only when an object of exactly that class was
*     created, so the comparison read a table full of zeros, and the
*     library crashed, if only objects of derived classes existed:
*
*       - A SpecFluxFrame is a CmpFrame, and could not be compared with
*         another SpecFluxFrame of a different Domain until some plain
*         CmpFrame had been created.
*       - A Region of any class could never be compared in this way, since
*         there is no such thing as a plain Region.  A CmpFrame that holds
*         a Region as one of its components was enough to reach it.
*
*     Each conversion here has no solution, so the test requires a NULL
*     pointer and no error.

*  Notes:
*     - This is a native C test with no Fortran original.
*     - Public C API only.
*     - The SpecFluxFrame case must run before anything in this program
*       creates a plain CmpFrame, so it comes first.
*/

#include <stdio.h>
#include "sae_par.h"
#include "ast.h"

/* Check that converting "from" to "to" finds no solution and reports no
   error.  Returns 0 if so, and 1 otherwise. */
static int check_no_conversion( const char *what, void *from, void *to ) {
   AstFrameSet *fs = astConvert( from, to, " " );
   int code = astStatus;
   if ( !astOK ) astClearStatus;
   if ( fs || code != SAI__OK ) {
      fprintf( stderr, "testcast: converting %s gave %s and status %d\n",
               what, fs ? "a FrameSet" : "no FrameSet", code );
      if ( fs ) fs = astAnnul( fs );
      return 1;
   }
   return 0;
}

int main( void ) {
   AstCircle *circle;
   AstFrame *cmp1;
   AstFrame *cmp2;
   AstSpecFluxFrame *sff1;
   AstSpecFluxFrame *sff2;
   AstSpecFrame *spec;
   double centre[ 2 ] = { 0.0, 0.0 };
   double radius = 1.0;
   int nfail = 0;
   int status = SAI__OK;

   astWatch( &status );
   astBegin;

/* Two SpecFluxFrames with different Domains.  No plain CmpFrame has been
   created at this point. */
   spec = astSpecFrame( "Domain=FOO" );
   sff1 = astSpecFluxFrame( spec, astFluxFrame( AST__BAD, astSpecFrame( " " ),
                                                " " ), " " );
   sff2 = astSpecFluxFrame( astSpecFrame( " " ),
                            astFluxFrame( AST__BAD, astSpecFrame( " " ), " " ),
                            " " );
   nfail += check_no_conversion( "between SpecFluxFrames with different "
                                 "Domains", sff1, sff2 );
   nfail += check_no_conversion( "between SpecFluxFrames with different "
                                 "Domains, reversed", sff2, sff1 );

/* A CmpFrame that holds a Region, and one that holds a Frame with another
   Domain in its place. */
   circle = astCircle( astFrame( 2, "Domain=A" ), 1, centre, &radius, NULL,
                       " " );
   cmp1 = (AstFrame *) astCmpFrame( circle, astSpecFrame( " " ), " " );
   cmp2 = (AstFrame *) astCmpFrame( astFrame( 2, "Domain=B" ),
                                    astSpecFrame( " " ), " " );
   nfail += check_no_conversion( "from a CmpFrame that holds a Region", cmp1,
                                 cmp2 );
   nfail += check_no_conversion( "to a CmpFrame that holds a Region", cmp2,
                                 cmp1 );

/* The same, for a Region in a sky coordinate system. */
   circle = astCircle( astSkyFrame( " " ), 1, centre, &radius, NULL, " " );
   cmp1 = (AstFrame *) astCmpFrame( circle, astSpecFrame( " " ), " " );
   nfail += check_no_conversion( "from a CmpFrame that holds a sky Region",
                                 cmp1, cmp2 );
   nfail += check_no_conversion( "to a CmpFrame that holds a sky Region",
                                 cmp2, cmp1 );

   astEnd;
   if ( !astOK ) astClearStatus;

   if ( nfail == 0 ) {
      printf( " All cast tests passed\n" );
      return 0;
   }
   printf( "Cast tests failed (%d)\n", nfail );
   return 1;
}
