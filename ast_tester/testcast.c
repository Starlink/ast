/*
*  Name:
*     testcast

*  Purpose:
*     Check the casting of a Frame into the class of another Frame.

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
*
*     The same three classes cast by copying the Frame as if it were of
*     the class asked for, which is only possible if it is an instance of
*     that class.  They used to do so for any class derived from their
*     own, so that a Frame was copied as if it were of a class it did not
*     belong to: an StcSearchLocation as the larger StcObsDataLocation, a
*     FrameSet as a Plot, a plain CmpFrame as a SpecFluxFrame.  The copy
*     read beyond the end of a Frame that is smaller than the class asked
*     for, and used whatever it found there as the components of that
*     class.  astConvert reaches the CmpFrame case, and read beyond the end
*     of an array when it used a CmpFrame of three axes as a SpecFluxFrame.
*     Such a cast must give no result, and the casts that are possible must
*     give a Frame of the class asked for.

*  Notes:
*     - This is a native C test with no Fortran original.
*     - Public C API only, except for astCast, which is a protected
*       function with no public equivalent.  It is declared here, and is
*       given pointers made by astMakePointer.
*     - The SpecFluxFrame case must run before anything in this program
*       creates a plain CmpFrame, so it comes first.
*/

#include <stdio.h>
#include <string.h>
#include "sae_par.h"
#include "ast.h"

/* The protected function that casts an Object into the class of another. */
AstObject *astCast_( AstObject *, AstObject *, int * );

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

/* Cast "this" into the class of "obj", and check that the result is of the
   class "wanted", or that there is no result if "wanted" is NULL.  Returns
   0 if so, and 1 otherwise. */
static int check_cast( const char *what, void *this, void *obj,
                       const char *wanted, int *status ) {
   AstObject *new;
   const char *class = NULL;
   int code;

   new = astCast_( astMakePointer( this ), astMakePointer( obj ), status );
   if ( astOK && new ) class = astGetC( astMakeId( new ), "Class" );
   code = astStatus;
   if ( !astOK ) astClearStatus;
   if ( code != SAI__OK || ( wanted ? !class || strcmp( class, wanted )
                                    : class != NULL ) ) {
      fprintf( stderr, "testcast: casting %s gave %s and status %d\n", what,
               class ? class : "nothing", code );
      return 1;
   }
   return 0;
}

int main( void ) {
   AstBox *box;
   AstCircle *circle;
   AstCircle *skycircle;
   AstFrame *cmp1;
   AstFrame *cmp2;
   AstFrame *frame;
   AstFrame *nested;
   AstFrameSet *fset;
   AstPlot *plot;
   AstPlot3D *plot3d;
   AstSpecFluxFrame *sff1;
   AstSpecFluxFrame *sff2;
   AstSpecFrame *spec;
   AstStcObsDataLocation *obsdata;
   AstStcSearchLocation *search;
   double bbox[ 4 ] = { 0.0, 0.0, 1.0, 1.0 };
   double bbox3[ 6 ] = { 0.0, 0.0, 0.0, 1.0, 1.0, 1.0 };
   double centre[ 2 ] = { 0.0, 0.0 };
   double corner[ 2 ] = { 1.0, 1.0 };
   double radius = 1.0;
   float gbox[ 4 ] = { 0.0f, 0.0f, 1.0f, 1.0f };
   float gbox3[ 6 ] = { 0.0f, 0.0f, 0.0f, 1.0f, 1.0f, 1.0f };
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
   skycircle = astCircle( astSkyFrame( " " ), 1, centre, &radius, NULL,
                          " " );
   cmp1 = (AstFrame *) astCmpFrame( skycircle, astSpecFrame( " " ), " " );
   nfail += check_no_conversion( "from a CmpFrame that holds a sky Region",
                                 cmp1, cmp2 );
   nfail += check_no_conversion( "to a CmpFrame that holds a sky Region",
                                 cmp2, cmp1 );

/* A plain CmpFrame that has a Domain, and components that have Domains,
   does not match a SpecFluxFrame, and is then cast into the class of the
   SpecFluxFrame. */
   cmp1 = (AstFrame *) astCmpFrame( astFrame( 1, "Domain=A" ),
                                    astFrame( 1, "Domain=B" ), "Domain=C" );
   nfail += check_no_conversion( "from a CmpFrame to a SpecFluxFrame", cmp1,
                                 sff2 );
   nfail += check_no_conversion( "from a SpecFluxFrame to a CmpFrame", sff2,
                                 cmp1 );

/* So is a CmpFrame that is a component of another CmpFrame.  This one has
   three axes, and was then matched as if it had the two axes of a
   SpecFluxFrame, so that the lists of matching axes were read beyond their
   ends. */
   nested = (AstFrame *) astCmpFrame( astCmpFrame( astSkyFrame( " " ),
                                                   astSpecFrame( " " ),
                                                   " " ),
                                      astFrame( 1, "Domain=B" ), " " );
   nfail += check_no_conversion( "from a SpecFluxFrame to a CmpFrame that "
                                 "holds a CmpFrame", sff2, nested );
   nfail += check_no_conversion( "from a CmpFrame that holds a CmpFrame to "
                                 "a SpecFluxFrame", nested, sff2 );

/* Objects of each class that does its own casting, and of classes derived
   from them.  The two Stc classes have the same parent, and an
   StcObsDataLocation is the larger. */
   frame = astFrame( 2, "Domain=A" );
   circle = astCircle( frame, 1, centre, &radius, NULL, " " );
   box = astBox( frame, 1, centre, corner, NULL, " " );
   search = astStcSearchLocation( circle, 0, NULL, " " );
   obsdata = astStcObsDataLocation( circle, 0, NULL, " " );
   fset = astFrameSet( frame, " " );
   plot = astPlot( frame, gbox, bbox, " " );
   plot3d = astPlot3D( astFrame( 3, " " ), gbox3, bbox3, " " );

/* A Frame cannot be cast into a class that it is not an instance of. */
   nfail += check_cast( "an StcSearchLocation into an StcObsDataLocation",
                        search, obsdata, NULL, &status );
   nfail += check_cast( "an StcObsDataLocation into an StcSearchLocation",
                        obsdata, search, NULL, &status );
   nfail += check_cast( "a Circle into a Box", circle, box, NULL, &status );
   nfail += check_cast( "a Circle into an StcSearchLocation", circle, search,
                        NULL, &status );
   nfail += check_cast( "a CmpFrame into a SpecFluxFrame", cmp1, sff2, NULL,
                        &status );
   nfail += check_cast( "a FrameSet into a Plot", fset, plot, NULL, &status );
   nfail += check_cast( "a Plot into a Plot3D", plot, plot3d, NULL,
                        &status );
   nfail += check_cast( "a Circle into a SkyFrame", circle,
                        astSkyFrame( " " ), NULL, &status );

/* It can be cast into its own class and into the classes that its own is
   derived from.  A Region or FrameSet is cast into a class of Frame that
   is not a Region or FrameSet by casting its current Frame. */
   nfail += check_cast( "an StcObsDataLocation into its own class", obsdata,
                        astCopy( obsdata ), "StcObsDataLocation", &status );
   nfail += check_cast( "a Circle into its own class", circle, skycircle,
                        "Circle", &status );
   nfail += check_cast( "a SpecFluxFrame into a CmpFrame", sff2, cmp1,
                        "CmpFrame", &status );
   nfail += check_cast( "a Plot into a FrameSet", plot, fset, "FrameSet",
                        &status );
   nfail += check_cast( "a Plot3D into a FrameSet", plot3d, fset, "FrameSet",
                        &status );
   nfail += check_cast( "a Circle into a Frame", circle, frame, "Frame",
                        &status );
   nfail += check_cast( "a sky Circle into a SkyFrame", skycircle,
                        astSkyFrame( " " ), "SkyFrame", &status );
   nfail += check_cast( "a CmpFrame into a Frame", cmp1, frame, "Frame",
                        &status );
   nfail += check_cast( "a FrameSet into a Frame", fset, frame, "Frame",
                        &status );

   astEnd;
   if ( !astOK ) astClearStatus;

   if ( nfail == 0 ) {
      printf( " All cast tests passed\n" );
      return 0;
   }
   printf( "Cast tests failed (%d)\n", nfail );
   return 1;
}
