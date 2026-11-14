/*
*  Name:
*     testplotcopy

*  Purpose:
*     Check that a copy of a Plot owns its graphics context.

*  Description:
*     A Plot holds a KeyMap, its graphics context, that it passes to the
*     drawing functions registered with astGrfSet.  The KeyMap is created
*     when it is first needed, for instance by astGetGrfContext.  A copy of
*     the Plot must have a graphics context of its own, which starts with
*     the same contents as the original, can be changed without changing
*     the original, and is released separately.

*     This test checks each of these for a Plot.  It also copies a Plot3D,
*     which gives each of its three encapsulated Plots a graphics context
*     that records the labelled corner, and checks that the copy can be
*     changed and annulled once the original has gone, and the reverse.

*  Notes:
*     - This is a native C test with no Fortran original.
*     - Public C API only.
*     - Releasing a graphics context twice is a use of freed memory, which
*       this test can only detect reliably when it is run under
*       AddressSanitizer.  The check that the two contexts change
*       independently does not depend on it.
*/

#include <stdio.h>
#include "sae_par.h"
#include "ast.h"

/* Return the integer stored with key "Value" in the graphics context of
   the supplied Plot, or -1 if there is none. */
static int context_value( AstPlot *plot ) {
   AstKeyMap *context = astGetGrfContext( plot );
   int value = -1;
   if ( !astMapGet0I( context, "Value", &value ) ) value = -1;
   context = astAnnul( context );
   return value;
}

/* Store an integer with key "Value" in the graphics context of the
   supplied Plot. */
static void set_context_value( AstPlot *plot, int value ) {
   AstKeyMap *context = astGetGrfContext( plot );
   astMapPut0I( context, "Value", value, " " );
   context = astAnnul( context );
}

/* Check that the graphics context of the supplied Plot holds "wanted".
   Returns 0 if it does, 1 otherwise. */
static int check_value( const char *what, AstPlot *plot, int wanted ) {
   int value = context_value( plot );
   if ( !astOK || value != wanted ) {
      fprintf( stderr, "testplotcopy: %s holds %d, expected %d\n", what,
               value, wanted );
      if ( !astOK ) astClearStatus;
      return 1;
   }
   return 0;
}

int main( void ) {
   int status = SAI__OK;
   int nfail = 0;
   float gbox[ 4 ] = { 0.0f, 0.0f, 1.0f, 1.0f };
   double bbox[ 4 ] = { 0.0, 0.0, 1.0, 1.0 };
   float gbox3d[ 6 ] = { -1.0f, -1.0f, -1.0f, 1.0f, 1.0f, 1.0f };
   double bbox3d[ 6 ] = { -1.0, -1.0, -1.0, 1.0, 1.0, 1.0 };
   AstPlot *plot;
   AstPlot *copy;
   AstPlot3D *plot3d;
   AstPlot3D *copy3d;

   astWatch( &status );

/* A Plot with nothing in its graphics context can be copied, and both
   Plots annulled. */
   plot = astPlot( AST__NULL, gbox, bbox, " " );
   copy = astCopy( plot );
   plot = astAnnul( plot );
   copy = astAnnul( copy );

/* The copy of a Plot starts with the contents of the original's graphics
   context. */
   plot = astPlot( AST__NULL, gbox, bbox, " " );
   set_context_value( plot, 1 );
   copy = astCopy( plot );
   nfail += check_value( "the context of a new copy", copy, 1 );

/* Changing the copy's graphics context leaves the original's alone. */
   set_context_value( copy, 2 );
   nfail += check_value( "the original context after changing the copy",
                         plot, 1 );

/* The copy's graphics context outlives the original Plot. */
   plot = astAnnul( plot );
   nfail += check_value( "the copy's context after annulling the original",
                         copy, 2 );
   copy = astAnnul( copy );

/* A copy of a Plot3D can be changed and annulled after the original has
   been annulled.  Changing the labelled corner updates the graphics
   context of each encapsulated Plot. */
   plot3d = astPlot3D( AST__NULL, gbox3d, bbox3d, "RootCorner=LLL" );
   copy3d = astCopy( plot3d );
   plot3d = astAnnul( plot3d );
   astSet( copy3d, "RootCorner=UUL" );
   copy3d = astAnnul( copy3d );

/* Likewise, the original can be changed and annulled after its copy has
   been annulled. */
   plot3d = astPlot3D( AST__NULL, gbox3d, bbox3d, "RootCorner=LLL" );
   copy3d = astCopy( plot3d );
   copy3d = astAnnul( copy3d );
   astSet( plot3d, "RootCorner=UUL" );
   plot3d = astAnnul( plot3d );

   if ( !astOK ) {
      fprintf( stderr, "testplotcopy: an error was reported\n" );
      astClearStatus;
      nfail++;
   }

   if ( nfail == 0 ) {
      printf( " All Plot copy tests passed\n" );
      return 0;
   }
   printf( "Plot copy tests failed (%d)\n", nfail );
   return 1;
}
