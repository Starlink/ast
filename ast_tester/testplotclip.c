/*
*  Name:
*     testplotclip

*  Purpose:
*     Check that a Plot survives an error raised while clipping.

*  Description:
*     A Plot with a clipping Frame set by astClip transforms each position
*     into that Frame, after it has transformed it between the base and
*     current Frames, so that it can be compared with the clipping bounds.
*     If this second step fails, the Plot has already filled the PointSet
*     in which its caller wants the results.  The caller owns that PointSet
*     and annuls it, so the Plot must report the error without annulling it
*     as well.

*     This test gives a Plot a clipping Frame that is reached through an
*     IntraMap whose transformation function can be told to fail.  It
*     checks that clipping works while the function succeeds, that
*     transforming positions and drawing with the Plot report the error
*     while the function fails, and that the Plot gives the original
*     results once the function succeeds again.

*  Notes:
*     - This is a native C test with no Fortran original.
*     - Public C API only.
*     - Annulling a PointSet twice is a use of freed memory, which this
*       test can only detect reliably when it is run under
*       AddressSanitizer.
*     - No graphics are drawn.  astBorder and astGrid are only called while
*       the transformation function is failing, and the error is raised
*       before either reaches the graphics system.
*/

#include <stdio.h>
#include "sae_par.h"
#include "ast.h"
#include "ast_err.h"

/* Non-zero if the IntraMap transformation function should fail. */
static int fail_now = 0;

/* An IntraMap transformation function that copies its inputs to its
   outputs, or reports an error if "fail_now" is set. */
static void copy_or_fail( AstMapping *this, int npoint, int ncoord_in,
                          const double *ptr_in[], int forward,
                          int ncoord_out, double *ptr_out[] ) {
   int coord;
   int point;
   (void) this;
   (void) ncoord_in;
   (void) forward;
   if ( fail_now ) {
      astSetStatus( AST__ITFER );
      return;
   }
   for ( coord = 0; coord < ncoord_out; coord++ ) {
      for ( point = 0; point < npoint; point++ ) {
         ptr_out[ coord ][ point ] = ptr_in[ coord ][ point ];
      }
   }
}

/* Transform two positions from graphics coordinates using the Plot, and
   check that the first, which is inside the clipping box, is unchanged and
   the second, which is outside it, is bad.  Returns the number of
   failures. */
static int check_clipping( const char *when, AstPlot *plot ) {
   double xin[ 2 ] = { 0.5, 0.95 };
   double yin[ 2 ] = { 0.5, 0.5 };
   double xout[ 2 ];
   double yout[ 2 ];

   astTran2( plot, 2, xin, yin, 1, xout, yout );
   if ( !astOK ) {
      fprintf( stderr, "testplotclip: could not transform positions %s\n",
               when );
      astClearStatus;
      return 1;
   }
   if ( xout[ 0 ] != 0.5 || yout[ 0 ] != 0.5 || xout[ 1 ] != AST__BAD ||
        yout[ 1 ] != AST__BAD ) {
      fprintf( stderr, "testplotclip: wrong positions %s: (%g,%g) and "
               "(%g,%g)\n", when, xout[ 0 ], yout[ 0 ], xout[ 1 ],
               yout[ 1 ] );
      return 1;
   }
   return 0;
}

/* Check that the operation just attempted reported the error raised by the
   transformation function, and clear it.  Returns 0 if it did, 1
   otherwise. */
static int check_failed( const char *what ) {
   int code = astStatus;
   if ( !astOK ) astClearStatus;
   if ( code != AST__ITFER ) {
      fprintf( stderr, "testplotclip: %s gave status %d, expected %d\n",
               what, code, AST__ITFER );
      return 1;
   }
   return 0;
}

int main( void ) {
   int status = SAI__OK;
   int nfail = 0;
   float gbox[ 4 ] = { 0.0f, 0.0f, 1.0f, 1.0f };
   double bbox[ 4 ] = { 0.0, 0.0, 1.0, 1.0 };
   double lbnd[ 2 ] = { 0.1, 0.1 };
   double ubnd[ 2 ] = { 0.9, 0.9 };
   double xin[ 1 ] = { 0.5 };
   double yin[ 1 ] = { 0.5 };
   double xout[ 1 ];
   double yout[ 1 ];
   AstFrameSet *fset;
   AstPlot *plot;

   astWatch( &status );
   astBegin;

/* Build a FrameSet in which the base Frame is joined to a clipping Frame
   by the IntraMap, and to the current Frame by a UnitMap.  The clipping
   Frame is then neither the base nor the current Frame, so a Plot has to
   use the IntraMap to find positions within it. */
   astIntraReg( "CopyOrFail", 2, 2, copy_or_fail, 0,
                "Copy positions, or fail on request", "AST test suite", "" );
   fset = astFrameSet( astFrame( 2, "Domain=GRID" ), " " );
   astAddFrame( fset, AST__BASE, astIntraMap( "CopyOrFail", 2, 2, " " ),
                astFrame( 2, "Domain=CLIP" ) );
   astAddFrame( fset, AST__BASE, astUnitMap( 2, " " ),
                astFrame( 2, "Domain=SKY" ) );

/* Create the Plot and clip it using a box in the clipping Frame.  The Plot
   puts its graphics Frame before the Frames of the FrameSet, which makes
   the clipping Frame the third Frame. */
   plot = astPlot( fset, gbox, bbox, " " );
   astClip( plot, 3, lbnd, ubnd );
   if ( !astOK ) {
      fprintf( stderr, "testplotclip: could not build the Plot\n" );
      astClearStatus;
      nfail++;

   } else {

/* Clipping works while the transformation function succeeds. */
      nfail += check_clipping( "before the error", plot );

/* While it fails, each operation that clips positions reports the error. */
      fail_now = 1;
      astTran2( plot, 1, xin, yin, 1, xout, yout );
      nfail += check_failed( "astTran2" );
      astBorder( plot );
      nfail += check_failed( "astBorder" );
      astGrid( plot );
      nfail += check_failed( "astGrid" );
      fail_now = 0;

/* The Plot is unharmed. */
      nfail += check_clipping( "after the error", plot );
   }

   astEnd;
   if ( !astOK ) {
      fprintf( stderr, "testplotclip: error on exit\n" );
      astClearStatus;
      nfail++;
   }

   if ( nfail == 0 ) {
      printf( " All Plot clipping tests passed\n" );
      return 0;
   }
   printf( "Plot clipping tests failed (%d)\n", nfail );
   return 1;
}
