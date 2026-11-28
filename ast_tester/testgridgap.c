/*
*  Name:
*     testgridgap
*  Purpose:
*     Check that a grid line is not drawn across a gap in the sky.
*  Description:
*     Draws Galactic and FK5 grids over a whole-sky HPX image, recording
*     each straight line through a Line function registered with
*     astGrfSet, and checks that the middle of every line, or a point just
*     round it, is on the sky (some lines run along the edge of a facet).
*     A line across the gap between two HPX facets has none of them on it.
*  Notes:
*     - This is a native C test with no Fortran original.
*     - Public C API only.  The other graphics functions are those of the
*       logging GRF module, grf_log.c, as for testgrid.
*/
#include <stdio.h>
#include "sae_par.h"
#include "ast.h"

#define MAXLINE 100000
static double midx[ MAXLINE ];
static double midy[ MAXLINE ];
static int nline = 0;

/* Record the middle of each straight line in a polyline. */
static int record_line( AstKeyMap *grfcon, int n, const float *x,
                        const float *y ) {
   int i;
   (void) grfcon;
   for ( i = 1; i < n && nline < MAXLINE; i++ ) {
      midx[ nline ] = 0.5 * ( x[ i - 1 ] + x[ i ] );
      midy[ nline ] = 0.5 * ( y[ i - 1 ] + y[ i ] );
      nline++;
   }
   return 1;
}

int main( void ) {
   static const char *cards[] = {
      "NAXIS   = 2", "NAXIS1  = 800", "NAXIS2  = 400",
      "CTYPE1  = 'GLON-HPX'", "CTYPE2  = 'GLAT-HPX'",
      "CRVAL1  = 0.0", "CRVAL2  = 0.0", "CRPIX1  = 400.5", "CRPIX2  = 200.5",
      "CDELT1  = -0.45", "CDELT2  = 0.45", NULL
   };
   static const char *systems[] = { "Galactic", "FK5", NULL };
   static const double dx[ 5 ] = { 0.0, 0.002, -0.002, 0.0, 0.0 };
   static const double dy[ 5 ] = { 0.0, 0.0, 0.0, 0.002, -0.002 };
   static double x[ MAXLINE ], y[ MAXLINE ], lon[ MAXLINE ], lat[ MAXLINE ];
   static int on_sky[ MAXLINE ];
   float gbox[ 4 ] = { 0.1f, 0.3f, 0.9f, 0.7f };
   double pbox[ 4 ] = { 0.5, 0.5, 800.5, 400.5 };
   int status = SAI__OK;
   int nfail = 0;
   int i, j, k, ngap;
   AstFitsChan *fc;
   AstFrameSet *fs;
   AstPlot *plot;

   astWatch( &status );

   fc = astFitsChan( NULL, NULL, " " );
   for ( i = 0; cards[ i ]; i++ ) astPutFits( fc, cards[ i ], 0 );
   astClear( fc, "Card" );
   fs = astRead( fc );

   for ( j = 0; systems[ j ] && astOK; j++ ) {
      plot = astPlot( fs, gbox, pbox, "System=%s, Border=0, NumLab=0, "
                      "TextLab=0, DrawTitle=0, MajTickLen=0, MinTickLen=0, "
                      "Grf=1", systems[ j ] );
      astGrfSet( plot, "Line", (AstGrfFun) record_line );
      nline = 0;
      astGrid( plot );
      for ( i = 0; i < nline; i++ ) on_sky[ i ] = 0;
      for ( k = 0; k < 5; k++ ) {
         for ( i = 0; i < nline; i++ ) {
            x[ i ] = midx[ i ] + dx[ k ];
            y[ i ] = midy[ i ] + dy[ k ];
         }
         astTran2( plot, nline, x, y, 1, lon, lat );
         for ( i = 0; i < nline; i++ ) {
            if ( lon[ i ] != AST__BAD ) on_sky[ i ] = 1;
         }
      }
      ngap = 0;
      for ( i = 0; i < nline; i++ ) ngap += !on_sky[ i ];
      if ( !astOK || nline == 0 || ngap > 0 ) {
         fprintf( stderr, "testgridgap: %d of the %d lines in the %s grid "
                  "cross a gap in the sky\n", ngap, nline, systems[ j ] );
         nfail++;
      }
      plot = astAnnul( plot );
   }
   fs = astAnnul( fs );
   fc = astAnnul( fc );

   if ( !astOK ) {
      fprintf( stderr, "testgridgap: an error was reported\n" );
      astClearStatus;
      nfail++;
   }
   if ( nfail == 0 ) {
      printf( " All grid gap tests passed\n" );
      return 0;
   }
   printf( "Grid gap tests failed (%d)\n", nfail );
   return 1;
}
