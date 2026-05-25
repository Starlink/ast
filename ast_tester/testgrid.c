/*
*  Name:
*     testgrid.c
*
*  Purpose:
*     Exercise the AST Plot/Grid system using the logging GRF module,
*     with no external graphics library dependency.
*
*  Usage:
*     testgrid [--no-simd] <fits file> <attr> <fattr> [<outfile>] [<xlo> <ylo> <xhi> <yhi>]
*
*     If <outfile> is "-" or omitted, no output is written (smoke-test
*     mode).  If a path ending in ".svg" is given, an SVG image is
*     produced.  Otherwise a GRF call log is written for regression
*     testing.
*
*     --no-simd (or -S) disables SIMD on any SphMap in the WCS mapping
*     before plotting.  The 1-ULP differences between libmvec vector math
*     and scalar libm shift a few integer pixel label positions by 1 px,
*     which makes the SVG output differ from the committed reference for
*     some projections (e.g. car3).  Disabling SIMD restores bit-for-bit
*     agreement with the scalar reference.
*
*  Description:
*     Reads a FITS header from <fits file>, builds a FrameSet via
*     FitsChan, creates an AST Plot against a virtual 0–1 viewport,
*     and calls astGrid.  The logging GRF module (grf_log.c) provides
*     the graphics backend — no PLplot or PGPLOT required.
*
*  Differences from testplotter.c:
*     - No PLplot dependency; uses grf_log instead.
*     - Graphics box computed directly (0.1–0.9 inset) rather than
*       queried from a PLplot viewport.
*     - No output image produced; this is for smoke-testing and
*       call-logging only.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ast.h"
#include "grf_log.h"

/* Recursively walk a (possibly compound) Mapping and set UseSIMD=0 on every
   component that has a SIMD-vectorised transform (SphMap, PolyMap, MatrixMap),
   so the whole pipeline matches the scalar reference bit-for-bit.

   Too bad astMapList isn't part of the public API; it might be nice to have
   for cases like this. */
static void disableSimd( AstMapping *map ) {
   AstMapping *map1 = NULL;
   AstMapping *map2 = NULL;
   int series, inv1, inv2;

   if( !astOK || map == AST__NULL )
      return;

   if( astIsASphMap( map ) || astIsAPolyMap( map ) || astIsAMatrixMap( map ) ) {
      astSet( map, "UseSIMD=0" );
      return;
   }

   astDecompose( map, &map1, &map2, &series, &inv1, &inv2 );

   if( map2 != AST__NULL ) {
      disableSimd( map1 );
      disableSimd( map2 );
   }

   if( map1 != AST__NULL )
      map1 = astAnnul( map1 );

   if( map2 != AST__NULL )
      map2 = astAnnul( map2 );
}

int main( int argc, char **argv ) {
   int status = 0;
   AstFitsChan *fc;
   AstFrameSet *fs;
   AstPlot *pl;
   char *file, *attr1, *attr2;
   FILE *fp;
   FILE *logfp = NULL;
   FILE *svgfp = NULL;
   char card[256];
   double pbox[4];
   float gbox[4];
   float delta, asp;
   int naxis1 = 100, naxis2 = 100;
   int arg_offset;
   int no_simd = 0;

   astWatch( &status );

   /* Optional leading --no-simd / -S flag.  Shift argv so the remaining
      positional argument handling is unchanged. */
   if( argc > 1 && ( strcmp( argv[1], "--no-simd" ) == 0 ||
                     strcmp( argv[1], "-S" ) == 0 ) ) {
      no_simd = 1;
      argv++;
      argc--;
   }

   if( argc < 4 ) {
      printf( "Usage: testgrid <fits file> <attrs> <fattrs> "
              "[<logfile>] [<xlo> <ylo> <xhi> <yhi>]\n" );
      return 1;
   }

   file  = argv[1];
   attr1 = argv[2];
   attr2 = argv[3];

   /* Determine whether a log/SVG file or BOX arguments follow. */
   arg_offset = 4;
   if( argc > 4 ) {
      char *endp;
      (void)strtod( argv[4], &endp );
      if( endp != argv[4] && *endp == '\0' ) {
         arg_offset = 4;
      } else if( strcmp( argv[4], "-" ) != 0 ) {
         const char *ext = strrchr( argv[4], '.' );
         if( ext && strcmp( ext, ".svg" ) == 0 ) {
            svgfp = fopen( argv[4], "w" );
            if( !svgfp ) {
               printf( "Failed to open SVG file %s\n", argv[4] );
               return 1;
            }
         } else {
            logfp = fopen( argv[4], "w" );
            if( !logfp ) {
               printf( "Failed to open log file %s\n", argv[4] );
               return 1;
            }
         }
         arg_offset = 5;
      } else {
         arg_offset = 5;
      }
   }

   /* Initialise the logging GRF module. */
   astGrfLogInit( logfp );
   if( svgfp ) astGrfLogSetSvg( svgfp, 720, 540 );

   fc = astFitsChan( NULL, NULL, "%s", attr2 );

   fp = fopen( file, "r" );
   if( !fp ) {
      printf( "Failed to open file %s\n", file );
      astGrfLogClose();
      if( logfp ) fclose( logfp );
      return 1;
   }

   while( fgets( card, sizeof(card), fp ) != NULL ) {
      size_t len = strlen( card );
      while( len > 0 && ( card[len-1] == '\n' || card[len-1] == '\r' ) ) {
         card[--len] = '\0';
      }
      astPutFits( fc, card, 0 );
   }
   fclose( fp );

   if( argc - arg_offset >= 4 ) {
      pbox[0] = atof( argv[arg_offset]     );
      pbox[1] = atof( argv[arg_offset + 1] );
      pbox[2] = atof( argv[arg_offset + 2] );
      pbox[3] = atof( argv[arg_offset + 3] );
   } else {
      astClear( fc, "Card" );
      if( astFindFits( fc, "NAXIS1", card, 1 ) ) {
         sscanf( card + 10, "%d", &naxis1 );
      } else {
         naxis1 = 100;
      }

      astClear( fc, "Card" );
      if( astFindFits( fc, "NAXIS2", card, 1 ) ) {
         sscanf( card + 10, "%d", &naxis2 );
      } else {
         naxis2 = 100;
      }

      pbox[0] = 0.5;
      pbox[1] = 0.5;
      pbox[2] = (double)naxis1 + 0.5;
      pbox[3] = (double)naxis2 + 0.5;
   }

   astClear( fc, "Card" );
   fs = astRead( fc );

   if( fs == AST__NULL ) {
      printf( "!!! No object read from FitsChan!!!\n" );
      astAnnul( fc );
      astGrfLogClose();
      if( logfp ) fclose( logfp );
      return 1;
   }

   /* Optionally disable SIMD on the base->current mapping.  astGetMapping
      returns a copy of the mapping, so we set UseSIMD=0 on the SIMD-capable
      components in that copy and rebuild the FrameSet around it; the value is
      preserved when the Plot subsequently copies the FrameSet. */
   if( no_simd && astOK ) {
      AstFrame *bfrm = astGetFrame( fs, AST__BASE );
      AstFrame *cfrm = astGetFrame( fs, AST__CURRENT );
      AstMapping *map = astGetMapping( fs, AST__BASE, AST__CURRENT );
      AstFrameSet *nfs;

      disableSimd( map );

      nfs = astFrameSet( bfrm, " " );
      astAddFrame( nfs, AST__BASE, map, cfrm );
      (void) astAnnul( fs );
      fs = nfs;

      bfrm = astAnnul( bfrm );
      cfrm = astAnnul( cfrm );
      map = astAnnul( map );
   }

   if( astOK ) {
      /* Graphics box: inset from the 0–1 virtual viewport, matching
         the ~5% shrink that testplotter applies to the PLplot viewport. */
      gbox[0] = 0.1f;
      gbox[1] = 0.1f;
      gbox[2] = 0.9f;
      gbox[3] = 0.9f;

      /* Aspect-ratio correction (same logic as testplotter.c). */
      asp = (float)( ( pbox[3] - pbox[1] ) / ( pbox[2] - pbox[0] ) );
      if( asp < 0.05f || asp > 20.0f ) asp = 1.0f;

      if( asp > 1.0f ) {
         delta = 0.5f * ( ( gbox[2] - gbox[0] )
                          - ( gbox[3] - gbox[1] ) / asp );
         gbox[2] -= delta;
         gbox[0] += delta;
      } else {
         delta = 0.5f * ( ( gbox[3] - gbox[1] )
                          - asp * ( gbox[2] - gbox[0] ) );
         gbox[3] -= delta;
         gbox[1] += delta;
      }

      pl = astPlot( fs, gbox, pbox, "title=A FITS test" );
      astSet( pl, "%s", attr1 );
      astGrid( pl );
      astAnnul( pl );
   }

   if( fs != AST__NULL ) astAnnul( fs );
   if( fc != AST__NULL ) astAnnul( fc );

   astGrfLogClose();
   if( logfp ) fclose( logfp );
   if( svgfp ) fclose( svgfp );

   return status == 0 ? 0 : 1;
}
