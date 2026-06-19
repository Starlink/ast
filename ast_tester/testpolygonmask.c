/*
*  Name:
*     testpolygonmask

*  Purpose:
*     Regression test for astMask on a SkyFrame Polygon whose interior is
*     the larger of the two lobes into which its boundary divides the sphere.

*  Description:
*     A Polygon defined within a SkyFrame divides the celestial sphere into
*     two lobes.  Which lobe is the interior is determined by the order in
*     which the vertices are supplied (see the Polygon documentation: the
*     interior is to the left of the boundary as the vertices are traversed).
*     Supplying the vertices in the "reverse" sense therefore selects the
*     LARGER lobe as the (un-negated) interior.  This is a perfectly legal
*     region: a finite region on a sphere has a finite negation, so
*     astGetBounded returns non-zero for it (see polygon.c GetBounded).
*
*     For such a Polygon, astTranN / astPointInRegion correctly report
*     points far from the vertices (but inside the large lobe) as being
*     inside the region.  However astMask<X> does NOT mask those points:
*     it masks zero pixels of a grid that lies entirely inside the region.
*
*     The cause is that astMask uses the region's boundary mesh / bounding
*     box (which encloses only the small vertex outline) and assumes that
*     pixels away from that outline lie outside the region.  That assumption
*     is valid only when the interior is the lobe enclosed by the vertices.
*     When the interior is the larger lobe, astMask and astTranN disagree.
*
*     Note that obtaining the same logical region by negating a
*     small-lobe Polygon (astSetI negated=1) masks correctly; only the
*     winding-selected large lobe triggers the bug.
*
*     This test asserts the CORRECT behaviour (every in-region pixel is
*     masked) and so fails against the buggy implementation.

*  Authors:
*     Demonstration test for the Starlink AST developers.
*/

#include <stdio.h>
#include <math.h>
#include "ast.h"

static void stopit( int *status, const char *text );
static void checkPolygonMaskLargeLobe( int *status );

int main( void ) {
   int status_value = 0;
   int *status = &status_value;

   astWatch( status );
   astBegin;

   checkPolygonMaskLargeLobe( status );

   astEnd;

   if( *status == 0 ) {
      printf( "All polygon-mask tests passed\n" );
   } else {
      printf( "Polygon-mask tests failed\n" );
   }
   return ( *status == 0 ) ? 0 : 1;
}

static void checkPolygonMaskLargeLobe( int *status ) {

/* A small square near the equator; geodesic (great-circle) edges are then
   almost straight, keeping the geometry easy to reason about. */
   const double ra0 = 0.5, dec0 = 0.0, half = 0.001;   /* radians */

/* Grid offset ~0.1 rad (~5.7 deg) from the square: far outside the vertex
   outline, but well inside the large lobe. */
   const double fra = ra0 + 0.1, fdec = dec0;
   const int g = 20;                 /* grid is (g+1) x (g+1) pixels */
   const double span = 0.002;        /* half-extent of the grid on the sky */

   AstSkyFrame *frm;
   AstPolygon *poly;
   AstWinMap *sky2pix;               /* maps sky (region) -> pixel (grid) */
   double verts[2][4];               /* verts[axis][vertex] */
   double pin[2], pout[2];
   double sky_lo[2], sky_hi[2], pix_lo[2], pix_hi[2];
   int lbnd[2], ubnd[2];
   int *data;
   int npix, i, nmasked;

   if( *status != 0 ) return;
   astBegin;

   frm = astSkyFrame( "System=ICRS" );

/* Vertices wound so that the UN-NEGATED interior is the large lobe.  (The
   reverse order would select the small square as the interior.) */
   verts[0][0] = ra0 - half;  verts[1][0] = dec0 - half;
   verts[0][1] = ra0 + half;  verts[1][1] = dec0 - half;
   verts[0][2] = ra0 + half;  verts[1][2] = dec0 + half;
   verts[0][3] = ra0 - half;  verts[1][3] = dec0 + half;
   poly = astPolygon( frm, 4, 4, (const double *)verts, AST__NULL, " " );

/* Guard: confirm the setup really does have the large lobe as its interior,
   using astTranN (interior points are returned unchanged; exterior points
   are set to AST__BAD).  The square's centre must be OUTSIDE, and a far
   point must be INSIDE. */
   pin[0] = ra0;  pin[1] = dec0;
   astTranN( poly, 1, 2, 1, (const double *)pin, 1, 2, 1, (double *)pout );
   if( pout[0] != AST__BAD ) {
      stopit( status, "setup: square centre is not outside the region "
                      "(vertex winding selected the small lobe)" );
   }

   pin[0] = fra;  pin[1] = fdec;
   astTranN( poly, 1, 2, 1, (const double *)pin, 1, 2, 1, (double *)pout );
   if( pout[0] == AST__BAD ) {
      stopit( status, "setup: far point is not inside the large lobe" );
   }
   if( *status != 0 ) {
      astEnd;
      return;
   }

/* Mapping from sky coordinates (the region's Frame) to the pixel grid, as
   required by astMask.  The sky window [sky_lo,sky_hi] maps to the pixel
   window [0,g] on each axis. */
   sky_lo[0] = fra - span;   sky_lo[1] = fdec - span;
   sky_hi[0] = fra + span;   sky_hi[1] = fdec + span;
   pix_lo[0] = 0.0;          pix_lo[1] = 0.0;
   pix_hi[0] = (double) g;   pix_hi[1] = (double) g;
   sky2pix = astWinMap( 2, sky_lo, sky_hi, pix_lo, pix_hi, " " );

/* Every pixel centre in this grid is inside the region (it is a small patch
   far from the boundary, entirely within the large lobe).  So masking with
   inside=1 should set every pixel. */
   lbnd[0] = 0;  lbnd[1] = 0;
   ubnd[0] = g;  ubnd[1] = g;
   npix = ( g + 1 ) * ( g + 1 );
   data = astMalloc( sizeof( int ) * (size_t) npix );
   if( astOK ) {
      for( i = 0; i < npix; i++ ) data[ i ] = 0;

      (void) astMaskI( poly, (AstMapping *) sky2pix, 1, 2, lbnd, ubnd,
                       data, 1 );

      nmasked = 0;
      for( i = 0; i < npix; i++ ) if( data[ i ] != 0 ) nmasked++;

      if( nmasked != npix ) {
         char buf[ 200 ];
         sprintf( buf, "astMaskI masked %d of %d in-region pixels "
                       "(astTranN reports all of them inside)",
                  nmasked, npix );
         stopit( status, buf );
      }
   }
   data = astFree( data );

   astEnd;
}

static void stopit( int *status, const char *text ) {
   if( *status != 0 ) return;
   *status = 1;
   printf( "Error: %s\n", text );
}
