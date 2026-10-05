/*
 * Benchmark astUnlock/astLock of a FITS-WCS FrameSet, which runs
 * ChangeThreadVtab on every component Object each time it is locked.
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "ast.h"

static const char *header[] = {
   "NAXIS   = 2",
   "NAXIS1  = 100",
   "NAXIS2  = 100",
   "CTYPE1  = 'RA---TAN'",
   "CTYPE2  = 'DEC--TAN'",
   "CRPIX1  = 50.0",
   "CRPIX2  = 50.0",
   "CRVAL1  = 10.0",
   "CRVAL2  = 20.0",
   "CDELT1  = -0.001",
   "CDELT2  = 0.001",
   "PC1_1   = 0.9",
   "PC1_2   = 0.1",
   "PC2_1   = -0.1",
   "PC2_2   = 0.9",
   NULL
};

static double now( void ) {
   struct timespec ts;

   clock_gettime( CLOCK_MONOTONIC, &ts );
   return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

/* Create one object of many classes so that the thread's list of known
   vtabs is about as long as in a real application. */
static void init_classes( void ) {
   double lut[ 3 ] = { 0.0, 1.0, 2.0 };
   double mat[ 4 ] = { 1.0, 0.0, 0.0, 1.0 };
   double lbnd[ 2 ] = { 0.0, 0.0 };
   double ubnd[ 2 ] = { 1.0, 1.0 };
   double centre[ 2 ] = { 0.0, 0.0 };
   double radius[ 1 ] = { 1.0 };
   double shift[ 2 ] = { 1.0, 1.0 };
   int perm[ 2 ] = { 2, 1 };
   const char *fwd[ 1 ] = { "y=x*2" };
   const char *inv[ 1 ] = { "x=y/2" };
   AstFrame *frame;

   astBegin;
   frame = astFrame( 2, " " );
   (void) astSkyFrame( " " );
   (void) astSpecFrame( " " );
   (void) astTimeFrame( " " );
   (void) astFluxFrame( 1.0, NULL, " " );
   (void) astDSBSpecFrame( " " );
   (void) astCmpFrame( frame, astFrame( 1, " " ), " " );
   (void) astUnitMap( 2, " " );
   (void) astZoomMap( 2, 2.0, " " );
   (void) astShiftMap( 2, shift, " " );
   (void) astWinMap( 2, lbnd, ubnd, ubnd, lbnd, " " );
   (void) astMatrixMap( 2, 2, 0, mat, " " );
   (void) astPermMap( 2, perm, 2, perm, NULL, " " );
   (void) astLutMap( 3, lut, 0.0, 1.0, " " );
   (void) astMathMap( 1, 1, 1, fwd, 1, inv, " " );
   (void) astSphMap( " " );
   (void) astNormMap( frame, " " );
   (void) astUnitNormMap( 2, centre, " " );
   (void) astKeyMap( " " );
   (void) astBox( frame, 1, lbnd, ubnd, NULL, " " );
   (void) astCircle( frame, 1, centre, radius, NULL, " " );
   (void) astInterval( frame, lbnd, ubnd, NULL, " " );
   (void) astNullRegion( frame, NULL, " " );
   (void) astTable( " " );
   astEnd;
}

int main( int argc, char **argv ) {
   AstFitsChan *fc;
   AstFrameSet *fs;
   double t0;
   double t1;
   int i;
   int n = argc > 1 ? atoi( argv[ 1 ] ) : 200000;
   int status = 0;

   astWatch( &status );
   init_classes();

   fc = astFitsChan( NULL, NULL, " " );
   for( i = 0; header[ i ]; i++ ) {
      astPutFits( fc, header[ i ], 0 );
   }
   astClear( fc, "Card" );
   fs = astRead( fc );
   fc = astAnnul( fc );

   /* Warm up. */
   for( i = 0; i < 1000; i++ ) {
      astUnlock( fs, 0 );
      astLock( fs, 0 );
   }

   t0 = now();
   for( i = 0; i < n; i++ ) {
      astUnlock( fs, 0 );
      astLock( fs, 0 );
   }
   t1 = now();

   printf( "%d unlock/lock pairs: %.1f ns per pair\n", n,
           1e9 * ( t1 - t0 ) / n );
   fs = astAnnul( fs );
   return status != 0;
}
