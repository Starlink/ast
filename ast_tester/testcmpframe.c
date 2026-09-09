/*
 *  Test the CmpFrame class (NormPoints).
 *  Converted from the Fortran test testcmpframe.f.
 *  Direct conversion; no material differences from the Fortran original.
 */
#include "ast.h"
#include <math.h>
#include <stdio.h>

static void stopit( int *status, const char *text ) {
   if( *status != 0 ) return;
   *status = 1;
   printf( "%s\n", text );
}

static void fill_input( double in[ static 12 ], const double ra[ static 4 ],
                        const double px[ static 4 ], const double dec[ static 4 ] ) {
   int i;
   for( i = 0; i < 4; i++ ) {
      in[ i ] = ra[ i ];
      in[ i + 4 ] = px[ i ];
      in[ i + 8 ] = dec[ i ];
   }
}

static int check_normed_output( const double out[ static 12 ],
                                const double ra[ static 4 ],
                                const double px[ static 4 ],
                                const double dec[ static 4 ] ) {
   int i;
   if( out[ 0 ] != ra[ 0 ] - 2*AST__DPI ) return 1;
   if( out[ 1 ] != ra[ 1 ] - 2*AST__DPI ) return 2;
   if( out[ 2 ] != ra[ 2 ] ) return 3;
   if( out[ 3 ] != ra[ 3 ] ) return 4;
   for( i = 0; i < 4; i++ ) {
      if( out[ i + 4 ] != px[ i ] ) return 5;
      if( out[ i + 8 ] != dec[ i ] ) return 6;
   }
   return 0;
}

/* Region bounds on a CmpFrame that holds a SkyFrame go through
 * CmpFrame::NormBox, which asks each component Frame to normalise its own
 * part of the box. The bounds must not depend on which component the
 * SkyFrame is: a Circle near the north pole on (SkyFrame, 1-D) and the same
 * Circle on (1-D, SkyFrame) must report the same limits, with the axes
 * permuted. Returns 0 on success, else the index of the first failed check. */
static int check_normbox_component_order( void ) {
   AstSkyFrame *sky;
   AstFrame *pfrm;
   AstCmpFrame *sky_first, *sky_last;
   AstCircle *c1, *c2;
   double centre1[ 3 ] = { 0.0, 1.2707963267948966, 50.0 };
   double centre2[ 3 ] = { 50.0, 0.0, 1.2707963267948966 };
   double radius = 0.5;
   double lb1[ 3 ], ub1[ 3 ], lb2[ 3 ], ub2[ 3 ];
   int i;

   sky = astSkyFrame( " " );
   pfrm = astFrame( 1, "Domain=FPLANE" );
   sky_first = astCmpFrame( sky, pfrm, " " );
   sky_last = astCmpFrame( pfrm, sky, " " );

   c1 = astCircle( sky_first, 1, centre1, &radius, NULL, " " );
   c2 = astCircle( sky_last, 1, centre2, &radius, NULL, " " );

   astGetRegionBounds( c1, lb1, ub1 );
   if( !astOK ) return 1;
   astGetRegionBounds( c2, lb2, ub2 );
   if( !astOK ) return 2;

   for( i = 0; i < 3; i++ ) {
      if( lb1[ i ] == AST__BAD || ub1[ i ] == AST__BAD ) return 3;
   }

/* Axis 3 of sky_first is axis 1 of sky_last, and axes 1,2 are axes 2,3. */
   if( fabs( lb1[ 2 ] - lb2[ 0 ] ) > 1.0e-10 ||
       fabs( ub1[ 2 ] - ub2[ 0 ] ) > 1.0e-10 ) return 4;
   for( i = 0; i < 2; i++ ) {
      if( fabs( lb1[ i ] - lb2[ i + 1 ] ) > 1.0e-10 ||
          fabs( ub1[ i ] - ub2[ i + 1 ] ) > 1.0e-10 ) return 5;
   }
   return 0;
}

static int check_passthrough_output( const double out[ static 12 ],
                                     const double ra[ static 4 ],
                                     const double px[ static 4 ],
                                     const double dec[ static 4 ] ) {
   int i;
   for( i = 0; i < 4; i++ ) {
      if( out[ i ] != ra[ i ] ) return 7;
      if( out[ i + 4 ] != px[ i ] ) return 8;
      if( out[ i + 8 ] != dec[ i ] ) return 9;
   }
   return 0;
}

int main( void ) {
   int status_value = 0;
   int *status = &status_value;
   int perm[3];
   AstSkyFrame *sfrm;
   AstFrame *pfrm;
   AstCmpFrame *cfrm;

   double ra[]  = { 6.1, 6.1, 0.04, 0.04 };
   double dec[] = { 0.2, -0.2, -0.2, 0.2 };
   double px[]  = { -100.0, -10.0, 10.0, 100.0 };
   double in[4*3], out[4*3];

   astWatch( status );
   astBegin;

   sfrm = astSkyFrame( " " );
   perm[0] = 2;
   perm[1] = 1;
   astPermAxes( sfrm, perm );

   pfrm = astFrame( 1, "Domain=FPLANE" );

   cfrm = astCmpFrame( pfrm, sfrm, " " );
   perm[0] = 3;
   perm[1] = 1;
   perm[2] = 2;
   astPermAxes( cfrm, perm );

   /* Fill input array: columns are RA, PX, Dec (after permutation).
      Fortran layout in(4,3) with column-major ordering becomes
      a flat array in[npoint*ncoord] with npoint=4, ncoord=3.
      in[i + npoint*coord] = value for point i, coordinate coord. */
   fill_input( in, ra, px, dec );

   /* Normalise with `contig=1` (coordinates are contiguous per axis). */
   astNormPoints( cfrm, 4, 3, 4, in, 1, 3, 4, out );

   switch( check_normed_output( out, ra, px, dec ) ) {
   case 1: stopit( status, "Error 1" ); break;
   case 2: stopit( status, "Error 2" ); break;
   case 3: stopit( status, "Error 3" ); break;
   case 4: stopit( status, "Error 4" ); break;
   case 5: stopit( status, "Error 5" ); break;
   case 6: stopit( status, "Error 6" ); break;
   }

   /* Normalise with `contig=0` (do not normalise). */
   astNormPoints( cfrm, 4, 3, 4, in, 0, 3, 4, out );

   switch( check_passthrough_output( out, ra, px, dec ) ) {
   case 7: stopit( status, "Error 7" ); break;
   case 8: stopit( status, "Error 8" ); break;
   case 9: stopit( status, "Error 9" ); break;
   }

   switch( check_normbox_component_order() ) {
   case 1: stopit( status, "NormBox: error getting bounds with the SkyFrame first" ); break;
   case 2: stopit( status, "NormBox: error getting bounds with the SkyFrame last" ); break;
   case 3: stopit( status, "NormBox: bad bound with the SkyFrame first" ); break;
   case 4: stopit( status, "NormBox: the 1-D axis bounds depend on component order" ); break;
   case 5: stopit( status, "NormBox: the sky axis bounds depend on component order" ); break;
   }

   astEnd;
   astFlushMemory( 1 );

   if( *status == 0 ) {
      printf( " All CmpFrame tests passed\n" );
   } else {
      printf( "CmpFrame tests failed\n" );
   }
   return *status;
}
