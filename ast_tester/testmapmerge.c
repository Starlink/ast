/*
*  Name:
*     testmapmerge

*  Purpose:
*     Check that simplifying a Mapping does not change what it does.

*  Description:
*     astSimplify replaces a compound Mapping by a simpler one that should
*     be equivalent to it.  Each class merges itself with its neighbours,
*     or swaps places with them, and each case here is a small Mapping for
*     which one of those steps used to give a Mapping that was not
*     equivalent, report an error, or read or write memory it should not
*     have.  The test simplifies every Mapping and requires no error, the
*     same numbers of inputs and outputs, and the same results from the
*     forward and inverse transformations.  It is intended to be run with
*     the address and undefined-behavior sanitizers enabled, which is how
*     the memory faults are detected.  The cases cover:
*
*       - PermMaps: a PermMap created without one of its permutation
*         arrays, on its own and beside another PermMap, and an inverted
*         PermMap with more inputs than outputs.
*       - A PermMap followed by two Mappings in parallel, where the PermMap
*         supplies constants to both of them, or its inverse copies one of
*         its outputs to an input that is otherwise unused.
*       - WcsMaps: two projections that differ only in the number of
*         projection parameters, and an inverted PermMap beside a WcsMap.
*       - A PolyMap beside its own inverse, which must be left with the
*         Invert value its owner gave it.
*       - A TimeMap in which a conversion is followed by its inverse.
*       - MatrixMaps: a full matrix that is really diagonal and is used
*         inverted, a matrix with no inverse that is used inverted, and
*         diagonal matrices that have different numbers of inputs and
*         outputs beside a WinMap.
*       - WinMaps: beside a unit MatrixMap, beside two Mappings in
*         parallel, with a scale of zero, and followed by a MatrixMap
*         that discards an axis the WinMap shifts.
*       - MathMaps: astEqual, which decides whether a Mapping and the
*         inverse of its neighbour cancel, for MathMaps whose functions
*         differ only in their last operation or are a single constant,
*         and a SwitchMap followed by the inverse of a SwitchMap with a
*         different route.

*  Notes:
*     - This is a native C test with no Fortran original.
*     - Public C API only.
*/

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "sae_par.h"
#include "ast.h"
#include "ast_err.h"

/* The largest number of coordinates, and the number of positions, used
   to compare two Mappings. */
#define MXCOORD 8
#define NPOINT 3

/* Fill "in" with NPOINT positions of "ncoord" coordinates, as astTranN
   expects them.  The values are small and positive so that they are sky
   positions that every projection used here can project, and projection
   plane positions that it can deproject. */
static void positions( int ncoord, double *in ) {
   int coord;
   int point;
   for ( coord = 0; coord < ncoord; coord++ ) {
      for ( point = 0; point < NPOINT; point++ ) {
         in[ coord*NPOINT + point ] = 0.11 + 0.07*coord + 0.05*point;
      }
   }
}

/* See if two transformed values agree: both bad, or equal to within a
   small fraction of their size. */
static int agree( double a, double b ) {
   if ( a == AST__BAD || b == AST__BAD ) return a == b;
   return fabs( a - b ) <= 1.0E-9*fmax( 1.0, fmax( fabs( a ), fabs( b ) ) );
}

/* Compare the forward ("forward" non-zero) or inverse transformations of
   "map" and "simp".  Nothing is compared if "map" does not have the
   transformation.  "in" holds the positions to transform, or is NULL to
   use those given by "positions".  Returns the number of failures. */
static int compare( const char *what, AstMapping *map, AstMapping *simp,
                    int forward, const double *in ) {
   const char *dir = forward ? "forward" : "inverse";
   const char *attr = forward ? "TranForward" : "TranInverse";
   double def[ MXCOORD*NPOINT ];
   double out1[ MXCOORD*NPOINT ];
   double out2[ MXCOORD*NPOINT ];
   int i;
   int nin = astGetI( map, forward ? "Nin" : "Nout" );
   int nout = astGetI( map, forward ? "Nout" : "Nin" );

   if ( !astGetL( map, attr ) ) return 0;
   if ( !astGetL( simp, attr ) ) {
      fprintf( stderr, "testmapmerge: %s: the %s transformation was lost\n",
               what, dir );
      return 1;
   }
   if ( !in ) {
      positions( nin, def );
      in = def;
   }
   astTranN( map, NPOINT, nin, NPOINT, in, forward, nout, NPOINT, out1 );
   astTranN( simp, NPOINT, nin, NPOINT, in, forward, nout, NPOINT, out2 );
   if ( !astOK ) {
      fprintf( stderr, "testmapmerge: %s: the %s transformation gave "
               "status %d\n", what, dir, astStatus );
      astClearStatus;
      return 1;
   }
   for ( i = 0; i < nout*NPOINT; i++ ) {
      if ( !agree( out1[ i ], out2[ i ] ) ) {
         fprintf( stderr, "testmapmerge: %s: %s coordinate %d of position "
                  "%d was %.17g and is %.17g after simplifying\n", what, dir,
                  i/NPOINT + 1, i%NPOINT + 1, out1[ i ], out2[ i ] );
         return 1;
      }
   }
   return 0;
}

/* Simplify "map" and check that the result is equivalent to it.  "fwd"
   and "inv" hold the positions with which to compare the forward and
   inverse transformations, or are NULL to use those given by
   "positions".  If "simplified" is not NULL the simplified Mapping is
   returned there (NULL if there is none), and is otherwise annulled.
   "map" is annulled.  Returns the number of failures. */
static int check( const char *what, void *map, const double *fwd,
                  const double *inv, AstMapping **simplified ) {
   AstMapping *simp;
   int nfail = 0;

   if ( simplified ) *simplified = NULL;
   if ( !astOK || !map ) {
      fprintf( stderr, "testmapmerge: %s: the Mapping could not be made "
               "(status %d)\n", what, astStatus );
      if ( !astOK ) astClearStatus;
      return 1;
   }
   simp = astSimplify( map );
   if ( !astOK || !simp ) {
      fprintf( stderr, "testmapmerge: %s: astSimplify gave status %d\n",
               what, astStatus );
      if ( !astOK ) astClearStatus;
      return 1;
   }
   if ( astGetI( simp, "Nin" ) != astGetI( map, "Nin" ) ||
        astGetI( simp, "Nout" ) != astGetI( map, "Nout" ) ) {
      fprintf( stderr, "testmapmerge: %s: had %d inputs and %d outputs, "
               "and has %d and %d after simplifying\n", what,
               astGetI( map, "Nin" ), astGetI( map, "Nout" ),
               astGetI( simp, "Nin" ), astGetI( simp, "Nout" ) );
      nfail++;
   } else {
      nfail += compare( what, map, simp, 1, fwd );
      nfail += compare( what, map, simp, 0, inv );
   }
   if ( simplified ) {
      *simplified = simp;
   } else {
      simp = astAnnul( simp );
   }
   map = astAnnul( map );
   if ( !astOK ) astClearStatus;
   return nfail;
}

/* Return a Mapping with one input that no neighbour can absorb, and that
   has an inverse. */
static AstMapping *nonlinear( void ) {
   const char *fwd[] = { "y=exp(x)" };
   const char *inv[] = { "x=log(y)" };
   return (AstMapping *) astMathMap( 1, 1, 1, fwd, 1, inv, " " );
}

/* PermMaps. */
static int test_permmap( void ) {
   AstMapping *map;
   AstMapping *perm;
   double con[] = { 9.0 };
   double in = 5.0;
   double out[ 2 ];
   int constant[] = { -1 };
   int first2[] = { 1, 2 };
   int in3[] = { 1, 2, 0 };
   int nfail = 0;
   int swap[] = { 2, 1 };

/* A PermMap with no "outperm" array copies each input to the
   corresponding output, and gives a bad value to an output that has no
   corresponding input.  It used to give the first constant instead. */
   perm = (AstMapping *) astPermMap( 1, constant, 2, NULL, con, " " );
   astTranN( perm, 1, 1, 1, &in, 1, 2, 1, out );
   if ( !astOK || out[ 0 ] != 5.0 || out[ 1 ] != AST__BAD ) {
      fprintf( stderr, "testmapmerge: a PermMap with no outperm array "
               "gave (%g,%g) for an input of 5, not (5,<bad>)\n", out[ 0 ],
               out[ 1 ] );
      if ( !astOK ) astClearStatus;
      nfail++;
   }

/* Merging it with a second PermMap must give the same. */
   map = (AstMapping *) astCmpMap( perm, astPermMap( 2, swap, 2, swap, NULL,
                                                     " " ), 1, " " );
   nfail += check( "PermMap with no outperm array, then a PermMap", map,
                   NULL, NULL, NULL );

/* A PermMap with no "inperm" array, simplified on its own.  The new
   "inperm" array was compared with the missing one. */
   nfail += check( "PermMap with no inperm array",
                   astPermMap( 2, NULL, 1, constant, con, " " ), NULL, NULL,
                   NULL );

/* An inverted PermMap with more inputs than outputs, in an inverted
   CmpMap.  The "outperm" array was read as if it had an element for each
   input. */
   perm = (AstMapping *) astPermMap( 3, in3, 2, first2, NULL, " " );
   astInvert( perm );
   map = (AstMapping *) astCmpMap( perm, astZoomMap( 3, 2.0, " " ), 1, " " );
   astInvert( map );
   nfail += check( "inverted CmpMap of an inverted PermMap and a ZoomMap",
                   map, NULL, NULL, NULL );
   return nfail;
}

/* A PermMap followed by two Mappings in parallel. */
static int test_permmap_parallel( void ) {
   AstMapping *map;
   AstMapping *par;
   double con3[] = { 0.1, 0.2, 5.0, 0.3, 0.4, 7.0 };
   double con4[] = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8 };
   double con1[] = { 5.0 };
   int bad3[] = { 0, 0, 0 };
   int in3[] = { -1, -2, -3 };
   int in4[] = { -1, -2, -3, -4 };
   int inlink[] = { 2, 1 };
   int nfail = 0;
   int out3[] = { -4, -5, -6 };
   int out4[] = { -5, -6, -7, -8 };
   int outlink[] = { -1, 1 };

/* Every PermMap output is a constant.  The PermMap that replaced it was
   made from arrays that were too small. */
   par = (AstMapping *) astCmpMap( astWcsMap( 2, AST__TAN, 1, 2, " " ),
                                   astZoomMap( 1, 3.0, " " ), 0, " " );
   map = (AstMapping *) astCmpMap( astPermMap( 3, in3, 3, out3, con3, " " ),
                                   par, 1, " " );
   nfail += check( "constant PermMap, then a WcsMap and ZoomMap in parallel",
                   map, NULL, NULL, NULL );
   map = (AstMapping *) astCmpMap( astPermMap( 3, bad3, 3, out3, con3, " " ),
                                   par, 1, " " );
   nfail += check( "constant PermMap with bad inputs, then a WcsMap and "
                   "ZoomMap in parallel", map, NULL, NULL, NULL );

/* The same with the parallel Mappings the other way round, where parts
   of the new permutation array were never set. */
   par = (AstMapping *) astCmpMap( astZoomMap( 1, 3.0, " " ),
                                   astWcsMap( 2, AST__TAN, 1, 2, " " ), 0,
                                   " " );
   map = (AstMapping *) astCmpMap( astPermMap( 3, in3, 3, out3, con3, " " ),
                                   par, 1, " " );
   nfail += check( "constant PermMap, then a ZoomMap and WcsMap in parallel",
                   map, NULL, NULL, NULL );

/* Two projections in parallel, where the inverse transformation was
   changed. */
   par = (AstMapping *) astCmpMap( astWcsMap( 2, AST__TAN, 1, 2, " " ),
                                   astWcsMap( 2, AST__SIN, 1, 2, " " ), 0,
                                   " " );
   map = (AstMapping *) astCmpMap( astPermMap( 4, in4, 4, out4, con4, " " ),
                                   par, 1, " " );
   nfail += check( "constant PermMap, then two WcsMaps in parallel", map,
                   NULL, NULL, NULL );

/* The first PermMap output is a constant and the second is the first
   input.  The second input is not used by the forward transformation,
   and the inverse transformation copies the constant output to it.  This
   was reported as an internal programming error. */
   par = (AstMapping *) astCmpMap( nonlinear(), nonlinear(), 0, " " );
   map = (AstMapping *) astCmpMap( astPermMap( 2, inlink, 2, outlink, con1,
                                               " " ), par, 1, " " );
   nfail += check( "PermMap whose inverse copies a constant output, then "
                   "two MathMaps in parallel", map, NULL, NULL, NULL );
   return nfail;
}

/* WcsMaps. */
static int test_wcsmap( void ) {
   AstMapping *map;
   AstMapping *simp;
   AstMapping *wcs1;
   AstMapping *wcs2;
   double con[] = { 5.0 };
   double plane[ 2*NPOINT ] = { 0.1, 0.2, 0.3, 1.2, 1.1, 1.0 };
   int in2[] = { 1, 2 };
   int nfail = 0;
   int out3[] = { 1, 2, -1 };

/* Two zenithal perspective projections, the second of which has a
   projection parameter that the first does not, do not cancel.  The
   forward transformation of a WcsMap goes from the projection plane to
   the sky, so this goes from one projection plane to the other. */
   wcs1 = (AstMapping *) astWcsMap( 2, AST__AZP, 1, 2, " " );
   wcs2 = (AstMapping *) astWcsMap( 2, AST__AZP, 1, 2, "PV2_1=0.5" );
   astInvert( wcs2 );
   map = (AstMapping *) astCmpMap( wcs1, wcs2, 1, " " );
   nfail += check( "AZP projection, then the inverse of a different AZP "
                   "projection", map, plane, NULL, &simp );
   if ( simp ) {
      if ( astIsAUnitMap( simp ) ) {
         fprintf( stderr, "testmapmerge: two different AZP projections "
                  "cancelled each other\n" );
         nfail++;
      }
      simp = astAnnul( simp );
   }

/* The same projection does cancel. */
   wcs1 = (AstMapping *) astWcsMap( 2, AST__AZP, 1, 2, "PV2_1=0.5" );
   wcs2 = (AstMapping *) astWcsMap( 2, AST__AZP, 1, 2, "PV2_1=0.5" );
   astInvert( wcs1 );
   map = (AstMapping *) astCmpMap( wcs1, wcs2, 1, " " );
   nfail += check( "inverse AZP projection, then the same projection", map,
                   NULL, NULL, &simp );
   if ( simp ) {
      if ( !astIsAUnitMap( simp ) ) {
         fprintf( stderr, "testmapmerge: an AZP projection did not cancel "
                  "its own inverse\n" );
         nfail++;
      }
      simp = astAnnul( simp );
   }

/* A PermMap that adds an axis, then a WcsMap, both used in the inverse
   direction.  The PermMap came out of the swap the wrong way round. */
   map = (AstMapping *) astCmpMap( astPermMap( 2, in2, 3, out3, con, " " ),
                                   astWcsMap( 3, AST__TAN, 1, 2, " " ), 1,
                                   " " );
   nfail += check( "PermMap then a WcsMap", astCopy( map ), NULL, NULL,
                   NULL );
   astInvert( map );
   nfail += check( "inverted CmpMap of a PermMap and a WcsMap", map, NULL,
                   NULL, NULL );
   return nfail;
}

/* A PolyMap and its inverse. */
static int test_polymap( void ) {
   AstMapping *map;
   AstMapping *simp;
   AstPolyMap *p;
   AstPolyMap *q;
   double fwd[] = { 1.0, 1.0, 2.0,  3.0, 1.0, 1.0 };
   double inv[] = { 0.5, 1.0, 1.0 };
   int nfail = 0;

/* A PolyMap is cancelled by a copy of itself that is used in the other
   direction (whether or not its two transformations are really inverses
   of each other, as they are not here).  The PolyMaps belong to the
   caller and must be left as they were.  The second was left with the
   Invert value it was being used with. */
   p = astPolyMap( 1, 1, 2, fwd, 1, inv, " " );
   q = astCopy( p );
   astInvert( q );
   map = (AstMapping *) astCmpMap( p, q, 1, " " );
   astInvert( map );
   simp = astSimplify( map );
   if ( !astOK || !simp || !astIsAUnitMap( simp ) ) {
      fprintf( stderr, "testmapmerge: a PolyMap did not cancel its own "
               "inverse (status %d)\n", astStatus );
      nfail++;
   } else if ( astGetI( p, "Invert" ) != 0 || astGetI( q, "Invert" ) != 1 ) {
      fprintf( stderr, "testmapmerge: simplifying a CmpMap changed the "
               "Invert attributes of its PolyMaps to %d and %d\n",
               astGetI( p, "Invert" ), astGetI( q, "Invert" ) );
      nfail++;
   }
   if ( !astOK ) astClearStatus;
   return nfail;
}

/* A TimeMap. */
static int test_timemap( void ) {
   AstMapping *simp;
   AstTimeMap *map;
   double args[ 3 ] = { 2000.0, 0.0, 0.0 };
   double mjd[ NPOINT ] = { 55000.0, 55000.25, 56000.5 };
   int nfail = 0;

/* A conversion from TDB to TCB and back again has no effect. */
   map = astTimeMap( 0, " " );
   astTimeAdd( map, "TDBTOTCB", 1, args );
   astTimeAdd( map, "TCBTOTDB", 1, args );
   nfail += check( "TimeMap from TDB to TCB and back", map, mjd, mjd, &simp );
   if ( simp ) {
      if ( !astIsAUnitMap( simp ) ) {
         fprintf( stderr, "testmapmerge: a TimeMap from TDB to TCB and "
                  "back was not simplified to a UnitMap\n" );
         nfail++;
      }
      simp = astAnnul( simp );
   }
   return nfail;
}

/* MatrixMaps. */
static int test_matrixmap( void ) {
   AstMapping *map;
   AstMapping *mat;
   const char *full = " Begin MatrixMap\n Nin = 2\n Invert = 1\n"
                      " IsA Mapping\n M0 = 2\n M1 = 0\n M2 = 0\n M3 = 4\n"
                      " IM0 = 0.5\n IM1 = 0\n IM2 = 0\n IM3 = 0.25\n"
                      " Form = \"Full\"\n End MatrixMap\n";
   double diag2[] = { 2.0, 3.0 };
   double general[] = { 1.0, 2.0, 3.0, 4.0 };
   double singular[] = { 1.0, 2.0, 2.0, 4.0 };
   double zero[] = { 0.0, 0.0 };
   int nfail = 0;

/* A full matrix with zeros off the diagonal, used inverted.  It was
   replaced by a diagonal matrix that was not inverted. */
   nfail += check( "inverted full MatrixMap that is diagonal",
                   astFromString( full ), NULL, NULL, NULL );

/* A matrix with no inverse, used inverted, has no forward transformation.
   Merging it with a neighbour reported an error. */
   mat = (AstMapping *) astMatrixMap( 2, 2, 0, singular, " " );
   astInvert( mat );
   map = (AstMapping *) astCmpMap( mat, astMatrixMap( 2, 2, 0, general, " " ),
                                   1, " " );
   nfail += check( "inverted singular MatrixMap, then a MatrixMap", map,
                   NULL, NULL, NULL );
   map = (AstMapping *) astCmpMap( astZoomMap( 2, 2.0, " " ), mat, 1, " " );
   nfail += check( "ZoomMap, then an inverted singular MatrixMap", map, NULL,
                   NULL, NULL );
   map = (AstMapping *) astCmpMap( astShiftMap( 2, diag2, " " ), mat, 1,
                                   " " );
   nfail += check( "ShiftMap, then an inverted singular MatrixMap", map,
                   NULL, NULL, NULL );

/* The inverse of a diagonal matrix of zeros has no usable elements. */
   mat = (AstMapping *) astMatrixMap( 2, 2, 1, zero, " " );
   astInvert( mat );
   map = (AstMapping *) astCmpMap( mat, astMatrixMap( 2, 2, 0, general, " " ),
                                   1, " " );
   nfail += check( "inverted MatrixMap of zeros, then a MatrixMap", map,
                   NULL, NULL, NULL );
   map = (AstMapping *) astCmpMap( mat, astShiftMap( 2, diag2, " " ), 1,
                                   " " );
   nfail += check( "inverted MatrixMap of zeros, then a ShiftMap", map, NULL,
                   NULL, NULL );
   return nfail;
}

/* WinMaps beside MatrixMaps. */
static int test_winmap( void ) {
   AstMapping *map;
   AstMapping *pair;
   AstMapping *win2;
   AstMapping *win3;
   double diag2[] = { 2.0, 3.0 };
   double ina[] = { 0.0, 0.0, 0.0 };
   double inb[] = { 1.0, 2.0, 3.0 };
   double outa[] = { 5.0, 6.0, 7.0 };
   double outb[] = { 6.0, 9.0, 12.0 };
   double same[] = { 0.0 };
   double scaled[] = { 4.0, 6.0 };
   double shift1[] = { 3.0 };
   double shift2[] = { 3.0, 4.0 };
   double shift3[] = { 1.0, 2.0, 3.0 };
   double unit[] = { 1.0, 1.0 };
   int nfail = 0;

   win2 = (AstMapping *) astWinMap( 2, ina, inb, outa, outb, " " );
   win3 = (AstMapping *) astWinMap( 3, ina, inb, outa, outb, " " );

/* A unit MatrixMap has no elements to merge with a WinMap, and the
   WinMap was returned as it stood: with too few axes in parallel, and
   ignoring the direction in which it was to be used in series. */
   map = (AstMapping *) astCmpMap( astShiftMap( 1, shift1, " " ),
                                   astMatrixMap( 2, 2, 2, NULL, " " ), 0,
                                   " " );
   nfail += check( "ShiftMap and a unit MatrixMap in parallel", map, NULL,
                   NULL, NULL );
   map = (AstMapping *) astCmpMap( astMatrixMap( 2, 2, 2, NULL, " " ), win2,
                                   1, " " );
   astInvert( map );
   nfail += check( "inverted CmpMap of a unit MatrixMap and a WinMap", map,
                   NULL, NULL, NULL );

/* A diagonal MatrixMap with different numbers of inputs and outputs
   cannot be merged with a WinMap, which has the same number of each. */
   map = (AstMapping *) astCmpMap( win3, astMatrixMap( 3, 2, 1, diag2, " " ),
                                   1, " " );
   nfail += check( "WinMap, then a diagonal MatrixMap with fewer outputs",
                   map, NULL, NULL, NULL );
   map = (AstMapping *) astCmpMap( astMatrixMap( 2, 3, 1, diag2, " " ), win3,
                                   1, " " );
   nfail += check( "diagonal MatrixMap with more outputs, then a WinMap",
                   map, NULL, NULL, NULL );
   map = (AstMapping *) astCmpMap( astMatrixMap( 3, 2, 1, diag2, " " ), win2,
                                   1, " " );
   nfail += check( "diagonal MatrixMap with fewer outputs, then a WinMap",
                   map, NULL, NULL, NULL );
   map = (AstMapping *) astCmpMap( win2, astMatrixMap( 2, 3, 1, diag2, " " ),
                                   1, " " );
   nfail += check( "WinMap, then a diagonal MatrixMap with more outputs",
                   map, NULL, NULL, NULL );

/* A WinMap that shifts an axis which the following MatrixMap discards.
   The inverse of the MatrixMap gives zero for that axis, to which the
   inverse of the WinMap is then applied.  Moving the WinMap after the
   MatrixMap, to bring it next to the ShiftMap, lost that. */
   pair = (AstMapping *) astCmpMap( astShiftMap( 3, shift3, " " ),
                                    astMatrixMap( 3, 2, 1, diag2, " " ), 1,
                                    " " );
   map = (AstMapping *) astCmpMap( pair, astShiftMap( 2, shift2, " " ), 1,
                                   " " );
   nfail += check( "ShiftMap, a MatrixMap that discards an axis, and a "
                   "ShiftMap", map, NULL, NULL, NULL );

/* A WinMap beside two Mappings in parallel that it cannot merge with.
   The WinMap is split into one WinMap for each of them, which must not
   change the number of axes that the WinMap itself is taken to have. */
   pair = (AstMapping *) astCmpMap( nonlinear(), nonlinear(), 0, " " );
   map = (AstMapping *) astCmpMap( astWinMap( 2, ina, unit, shift2, scaled,
                                              " " ), pair, 1, " " );
   nfail += check( "WinMap, then two MathMaps in parallel", map, NULL, NULL,
                   NULL );
   map = (AstMapping *) astCmpMap( pair, astWinMap( 2, ina, unit, shift2,
                                                    scaled, " " ), 1, " " );
   nfail += check( "two MathMaps in parallel, then a WinMap", map, NULL, NULL,
                   NULL );

/* A WinMap that maps every input to zero has no inverse scale, so no
   diagonal MatrixMap is equivalent to its inverse. */
   map = (AstMapping *) astWinMap( 1, ina, inb, same, same, " " );
   astInvert( map );
   nfail += check( "inverted WinMap with a scale of zero", map, NULL, NULL,
                   NULL );
   return nfail;
}

/* Return a MathMap with one input and one output that has the forward
   function "fwd" and the inverse function "inv". */
static AstMapping *mathmap( const char *fwd, const char *inv ) {
   const char *fwdfun[ 1 ];
   const char *invfun[ 1 ];
   fwdfun[ 0 ] = fwd;
   invfun[ 0 ] = inv;
   return (AstMapping *) astMathMap( 1, 1, 1, fwdfun, 1, invfun, " " );
}

/* Return a SwitchMap with one input and one output that always uses
   "route", its only route Mapping. */
static AstMapping *switchmap( AstMapping *route ) {
   const char *pick[] = { "r=1" };
   const char *none[] = { "a" };
   AstMapping *routes[ 1 ];
   routes[ 0 ] = route;
   return (AstMapping *) astSwitchMap( astMathMap( 1, 1, 1, pick, 1, none,
                                                   " " ),
                                       astMathMap( 1, 1, 1, none, 1, pick,
                                                   " " ),
                                       1, (void **) routes, " " );
}

/* Check that astEqual gives "expected" for the Mappings "a" and "b", in
   either order.  Both are annulled.  Returns the number of failures. */
static int check_equal( const char *what, AstMapping *a, AstMapping *b,
                        int expected ) {
   int nfail = 0;
   int ab;
   int ba;

   if ( !astOK || !a || !b ) {
      fprintf( stderr, "testmapmerge: %s: the Mappings could not be made "
               "(status %d)\n", what, astStatus );
      if ( !astOK ) astClearStatus;
      return 1;
   }
   ab = astEqual( a, b );
   ba = astEqual( b, a );
   if ( !astOK ) {
      fprintf( stderr, "testmapmerge: %s: astEqual gave status %d\n", what,
               astStatus );
      astClearStatus;
      nfail++;
   } else if ( ab != expected || ba != expected ) {
      fprintf( stderr, "testmapmerge: %s: astEqual gave %d and %d but %d "
               "was expected\n", what, ab, ba, expected );
      nfail++;
   }
   a = astAnnul( a );
   b = astAnnul( b );
   return nfail;
}

/* MathMaps that are nearly the same. */
static int test_mathmap( void ) {
   AstMapping *first;
   AstMapping *second;
   int nfail = 0;

/* Two MathMaps are equal if their functions compile to the same operations
   and constants.  The last operation of each function used not to be
   compared, so MathMaps that ended in different operations were equal. */
   nfail += check_equal( "the same MathMap twice",
                         mathmap( "y=x+1", "x=y-1" ),
                         mathmap( "y=x+1", "x=y-1" ), 1 );
   nfail += check_equal( "MathMaps that add and subtract",
                         mathmap( "y=x+1", "x=y-1" ),
                         mathmap( "y=x-1", "x=y+1" ), 0 );
   nfail += check_equal( "MathMaps that multiply and divide",
                         mathmap( "y=x*2", "x=y/2" ),
                         mathmap( "y=x/2", "x=y*2" ), 0 );
   nfail += check_equal( "MathMaps that take a sine and a cosine",
                         mathmap( "y=sin(x)", "x" ),
                         mathmap( "y=cos(x)", "x" ), 0 );
   nfail += check_equal( "MathMaps that add different constants",
                         mathmap( "y=x+1", "x=y-1" ),
                         mathmap( "y=x+2", "x=y-2" ), 0 );
   nfail += check_equal( "MathMaps that take a maximum and a minimum",
                         mathmap( "y=max(x,1,2)", "x" ),
                         mathmap( "y=min(x,1,2)", "x" ), 0 );
   nfail += check_equal( "MathMaps that take different maxima",
                         mathmap( "y=max(x,1,2)", "x" ),
                         mathmap( "y=max(x,1,3)", "x" ), 0 );
   nfail += check_equal( "MathMaps that take the same maximum",
                         mathmap( "y=max(x,1,2)", "x" ),
                         mathmap( "y=max(x,1,2)", "x" ), 1 );

/* An inverted MathMap is compared through the functions it then uses in
   each direction. */
   first = mathmap( "y=x+1", "x=y-1" );
   second = mathmap( "y=x-1", "x=y+1" );
   astInvert( second );
   nfail += check_equal( "a MathMap and the inverse of its opposite", first,
                         second, 1 );

/* A function that is a single constant has one operation and no list of
   constants.  A list was read for it all the same. */
   nfail += check_equal( "MathMaps that give the same constant",
                         mathmap( "y=<pi>", "x" ), mathmap( "y=<pi>", "x" ),
                         1 );
   nfail += check_equal( "MathMaps that give different constants",
                         mathmap( "y=<pi>", "x" ), mathmap( "y=<e>", "x" ),
                         0 );
   nfail += check_equal( "a MathMap that copies its input and one that "
                         "gives a constant", mathmap( "y=x", "x=y" ),
                         mathmap( "y=<pi>", "x=y" ), 0 );

/* A SwitchMap followed by the inverse of an equal SwitchMap is removed
   when the pair is simplified.  SwitchMaps whose routes are MathMaps that
   differ only in their last operations were removed in the same way. */
   first = switchmap( mathmap( "y=x+1", "x=y-1" ) );
   second = switchmap( mathmap( "y=x-1", "x=y+1" ) );
   astInvert( second );
   nfail += check( "SwitchMap, then the inverse of a different SwitchMap",
                   astCmpMap( first, second, 1, " " ), NULL, NULL, NULL );
   return nfail;
}

int main( void ) {
   int status = SAI__OK;
   int nfail = 0;
   astWatch( &status );
   astBegin;

   nfail += test_permmap();
   nfail += test_permmap_parallel();
   nfail += test_wcsmap();
   nfail += test_polymap();
   nfail += test_timemap();
   nfail += test_matrixmap();
   nfail += test_winmap();
   nfail += test_mathmap();

   astEnd;
   if ( !astOK ) astClearStatus;

   if ( nfail == 0 ) {
      printf( " All Mapping merger tests passed\n" );
      return 0;
   }
   printf( "Mapping merger tests failed (%d)\n", nfail );
   return 1;
}
