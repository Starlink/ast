/*
*  Purpose:
*     Exercise the XphMap loader: the projection-type names it accepts, the
*     components that survive a dump/load round trip, and the transformation
*     of a loaded XphMap.
*
*     Every HEALPix pixel on the sky is one used cell in the grid of each
*     projection, so of the pixel centres in the grid of a projection, the
*     number that the forward transformation of an XphMap transforms must be
*     the number of pixels on the sky, and each must come back through the
*     inverse.  An HPX grid also holds an unused second copy of each pixel of
*     the facet that is split between its two ends, and the inverse
*     transformation, which takes positions in an HPX12 grid, must treat the
*     unused copy as it treats the used one.
*
*     A position is on the sky if it is within a used cell or on its boundary,
*     and the transformation within a cell is a shift and a rotation.  So a
*     position on the boundary of a used cell must be transformed to the
*     position that continues the transformation of that cell, or of another
*     used cell that it bounds, and a position that bounds no used cell must
*     give bad values.  That includes the poles: a pole is the centre or the
*     four corners of an XPH grid but is four separate corners of polar
*     facets in an HPX grid, and all of them must transform.  A position that
*     rounding has put marginally off the edge of the sky must be treated as
*     being on the edge.
*
*     Positions too far from the grid for the index of a grid cell to be held
*     in an integer, and positions that are not finite, must give bad values.
*     So must every position outside the grid of the projection being
*     transformed, however near.  An Order that the constructor would not
*     accept must be rejected by the loader, and the largest Order must work.
*
*  Notes:
*     - Public C API only, with astWatch/astOK for error checking.  XphMap has
*       no public constructor, so every XphMap here is read from a dump; that
*       is also the only route to astLoadXphMap, which is what this test is
*       about.
*     - Order and Type are dump components rather than readable attributes, so
*       they are checked by re-dumping the loaded object.
*     - astLoadXphMap used to leak the string it read for the Type component.
*       No output value can reveal a leak, so nothing here asserts it is gone;
*       the check is an AddressSanitizer run of this test with detect_leaks=1,
*       which is not wired into either build system because LeakSanitizer does
*       not exist on Darwin and every other test in this directory suppresses
*       leak reports (AST's memory system caches freed blocks).
*     - The transformation used to convert the index of a grid cell to an
*       integer whether or not an integer could hold it, which is undefined.
*       On most machines the far positions here gave bad values regardless,
*       so most of that check can fail only under UndefinedBehaviorSanitizer.
*/
#include "ast.h"
#include "ast_err.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

static int status = 0;

static void fail( const char *text ) {
   if( status == 0 ) status = 1;
   printf( "%s\n", text );
}

/* Build a dump for an XphMap of the given order and projection type. */
static void make_dump( char *buf, size_t len, int order, const char *type ) {
   snprintf( buf, len,
             " Begin XphMap\n"
             "    Nin = 2\n"
             " IsA Mapping\n"
             "    Order = %d\n"
             "    Type = \"%s\"\n"
             " End XphMap\n", order, type );
}

/* Load an XphMap of the named type, then check its class, its arity, that a
   re-dump still carries the order and type, and that two independent loads of
   the same text compare equal. */
static void check_type( const char *type, int order ) {
   AstMapping *map, *again;
   char dump[ 256 ];
   char *text;

   if( status != 0 ) return;
   astBegin;

   make_dump( dump, sizeof dump, order, type );
   map = astFromString( dump );
   if( !astOK || !map ) {
      printf( "testxphmap: could not load an XphMap of type %s\n", type );
      fail( "testxphmap: astFromString failed" );
      astEnd;
      return;
   }

   if( strcmp( astGetC( map, "Class" ), "XphMap" ) ) {
      fail( "testxphmap: loaded object is not an XphMap" );
   } else if( astGetI( map, "Nin" ) != 2 || astGetI( map, "Nout" ) != 2 ) {
      fail( "testxphmap: loaded XphMap is not 2-in, 2-out" );
   }

/* Re-dump and look for the two components back again.  They are written as
   class data rather than exposed as attributes, so the text is the only
   place to read them. */
   text = astToString( map );
   if( !astOK || !text ) {
      fail( "testxphmap: astToString failed" );
   } else {
      char want[ 64 ];
      snprintf( want, sizeof want, "Type = \"%s\"", type );
      if( !strstr( text, want ) ) {
         printf( "testxphmap: re-dump of type %s lost its Type card\n", type );
         fail( "testxphmap: Type did not survive the round trip" );
      }
      snprintf( want, sizeof want, "Order = %d", order );
      if( !strstr( text, want ) ) {
         printf( "testxphmap: re-dump of type %s lost its Order card\n", type );
         fail( "testxphmap: Order did not survive the round trip" );
      }
      text = astFree( text );
   }

   again = astFromString( dump );
   if( astOK && !astEqual( again, map ) ) {
      printf( "testxphmap: two loads of the same type %s dump differ\n", type );
      fail( "testxphmap: astEqual failed on identical XphMaps" );
   }

   astEnd;
}

/* The forward transformation of a loaded XphMap must be invertible wherever
   it is defined.  The inputs are grid coordinates within the projection, and
   only a diamond-shaped part of the bounding box is on the sphere, so the
   test scans a grid, keeps the points the forward transformation accepts, and
   requires every one of those to come back through the inverse. */
static void check_transform( void ) {
#define XPH_ORDER 8
#define XPH_STEP  37
   AstMapping *map;
   char dump[ 256 ];
   double xin[ 1024 ], yin[ 1024 ];
   double xm[ 1024 ], ym[ 1024 ], xout[ 1024 ], yout[ 1024 ];
   double lim;
   double x, y;
   int i, n, ngood;

   if( status != 0 ) return;
   astBegin;

   make_dump( dump, sizeof dump, XPH_ORDER, "HPX12" );
   map = astFromString( dump );
   if( !astOK || !map ) {
      fail( "testxphmap: could not load the XphMap to transform" );
      astEnd;
      return;
   }

/* An order N HPX grid has 2^N pixels along a facet edge and spans about five
   facets, so this covers the whole projection with room to spare. */
   lim = 5.0 * ( 1 << XPH_ORDER );

   n = 0;
   for( x = 0.0; x <= lim && n < 1024; x += XPH_STEP ) {
      for( y = 0.0; y <= lim && n < 1024; y += XPH_STEP ) {
         xin[ n ] = x;
         yin[ n ] = y;
         n++;
      }
   }

   astTran2( map, n, xin, yin, 1, xm, ym );
   astTran2( map, n, xm, ym, 0, xout, yout );
   if( !astOK ) {
      fail( "testxphmap: astTran2 failed" );
      astEnd;
      return;
   }

   ngood = 0;
   for( i = 0; i < n; i++ ) {
      if( xm[ i ] == AST__BAD || ym[ i ] == AST__BAD ) continue;
      ngood++;
      if( xout[ i ] == AST__BAD || yout[ i ] == AST__BAD ) {
         printf( "testxphmap: (%g,%g) transformed but did not come back\n",
                 xin[ i ], yin[ i ] );
         fail( "testxphmap: inverse returned a bad value for a good point" );
         break;
      }
      if( fabs( xout[ i ] - xin[ i ] ) > 1.0E-8 ||
          fabs( yout[ i ] - yin[ i ] ) > 1.0E-8 ) {
         printf( "testxphmap: (%g,%g) round tripped to (%g,%g)\n",
                 xin[ i ], yin[ i ], xout[ i ], yout[ i ] );
         fail( "testxphmap: round trip did not recover the input" );
         break;
      }
   }

/* Guard against the scan silently finding nothing, which would make the loop
   above vacuous. */
   if( status == 0 && ngood < 20 ) {
      printf( "testxphmap: only %d of %d scanned points were on the sphere\n",
              ngood, n );
      fail( "testxphmap: the scan found too few good points to be a test" );
   }

   astEnd;
#undef XPH_ORDER
#undef XPH_STEP
}

/* A dump whose Type is not a known projection name must be rejected.  This is
   the error path of the loader, where the type string is also read. */
static void check_bad_type( void ) {
   AstObject *obj;
   char dump[ 256 ];

   if( status != 0 ) return;
   astBegin;

   make_dump( dump, sizeof dump, 8, "NOSUCHPROJECTION" );
   obj = astFromString( dump );
   if( astOK ) {
      fail( "testxphmap: an unknown Type was accepted" );
   } else {
      astClearStatus;
   }
   (void) obj;

   astEnd;
}

/* Count the pixel centres in the grid of a projection that an XphMap of the
   given type and order transforms in the given direction, and send each one
   back.  The forward transformation takes positions in the grid of the named
   projection, which is four facets wide for an XPH projection and five for
   an HPX projection.  The inverse takes positions in an HPX grid. */
static void check_coverage( const char *type, int order, int forward ) {
#define XPH_MAXSCAN 1600
   static double xin[ XPH_MAXSCAN ], yin[ XPH_MAXSCAN ];
   static double xm[ XPH_MAXSCAN ], ym[ XPH_MAXSCAN ];
   static double xout[ XPH_MAXSCAN ], yout[ XPH_MAXSCAN ];
   AstMapping *map;
   char dump[ 256 ];
   const char *dir = forward ? "forward" : "inverse";
   int nppf = 1 << order;
   int ncell = ( forward && type[ 0 ] == 'X' ) ? 4*nppf : 5*nppf;
   int i, ix, iy, n, ncopy, ngood;

   if( status != 0 ) return;
   astBegin;

   make_dump( dump, sizeof dump, order, type );
   map = astFromString( dump );

   n = 0;
   for( iy = 1; iy <= ncell && n < XPH_MAXSCAN; iy++ ) {
      for( ix = 1; ix <= ncell && n < XPH_MAXSCAN; ix++ ) {
         xin[ n ] = ix;
         yin[ n ] = iy;
         n++;
      }
   }

   astTran2( map, n, xin, yin, forward, xm, ym );
   astTran2( map, n, xm, ym, !forward, xout, yout );
   if( !astOK ) {
      printf( "testxphmap: could not scan type %s order %d\n", type, order );
      fail( "testxphmap: astTran2 failed in the coverage scan" );
      astEnd;
      return;
   }

/* A pixel centre must come back as itself.  The exception is the unused
   copy of a pixel in an HPX12 grid, which the inverse transformation accepts
   and which must come back as the used copy, four facets away on both
   axes. */
   ngood = 0;
   ncopy = 0;
   for( i = 0; i < n && status == 0; i++ ) {
      if( xm[ i ] == AST__BAD || ym[ i ] == AST__BAD ) continue;
      ngood++;
      if( xout[ i ] == xin[ i ] && yout[ i ] == yin[ i ] ) continue;
      if( !forward && fabs( xout[ i ] - xin[ i ] ) == 4*nppf &&
          yout[ i ] - yin[ i ] == xout[ i ] - xin[ i ] ) {
         ncopy++;
         continue;
      }
      printf( "testxphmap: type %s order %d %s: (%g,%g) came back as "
              "(%g,%g)\n", type, order, dir, xin[ i ], yin[ i ],
              xout[ i ], yout[ i ] );
      fail( "testxphmap: a pixel centre did not make the return trip" );
   }

   if( status == 0 && ngood - ncopy != 12*nppf*nppf ) {
      printf( "testxphmap: type %s order %d %s: %d pixel centres were "
              "transformed, but the sky has %d pixels\n", type, order, dir,
              ngood - ncopy, 12*nppf*nppf );
      fail( "testxphmap: wrong number of pixels in the projection" );
   }

/* The split facet has as many pixels as any other, and each has one unused
   copy. */
   if( status == 0 && ncopy != ( forward ? 0 : nppf*nppf ) ) {
      printf( "testxphmap: type %s order %d %s: %d unused copies of a pixel "
              "were transformed\n", type, order, dir, ncopy );
      fail( "testxphmap: wrong number of unused copies of a pixel" );
   }

   astEnd;
#undef XPH_MAXSCAN
}

/* Transform positions that no grid index held in an integer can describe,
   and positions that are not finite.  Each far value is tried on both axes,
   with every other far value and with a position inside the grid on the
   other axis, since the quadrant of an XPH grid that a position falls in
   selects the arithmetic used.  Every one must give bad values. */
static void check_far_positions( const char *type, int order ) {
   AstMapping *map;
   char dump[ 256 ];
   double far[ 8 ];
   double xin[ 80 ], yin[ 80 ], xout[ 80 ], yout[ 80 ];
   int forward, i, j, n;
   int nfar = 0;

   if( status != 0 ) return;
   astBegin;

   far[ nfar++ ] = 1.0E10;
   far[ nfar++ ] = -1.0E10;
   far[ nfar++ ] = 1.0E300;
   far[ nfar++ ] = -1.0E300;
   far[ nfar++ ] = HUGE_VAL;
   far[ nfar++ ] = -HUGE_VAL;
   far[ nfar++ ] = nan( "" );

   n = 0;
   for( j = 0; j < nfar; j++ ) {
      for( i = 0; i < nfar; i++ ) {
         xin[ n ] = far[ i ];
         yin[ n ] = far[ j ];
         n++;
      }
      xin[ n ] = far[ j ];
      yin[ n ] = 3.0;
      n++;
      xin[ n ] = 3.0;
      yin[ n ] = far[ j ];
      n++;
   }

   make_dump( dump, sizeof dump, order, type );
   map = astFromString( dump );

   for( forward = 1; forward >= 0 && status == 0; forward-- ) {
      astTran2( map, n, xin, yin, forward, xout, yout );
      if( !astOK ) {
         fail( "testxphmap: astTran2 failed for far positions" );
         break;
      }
      for( i = 0; i < n; i++ ) {
         if( xout[ i ] != AST__BAD || yout[ i ] != AST__BAD ) {
            printf( "testxphmap: type %s order %d %s: (%g,%g) gave (%g,%g)\n",
                    type, order, forward ? "forward" : "inverse", xin[ i ],
                    yin[ i ], xout[ i ], yout[ i ] );
            fail( "testxphmap: wrong result for a far position" );
            break;
         }
      }
   }

   astEnd;
}

/* The grid of an XPH projection is four facets wide and the grid of an HPX
   projection is five.  The forward transformation takes positions in the
   grid of the named projection and the inverse takes positions in an HPX12
   grid.  Each gore of an XPH grid runs diagonally from the centre of the
   grid to a corner, where it ends.  Scan every quarter of a pixel over the
   grid and a margin wider than the grid on every side, so that the scan
   crosses the line of each gore beyond its corner.  The transformation must
   give bad values for every position outside the grid.  The edges of the
   grid are part of it. */
static void check_outside_grid( const char *type, int order, int forward ) {
#define XPH_MAXROW 800
   static double xin[ XPH_MAXROW ], yin[ XPH_MAXROW ];
   static double xout[ XPH_MAXROW ], yout[ XPH_MAXROW ];
   AstMapping *map;
   char dump[ 256 ];
   const char *dir = forward ? "forward" : "inverse";
   int nfacet = ( forward && type[ 0 ] == 'X' ) ? 4 : 5;
   double hi = nfacet*( 1 << order ) + 0.5;
   double lo = 0.5;
   double margin = 5.0*( 1 << order );
   int i, j, n, ninside, noutside;

   if( status != 0 ) return;
   astBegin;

   make_dump( dump, sizeof dump, order, type );
   map = astFromString( dump );

   n = (int)( ( hi - lo + 2.0*margin )*4.0 ) + 1;
   if( n > XPH_MAXROW ) {
      fail( "testxphmap: the scan of the surroundings of a grid is too big" );
      n = 0;
   }

   ninside = 0;
   noutside = 0;
   for( j = 0; j < n && status == 0; j++ ) {
      for( i = 0; i < n; i++ ) {
         xin[ i ] = lo - margin + 0.25*i;
         yin[ i ] = lo - margin + 0.25*j;
      }
      astTran2( map, n, xin, yin, forward, xout, yout );
      if( !astOK ) {
         fail( "testxphmap: astTran2 failed around the grid" );
         break;
      }
      for( i = 0; i < n; i++ ) {
         if( xin[ i ] >= lo && xin[ i ] <= hi &&
             yin[ i ] >= lo && yin[ i ] <= hi ) {
            if( xout[ i ] != AST__BAD && yout[ i ] != AST__BAD ) ninside++;
         } else {
            noutside++;
            if( xout[ i ] != AST__BAD || yout[ i ] != AST__BAD ) {
               printf( "testxphmap: type %s order %d %s: (%g,%g) is outside "
                       "the grid but gave (%g,%g)\n", type, order, dir,
                       xin[ i ], yin[ i ], xout[ i ], yout[ i ] );
               fail( "testxphmap: a position outside the grid was "
                     "transformed" );
               break;
            }
         }
      }
   }

/* Guard against a scan that missed the grid or its surroundings, which
   would make the loop above vacuous. */
   if( status == 0 && ( ninside < 12 || noutside < 12 ) ) {
      printf( "testxphmap: type %s order %d %s: the scan transformed %d "
              "positions inside the grid and tried %d outside it\n", type,
              order, dir, ninside, noutside );
      fail( "testxphmap: the scan around the grid is not a test" );
   }

   astEnd;
#undef XPH_MAXROW
}

/* Check the boundary of every cell.  The scan is a lattice of positions a
   quarter of a pixel apart, over the grid and a margin of one pixel, so it
   holds the centre of each cell, positions within the cell a quarter of a
   pixel from the centre along each axis, and the corners and edge centres
   of the cell.  A cell is used if its centre is transformed, and the
   positions within it then give the shift and rotation that the cell is
   transformed by.  Each corner and edge centre must be transformed to the
   position that continues the transformation of one of the used cells it
   bounds, or to bad values if it bounds none.  All the values involved are
   multiples of a quarter, so the comparison is exact. */
static void check_boundaries( const char *type, int order, int forward ) {
#define XPH_MAXLAT 89
   static double xin[ XPH_MAXLAT*XPH_MAXLAT ], yin[ XPH_MAXLAT*XPH_MAXLAT ];
   static double xout[ XPH_MAXLAT*XPH_MAXLAT ], yout[ XPH_MAXLAT*XPH_MAXLAT ];
   AstMapping *map;
   char dump[ 256 ];
   const char *dir = forward ? "forward" : "inverse";
   int nppf = 1 << order;
   int ncell = ( forward && type[ 0 ] == 'X' ) ? 4*nppf : 5*nppf;
   int nlat = 4*( ncell + 2 ) + 1;
   int a, b, c, ca, cb, ia, ib, k, kx[ 2 ], ky[ 2 ], match, nkx, nky;
   int nbound, nedge, nused;
   double wantx, wanty;

   if( status != 0 ) return;
   if( nlat > XPH_MAXLAT ) {
      fail( "testxphmap: the scan of the cell boundaries is too big" );
      return;
   }
   astBegin;

   make_dump( dump, sizeof dump, order, type );
   map = astFromString( dump );

/* Lattice index "a" is at grid coordinate ( 0.25*a - 0.5 ), so the centre
   of the cell with one-based index "k" is at lattice index ( 4*k + 2 ) and
   its edges are two lattice steps either side. */
   for( b = 0; b < nlat; b++ ) {
      for( a = 0; a < nlat; a++ ) {
         xin[ b*nlat + a ] = 0.25*a - 0.5;
         yin[ b*nlat + a ] = 0.25*b - 0.5;
      }
   }
   astTran2( map, nlat*nlat, xin, yin, forward, xout, yout );
   if( !astOK ) {
      fail( "testxphmap: astTran2 failed in the scan of cell boundaries" );
      astEnd;
      return;
   }

/* Visit each corner, edge centre and cell centre. */
   nedge = 0;
   nbound = 0;
   for( b = 0; b < nlat && status == 0; b += 2 ) {
      for( a = 0; a < nlat && status == 0; a += 2 ) {
         c = b*nlat + a;

/* The one-based indices of the cells that hold this position on each
   axis: one cell if the position is at the centre of a cell on that axis,
   and two if it is on an edge. */
         nkx = 0;
         if( a % 4 == 2 ) {
            kx[ nkx++ ] = ( a - 2 )/4;
         } else {
            kx[ nkx++ ] = a/4 - 1;
            kx[ nkx++ ] = a/4;
         }
         nky = 0;
         if( b % 4 == 2 ) {
            ky[ nky++ ] = ( b - 2 )/4;
         } else {
            ky[ nky++ ] = b/4 - 1;
            ky[ nky++ ] = b/4;
         }
         if( nkx*nky == 1 ) continue;
         nedge++;

/* See if the position continues the transformation of any of the used
   cells among them. */
         nused = 0;
         match = 0;
         for( ib = 0; ib < nky; ib++ ) {
            for( ia = 0; ia < nkx; ia++ ) {
               if( kx[ ia ] < 1 || kx[ ia ] > ncell ||
                   ky[ ib ] < 1 || ky[ ib ] > ncell ) continue;
               ca = 4*kx[ ia ] + 2;
               cb = 4*ky[ ib ] + 2;
               k = cb*nlat + ca;
               if( xout[ k ] == AST__BAD ) continue;
               nused++;
               if( xout[ k + 1 ] == AST__BAD || xout[ k + nlat ] == AST__BAD ) {
                  printf( "testxphmap: type %s order %d %s: the cell at "
                          "(%g,%g) is not transformed as a whole\n", type,
                          order, dir, xin[ k ], yin[ k ] );
                  fail( "testxphmap: a used cell is not transformed" );
                  continue;
               }
               wantx = xout[ k ] + ( a - ca )*( xout[ k + 1 ] - xout[ k ] )
                                 + ( b - cb )*( xout[ k + nlat ] - xout[ k ] );
               wanty = yout[ k ] + ( a - ca )*( yout[ k + 1 ] - yout[ k ] )
                                 + ( b - cb )*( yout[ k + nlat ] - yout[ k ] );
               if( xout[ c ] == wantx && yout[ c ] == wanty ) match = 1;
            }
         }

         if( nused ) nbound++;
         if( nused && !match ) {
            printf( "testxphmap: type %s order %d %s: (%g,%g) is on the "
                    "boundary of a used cell but gave (%g,%g)\n", type,
                    order, dir, xin[ c ], yin[ c ], xout[ c ], yout[ c ] );
            fail( "testxphmap: wrong result on the boundary of a cell" );
         } else if( !nused && ( xout[ c ] != AST__BAD ||
                                yout[ c ] != AST__BAD ) ) {
            printf( "testxphmap: type %s order %d %s: (%g,%g) bounds no used "
                    "cell but gave (%g,%g)\n", type, order, dir, xin[ c ],
                    yin[ c ], xout[ c ], yout[ c ] );
            fail( "testxphmap: a position off the sky was transformed" );
         }
      }
   }

/* Guard against a scan that found no boundaries, or none off the sky. */
   if( status == 0 && ( nbound < 12 || nedge - nbound < 12 ) ) {
      printf( "testxphmap: type %s order %d %s: %d of %d boundary positions "
              "bound a used cell\n", type, order, dir, nbound, nedge );
      fail( "testxphmap: the scan of cell boundaries is not a test" );
   }

   astEnd;
#undef XPH_MAXLAT
}

/* The transformation uses the order as a shift count, so the loader must
   reject an order outside the range the constructor accepts, 0 to 27. */
static void check_bad_order( int order ) {
   AstObject *obj;
   char dump[ 256 ];

   if( status != 0 ) return;
   astBegin;

   make_dump( dump, sizeof dump, order, "HPX12" );
   obj = astFromString( dump );
   if( astOK ) {
      printf( "testxphmap: an Order of %d was accepted\n", order );
      fail( "testxphmap: an Order out of range was accepted" );
   } else if( status != AST__LDERR ) {
      printf( "testxphmap: an Order of %d gave error %d, not AST__LDERR\n",
              order, status );
      astClearStatus;
      fail( "testxphmap: wrong error for an Order out of range" );
   } else {
      astClearStatus;
   }
   (void) obj;

   astEnd;
}

/* Transform one position and compare with the wanted result, where a wanted
   X of AST__BAD means that both outputs must be bad. */
static void check_point( AstMapping *map, int forward, const char *what,
                         double x, double y, double xwant, double ywant ) {
   double xout, yout;

   if( status != 0 ) return;
   astTran2( map, 1, &x, &y, forward, &xout, &yout );
   if( !astOK ) {
      fail( "testxphmap: astTran2 failed for a single position" );
   } else if( xwant == AST__BAD ? ( xout != AST__BAD || yout != AST__BAD )
                                : ( xout != xwant || yout != ywant ) ) {
      printf( "testxphmap: %s: (%.17g,%.17g) gave (%.17g,%.17g)\n", what,
              x, y, xout, yout );
      fail( "testxphmap: wrong result for a single position" );
   }
}

/* Transform one position and require the result to be one of four given
   positions. */
static void check_one_of( AstMapping *map, int forward, const char *what,
                          double x, double y, const double xset[ 4 ],
                          const double yset[ 4 ] ) {
   double xout, yout;
   int i, found;

   if( status != 0 ) return;
   astTran2( map, 1, &x, &y, forward, &xout, &yout );
   found = 0;
   for( i = 0; i < 4; i++ ) {
      if( xout == xset[ i ] && yout == yset[ i ] ) found = 1;
   }
   if( !astOK ) {
      fail( "testxphmap: astTran2 failed at a pole" );
   } else if( !found ) {
      printf( "testxphmap: %s: (%.17g,%.17g) gave (%.17g,%.17g)\n", what,
              x, y, xout, yout );
      fail( "testxphmap: a pole was not transformed to a pole" );
   }
}

/* A pole is one point on the sky.  In an HPX grid the north pole is the top
   left corner of each of the four north polar facets and the south pole is
   the bottom right corner of each of the four south polar facets.  In an XPH
   grid the pole the projection is centred on is the centre of the grid and
   the other pole is its four corners.  Every one of these positions must be
   transformed, and to a position of the same pole. */
static void check_poles( int order ) {
   static const char *types[ 4 ] = { "HPX0", "HPX12", "XPHN", "XPHS" };
   AstMapping *map;
   char dump[ 256 ];
   char what[ 80 ];
   double n = (double)( 1 << order );
   double ntip_x[ 4 ], ntip_y[ 4 ], stip_x[ 4 ], stip_y[ 4 ];
   double corner_x[ 4 ], corner_y[ 4 ], centre_x[ 4 ], centre_y[ 4 ];
   const double *nx, *ny, *sx, *sy;
   int itype, k;

   if( status != 0 ) return;
   astBegin;

   for( k = 0; k < 4; k++ ) {
      ntip_x[ k ] = k*n + 0.5;
      ntip_y[ k ] = ( k + 2 )*n + 0.5;
      stip_x[ k ] = ( k + 2 )*n + 0.5;
      stip_y[ k ] = k*n + 0.5;
      corner_x[ k ] = ( k % 2 ) ? 4*n + 0.5 : 0.5;
      corner_y[ k ] = ( k / 2 ) ? 4*n + 0.5 : 0.5;
      centre_x[ k ] = 2*n + 0.5;
      centre_y[ k ] = 2*n + 0.5;
   }

   for( itype = 0; itype < 4 && status == 0; itype++ ) {
      make_dump( dump, sizeof dump, order, types[ itype ] );
      map = astFromString( dump );

/* The positions of the north and south poles in the grid of this type. */
      if( itype == 2 ) {
         nx = centre_x;
         ny = centre_y;
         sx = corner_x;
         sy = corner_y;
      } else if( itype == 3 ) {
         nx = corner_x;
         ny = corner_y;
         sx = centre_x;
         sy = centre_y;
      } else {
         nx = ntip_x;
         ny = ntip_y;
         sx = stip_x;
         sy = stip_y;
      }

/* The forward transformation takes each to a pole of an HPX12 grid, and
   the inverse takes each pole of an HPX12 grid to one of them. */
      for( k = 0; k < 4; k++ ) {
         snprintf( what, sizeof what, "type %s order %d forward north pole %d",
                   types[ itype ], order, k );
         check_one_of( map, 1, what, nx[ k ], ny[ k ], ntip_x, ntip_y );
         snprintf( what, sizeof what, "type %s order %d forward south pole %d",
                   types[ itype ], order, k );
         check_one_of( map, 1, what, sx[ k ], sy[ k ], stip_x, stip_y );
         snprintf( what, sizeof what, "type %s order %d inverse north pole %d",
                   types[ itype ], order, k );
         check_one_of( map, 0, what, ntip_x[ k ], ntip_y[ k ], nx, ny );
         snprintf( what, sizeof what, "type %s order %d inverse south pole %d",
                   types[ itype ], order, k );
         check_one_of( map, 0, what, stip_x[ k ], stip_y[ k ], sx, sy );
      }
   }

   astEnd;
}

/* The grid coordinates of a position on the edge of the sky are calculated
   from its sky coordinates, and rounding can leave it marginally off the
   edge.  The left edge of a north polar facet of an HPX12 grid is such an
   edge: the cells to its left are not used.  A position displaced to the
   left of it by less than the margin the HPX projection allows for rounding,
   1E-12 degrees, where a pixel is 45/2^order degrees, must be transformed as
   the position on the edge is.  One displaced by much more must give bad
   values. */
static void check_rounding( const char *type, int order ) {
   AstMapping *map;
   char dump[ 256 ];
   char what[ 80 ];
   double n = (double)( 1 << order );
   double tol = 1.0E-12*n/45.0;
   double x = n + 0.5;
   double y = 2*n + 0.75;
   double xedge, yedge;

   if( status != 0 ) return;
   astBegin;

   make_dump( dump, sizeof dump, order, type );
   map = astFromString( dump );

   astTran2( map, 1, &x, &y, 0, &xedge, &yedge );
   if( !astOK || xedge == AST__BAD || yedge == AST__BAD ) {
      printf( "testxphmap: type %s order %d: (%.17g,%.17g) on the edge of a "
              "polar facet was not transformed\n", type, order, x, y );
      fail( "testxphmap: a position on the edge of the sky gave bad values" );
   }

   snprintf( what, sizeof what, "type %s order %d just off the edge", type,
             order );
   check_point( map, 0, what, x - 0.5*tol, y, xedge, yedge );
   snprintf( what, sizeof what, "type %s order %d well off the edge", type,
             order );
   check_point( map, 0, what, x - 100.0*tol, y, AST__BAD, 0.0 );

   astEnd;
}

/* At the largest order an HPX grid is 5*2^27 pixels wide and an XPH grid is
   4*2^27.  Check pixels at both ends of the grid, positions beyond them, the
   unused copy of a pixel, and a pole. */
static void check_max_order( void ) {
   AstMapping *map;
   char dump[ 256 ];
   double n = 134217728.0;

   if( status != 0 ) return;
   astBegin;

   make_dump( dump, sizeof dump, 27, "HPX12" );
   map = astFromString( dump );
   check_point( map, 1, "first HPX12 column", 1.0, n + 1.0, 1.0, n + 1.0 );
   check_point( map, 1, "before the HPX12 grid", -n, n + 1.0, AST__BAD, 0.0 );
   check_point( map, 1, "just before the HPX12 grid", 0.0, 1.5*n + 1.0,
                AST__BAD, 0.0 );
   check_point( map, 1, "last HPX12 column", 5.0*n, 4.0*n + 1.0,
                5.0*n, 4.0*n + 1.0 );
   check_point( map, 1, "after the HPX12 grid", 6.0*n, 4.0*n + 1.0,
                AST__BAD, 0.0 );
   check_point( map, 1, "unused HPX12 corner", 5.0*n, 5.0*n, AST__BAD, 0.0 );
   check_point( map, 0, "unused HPX12 corner, inverse", 5.0*n, 5.0*n, n, n );

   make_dump( dump, sizeof dump, 27, "XPHN" );
   map = astFromString( dump );
   check_point( map, 1, "last XPHN column", 4.0*n, 1.0, 2.0*n, 1.0 );
   check_point( map, 1, "after the XPHN grid", 4.0*n + 1.0, 1.0,
                AST__BAD, 0.0 );
   check_point( map, 0, "north pole to XPHN", n + 0.5, 3.0*n + 0.5,
                2.0*n + 0.5, 2.0*n + 0.5 );
   check_point( map, 1, "north pole from XPHN", 2.0*n + 0.5, 2.0*n + 0.5,
                3.0*n + 0.5, 5.0*n + 0.5 );

   astEnd;
}

int main( void ) {
   static const char *types[ 4 ] = { "HPX0", "HPX12", "XPHN", "XPHS" };
   int forward, itype, order;

   astWatch( &status );

/* Every name in the loader's table, so each arm of its lookup is taken. */
   check_type( "HPX0", 14 );
   check_type( "HPX12", 8 );
   check_type( "XPHN", 10 );
   check_type( "XPHS", 10 );
   check_transform();
   check_bad_type();

/* Every type, in both directions, from the order at which a facet is one
   pixel to one at which every facet has an interior. */
   for( itype = 0; itype < 4; itype++ ) {
      for( order = 0; order <= 3; order++ ) {
         for( forward = 1; forward >= 0; forward-- ) {
            check_coverage( types[ itype ], order, forward );
         }
      }
      check_far_positions( types[ itype ], 0 );
      check_far_positions( types[ itype ], 3 );
      check_far_positions( types[ itype ], 27 );
   }

/* Every type, in both directions: the surroundings of the grid, and the
   boundary of every cell. */
   for( itype = 0; itype < 4; itype++ ) {
      for( forward = 1; forward >= 0; forward-- ) {
         for( order = 0; order <= 3; order++ ) {
            check_outside_grid( types[ itype ], order, forward );
         }
         for( order = 0; order <= 2; order++ ) {
            check_boundaries( types[ itype ], order, forward );
         }
      }
      check_rounding( types[ itype ], 0 );
      check_rounding( types[ itype ], 3 );
      check_rounding( types[ itype ], 27 );
   }

   check_poles( 0 );
   check_poles( 1 );
   check_poles( 3 );
   check_poles( 27 );

   check_bad_order( -1 );
   check_bad_order( 28 );
   check_bad_order( 31 );
   check_bad_order( 32 );
   check_max_order();

   if( status == 0 ) {
      printf( "All XphMap tests passed\n" );
   } else {
      printf( "XphMap tests failed\n" );
   }
   return status != 0;
}
