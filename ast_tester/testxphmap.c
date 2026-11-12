/*
*  Purpose:
*     Exercise the XphMap loader: the projection-type names it accepts, the
*     components that survive a dump/load round trip, and the transformation
*     of a loaded XphMap.
*
*     Every HEALPix pixel on the sky appears exactly once in each projection,
*     so of the pixel centres in the grid of a projection, the number that an
*     XphMap transforms must be the number of pixels on the sky, and each
*     must come back through the opposite transformation.  Positions too far
*     from the grid for the index of a grid cell to be held in an integer,
*     and positions that are not finite, must not give values that are not
*     finite, and must give bad values wherever the transformation tests a
*     position against the extent of the grid.  An Order that the constructor
*     would not accept must be rejected by the loader, and the largest Order
*     must work.
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
   int i, ix, iy, n, ngood;

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

   ngood = 0;
   for( i = 0; i < n && status == 0; i++ ) {
      if( xm[ i ] == AST__BAD || ym[ i ] == AST__BAD ) continue;
      ngood++;
      if( xout[ i ] != xin[ i ] || yout[ i ] != yin[ i ] ) {
         printf( "testxphmap: type %s order %d %s: (%g,%g) came back as "
                 "(%g,%g)\n", type, order, dir, xin[ i ], yin[ i ],
                 xout[ i ], yout[ i ] );
         fail( "testxphmap: a pixel centre did not make the return trip" );
      }
   }

   if( status == 0 && ngood != 12*nppf*nppf ) {
      printf( "testxphmap: type %s order %d %s: %d pixel centres were "
              "transformed, but the sky has %d pixels\n", type, order, dir,
              ngood, 12*nppf*nppf );
      fail( "testxphmap: wrong number of pixels in the projection" );
   }

   astEnd;
#undef XPH_MAXSCAN
}

/* Transform positions that no grid index held in an integer can describe,
   and positions that are not finite.  Each far value is tried on both axes,
   with every other far value and with a position inside the grid on the
   other axis, since the quadrant of an XPH grid that a position falls in
   selects the arithmetic used.

   A position that is not finite must give bad values.  So must any far
   position, wherever the transformation tests a position against the extent
   of the grid.  The forward transformation of an XPH projection does not:
   it tests a position against the sides of a gore but not against the
   outer edges of the grid, so a far position in line with a gore gives
   good values, and those are only required to be finite. */
static void check_far_positions( const char *type, int order ) {
   AstMapping *map;
   char dump[ 256 ];
   double far[ 8 ];
   double xin[ 80 ], yin[ 80 ], xout[ 80 ], yout[ 80 ];
   int forward, i, j, n, bounded, wantbad;
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
      bounded = !( forward && type[ 0 ] == 'X' );
      astTran2( map, n, xin, yin, forward, xout, yout );
      if( !astOK ) {
         fail( "testxphmap: astTran2 failed for far positions" );
         break;
      }
      for( i = 0; i < n; i++ ) {
         wantbad = bounded || !isfinite( xin[ i ] ) || !isfinite( yin[ i ] );
         if( wantbad ? ( xout[ i ] != AST__BAD || yout[ i ] != AST__BAD )
                     : ( !isfinite( xout[ i ] ) || !isfinite( yout[ i ] ) ) ) {
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
static void check_point( AstMapping *map, const char *what, double x,
                         double y, double xwant, double ywant ) {
   double xout, yout;

   if( status != 0 ) return;
   astTran2( map, 1, &x, &y, 1, &xout, &yout );
   if( !astOK ) {
      fail( "testxphmap: astTran2 failed at the largest order" );
   } else if( xwant == AST__BAD ? ( xout != AST__BAD || yout != AST__BAD )
                                : ( xout != xwant || yout != ywant ) ) {
      printf( "testxphmap: %s: (%.17g,%.17g) gave (%.17g,%.17g)\n", what,
              x, y, xout, yout );
      fail( "testxphmap: wrong result at the largest order" );
   }
}

/* At the largest order an HPX grid is 5*2^27 pixels wide.  Check pixels at
   both ends of the grid, and positions a facet or more beyond them. */
static void check_max_order( void ) {
   AstMapping *map;
   char dump[ 256 ];
   double n = 134217728.0;

   if( status != 0 ) return;
   astBegin;

   make_dump( dump, sizeof dump, 27, "HPX12" );
   map = astFromString( dump );
   check_point( map, "first HPX12 column", 1.0, n + 1.0, 1.0, n + 1.0 );
   check_point( map, "before the HPX12 grid", -n, n + 1.0, AST__BAD, 0.0 );
   check_point( map, "last HPX12 column", 5.0*n, 4.0*n + 1.0,
                5.0*n, 4.0*n + 1.0 );
   check_point( map, "after the HPX12 grid", 6.0*n, 4.0*n + 1.0,
                AST__BAD, 0.0 );
   check_point( map, "unused HPX12 corner", 5.0*n, 5.0*n, AST__BAD, 0.0 );

   make_dump( dump, sizeof dump, 27, "XPHN" );
   map = astFromString( dump );
   check_point( map, "last XPHN column", 4.0*n, 1.0, 2.0*n, 1.0 );

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
