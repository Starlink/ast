/*
*  Purpose:
*     Exercise the XphMap loader: the projection-type names it accepts, the
*     components that survive a dump/load round trip, and the transformation
*     of a loaded XphMap.
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
*/
#include "ast.h"

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

int main( void ) {
   astWatch( &status );

/* Every name in the loader's table, so each arm of its lookup is taken. */
   check_type( "HPX0", 14 );
   check_type( "HPX12", 8 );
   check_type( "XPHN", 10 );
   check_type( "XPHS", 10 );
   check_transform();
   check_bad_type();

   if( status == 0 ) {
      printf( "All XphMap tests passed\n" );
   } else {
      printf( "XphMap tests failed\n" );
   }
   return status != 0;
}
