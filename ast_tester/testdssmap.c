/*
*  Name:
*     testdssmap

*  Purpose:
*     Check the DssMap class, which can only be reached by loading a dump.

*  Description:
*     A DssMap describes a Digitised Sky Survey plate solution.  AST no
*     longer creates DssMaps, since it now reads a DSS header as a
*     combination of other Mappings, and the class has no public
*     constructor.  It remains so that dumps written by older versions of
*     AST can still be read, which makes its loader the only way in.

*     This test reads the dump of a DssMap from a fixture file and checks:
*
*       - that pixel positions on the plate are transformed to the expected
*         sky positions, and back again;
*       - that sky positions on the far side of the sky from the plate
*         centre, which a plate cannot show, are given bad pixel positions;
*       - that a position too large for the plate solution to evaluate is
*         given a bad sky position rather than one that is not a number;
*       - that a dump is rejected if a required value is missing, if the
*         numbers of inputs and outputs are not two, or if the plate centre
*         or pixel offsets are too large to be described by the DSS FITS
*         keywords from which a DssMap is built;
*       - that a DssMap with a pixel size of zero loads, but gives bad
*         pixel positions;
*       - that, when it is simplified, a DssMap absorbs a neighbouring
*         WinMap that leaves it with a whole number of pixels as its pixel
*         offsets, and leaves alone one that would need a fraction of a
*         pixel or more pixels than the DSS FITS keywords can record, with
*         the same transformation resulting in every case.

*  Notes:
*     - This is a native C test with no Fortran original.
*     - Public C API only.
*     - The fixture is fixtures/serialisation/c_dssmap.ast, found through
*       the "srcdir" environment variable.
*     - The oversized WinMap shift used to take an out-of-range conversion
*       of a floating point value to an integer, which this test can only
*       detect when it is run under UndefinedBehaviorSanitizer.
*/

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sae_par.h"
#include "ast.h"
#include "ast_err.h"

/* The plate centre of the DssMap in the fixture, in radians. */
#define PLATE_RA 2.8272374655684112
#define PLATE_DEC ( -1.0518122540502668 )

/* Read the fixture file into memory allocated by astMalloc.  Returns NULL
   if the file cannot be read. */
static char *read_fixture( void ) {
   const char *srcdir = getenv( "srcdir" );
   char fname[ PATH_MAX ];
   char *dump = NULL;
   FILE *in;
   long size;

   if ( !srcdir || !*srcdir ) srcdir = ".";
   (void) snprintf( fname, sizeof( fname ), "%s/fixtures/serialisation/"
                    "c_dssmap.ast", srcdir );
   in = fopen( fname, "r" );
   if ( !in ) {
      fprintf( stderr, "testdssmap: could not open '%s'\n", fname );
      return NULL;
   }
   if ( fseek( in, 0L, SEEK_END ) == 0 && ( size = ftell( in ) ) > 0 &&
        fseek( in, 0L, SEEK_SET ) == 0 ) {
      dump = astMalloc( (size_t) size + 1 );
      if ( dump ) {
         if ( fread( dump, 1, (size_t) size, in ) == (size_t) size ) {
            dump[ size ] = '\0';
         } else {
            dump = astFree( dump );
         }
      }
   }
   (void) fclose( in );
   if ( !dump ) fprintf( stderr, "testdssmap: could not read '%s'\n", fname );
   return dump;
}

/* Return a copy of "dump" in which the line holding the item called "key"
   is replaced by "text", which may be empty to remove the item.  The copy
   is allocated with astMalloc.  Returns NULL if the item is not found. */
static char *set_line( const char *dump, const char *key, const char *text ) {
   char needle[ 64 ];
   const char *pos;
   const char *bol;
   const char *eol;
   char *out;
   size_t pre;
   size_t tlen = strlen( text );

   (void) snprintf( needle, sizeof( needle ), " %s = ", key );
   pos = strstr( dump, needle );
   if ( !pos ) return NULL;
   bol = pos;
   while ( bol > dump && bol[ -1 ] != '\n' ) bol--;
   eol = strchr( pos, '\n' );
   if ( !eol ) return NULL;

   pre = (size_t) ( bol - dump );
   out = astMalloc( pre + tlen + strlen( eol + 1 ) + 1 );
   if ( out ) {
      memcpy( out, dump, pre );
      memcpy( out + pre, text, tlen );
      strcpy( out + pre + tlen, eol + 1 );
   }
   return out;
}

/* Attempt to load a dump that should be rejected with the error code
   "expected", and free it.  Returns 0 if no object was returned and that
   error was reported, 1 otherwise. */
static int check_rejected( const char *what, char *dump, int expected ) {
   AstObject *obj;
   int code;
   int loaded;

   if ( !dump ) {
      fprintf( stderr, "testdssmap: could not build the dump for %s\n",
               what );
      if ( !astOK ) astClearStatus;
      return 1;
   }
   obj = astFromString( dump );
   code = astStatus;
   loaded = ( obj != NULL );
   if ( !astOK ) astClearStatus;
   if ( obj ) obj = astAnnul( obj );
   dump = astFree( dump );
   if ( loaded ) {
      fprintf( stderr, "testdssmap: %s was accepted but should have been "
               "rejected\n", what );
      return 1;
   }
   if ( code != expected ) {
      fprintf( stderr, "testdssmap: %s was rejected with status %d but "
               "status %d was expected\n", what, code, expected );
      return 1;
   }
   return 0;
}

/* Check that transforming the position (xin,yin) with the Mapping, in the
   direction given by "forward", gives (xwant,ywant) to within "tol".  Bad
   wanted values must be matched exactly.  Returns 0 if so, 1 otherwise. */
static int check_tran( const char *what, AstMapping *map, int forward,
                       double xin, double yin, double xwant, double ywant,
                       double tol ) {
   double xout = 0.0;
   double yout = 0.0;
   int ok;

   astTran2( map, 1, &xin, &yin, forward, &xout, &yout );
   if ( !astOK ) {
      fprintf( stderr, "testdssmap: %s reported an error\n", what );
      astClearStatus;
      return 1;
   }
   if ( xwant == AST__BAD ) {
      ok = ( xout == AST__BAD && yout == AST__BAD );
   } else {
      ok = ( xout != AST__BAD && yout != AST__BAD &&
             fabs( xout - xwant ) <= tol && fabs( yout - ywant ) <= tol );
   }
   if ( !ok ) {
      fprintf( stderr, "testdssmap: %s gave (%.17g,%.17g), expected "
               "(%.17g,%.17g)\n", what, xout, yout, xwant, ywant );
      return 1;
   }
   return 0;
}

/* Check that two Mappings give the same sky position for a set of pixel
   positions.  Returns the number of positions at which they differ. */
static int check_same( const char *what, AstMapping *map1, AstMapping *map2 ) {
   double xin[ 4 ] = { 1.0, 250.0, -300.0, 1000.0 };
   double yin[ 4 ] = { 1.0, -125.0, 750.0, 1000.0 };
   double x1[ 4 ];
   double y1[ 4 ];
   double x2[ 4 ];
   double y2[ 4 ];
   int i;
   int nfail = 0;

   astTran2( map1, 4, xin, yin, 1, x1, y1 );
   astTran2( map2, 4, xin, yin, 1, x2, y2 );
   if ( !astOK ) {
      fprintf( stderr, "testdssmap: %s reported an error\n", what );
      astClearStatus;
      return 1;
   }
   for ( i = 0; i < 4; i++ ) {
      if ( x1[ i ] == AST__BAD || y1[ i ] == AST__BAD ||
           fabs( x1[ i ] - x2[ i ] ) > 1.0e-9 ||
           fabs( y1[ i ] - y2[ i ] ) > 1.0e-9 ) {
         fprintf( stderr, "testdssmap: %s differ at (%g,%g): (%.17g,%.17g) "
                  "and (%.17g,%.17g)\n", what, xin[ i ], yin[ i ], x1[ i ],
                  y1[ i ], x2[ i ], y2[ i ] );
         nfail++;
      }
   }
   return nfail;
}

/* Put a WinMap that doubles both pixel axes and then shifts them by
   "shift" in front of the DssMap, simplify the combination, and check that
   the result has the class "class" and transforms positions as the
   combination does.  Returns the number of failures. */
static int check_merge( const char *what, AstMapping *dss, double shift,
                        const char *class ) {
   double ina[ 2 ] = { 0.0, 0.0 };
   double inb[ 2 ] = { 1.0, 1.0 };
   double outa[ 2 ];
   double outb[ 2 ];
   AstWinMap *winmap;
   AstCmpMap *cmpmap;
   AstMapping *simple;
   int nfail = 0;

   outa[ 0 ] = shift;
   outa[ 1 ] = shift;
   outb[ 0 ] = shift + 2.0;
   outb[ 1 ] = shift + 2.0;
   winmap = astWinMap( 2, ina, inb, outa, outb, " " );
   cmpmap = astCmpMap( winmap, dss, 1, " " );
   simple = astSimplify( cmpmap );
   if ( !astOK ) {
      fprintf( stderr, "testdssmap: simplifying %s reported an error\n",
               what );
      astClearStatus;
      return 1;
   }
   if ( strcmp( astGetC( simple, "Class" ), class ) != 0 ) {
      fprintf( stderr, "testdssmap: %s simplified to a %s, expected a %s\n",
               what, astGetC( simple, "Class" ), class );
      nfail++;
   }

/* The largest shift moves the positions used by check_same far off the
   plate, where the two Mappings need not agree to the same accuracy, so
   only compare them for the other shifts. */
   if ( fabs( shift ) < 1.0e6 ) {
      nfail += check_same( what, (AstMapping *) cmpmap, simple );
   }
   simple = astAnnul( simple );
   cmpmap = astAnnul( cmpmap );
   winmap = astAnnul( winmap );
   return nfail;
}

int main( void ) {
   int status = SAI__OK;
   int nfail = 0;
   char *dump;
   AstObject *obj;
   AstMapping *dss;
   AstMapping *copy;

   astWatch( &status );

   dump = read_fixture();
   if ( !dump ) {
      printf( "DssMap tests failed (1)\n" );
      return 1;
   }

/* The dump in the fixture loads as a DssMap with two inputs and two
   outputs. */
   obj = astFromString( dump );
   if ( !astOK || !obj || strcmp( astGetC( obj, "Class" ), "DssMap" ) != 0 ||
        astGetI( obj, "Nin" ) != 2 || astGetI( obj, "Nout" ) != 2 ) {
      fprintf( stderr, "testdssmap: the fixture did not load as a DssMap\n" );
      if ( !astOK ) astClearStatus;
      printf( "DssMap tests failed (1)\n" );
      return 1;
   }
   dss = (AstMapping *) obj;

/* Pixel positions on the plate transform to the expected sky positions,
   and back again. */
   nfail += check_tran( "the sky position of pixel (1,1)", dss, 1, 1.0, 1.0,
                        2.8243236575347148, -1.050832220924965, 1.0e-12 );
   nfail += check_tran( "the sky position of pixel (1000,1000)", dss, 1,
                        1000.0, 1000.0, 2.8080654084049983,
                        -1.0425371444459393, 1.0e-12 );
   nfail += check_tran( "the pixel position of the first sky position", dss,
                        0, 2.8243236575347148, -1.050832220924965, 1.0, 1.0,
                        1.0e-8 );
   nfail += check_tran( "the pixel position of the second sky position", dss,
                        0, 2.8080654084049983, -1.0425371444459393, 1000.0,
                        1000.0, 1.0e-8 );

/* A plate shows only the half of the sky on the same side as its centre.
   Sky positions in the other half have no pixel position: the point
   opposite the plate centre, the pole in the other hemisphere, and a
   point on the equator that is more than 90 degrees from the centre. */
   nfail += check_tran( "the pixel position of the point opposite the plate "
                        "centre", dss, 0, PLATE_RA + AST__DPI, -PLATE_DEC,
                        AST__BAD, AST__BAD, 0.0 );
   nfail += check_tran( "the pixel position of the far pole", dss, 0, 0.0,
                        AST__DPIBY2, AST__BAD, AST__BAD, 0.0 );
   nfail += check_tran( "the pixel position of a point on the far side of "
                        "the sky", dss, 0, 0.0, 0.0, AST__BAD, AST__BAD,
                        0.0 );

/* A pixel position so large that the plate solution overflows has no sky
   position. */
   nfail += check_tran( "the sky position of an enormous pixel position",
                        dss, 1, 1.0e300, 1.0e300, AST__BAD, AST__BAD, 0.0 );

/* A copy is equal to the original. */
   copy = astCopy( dss );
   if ( !astOK || !astEqual( copy, dss ) ) {
      fprintf( stderr, "testdssmap: a copy of the DssMap is not equal to "
               "it\n" );
      if ( !astOK ) astClearStatus;
      nfail++;
   }
   if ( copy ) copy = astAnnul( copy );

/* A dump is rejected if a required value is missing. */
   nfail += check_rejected( "DssMap with no plate centre RA",
                            set_line( dump, "PltRA", "" ), AST__RDERR );
   nfail += check_rejected( "DssMap with no X pixel size",
                            set_line( dump, "XPixSz", "" ), AST__RDERR );
   nfail += check_rejected( "DssMap with a missing plate solution "
                            "coefficient", set_line( dump, "AMDY5", "" ),
                            AST__RDERR );

/* A DssMap has two inputs and two outputs. */
   nfail += check_rejected( "DssMap with three inputs",
                            set_line( dump, "Nin", " Nin = 3\n Nout = 2\n" ),
                            AST__LDERR );
   nfail += check_rejected( "DssMap with one output",
                            set_line( dump, "Nin", " Nin = 2\n Nout = 1\n" ),
                            AST__LDERR );

/* A DssMap is built from DSS FITS keywords, which hold the pixel offsets
   and the hours and degrees of the plate centre as integers.  Values that
   are too large for them cannot have come from such keywords. */
   nfail += check_rejected( "DssMap with an enormous pixel offset",
                            set_line( dump, "CNPix1", " CNPix1 = 1e300\n" ),
                            AST__LDERR );
   nfail += check_rejected( "DssMap with an enormous negative pixel offset",
                            set_line( dump, "CNPix2", " CNPix2 = -3e9\n" ),
                            AST__LDERR );
   nfail += check_rejected( "DssMap with an enormous plate centre RA",
                            set_line( dump, "PltRA", " PltRA = 1e300\n" ),
                            AST__LDERR );
   nfail += check_rejected( "DssMap with an enormous plate centre Dec",
                            set_line( dump, "PltDec", " PltDec = -1e300\n" ),
                            AST__LDERR );

/* A pixel size of zero cannot be divided into a distance on the plate, so
   a DssMap that has one gives bad pixel positions. */
   {
      char *zero = set_line( dump, "XPixSz", " XPixSz = 0\n" );
      AstObject *zobj = zero ? astFromString( zero ) : NULL;
      if ( !astOK || !zobj ) {
         fprintf( stderr, "testdssmap: a DssMap with a pixel size of zero "
                  "failed to load\n" );
         if ( !astOK ) astClearStatus;
         nfail++;
      } else {
         nfail += check_tran( "the pixel position given by a DssMap with a "
                              "pixel size of zero", (AstMapping *) zobj, 0,
                              PLATE_RA, PLATE_DEC, AST__BAD, AST__BAD, 0.0 );
      }
      if ( zobj ) zobj = astAnnul( zobj );
      if ( zero ) zero = astFree( zero );
   }

/* When simplified, a DssMap absorbs a preceding WinMap by changing its
   pixel sizes and offsets, but only if the new offsets are whole numbers
   of pixels, since the DSS FITS keywords hold them as integers.  The
   offsets in the fixture are 4201 pixels, so doubling the pixel
   coordinates and adding half a pixel gives new offsets of 2101 pixels.
   Adding three quarters of a pixel instead gives a fraction of a pixel,
   and adding 2e10 pixels gives more pixels than an integer can hold. */
   nfail += check_merge( "a DssMap and a WinMap giving whole pixel offsets",
                         dss, 0.5, "DssMap" );
   nfail += check_merge( "a DssMap and a WinMap giving fractional pixel "
                         "offsets", dss, 0.75, "CmpMap" );
   nfail += check_merge( "a DssMap and a WinMap giving enormous pixel "
                         "offsets", dss, 2.0e10 + 0.5, "CmpMap" );

   dss = astAnnul( dss );
   dump = astFree( dump );

   if ( !astOK ) {
      fprintf( stderr, "testdssmap: error on exit\n" );
      astClearStatus;
      nfail++;
   }

   if ( nfail == 0 ) {
      printf( " All DssMap tests passed\n" );
      return 0;
   }
   printf( "DssMap tests failed (%d)\n", nfail );
   return 1;
}
