/*
*  Name:
*     testloaderhardening

*  Purpose:
*     Check that object loaders validate structural data read from a dump.

*  Description:
*     A class loader reads counts, sizes, indices and component objects from
*     a dump.  The constructor of the same class checks the equivalent
*     values, so every object it creates satisfies invariants that the rest
*     of the class relies on to stay within its arrays.  This test confirms
*     that each loader enforces those invariants too.  For each class it
*     loads a hand-crafted dump that is valid, then loads variants of it
*     that each break one invariant and checks that every variant is
*     rejected at load time with the expected error code.  It also checks
*     that an object built by the constructor survives a dump and load.
*     The crafted cases cover:
*
*       - PolyMap: an axis power larger than the stated maximum power, a
*         negative number of coefficients and a negative maximum power.
*       - LutMap: fewer than two lookup table elements, a zero increment,
*         and an inverse transformation claimed for a table with fewer than
*         two good values (which is loaded, but without the inverse).
*       - WcsMap: longitude and latitude axis indices that are out of range
*         or equal, and differing numbers of inputs and outputs.
*       - SplineMap: too few coefficients for the spline order, knots that
*         are out of order, coefficient counts whose product overflows, and
*         a number of inputs other than two.
*       - Plot: a base Frame with more axes than the fixed per-axis arrays,
*         a negative number of tick mark values, and clipping information
*         that matches no Frame in the Plot.
*       - Plot3D: an axis labelled by a non-existent Plot or Plot axis, and
*         an encapsulated Plot that is missing or is not a Plot.
*       - Moc: ranges without a MaxOrder, a negative number of ranges, a
*         negative lower bound, and more ranges claimed than are present.
*
*     Two things that the loaders depend on are checked as well.  The memory
*     allocation functions must report an error for a size that cannot be
*     represented, rather than return a smaller block.  A SplineMap whose
*     end knots are not coincident, which passes every check made by the
*     loader, must give bad values outside the interval over which the
*     spline is defined.

*  Notes:
*     - This is a native C test with no Fortran original.
*     - Public C API only.  Errors expected from a deliberately invalid dump
*       are caught through astStatus and cleared with astClearStatus so the
*       remaining cases can run; failures of the test itself are counted and
*       reported at the end.
*/

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "sae_par.h"
#include "ast.h"
#include "ast_err.h"

/* Attempt to load a dump that should be rejected with the error code
   "expected".  Returns 0 if no object was returned and that error was
   reported, 1 otherwise.  Any raised status is cleared before returning so
   the caller can continue. */
static int check_rejected( const char *what, const char *dump, int expected ) {
   AstObject *obj = astFromString( dump );
   int code = astStatus;
   int loaded = ( obj != NULL );
   if ( !astOK ) astClearStatus;
   if ( obj ) obj = astAnnul( obj );
   if ( loaded ) {
      fprintf( stderr, "testloaderhardening: %s was accepted but should "
               "have been rejected\n", what );
      return 1;
   }
   if ( code != expected ) {
      fprintf( stderr, "testloaderhardening: %s was rejected with status %d "
               "but status %d was expected\n", what, code, expected );
      return 1;
   }
   return 0;
}

/* As check_rejected, but for a dump held in memory allocated by astMalloc,
   which is freed.  A NULL dump is reported as a failure. */
static int check_rejected_free( const char *what, char *dump, int expected ) {
   int nfail;
   if ( !dump ) {
      fprintf( stderr, "testloaderhardening: could not build the dump for "
               "%s\n", what );
      if ( !astOK ) astClearStatus;
      return 1;
   }
   nfail = check_rejected( what, dump, expected );
   dump = astFree( dump );
   return nfail;
}

/* Load a dump that should be accepted and return the object, which the
   caller should annul.  Returns NULL, after reporting the failure and
   incrementing "*nfail", if the dump does not load as class "class". */
static AstObject *load_accepted( const char *what, const char *dump,
                                 const char *class, int *nfail ) {
   AstObject *obj = astFromString( dump );
   if ( !astOK || !obj ) {
      fprintf( stderr, "testloaderhardening: %s failed to load\n", what );
      if ( !astOK ) astClearStatus;
      if ( obj ) obj = astAnnul( obj );
      ( *nfail )++;
      return NULL;
   }
   if ( strcmp( astGetC( obj, "Class" ), class ) != 0 ) {
      fprintf( stderr, "testloaderhardening: %s loaded as %s\n", what,
               astGetC( obj, "Class" ) );
      obj = astAnnul( obj );
      ( *nfail )++;
      return NULL;
   }
   return obj;
}

/* Check that a dump that should be accepted loads as class "class".
   Returns 0 on success, 1 on failure. */
static int check_accepted( const char *what, const char *dump,
                           const char *class ) {
   int nfail = 0;
   AstObject *obj = load_accepted( what, dump, class, &nfail );
   if ( obj ) obj = astAnnul( obj );
   return nfail;
}

/* Write a valid object to a dump and confirm that the dump loads again as
   the expected class.  The supplied object is not annulled.  Returns 0 on
   success, 1 on failure. */
static int check_valid( const char *what, AstObject *obj, const char *class ) {
   char *dump;
   int nfail;

   if ( !astOK || !obj ) {
      fprintf( stderr, "testloaderhardening: could not build %s\n", what );
      if ( !astOK ) astClearStatus;
      return 1;
   }

   dump = astToString( obj );
   nfail = check_accepted( what, dump, class );
   dump = astFree( dump );
   return nfail;
}

/* Return a copy of "dump" in which the text between "start" and "end"
   (both of which point into "dump") is replaced by "text".  The copy is
   allocated with astMalloc and should be freed by the caller. */
static char *splice( const char *dump, const char *start, const char *end,
                     const char *text ) {
   size_t pre = (size_t) ( start - dump );
   size_t tlen = strlen( text );
   char *out = astMalloc( pre + tlen + strlen( end ) + 1 );
   if ( out ) {
      memcpy( out, dump, pre );
      memcpy( out + pre, text, tlen );
      strcpy( out + pre + tlen, end );
   }
   return out;
}

/* Return a copy of "dump" in which the value of the first item called
   "key" is replaced by "newval".  Returns NULL if the item is not found. */
static char *set_value( const char *dump, const char *key,
                        const char *newval ) {
   char needle[ 64 ];
   const char *pos;
   const char *eol;

   (void) snprintf( needle, sizeof( needle ), " %s = ", key );
   pos = dump ? strstr( dump, needle ) : NULL;
   if ( !pos ) return NULL;
   pos += strlen( needle );
   eol = strchr( pos, '\n' );
   return eol ? splice( dump, pos, eol, newval ) : NULL;
}

/* Return a copy of "dump" in which "lines" is inserted before the final
   "End" line, so that it becomes part of the data of the class to which
   the outermost object belongs.  Returns NULL if there is no such line. */
static char *add_lines( const char *dump, const char *lines ) {
   const char *pos = NULL;
   const char *next = dump ? strstr( dump, " End " ) : NULL;
   while ( next ) {
      pos = next;
      next = strstr( next + 1, " End " );
   }
   return pos ? splice( dump, pos, pos, lines ) : NULL;
}

/* Return a copy of "dump" in which the component object called "key",
   which extends to the first following occurrence of "endline", is
   replaced by "text".  Returns NULL if the object is not found. */
static char *set_object( const char *dump, const char *key,
                         const char *endline, const char *text ) {
   char needle[ 64 ];
   const char *start;
   const char *end;

   (void) snprintf( needle, sizeof( needle ), " %s =\n", key );
   start = dump ? strstr( dump, needle ) : NULL;
   end = start ? strstr( start, endline ) : NULL;
   if ( !end ) return NULL;
   return splice( dump, start, end + strlen( endline ), text );
}

/* Check that a memory allocation function reported AST__NOMEM and returned
   "expected".  Clears the status.  Returns 0 on success, 1 on failure. */
static int check_nomem( const char *what, void *result, void *expected ) {
   int code = astStatus;
   if ( !astOK ) astClearStatus;
   if ( code != AST__NOMEM || result != expected ) {
      fprintf( stderr, "testloaderhardening: %s gave status %d and %s "
               "pointer\n", what, code,
               ( result == expected ) ? "the expected" : "an unexpected" );
      return 1;
   }
   return 0;
}

/* Transform positions with a SplineMap of order 1 in X and order 2 in Y
   that has a single X coefficient, Y knots of 0, 0.5, 1, 1.5 and 2, and
   coefficients of 10, 20 and 30 (first output) and 1, 2 and 3 (second
   output).  The Y knots are all distinct, so the spline is defined only
   between the second and fourth of them, where it interpolates linearly
   between the coefficients.  Positions outside that interval, but inside
   the range of the knots, must give bad values.  Returns the number of
   failures. */
static int check_spline_values( const char *what, AstSplineMap *map ) {
   double xin[ 5 ] = { 0.5, 0.5, 0.5, 0.5, 0.5 };
   double yin[ 5 ] = { 0.25, 0.5, 0.75, 1.5, 1.75 };
   double uexp[ 5 ] = { AST__BAD, 10.0, 15.0, 30.0, AST__BAD };
   double vexp[ 5 ] = { AST__BAD, 1.0, 1.5, 3.0, AST__BAD };
   double uout[ 5 ];
   double vout[ 5 ];
   int i;
   int nfail = 0;

   astTran2( map, 5, xin, yin, 1, uout, vout );
   if ( !astOK ) {
      fprintf( stderr, "testloaderhardening: %s could not be transformed\n",
               what );
      astClearStatus;
      return 1;
   }
   for ( i = 0; i < 5; i++ ) {
      if ( ( uexp[ i ] == AST__BAD ) ? ( uout[ i ] != AST__BAD ||
                                         vout[ i ] != AST__BAD )
                                     : ( fabs( uout[ i ] - uexp[ i ] ) > 1e-12 ||
                                         fabs( vout[ i ] - vexp[ i ] ) > 1e-12 ) ) {
         fprintf( stderr, "testloaderhardening: %s at Y=%g gave (%g,%g), "
                  "expected (%g,%g)\n", what, yin[ i ], uout[ i ], vout[ i ],
                  uexp[ i ], vexp[ i ] );
         nfail++;
      }
   }
   return nfail;
}

/* The start of a Moc dump, up to the data for the Moc class itself. */
#define MOC_HEAD \
   " Begin Moc\n Adapt = 1\n" \
   " Frm =\n Begin SkyFrame\n Naxes = 2\n System = \"ICRS\"\n IsA Frame\n" \
   " End SkyFrame\n RegAxes = 2\n IsA Region\n"

/* The knots and coefficients of a SplineMap dump with two coefficients on
   each axis. */
#define SPLINE_TAIL \
   " TX1 = 0\n TX2 = 0.5\n TX3 = 1\n TY1 = 0\n TY2 = 0.5\n TY3 = 1\n" \
   " CU1 = 0\n CU2 = 1\n CU3 = 2\n CU4 = 3\n" \
   " CV1 = 0\n CV2 = 1\n CV3 = 2\n CV4 = 3\n End SplineMap\n"

int main( void ) {
   int status = SAI__OK;
   int nfail = 0;
   astWatch( &status );

/* PolyMap: a coefficient power must not exceed the maximum power read for
   that axis, or the transformation indexes its work array out of bounds.
   The number of coefficients and the maximum power size arrays, so neither
   can be negative. */
   nfail += check_accepted( "crafted PolyMap",
      " Begin PolyMap\n Nin = 1\n IsA Mapping\n"
      " MPF1 = 2\n NCF1 = 1\n CF1 = 1\n PF1 = 1\n IterInv = 0\n"
      " End PolyMap\n", "PolyMap" );
   nfail += check_rejected( "PolyMap power exceeding max power",
      " Begin PolyMap\n Nin = 1\n IsA Mapping\n"
      " MPF1 = 2\n NCF1 = 1\n CF1 = 1\n PF1 = 1000\n IterInv = 0\n"
      " End PolyMap\n", AST__BADPW );
   nfail += check_rejected( "PolyMap with a negative coefficient count",
      " Begin PolyMap\n Nin = 1\n IsA Mapping\n"
      " MPF1 = 2\n NCF1 = -1\n IterInv = 0\n"
      " End PolyMap\n", AST__LDERR );
   nfail += check_rejected( "PolyMap with a negative max power",
      " Begin PolyMap\n Nin = 1\n IsA Mapping\n"
      " MPF1 = -1\n NCF1 = 0\n IterInv = 0\n"
      " End PolyMap\n", AST__LDERR );
   nfail += check_rejected( "PolyMap with a negative inverse coefficient "
                            "count",
      " Begin PolyMap\n Nin = 1\n IsA Mapping\n"
      " MPF1 = 2\n NCF1 = 1\n CF1 = 1\n PF1 = 1\n"
      " MPI1 = 2\n NCI1 = -1\n IterInv = 0\n"
      " End PolyMap\n", AST__LDERR );
   {
      double cf[] = { 1.0, 1.0, 2.0 };
      AstPolyMap *pm = astPolyMap( 1, 1, 1, cf, 0, NULL, "IterInverse=0" );
      nfail += check_valid( "PolyMap", (AstObject *) pm, "PolyMap" );
      if ( pm ) pm = astAnnul( pm );
   }

/* LutMap: at least two table elements are needed for interpolation, and
   the increment between them cannot be zero. */
   nfail += check_rejected( "LutMap with one element",
      " Begin LutMap\n Nin = 1\n IsA Mapping\n"
      " Nlut = 1\n LutInt = 0\n L1 = 1\n End LutMap\n", AST__LUTIN );
   nfail += check_rejected( "LutMap with zero increment",
      " Begin LutMap\n Nin = 1\n IsA Mapping\n"
      " Nlut = 2\n Incr = 0\n L1 = 1\n L2 = 2\n End LutMap\n", AST__LUTII );

/* The inverse transformation interpolates between two good table values,
   so a dump cannot claim an inverse for a table that has fewer than two.
   Such a table is loaded without an inverse, which must then be refused.
   A table that is not monotonic but has two or more good values keeps the
   inverse that its dump records, since dumps written by old versions of
   AST record one, and that inverse must stay within the table. */
   {
      static const char *dumps[ 4 ] = {
         " Begin LutMap\n Nin = 1\n IsA Mapping\n"
         " Nlut = 2\n L1 = 1\n L2 = 2\n End LutMap\n",
         " Begin LutMap\n Nin = 1\n IsA Mapping\n"
         " Nlut = 3\n L1 = 1\n L2 = 3\n L3 = 2\n End LutMap\n",
         " Begin LutMap\n Nin = 1\n IsA Mapping\n"
         " Nlut = 2\n L1 = 1\n End LutMap\n",
         " Begin LutMap\n Nin = 1\n IsA Mapping\n"
         " Nlut = 2\n End LutMap\n" };
      static const char *names[ 4 ] = { "monotonic LutMap",
                                        "non-monotonic LutMap",
                                        "LutMap with one good value",
                                        "LutMap with no good values" };
      int wanted[ 4 ] = { 1, 1, 0, 0 };
      double xin[ 3 ] = { 0.5, 1.0, 2.5 };
      double xout[ 3 ];
      int i;
      for ( i = 0; i < 4; i++ ) {
         AstObject *m = load_accepted( names[ i ], dumps[ i ], "LutMap",
                                       &nfail );
         if ( m ) {
            if ( astGetI( m, "TranInverse" ) != wanted[ i ] ) {
               fprintf( stderr, "testloaderhardening: %s has TranInverse=%d, "
                        "expected %d\n", names[ i ],
                        astGetI( m, "TranInverse" ), wanted[ i ] );
               nfail++;
            }
            astTran1( m, 3, xin, 0, xout );
            if ( astStatus != ( wanted[ i ] ? SAI__OK : AST__TRNND ) ) {
               fprintf( stderr, "testloaderhardening: the inverse of the %s "
                        "gave status %d\n", names[ i ], astStatus );
               nfail++;
            }
            if ( !astOK ) astClearStatus;
            m = astAnnul( m );
         }
      }
   }
   {
      double lut[ 4 ] = { 1.0, 2.0, 4.0, 8.0 };
      AstLutMap *m = astLutMap( 4, lut, 0.0, 1.0, " " );
      nfail += check_valid( "LutMap", (AstObject *) m, "LutMap" );
      if ( m ) m = astAnnul( m );
   }

/* WcsMap: the longitude and latitude axis indices must be distinct axes of
   the WcsMap, or the transformation and the projection parameter code
   index their per-axis arrays out of bounds.  The numbers of inputs and
   outputs must be equal, since each input axis is copied to the output
   axis with the same index. */
   nfail += check_accepted( "crafted WcsMap",
      " Begin WcsMap\n Nin = 2\n IsA Mapping\n"
      " Type = \"TAN\"\n End WcsMap\n", "WcsMap" );
   nfail += check_rejected( "WcsMap with out-of-range latitude axis",
      " Begin WcsMap\n Nin = 2\n IsA Mapping\n"
      " Type = \"TAN\"\n WcsAx2 = 99\n End WcsMap\n", AST__WCSAX );
   nfail += check_rejected( "WcsMap with equal longitude and latitude axes",
      " Begin WcsMap\n Nin = 2\n IsA Mapping\n"
      " Type = \"TAN\"\n WcsAx1 = 2\n WcsAx2 = 2\n End WcsMap\n",
      AST__WCSAX );
   nfail += check_rejected( "WcsMap with more inputs than outputs",
      " Begin WcsMap\n Nin = 3\n Nout = 2\n IsA Mapping\n"
      " Type = \"TAN\"\n WcsAx1 = 3\n WcsAx2 = 1\n End WcsMap\n",
      AST__LDERR );
   nfail += check_rejected( "WcsMap with more outputs than inputs",
      " Begin WcsMap\n Nin = 2\n Nout = 9\n IsA Mapping\n"
      " Type = \"TAN\"\n End WcsMap\n", AST__LDERR );
   {
      AstWcsMap *m = astWcsMap( 2, AST__SIN, 1, 2, " " );
      nfail += check_valid( "WcsMap", (AstObject *) m, "WcsMap" );
      if ( m ) m = astAnnul( m );
   }

/* The constructor does not check the axis indices of a WcsMap that has no
   projection (type AST__WCSBAD), since it never uses them, so neither can
   the loader. */
   {
      AstWcsMap *m = astWcsMap( 1, AST__WCSBAD, 1, 1, " " );
      nfail += check_valid( "one-axis WcsMap with no projection",
                            (AstObject *) m, "WcsMap" );
      if ( m ) m = astAnnul( m );
      m = astWcsMap( 2, AST__WCSBAD, 1, 1, " " );
      nfail += check_valid( "WcsMap with no projection and equal axes",
                            (AstObject *) m, "WcsMap" );
      if ( m ) m = astAnnul( m );
      m = astWcsMap( 3, AST__WCSBAD, 0, 0, " " );
      nfail += check_valid( "WcsMap with no projection and no axes",
                            (AstObject *) m, "WcsMap" );
      if ( m ) m = astAnnul( m );
   }

/* SplineMap: the number of coefficients must be at least the spline order,
   the knots must be in order, the number of coefficients must be small
   enough to index, and there must be two inputs and two outputs. */
   nfail += check_accepted( "crafted SplineMap",
      " Begin SplineMap\n Nin = 2\n IsA Mapping\n"
      " KX = 1\n KY = 1\n NX = 2\n NY = 2\n" SPLINE_TAIL, "SplineMap" );
   nfail += check_rejected( "SplineMap with too few X coefficients",
      " Begin SplineMap\n Nin = 2\n IsA Mapping\n"
      " KX = 3\n KY = 1\n NX = 2\n NY = 2\n"
      " TX1 = 0\n TX2 = 0.5\n TX3 = 1\n TX4 = 1\n TX5 = 1\n"
      " TY1 = 0\n TY2 = 0.5\n TY3 = 1\n"
      " CU1 = 0\n CU2 = 1\n CU3 = 2\n CU4 = 3\n"
      " CV1 = 0\n CV2 = 1\n CV3 = 2\n CV4 = 3\n End SplineMap\n",
      AST__ATTIN );
   nfail += check_rejected( "SplineMap with non-monotonic X knots",
      " Begin SplineMap\n Nin = 2\n IsA Mapping\n"
      " KX = 1\n KY = 1\n NX = 3\n NY = 3\n"
      " TX1 = 0\n TX2 = 0.9\n TX3 = 0.1\n TX4 = 1\n"
      " TY1 = 0\n TY2 = 0.4\n TY3 = 0.6\n TY4 = 1\n"
      " CU1=0\n CU2=1\n CU3=2\n CU4=3\n CU5=4\n CU6=5\n CU7=6\n CU8=7\n CU9=8\n"
      " CV1=0\n CV2=1\n CV3=2\n CV4=3\n CV5=4\n CV6=5\n CV7=6\n CV8=7\n CV9=8\n"
      " End SplineMap\n", AST__ATTIN );
   nfail += check_rejected( "SplineMap with too many coefficients",
      " Begin SplineMap\n Nin = 2\n IsA Mapping\n"
      " KX = 1\n KY = 1\n NX = 65536\n NY = 65537\n End SplineMap\n",
      AST__ATTIN );
   nfail += check_rejected( "SplineMap with one input",
      " Begin SplineMap\n Nin = 1\n IsA Mapping\n"
      " KX = 1\n KY = 1\n NX = 2\n NY = 2\n" SPLINE_TAIL, AST__LDERR );
   nfail += check_rejected( "SplineMap with one output",
      " Begin SplineMap\n Nin = 2\n Nout = 1\n IsA Mapping\n"
      " KX = 1\n KY = 1\n NX = 2\n NY = 2\n" SPLINE_TAIL, AST__LDERR );
   {
      double tx[ 3 ] = { 0.0, 0.5, 1.0 }, ty[ 3 ] = { 0.0, 0.5, 1.0 };
      double cu[ 4 ] = { 0.0, 1.0, 2.0, 3.0 }, cv[ 4 ] = { 0.0, 1.0, 2.0, 3.0 };
      AstSplineMap *m = astSplineMap( 1, 1, 2, 2, tx, ty, cu, cv, " " );
      nfail += check_valid( "SplineMap", (AstObject *) m, "SplineMap" );
      if ( m ) m = astAnnul( m );
   }

/* A SplineMap in which the first and last knots are not repeated satisfies
   every invariant, but its spline is not defined over the whole range of
   its knots.  Check that it is evaluated only where it is defined, both as
   constructed and as loaded from its dump. */
   {
      double tx[ 2 ] = { 0.0, 1.0 };
      double ty[ 5 ] = { 0.0, 0.5, 1.0, 1.5, 2.0 };
      double cu[ 3 ] = { 10.0, 20.0, 30.0 }, cv[ 3 ] = { 1.0, 2.0, 3.0 };
      AstSplineMap *m = astSplineMap( 1, 2, 1, 3, tx, ty, cu, cv, " " );
      if ( astOK && m ) {
         char *dump = astToString( m );
         AstObject *m2 = load_accepted( "SplineMap with distinct knots",
                                        dump, "SplineMap", &nfail );
         nfail += check_spline_values( "SplineMap with distinct knots", m );
         if ( m2 ) {
            nfail += check_spline_values( "loaded SplineMap with distinct "
                                          "knots", (AstSplineMap *) m2 );
            m2 = astAnnul( m2 );
         }
         dump = astFree( dump );
      } else {
         fprintf( stderr, "testloaderhardening: could not build SplineMap "
                  "with distinct knots\n" );
         if ( !astOK ) astClearStatus;
         nfail++;
      }
      if ( m ) m = astAnnul( m );
   }

/* Plot: the base Frame must not have more axes than the fixed-length
   per-axis arrays in the Plot structure (sized for a Plot3D). */
   nfail += check_rejected( "Plot with a four-axis base Frame",
      " Begin Plot\n Nframe = 2\n Currnt = 2\n Lnk2 = 1\n"
      " Frm1 =\n Begin Frame\n Naxes = 4\n Domain = \"GRAPHICS\"\n End Frame\n"
      " Frm2 =\n Begin Frame\n Naxes = 4\n Domain = \"GRID\"\n End Frame\n"
      " Map2 =\n Begin WinMap\n Nin = 4\n IsA Mapping\n End WinMap\n"
      " IsA FrameSet\n End Plot\n", AST__NAXIN );

/* The numbers of tick mark values size arrays, so cannot be negative.  The
   clipping Frame must be a Frame of the Plot, with one pair of clipping
   bounds for each of its axes.  A full Plot dump is too large to write by
   hand, so build a real one, write it out and add the invalid items to the
   text. */
   {
      float gbox[ 4 ] = { 0.0f, 0.0f, 1.0f, 1.0f };
      double bbox[ 4 ] = { 0.0, 0.0, 1.0, 1.0 };
      double lbnd[ 2 ] = { 0.1, 0.1 };
      double ubnd[ 2 ] = { 0.9, 0.9 };
      AstPlot *p = astPlot( AST__NULL, gbox, bbox, " " );
      char *dump = ( astOK && p ) ? astToString( p ) : NULL;
      if ( dump ) {
         nfail += check_rejected_free( "Plot with a negative number of major "
            "tick values", add_lines( dump, " NMjTk1 = -1\n" ), AST__LDERR );
         nfail += check_rejected_free( "Plot with a negative number of minor "
            "tick values", add_lines( dump, " NMnTk2 = -3\n" ), AST__LDERR );
         nfail += check_rejected_free( "Plot with a non-existent clipping "
            "Frame", add_lines( dump, " ClpAxs = 2\n ClpFrm = 77\n" ),
            AST__LDERR );
         nfail += check_rejected_free( "Plot with too many clipping bounds",
            add_lines( dump, " ClpAxs = 50\n ClpFrm = 2\n" ), AST__LDERR );
         nfail += check_valid( "Plot", (AstObject *) p, "Plot" );
         dump = astFree( dump );

/* Clipping information stored by astClip must still load. */
         astClip( p, AST__CURRENT, lbnd, ubnd );
         nfail += check_valid( "Plot with clipping", (AstObject *) p, "Plot" );
      } else {
         fprintf( stderr, "testloaderhardening: could not build Plot\n" );
         if ( !astOK ) astClearStatus;
         nfail++;
      }
      if ( p ) p = astAnnul( p );
   }

/* Plot3D: each 3D axis is labelled using a nominated axis of one of the
   three encapsulated Plots.  The plane identifier must select one of those
   Plots and the axis index must be a valid Plot axis.  Each of the three
   Plots must be present, and must be a Plot.  Again, build a real Plot3D,
   write it out and corrupt the text. */
   {
      float gbox[ 6 ] = { -1.0f, -1.0f, -1.0f, 1.0f, 1.0f, 1.0f };
      double bbox[ 6 ] = { -1.0, -1.0, -1.0, 1.0, 1.0, 1.0 };
      AstPlot3D *p = astPlot3D( AST__NULL, gbox, bbox, "MinTickLen=0" );
      char *dump = ( astOK && p ) ? astToString( p ) : NULL;
      if ( dump ) {
         nfail += check_rejected_free( "Plot3D labelling an axis with a "
            "non-existent Plot", set_value( dump, "AxPlX1", "9" ),
            AST__LDERR );
         nfail += check_rejected_free( "Plot3D labelling an axis with an "
            "out-of-range Plot axis", set_value( dump, "AxInX1", "99" ),
            AST__LDERR );
         nfail += check_rejected_free( "Plot3D with an invalid base Plot",
            set_value( dump, "BasePl", "7" ), AST__LDERR );
         nfail += check_rejected_free( "Plot3D with no XY Plot",
            set_object( dump, "PlotXY", " End Plot\n", "" ), AST__LDERR );
         nfail += check_rejected_free( "Plot3D with a Frame as its XZ Plot",
            set_object( dump, "PlotXZ", " End Plot\n",
                        " PlotXZ =\n Begin Frame\n Naxes = 2\n End Frame\n" ),
            AST__LDERR );
         nfail += check_rejected_free( "Plot3D with a UnitMap as its YZ Plot",
            set_object( dump, "PlotYZ", " End Plot\n",
                        " PlotYZ =\n Begin UnitMap\n Nin = 2\n IsA Mapping\n"
                        " End UnitMap\n" ), AST__LDERR );
         nfail += check_valid( "Plot3D", (AstObject *) p, "Plot3D" );
         dump = astFree( dump );
      } else {
         fprintf( stderr, "testloaderhardening: could not build Plot3D\n" );
         if ( !astOK ) astClearStatus;
         nfail++;
      }
      if ( p ) p = astAnnul( p );
   }

/* Moc: a Moc that contains any ranges must carry a valid MaxOrder, since
   it bounds the cell indices.  The bounds used here are small enough to be
   accepted if the missing MaxOrder is not noticed.  The number of ranges
   and the lower bound of a range cannot be negative, and every range that
   is claimed must be present. */
   nfail += check_accepted( "crafted Moc",
      MOC_HEAD " MaxOrder = 4\n NRange = 1\n I0 = 80\n K0 = 95\n End Moc\n",
      "Moc" );
   nfail += check_rejected( "Moc with ranges but no MaxOrder",
      MOC_HEAD " NRange = 1\n I0 = 8\n K0 = 9\n End Moc\n", AST__LDERR );
   nfail += check_rejected( "Moc with a negative number of ranges",
      MOC_HEAD " MaxOrder = 4\n NRange = -1\n End Moc\n", AST__LDERR );
   nfail += check_rejected( "Moc with a negative lower bound",
      MOC_HEAD " MaxOrder = 4\n NRange = 1\n I0 = -5\n J0 = -1\n K0 = 95\n"
      " End Moc\n", AST__LDERR );
   nfail += check_rejected( "Moc with fewer ranges than it claims",
      MOC_HEAD " MaxOrder = 4\n NRange = 2147483647\n I0 = 80\n K0 = 95\n"
      " End Moc\n", AST__LDERR );
   {
      AstMoc *m = astMoc( "MaxOrder=4" );
      astAddCell( m, AST__OR, 2, 5 );
      nfail += check_valid( "Moc", (AstObject *) m, "Moc" );
      if ( m ) m = astAnnul( m );
   }

/* Memory allocation: a count read from a dump is turned into a size in
   bytes for the allocation functions.  If that size, or the size including
   the header that precedes each block, is too large to be represented, the
   functions must report an error rather than allocate the smaller block
   that the wrapped size describes. */
   {
      size_t toobig = ( (size_t) -1 )/4;
      void *ptr = astMalloc( 16 );
      nfail += check_nomem( "astCalloc with an overflowing size",
                            astCalloc( toobig, sizeof( double ) ), NULL );
      nfail += check_nomem( "astGrow with an overflowing size",
                            astGrow( NULL, toobig, sizeof( double ) ), NULL );
      nfail += check_nomem( "astMalloc with an overflowing size",
                            astMalloc( (size_t) -8 ), NULL );
      nfail += check_nomem( "astRealloc with an overflowing size",
                            astRealloc( ptr, (size_t) -8 ), ptr );
      ptr = astFree( ptr );
   }

   if ( !astOK ) astClearStatus;

   if ( nfail == 0 ) {
      printf( " All loader hardening tests passed\n" );
      return 0;
   }
   printf( "Loader hardening tests failed (%d)\n", nfail );
   return 1;
}
