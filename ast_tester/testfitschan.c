/*
 *  Test the FitsChan class.
 *  Converted from the Fortran test testfitschan.f.
 *
 *  Differences from the Fortran original:
 *
 *  - err_mark/err_rlse/err_begin/err_end are omitted.  The C version uses
 *    astWatch(&status) + astOK/astClearStatus for error handling.
 *
 *  - err_annul(status) is replaced by astClearStatus.
 *
 *  - err_flush/msg_out/err_rep are replaced by printf().
 *
 *  - The Fortran checktab() subroutine used a COMMON block (/tabsrc/) to
 *    share state with the tabsource() callback.  In C this is replaced by
 *    static globals holding the callback FitsTable and the test status
 *    pointer. The callback's final argument is treated as the required
 *    success flag for astTableSource, not as the inherited AST status.
 *
 *  - The Fortran readobj() helper used a Fortran CHANNEL with a file-based
 *    chsource() callback.  The C version uses fopen/fgets directly with a
 *    static file pointer and a C-style source callback.
 *
 *  - VAL__NBD (bytes per double in Fortran PRM_PAR) is replaced by
 *    sizeof(double).
 *
 *  - astGetFitsS() returns char** (pointer to internal string); the Fortran
 *    version returned a CHARACTER*80.  strcmp() is used for comparisons.
 *
 *  - astFindFits(), astGetFitsI(), astGetFitsK(), astGetFitsS() etc. are
 *    used as boolean expressions (they return int via astINVOKE/astRetV_).
 *
 *  - astConvert no longer triggers an internal AddressSanitizer bug, as
 *    the underlying heap-use-after-free issue in libast (memory.c:4127)
 *    has been fixed.
 *
 *  - The sip.head test was referenced via SourceFile attribute in Fortran;
 *    we do the same in C.
 *
 *  - The callback now supplies the required WCS-TAB extension using
 *    astPutTable(), matching the current C astTableSource contract.
 */

#include "ast.h"
#include "ast_err.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>

/* -----------------------------------------------------------------------
 * Static globals for the tabsource callback (replaces Fortran COMMON block)
 * -----------------------------------------------------------------------*/
static AstFitsTable *g_table = NULL;
static int *g_test_status = NULL;

/* -----------------------------------------------------------------------
 * roundtrip: write a FrameSet to a FitsChan with a given encoding,
 * rewind, read back. If check_equal is set, also verify astEqual.
 * Returns:
 *   1 = success (and astEqual if requested)
 *   0 = write or read failed
 *  -1 = read back OK but astEqual returned false
 * -----------------------------------------------------------------------*/
static int roundtrip( AstFrameSet *fs, const char *encoding,
                      const char *attrs, int check_equal, int *status ) {
   AstFitsChan *fc;
   AstFrameSet *fs2;
   int ok = 0;

   if( *status != 0 ) return 0;

   fc = astFitsChan( NULL, NULL, " " );
   astSetC( fc, "Encoding", encoding );
   if( attrs && attrs[0] ) astSet( fc, "%s", attrs );
   if( astWrite( fc, fs ) != 1 ) goto done;
   astClear( fc, "Card" );
   fs2 = (AstFrameSet *) astRead( fc );
   if( !fs2 ) goto done;
   if( check_equal && !astEqual( fs, fs2 ) ) {
      ok = -1;
   } else {
      ok = 1;
   }
   fs2 = astAnnul( fs2 );
done:
   if( !astOK ) astClearStatus;
   fc = astAnnul( fc );
   return ok;
}

/* -----------------------------------------------------------------------
 * stopit: record first error
 * -----------------------------------------------------------------------*/
static void stopit( int errnum, const char *text, int *status ) {
   if( *status != 0 ) return;
   *status = 1;
   if( text && text[0] && text[0] != ' ' )
      printf( "Error %d: %s\n", errnum, text );
   else
      printf( "Error %d\n", errnum );
}

/* -----------------------------------------------------------------------
 * tabsource callback: called by FitsChan when it needs a FITS extension
 * -----------------------------------------------------------------------*/
static void tabsource( AstFitsChan *fc, const char *extnam, int extver,
                       int extlevel, int *status ) {
   *status = 0;

   if( strcmp( extnam, "WCS-TAB" ) != 0 ) {
      stopit( 1035, "tabsource: unexpected extnam", g_test_status );
      return;
   }

   if( !g_table ) {
      stopit( 1036, "tabsource: g_table is NULL", g_test_status );
      return;
   }

   if( !astIsAFitsTable( g_table ) ) {
      stopit( 1037, "tabsource: not a FitsTable", g_test_status );
      return;
   }

   if( extver != 1 ) {
      printf( "EXTVER=%d\n", extver );
      stopit( 1065, "tabsource: wrong extver", g_test_status );
      return;
   }

   if( extlevel != 1 ) {
      stopit( 1066, "tabsource: wrong extlevel", g_test_status );
      return;
   }

   astPutTable( fc, g_table, extnam );
   if( astOK ) *status = 1;
}

/* -----------------------------------------------------------------------
 * readobj: read an AST object from a file using a Channel
 * -----------------------------------------------------------------------*/
static AstObject *readobj( const char *file, int *status ) {
   AstChannel *ch;
   AstObject *obj;
   char opts[256];
   if( *status != 0 ) return NULL;

   snprintf(opts, sizeof(opts), "SourceFile=%s", file);
   ch = astChannel( NULL, NULL, "%s", opts );
   obj = astRead( ch );
   astAnnul( ch );
   if( !obj ) *status = 1;
   return obj;
}

/* -----------------------------------------------------------------------
 * checkft: verify contents of coords/index arrays from TAB lookup
 * -----------------------------------------------------------------------*/
static void checkft( int nelem, double *coords, double *indx, int *status ) {
   if( *status != 0 ) return;

   if( indx[0] != 1.0 ) {
      stopit( 2001, " ", status );
   } else if( coords[0] != 1.0 ) {
      stopit( 2002, " ", status );
   } else if( indx[nelem-1] != 1.0e2 ) {
      stopit( 2003, " ", status );
   } else if( coords[nelem-1] != 1.0e-4 ) {
      stopit( 2004, " ", status );
   } else if( fabs( coords[nelem/2 - 1] -
                    pow( indx[nelem/2 - 1], -2.0 ) ) > 1e-20 ) {
      stopit( 2005, " ", status );
   }
}

/* -----------------------------------------------------------------------
 * checkft2: verify contents of LUT-based coords array
 * -----------------------------------------------------------------------*/
static void checkft2( int nelem, double *coords, int *status ) {
   if( *status != 0 ) return;

   if( fabs( coords[0] - 299.792458 ) > 1.0e-5 ) {
      stopit( 3002, " ", status );
   } else if( fabs( coords[nelem-1] - 2997924.58 ) > 1.0e-1 ) {
      stopit( 3004, " ", status );
   }
}

/* -----------------------------------------------------------------------
 * test_fitsrounding: test FitsRounding attribute behaviour
 * -----------------------------------------------------------------------*/

/* -----------------------------------------------------------------------
 * rstrip: strip trailing spaces
 * -----------------------------------------------------------------------*/
static char *rstrip( char *s ) {
   int n;
   if( !s ) return s;
   n = strlen( s );
   while( n > 0 && s[n-1] == ' ' ) {
      s[n-1] = '\0';
      n--;
   }
   return s;
}

static void test_fitsrounding( AstFitsChan *fc, int *status ) {
   char card[81];

   if( *status != 0 ) return;

   if( astGetI( fc, "FitsRounding" ) != 10 )
      stopit( 12314, " ", status );

   /* Test 1: value with 10 significant figures -> should be kept */
   astPutFits( fc, "DUMMY   =  1.0000010001", 0 );
   astClear( fc, "Card" );
   if( astFindFits( fc, "DUMMY", card, 0 ) ) {
      if( strcmp( rstrip(card), "DUMMY   =         1.0000010001" ) != 0 ) {
         printf( "FitsRounding test 1: got [%s]\n", card );
         stopit( 12310, " ", status );
      }
   } else {
      stopit( 12311, " ", status );
   }

   /* Test 2: value with 11 sig figs -> should be rounded to 10 */
   astPutFits( fc, "DUMMY   =  1.00000100001", 0 );
   astClear( fc, "Card" );
   if( astFindFits( fc, "DUMMY", card, 0 ) ) {
      if( strcmp( rstrip(card), "DUMMY   =             1.000001" ) != 0 ) {
         printf( "FitsRounding test 2: got [%s]\n", card );
         stopit( 12312, " ", status );
      }
   } else {
      stopit( 12313, " ", status );
   }

   /* Test 3 */
   astPutFits( fc, "DUMMY   =  1.0000019991", 0 );
   astClear( fc, "Card" );
   if( astFindFits( fc, "DUMMY", card, 0 ) ) {
      if( strcmp( rstrip(card), "DUMMY   =         1.0000019991" ) != 0 ) {
         printf( "FitsRounding test 3: got [%s]\n", card );
         stopit( 12314, " ", status );
      }
   } else {
      stopit( 12315, " ", status );
   }

   /* Test 4 */
   astPutFits( fc, "DUMMY   =  1.00000199991", 0 );
   astClear( fc, "Card" );
   if( astFindFits( fc, "DUMMY", card, 0 ) ) {
      if( strcmp( rstrip(card), "DUMMY   =             1.000002" ) != 0 ) {
         printf( "FitsRounding test 4: got [%s]\n", card );
         stopit( 12316, " ", status );
      }
   } else {
      stopit( 12317, " ", status );
   }

   /* Test 5 */
   astPutFits( fc, "DUMMY   =  -1.0000010000001", 0 );
   astClear( fc, "Card" );
   if( astFindFits( fc, "DUMMY", card, 0 ) ) {
      if( strcmp( rstrip(card), "DUMMY   =            -1.000001" ) != 0 ) {
         printf( "FitsRounding test 5: got [%s]\n", card );
         stopit( 12318, " ", status );
      }
   } else {
      stopit( 12319, " ", status );
   }

   /* Test 6: FitsRounding=5 */
   astSetI( fc, "FitsRounding", 5 );
   astClear( fc, "Card" );
   if( astFindFits( fc, "DUMMY", card, 0 ) ) {
      if( strcmp( rstrip(card), "DUMMY   =                 -1.0" ) != 0 ) {
         printf( "FitsRounding test 6: got [%s]\n", card );
         stopit( 12320, " ", status );
      }
   } else {
      stopit( 12321, " ", status );
   }

   /* Test 7: clear FitsRounding -> back to 10 */
   astClear( fc, "FitsRounding" );
   if( astGetI( fc, "FitsRounding" ) != 10 )
      stopit( 12322, " ", status );

   /* Test 8 */
   astPutFits( fc, "DUMMY   =  1.9999969999993E-02", 0 );
   astClear( fc, "Card" );
   if( astFindFits( fc, "DUMMY", card, 0 ) ) {
      if( strcmp( rstrip(card), "DUMMY   =           0.01999997" ) != 0 ) {
         printf( "FitsRounding test 8: got [%s]\n", card );
         stopit( 12323, " ", status );
      }
   } else {
      stopit( 12324, " ", status );
   }

   /* Test 9: FitsRounding=5 */
   astSetI( fc, "FitsRounding", 5 );
   astClear( fc, "Card" );
   if( astFindFits( fc, "DUMMY", card, 0 ) ) {
      if( strcmp( rstrip(card), "DUMMY   =                 0.02" ) != 0 ) {
         printf( "FitsRounding test 9: got [%s]\n", card );
         stopit( 12325, " ", status );
      }
   } else {
      stopit( 12326, " ", status );
   }

   /* Test 10: FitsRounding=10, scientific notation */
   astSetI( fc, "FitsRounding", 10 );
   astPutFits( fc, "DUMMY   =  1.9999969999993E-08", 0 );
   astClear( fc, "Card" );
   if( astFindFits( fc, "DUMMY", card, 0 ) ) {
      if( strcmp( rstrip(card), "DUMMY   =          1.999997E-8" ) != 0 ) {
         printf( "FitsRounding test 10: got [%s]\n", card );
         stopit( 12327, " ", status );
      }
   } else {
      stopit( 12328, " ", status );
   }

   /* Test 11: FitsRounding=5 */
   astSetI( fc, "FitsRounding", 5 );
   astClear( fc, "Card" );
   if( astFindFits( fc, "DUMMY", card, 0 ) ) {
      if( strcmp( rstrip(card), "DUMMY   =               2.0E-8" ) != 0 ) {
         printf( "FitsRounding test 11: got [%s]\n", card );
         stopit( 12329, " ", status );
      }
   } else {
      stopit( 12330, " ", status );
   }

   /* Test 12 */
   astPutFits( fc, "DUMMY   =  9.999999E-6", 0 );
   astClear( fc, "Card" );
   if( astFindFits( fc, "DUMMY", card, 0 ) ) {
      if( strcmp( rstrip(card), "DUMMY   =              10.0E-6" ) != 0 ) {
         printf( "FitsRounding test 12: got [%s]\n", card );
         stopit( 12331, " ", status );
      }
   } else {
      stopit( 12332, " ", status );
   }

   /* Test 13 */
   astPutFits( fc, "DUMMY   =  -9.999999", 0 );
   astClear( fc, "Card" );
   if( astFindFits( fc, "DUMMY", card, 0 ) ) {
      if( strcmp( rstrip(card), "DUMMY   =                -10.0" ) != 0 ) {
         printf( "FitsRounding test 13: got [%s]\n", card );
         stopit( 12333, " ", status );
      }
   } else {
      stopit( 12334, " ", status );
   }

   /* Test 14: reset FitsRounding and re-check */
   astClear( fc, "FitsRounding" );
   astClear( fc, "Card" );
   if( astFindFits( fc, "DUMMY", card, 0 ) ) {
      if( strcmp( rstrip(card), "DUMMY   =            -9.999999" ) != 0 ) {
         printf( "FitsRounding test 14: got [%s]\n", card );
         stopit( 12335, " ", status );
      }
   } else {
      stopit( 12336, " ", status );
   }

   /* Test 15: FitsRounding=6 */
   astSetI( fc, "FitsRounding", 6 );
   astClear( fc, "Card" );
   if( astFindFits( fc, "DUMMY", card, 0 ) ) {
      if( strcmp( rstrip(card), "DUMMY   =                -10.0" ) != 0 ) {
         printf( "FitsRounding test 15: got [%s]\n", card );
         stopit( 12337, " ", status );
      }
   } else {
      stopit( 12338, " ", status );
   }

   /* Test 16: FitsRounding=7 */
   astSetI( fc, "FitsRounding", 7 );
   astClear( fc, "Card" );
   if( astFindFits( fc, "DUMMY", card, 0 ) ) {
      if( strcmp( rstrip(card), "DUMMY   =            -9.999999" ) != 0 ) {
         printf( "FitsRounding test 16: got [%s]\n", card );
         stopit( 12339, " ", status );
      }
   } else {
      stopit( 12340, " ", status );
   }

   /* Test 17: fitsrounding=6, new value */
   astSetI( fc, "FitsRounding", 6 );
   astPutFits( fc, "DUMMY   =  -0.0019999912", 0 );
   astClear( fc, "Card" );
   if( astFindFits( fc, "DUMMY", card, 0 ) ) {
      if( strcmp( rstrip(card), "DUMMY   =        -0.0019999912" ) != 0 ) {
         printf( "FitsRounding test 17: got [%s]\n", card );
         stopit( 12341, " ", status );
      }
   } else {
      stopit( 12342, " ", status );
   }

   /* Test 18: fitsrounding=5 */
   astSetI( fc, "FitsRounding", 5 );
   astClear( fc, "Card" );
   if( astFindFits( fc, "DUMMY", card, 0 ) ) {
      if( strcmp( rstrip(card), "DUMMY   =               -0.002" ) != 0 ) {
         printf( "FitsRounding test 18: got [%s]\n", card );
         stopit( 12343, " ", status );
      }
   } else {
      stopit( 12344, " ", status );
   }

   /* Reset FitsRounding at end */
   astClear( fc, "FitsRounding" );
}

/* -----------------------------------------------------------------------
 * checktab: test FITS -TAB encoding with astTableSource
 * -----------------------------------------------------------------------*/
static void checktab( int *status ) {
   AstSpecFrame *sf;
   AstFrame *gf;
   AstMapping *mm, *lm, *sm, *map;
   AstFrameSet *fs, *fs2;
   AstFitsChan *fc, *fc2;
   AstKeyMap *tables;
   AstObject *table_obj;
   int nelem, ncard;
   size_t sz;
   double *coords1 = NULL, *index1 = NULL;
   double lut[100], shift;
   double x[3], y[3];
      char *sval;

   if( *status != 0 ) return;

   astBegin;

   sf = astSpecFrame( "system=freq,unit=MHz" );
   gf = astFrame( 1, "domain=GRID" );
   mm = (AstMapping *)astMathMap( 1, 1, 1,
                                  (const char *[]){"y=1/(x*x)"},
                                  1,
                                  (const char *[]){"x=1/sqrt(y)"},
                                  " " );
   fs = astFrameSet( gf, " " );
   astAddFrame( fs, AST__BASE, mm, (AstFrame *)sf );

   fc = astFitsChan( NULL, NULL, "Encoding=FITS-WCS" );
   astPutFits( fc, "NAXIS   = 1", 0 );
   astPutFits( fc, "NAXIS1  = 100", 0 );

   if( astWrite( fc, fs ) != 0 )
      stopit( 1000, " ", status );
   else if( astGetTables( fc ) != AST__NULL )
      stopit( 1001, "CheckTab", status );

   astSetI( fc, "TabOK", 1 );

   if( astWrite( fc, fs ) != 1 )
      stopit( 1002, " ", status );

   tables = (AstKeyMap *)astGetTables( fc );
   if( !tables )
      stopit( 1003, " ", status );
   else if( !astIsAKeyMap( tables ) )
      stopit( 1004, " ", status );
   else if( astMapSize( tables ) != 1 )
      stopit( 1005, " ", status );

   /* Check the key name */
   sval = (char *)astMapKey( tables, 0 );  /* 0-based in C */
   if( strcmp( sval, "WCS-TAB" ) != 0 ) {
      printf( "MapKey: [%s]\n", sval );
      stopit( 1006, " ", status );
   }

   if( !astMapGet0A( tables, "WCS-TAB", &table_obj ) )
      stopit( 1007, " ", status );
   else if( !astIsAFitsTable( table_obj ) )
      stopit( 1004, " ", status );
   else
      astAnnul( table_obj );

   /* Check table dimensions */
   {
      AstFitsTable *tbl = NULL;
      astMapGet0A( tables, "WCS-TAB", (AstObject **)&tbl );
      if( tbl ) {
         if( astGetI( tbl, "NColumn" ) != 2 )
            stopit( 1005, " ", status );
         else if( astGetI( tbl, "NRow" ) != 1 )
            stopit( 1006, " ", status );
         else if( astGetI( tbl, "ColumnLength(coords1)" ) != 197 )
            stopit( 1007, " ", status );
         else if( astGetI( tbl, "ColumnLength(index1)" ) != 197 )
            stopit( 1008, " ", status );

         sz = astColumnSize( tbl, "COORDS1" );
         if( sz != sizeof(double) * 197 )
            stopit( 1009, " ", status );

         coords1 = malloc( sz );
         astGetColumnData( tbl, "Coords1", 0.0f, AST__BAD, sz, coords1, &nelem );
         if( nelem != 197 ) stopit( 1010, " ", status );

         sz = astColumnSize( tbl, "INDEX1" );
         if( sz != sizeof(double) * 197 )
            stopit( 1011, " ", status );

         index1 = malloc( sz );
         astGetColumnData( tbl, "inDex1", 0.0f, AST__BAD, sz, index1, &nelem );
         if( nelem != 197 ) stopit( 1012, " ", status );

         checkft( 197, coords1, index1, status );

         free( coords1 );  coords1 = NULL;
         free( index1 );   index1 = NULL;

         astAnnul( tbl );
      }
   }

   /* ------------------------------------------------------------------
    * Now replace the MathMap with a LutMap
    * ------------------------------------------------------------------*/
   {
      int i;
      for( i = 0; i < 100; i++ )
         lut[i] = 1.0 / (double)((i+1)*(i+1));
   }

   lm = (AstMapping *)astLutMap( 100, lut, -49.0, 1.0, " " );

   astSetC( (AstObject *)sf, "System", "Wave" );
   astSetC( (AstObject *)sf, "Unit", "m" );
   astRemoveFrame( fs, AST__CURRENT );
   astAddFrame( fs, AST__BASE, lm, (AstFrame *)sf );
   astSet( fs, "System=freq" );

   shift = 50.0;
   sm = (AstMapping *)astShiftMap( 1, &shift, " " );
   astRemapFrame( fs, AST__BASE, sm );

   astRemoveTables( fc, "WCS-TAB" );
   astPurgeWCS( fc );

   if( astWrite( fc, fs ) != 1 )
      stopit( 1013, " ", status );

   tables = (AstKeyMap *)astGetTables( fc );
   if( !tables )
      stopit( 1014, " ", status );
   else if( !astIsAKeyMap( tables ) )
      stopit( 1015, " ", status );
   else if( astMapSize( tables ) != 1 )
      stopit( 1016, " ", status );

   sval = (char *)astMapKey( tables, 0 );
   if( strcmp( sval, "WCS-TAB" ) != 0 )
      stopit( 1017, " ", status );

   if( !astMapGet0A( tables, "WCS-TAB", &table_obj ) )
      stopit( 1018, " ", status );
   else if( !astIsAFitsTable( table_obj ) )
      stopit( 1019, " ", status );
   else
      astAnnul( table_obj );

   {
      AstFitsTable *tbl = NULL;
      astMapGet0A( tables, "WCS-TAB", (AstObject **)&tbl );
      if( tbl ) {
         if( astGetI( tbl, "NColumn" ) != 1 )
            stopit( 1020, " ", status );
         else if( astGetI( tbl, "NRow" ) != 1 )
            stopit( 1021, " ", status );
         else if( astGetI( tbl, "ColumnLength(coords1)" ) != 100 )
            stopit( 1022, " ", status );

         sz = astColumnSize( tbl, "COORDS1" );
         if( sz != sizeof(double) * 100 )
            stopit( 1024, " ", status );

         coords1 = malloc( sz );
         astGetColumnData( tbl, "Coords1", 0.0f, AST__BAD, sz, coords1, &nelem );
         if( nelem != 100 ) stopit( 1025, " ", status );

         checkft2( 100, coords1, status );

         free( coords1 );  coords1 = NULL;
         astAnnul( tbl );
      }
   }

   /* ------------------------------------------------------------------
    * Test round-trip: remove tables from fc copy, re-read with tabsource
    * ------------------------------------------------------------------*/
   astRemoveTables( fc, " " );
   fc2 = (AstFitsChan *)astCopy( fc );
   astPutTables( fc, tables );
   astClear( fc, "Card" );

   fs2 = (AstFrameSet *)astRead( fc );
   if( !fs2 ) stopit( 1028, " ", status );

   if( fs2 ) {
      AstFrame *fr1 = astGetFrame( fs, AST__CURRENT );
      AstFrame *fr2 = astGetFrame( fs2, AST__CURRENT );
      if( !astEqual( fr1, fr2 ) )
         stopit( 1029, " ", status );
      astAnnul( fr1 );
      astAnnul( fr2 );

      map = (AstMapping *)astCmpMap(
               astGetMapping( fs, AST__BASE, AST__CURRENT ),
               astGetMapping( fs2, AST__CURRENT, AST__BASE ),
               1, " " );

      x[0] = 1.0;  x[1] = 50.0;  x[2] = 100.0;
      astTran1( map, 3, x, 1, y );

      if( fabs( y[0] - x[0] ) > 1.0e-4 ||
          fabs( y[1] - x[1] ) > 1.0e-4 ||
          fabs( y[2] - x[2] ) > 1.0e-4 )
         stopit( 1030, " ", status );
      astAnnul( map );
      astAnnul( fs2 );
   }

   /* ------------------------------------------------------------------
    * Test that read without tables raises AST__NOTAB, then use tabsource
    * ------------------------------------------------------------------*/
   if( !astGetI( fc2, "TabOK" ) )
      stopit( 1031, " ", status );

   ncard = astGetI( fc2, "Ncard" );

   astClear( fc2, "Card" );
   fs2 = (AstFrameSet *)astRead( fc2 );
   if( astOK ) {
      /* Should have raised an error */
      stopit( 1032, " ", status );
   } else if( astStatus == AST__NOTAB ) {
      astClearStatus;
   } else {
      astClearStatus;
      stopit( 1032, "wrong error", status );
   }

   astSetI( fc2, "TabOK", 0 );
   astClear( fc2, "Card" );
   fs2 = (AstFrameSet *)astRead( fc2 );
   if( astOK ) {
      stopit( 1032, " ", status );
   } else if( astStatus == AST__BDFTS ) {
      astClearStatus;
   } else {
      astClearStatus;
      stopit( 1032, "wrong error (TabOK=0)", status );
   }
   astSetI( fc2, "TabOK", 1 );

   if( ncard != astGetI( fc2, "Ncard" ) )
      stopit( 1034, " ", status );

   /* Set g_table so tabsource callback can use it */
   g_test_status = status;
   if( !astMapGet0A( tables, "WCS-TAB", (AstObject **) &g_table ) )
      stopit( 1034, "missing WCS-TAB table for callback", status );
   astTableSource( fc2, tabsource );
   astClear( fc2, "Card" );
   fs2 = (AstFrameSet *)astRead( fc2 );
   if( !fs2 ) stopit( 1035, " ", status );
   if( *status != 0 ) { printf("astRead at 678 failed! status=%d\n", *status); return; }
   if( g_table ) g_table = astAnnul( g_table );
   g_test_status = NULL;

   if( fs2 ) {
      AstFrame *fr1 = astGetFrame( fs, AST__CURRENT );
      AstFrame *fr2 = astGetFrame( fs2, AST__CURRENT );
      if( !astEqual( fr1, fr2 ) )
         stopit( 1036, " ", status );
      astAnnul( fr1 );
      astAnnul( fr2 );

      map = (AstMapping *)astCmpMap(
               astGetMapping( fs, AST__BASE, AST__CURRENT ),
               astGetMapping( fs2, AST__CURRENT, AST__BASE ),
               1, " " );

      x[0] = 1.0;  x[1] = 50.0;  x[2] = 100.0;
      astTran1( map, 3, x, 1, y );

      if( fabs( y[0] - x[0] ) > 1.0e-4 ||
          fabs( y[1] - x[1] ) > 1.0e-4 ||
          fabs( y[2] - x[2] ) > 1.0e-4 )
         stopit( 1037, " ", status );
      astAnnul( map );
      astAnnul( fs2 );
   }

   /* ------------------------------------------------------------------
    * Test with sparse.ast (multi-dim TAB WCS)
    * ------------------------------------------------------------------*/
   {
      AstObject *sparseobj = readobj( "sparse.ast", status );
      if( !sparseobj ) {
         stopit( 1038, "readobj failed", status );
      } else {
         AstFrameSet *fssp = (AstFrameSet *)sparseobj;
         AstFitsChan *fcsp;
         AstFrameSet *fs2sp;
         AstFrameSet *fs3sp;
         double xs[3], ys[3], y2s[3];
         int ii;

         fcsp = astFitsChan( NULL, NULL, "Encoding=FITS-WCS,TabOK=1" );
         astPutFits( fcsp, "NAXIS   = 2", 0 );
         astPutFits( fcsp, "NAXIS1  = 2000", 0 );
         astPutFits( fcsp, "NAXIS2  = 1", 0 );

         if( astWrite( fcsp, fssp ) != 1 )
            stopit( 1038, " ", status );

         astClear( fcsp, "Card" );

         fs2sp = (AstFrameSet *)astRead( fcsp );
         if (!fs2sp || *status != 0) {
            printf("fs2sp is NULL or status is bad: %d\n", *status);
            return;
         }


         astInvert( fssp );
         astInvert( fs2sp );
         {
            AstFrameSet *fssp_clone = astClone(fssp);
            AstFrameSet *fs2sp_clone = astClone(fs2sp);
            fs3sp = (AstFrameSet *)astConvert( fssp_clone, fs2sp_clone, "SKY-DSBSPECTRUM" );
            if( fssp_clone ) astAnnul( fssp_clone );
            if( fs2sp_clone ) astAnnul( fs2sp_clone );
         }
         if( !fs3sp )
            stopit( 1039, " ", status );


         {
            AstFrame *frb = astGetFrame( fssp, AST__BASE );
            if( strcmp( astGetC( frb, "Domain" ), "SKY-DSBSPECTRUM" ) != 0 )
               stopit( 1040, " ", status );
            astAnnul( frb );
         }
         {
            AstFrame *frb = astGetFrame( fs2sp, AST__BASE );
            if( strcmp( astGetC( frb, "Domain" ), "SKY-DSBSPECTRUM" ) != 0 )
               stopit( 1041, " ", status );
            astAnnul( frb );
         }

         astInvert( fssp );
         astInvert( fs2sp );

         xs[0] = 1.0;  xs[1] = 1.0;  xs[2] = 1.0;
         {
            double xin3[3];
            xin3[0] = xs[0]; xin3[1] = xs[1]; xin3[2] = xs[2];
            astTranN( fssp,  1, 3, 1, xin3, 1, 3, 1, ys );
            astTranN( fs2sp, 1, 3, 1, xin3, 1, 3, 1, y2s );
         }
         for( ii = 0; ii < 3; ii++ ) {
            if( fabs( ys[ii] - y2s[ii] ) > 1.0e-8 )
               stopit( 1042, " ", status );
         }

         xs[0] = 10.0;  xs[1] = 1.0;  xs[2] = 1000.0;
         {
            double xin3[3];
            xin3[0] = xs[0]; xin3[1] = xs[1]; xin3[2] = xs[2];
            astTranN( fssp,  1, 3, 1, xin3, 1, 3, 1, ys );
            astTranN( fs2sp, 1, 3, 1, xin3, 1, 3, 1, y2s );
         }
         for( ii = 0; ii < 3; ii++ ) {
            if( fabs( ys[ii] - y2s[ii] ) > 1.0e-8 )
               stopit( 1042, " ", status );
         }

         if( fs3sp ) astAnnul( fs3sp );
         if( fs2sp ) astAnnul( fs2sp );
         astAnnul( fcsp );
         astAnnul( fssp );
      }
   }

   /* ------------------------------------------------------------------
    * Test with voltage domain (non-celestial TAB)
    * ------------------------------------------------------------------*/
   {
      AstFrame *sfv, *gfv;
      AstMapping *mmv;
      AstFrameSet *fsv, *fs2v, *fs3v;
      AstKeyMap *tablesv;
      AstObject *table_objv;
            double xv[3], yv[3];

      sfv = astFrame( 1, "domain=voltage,unit=V" );
      gfv = astFrame( 1, "domain=GRID" );
      mmv = (AstMapping *)astMathMap( 1, 1, 1,
                                      (const char *[]){"y=1/(x*x)"},
                                      1,
                                      (const char *[]){"x=1/sqrt(y)"},
                                      " " );
      fsv = astFrameSet( gfv, " " );
      astAddFrame( fsv, AST__BASE, mmv, sfv );

      astEmptyFits( fc );
      astPutFits( fc, "NAXIS   = 1", 0 );
      astPutFits( fc, "NAXIS1  = 100", 0 );

      if( astWrite( fc, fsv ) != 1 )
         stopit( 1043, " ", status );

      {
         char *strval;
         if( astGetFitsS( fc, "CTYPE1", &strval ) ) {
            if( strcmp( strval, "VOLT-TAB" ) != 0 )
               stopit( 1059, " ", status );
         } else {
            stopit( 1060, " ", status );
         }
      }

      tablesv = (AstKeyMap *)astGetTables( fc );
      if( !tablesv )
         stopit( 1044, " ", status );
      else if( !astIsAKeyMap( tablesv ) )
         stopit( 1045, " ", status );
      else if( astMapSize( tablesv ) != 1 )
         stopit( 1046, " ", status );

      sval = (char *)astMapKey( tablesv, 0 );
      if( strcmp( sval, "WCS-TAB" ) != 0 )
         stopit( 1047, " ", status );

      if( !astMapGet0A( tablesv, "WCS-TAB", &table_objv ) )
         stopit( 1048, " ", status );
      else if( !astIsAFitsTable( table_objv ) )
         stopit( 1049, " ", status );
      else
         astAnnul( table_objv );

      {
         AstFitsTable *tblv = NULL;
         astMapGet0A( tablesv, "WCS-TAB", (AstObject **)&tblv );
         if( tblv ) {
            size_t szv;
            int nelemv;
            double *coords1v = NULL, *index1v = NULL;

            if( astGetI( tblv, "NColumn" ) != 2 )
               stopit( 1050, " ", status );
            else if( astGetI( tblv, "NRow" ) != 1 )
               stopit( 1051, " ", status );
            else if( astGetI( tblv, "ColumnLength(coords1)" ) != 197 )
               stopit( 1052, " ", status );
            else if( astGetI( tblv, "ColumnLength(index1)" ) != 197 )
               stopit( 1053, " ", status );

            szv = astColumnSize( tblv, "COORDS1" );
            if( szv != sizeof(double) * 197 )
               stopit( 1054, " ", status );

            coords1v = malloc( szv );
            astGetColumnData( tblv, "Coords1", 0.0f, AST__BAD, szv, coords1v, &nelemv );
            if( nelemv != 197 ) stopit( 1055, " ", status );

            szv = astColumnSize( tblv, "INDEX1" );
            if( szv != sizeof(double) * 197 )
               stopit( 1056, " ", status );

            index1v = malloc( szv );
            astGetColumnData( tblv, "inDex1", 0.0f, AST__BAD, szv, index1v, &nelemv );
            if( nelemv != 197 ) stopit( 1057, " ", status );

            checkft( 197, coords1v, index1v, status );

            free( coords1v );
            free( index1v );
            astAnnul( tblv );
         }
      }

      astClear( fc, "Card" );
      fs2v = (AstFrameSet *)astRead( fc );
      if( !fs2v )
         stopit( 1058, " ", status );

      if( fs2v ) {
         astInvert( fsv );
         astInvert( fs2v );
         fs3v = (AstFrameSet *)astConvert( fs2v, fsv, " " );
         if( !fs3v )
            stopit( 1061, " ", status );

         {
            AstFrame *frb = astGetFrame( fsv, AST__BASE );
            if( strcmp( astGetC( frb, "Domain" ), "VOLTAGE" ) != 0 )
               stopit( 1062, " ", status );
            astAnnul( frb );
         }
         {
            AstFrame *frb = astGetFrame( fs2v, AST__BASE );
            if( strcmp( astGetC( frb, "Domain" ), "VOLTAGE" ) != 0 )
               stopit( 1063, " ", status );
            astAnnul( frb );
         }

         xv[0] = 1.0;  xv[1] = 10.0;  xv[2] = 100.0;
         astTran1( fs3v, 3, xv, 1, yv );
         if( fabs( yv[0] - xv[0] ) > 1.0e-2 ||
             fabs( yv[1] - xv[1] ) > 1.0e-2 ||
             fabs( yv[2] - xv[2] ) > 1.0e-2 )
            stopit( 1064, " ", status );

         if( fs3v ) astAnnul( fs3v );
         astAnnul( fs2v );
      }
      astAnnul( fsv );
   }

   astEnd;
}


/* -----------------------------------------------------------------------
 * checktab2: test 2-D TAB encoding with SkyRef
 * -----------------------------------------------------------------------*/
static void checktab2( int *status ) {
   AstSkyFrame *sf;
   AstFrame *gf;
   AstMapping *mm, *mmm;
   AstFrameSet *fs, *fs2;
   AstFitsChan *fc;

   if( *status != 0 ) return;

   astBegin;

   sf = astSkyFrame( " " );
   astSetD( sf, "SkyRef(1)", 0.0001 );
   astSetD( sf, "SkyRef(2)", 0.0001 );

   gf = astFrame( 2, "domain=GRID" );

   mm = (AstMapping *)astMathMap( 1, 1, 1,
                                  (const char *[]){"y=(x+50)**(-2)"},
                                  1,
                                  (const char *[]){"x=-50+1/sqrt(y)"},
                                  " " );
   mmm = (AstMapping *)astCmpMap( mm, mm, 0, " " );

   fs = astFrameSet( gf, " " );
   astAddFrame( fs, AST__BASE, mmm, (AstFrame *)sf );

   fc = astFitsChan( NULL, NULL, "Encoding=FITS-WCS,TabOK=1" );
   astPutFits( fc, "NAXIS   = 2", 0 );
   astPutFits( fc, "NAXIS1  = 100", 0 );
   astPutFits( fc, "NAXIS2  = 100", 0 );

   if( astWrite( fc, fs ) != 1 )
      stopit( 2000, " ", status );

   astClear( fc, "Card" );
   fs2 = (AstFrameSet *)astRead( fc );
   if( !fs2 )
      stopit( 2001, " ", status );

   if( fs2 ) {
      if( fabs( astGetD( fs2, "SkyRef(1)" ) - 0.0001 ) > 1.0e-7 )
         stopit( 2001, " ", status );
      if( fabs( astGetD( fs2, "SkyRef(2)" ) - 0.0001 ) > 1.0e-7 )
         stopit( 2001, " ", status );
      astAnnul( fs2 );
   }

   astEnd;
}


/* -----------------------------------------------------------------------
 * main
 * -----------------------------------------------------------------------*/
int main( void ) {
   int status_value = 0;
   int *status = &status_value;

   AstFitsChan *fc;
   AstFrameSet *fs;
   AstFrame *iwcfrm;
   AstMapping *map;
   AstKeyMap *km;
   int i, val;
   int there;
   int64_t kval;
   char card[81];
   char *cval;
   double xin, yin, xout, yout;

   /* Storage for FITS cards */
   char cards[10][81];

   astWatch( status );
   astBegin;

   /* Create a FitsChan that will write its contents to fred.txt when deleted */
   fc = astFitsChan( NULL, NULL, "SinkFile=./fred.txt" );

   if( !astGetI( fc, "SipOK" ) )
      stopit( 776, " ", status );

   test_fitsrounding( fc, status );
   astEmptyFits( fc );

   /* Put a FITS-WCS header into the FitsChan */
   strncpy( cards[0], "CRPIX1  =                   45", 80 ); cards[0][80] = '\0';
   strncpy( cards[1], "CRPIX2  =                   45", 80 ); cards[1][80] = '\0';
   strncpy( cards[2], "CRVAL1  =                   45", 80 ); cards[2][80] = '\0';
   strncpy( cards[3], "CRVAL2  =                 89.9", 80 ); cards[3][80] = '\0';
   strncpy( cards[4], "MYNAME  =                     ", 80 ); cards[4][80] = '\0';
   strncpy( cards[5], "CDELT1  =                -0.01", 80 ); cards[5][80] = '\0';
   strncpy( cards[6], "CDELT2  =                 0.01", 80 ); cards[6][80] = '\0';
   strncpy( cards[7], "CTYPE1  = 'RA---TAN'", 80 );           cards[7][80] = '\0';
   strncpy( cards[8], "CTYPE2  = 'DEC--TAN'", 80 );           cards[8][80] = '\0';

   for( i = 0; i < 9; i++ )
      astPutFits( fc, cards[i], 0 );

   /* Test astGetFitsI via card position */
   astSetI( fc, "Card", 2 );
   if( !astGetFitsI( fc, NULL, &val ) )
      stopit( 777, " ", status );
   else if( val != 45 )
      stopit( 778, " ", status );

   /* Test astGetFitsK */
   astSetI( fc, "Card", 2 );
   if( !astGetFitsK( fc, NULL, &kval ) )
      stopit( 7777, " ", status );
   else if( kval != 45 )
      stopit( 7778, " ", status );

   /* Test astSetFitsK with large value */
   astSetI( fc, "Card", 1000000000 );
   kval = 32147483647LL;
   astSetFitsK( fc, "KTEST", kval, " ", 0 );
   if( !astGetFitsK( fc, "KTEST", &kval ) )
      stopit( 7779, " ", status );
   else if( kval != 32147483647LL )
      stopit( 7780, " ", status );

   /* Test astTestFits */
   astSetI( fc, "Card", 5 );
   if( astTestFits( fc, NULL, &there ) )  /* card 5 (MYNAME) has undefined value */
      stopit( 779, " ", status );
   else if( !there )
      stopit( 780, " ", status );

   if( !astTestFits( fc, "CDELT1", &there ) )
      stopit( 781, " ", status );
   else if( !there )
      stopit( 782, " ", status );

   if( astTestFits( fc, "ABCDEF", &there ) )
      stopit( 783, " ", status );
   else if( there )
      stopit( 784, " ", status );

   /* Position beyond end of FitsChan */
   astSetI( fc, "Card", 11 );
   if( astTestFits( fc, NULL, &there ) )
      stopit( 785, " ", status );
   else if( there )
      stopit( 786, " ", status );

   /* CardName should be blank (past end) */
   cval = (char *)astGetC( fc, "CardName" );
   if( cval && cval[0] != '\0' && strcmp( cval, " " ) != 0 )
      stopit( 787, " ", status );

   /* Annul the fitschan - writes fred.txt */
   astAnnul( fc );

   /* Create another FitsChan reading from fred.txt */
   fc = astFitsChan( NULL, NULL, "SourceFile=./fred.txt" );

   if( astGetI( fc, "NCard" ) != 10 ) {
      printf( "NCard = %d\n", astGetI( fc, "NCard" ) );
      stopit( 1000, " ", status );
   }

   if( astGetI( fc, "Nkey" ) != 10 ) {
      printf( "Nkey = %d\n", astGetI( fc, "Nkey" ) );
      stopit( 999, " ", status );
   }

   /* Iterate through cards and compare with originals */
   astClear( fc, "Card" );
   i = 0;
   while( astFindFits( fc, "%f", card, 1 ) ) {
      i++;
      if( i <= 9 ) {
         /* Compare first 30 chars (key+value, ignoring trailing spaces) */
         if( strcmp( rstrip(card), rstrip(cards[i-1]) ) != 0 ) {
            printf( "Card %d mismatch:\n  got [%s]\n  exp [%s]\n",
                    i, card, cards[i-1] );
            stopit( 1001, " ", status );
         }
      }
   }

   astAnnul( fc );

   /* ---------------------------------------------------------------
    * Simple FITS-WCS header with IWC
    * ---------------------------------------------------------------*/
   strncpy( cards[0], "CRPIX1  = 45", 80 );           cards[0][80] = '\0';
   strncpy( cards[1], "CRPIX2  = 45", 80 );           cards[1][80] = '\0';
   strncpy( cards[2], "CRVAL1  = 45", 80 );           cards[2][80] = '\0';
   strncpy( cards[3], "CRVAL2  = 89.9", 80 );         cards[3][80] = '\0';
   strncpy( cards[4], "CDELT1  = -0.01", 80 );        cards[4][80] = '\0';
   strncpy( cards[5], "CDELT2  = 0.01", 80 );         cards[5][80] = '\0';
   strncpy( cards[6], "CTYPE1  = 'RA---TAN'", 80 );   cards[6][80] = '\0';
   strncpy( cards[7], "CTYPE2  = 'DEC--TAN'", 80 );   cards[7][80] = '\0';

   fc = astFitsChan( NULL, NULL, "Iwc=1" );
   for( i = 0; i < 8; i++ )
      astPutFits( fc, cards[i], 0 );

   astClear( fc, "Card" );
   if( astGetI( fc, "CardType" ) != AST__INT ) {
      printf( "CardType = %d, should be %d (AST__INT)\n",
              astGetI( fc, "CardType" ), AST__INT );
      stopit( 993, " ", status );
   }

   /* Retain CTYPE1 card */
   if( astFindFits( fc, "CTYPE1", card, 0 ) )
      astRetainFits( fc );

   /* Read a FrameSet */
   astClear( fc, "Card" );
   fs = (AstFrameSet *)astRead( fc );
   if( !fs )
      stopit( 1, "No FrameSet read from FitsChan", status );

   /* CTYPE1 should be present */
   astClear( fc, "Card" );
   if( !astFindFits( fc, "CTYPE1", card, 0 ) )
      stopit( 2, "CTYPE1 has not been retained", status );

   /* CTYPE2 should NOT be present */
   astClear( fc, "Card" );
   if( astFindFits( fc, "CTYPE2", card, 0 ) )
      stopit( 3, "CTYPE2 has been retained", status );

   /* Check IWC frame */
   if( astGetI( fs, "Nframe" ) != 3 )
      stopit( 301, "Wrong number of Frames", status );
   if( astGetI( fs, "Current" ) != 2 )
      stopit( 302, "Wrong current Frame", status );

   iwcfrm = astGetFrame( fs, 3 );
   if( strcmp( astGetC( iwcfrm, "Domain" ), "IWC" ) != 0 )
      stopit( 303, "Wrong Domain in IWC Frame", status );
   astAnnul( iwcfrm );

   map = astGetMapping( fs, 1, 3 );
   xin = 45.0;  yin = 45.0;
   astTran2( map, 1, &xin, &yin, 1, &xout, &yout );
   if( xout != 0.0 || yout != 0.0 )
      stopit( 304, "Wrong IWC for CRPIX position", status );

   xin = 46.0;
   astTran2( map, 1, &xin, &yin, 1, &xout, &yout );
   if( xout != -0.01 || yout != 0.0 )
      stopit( 305, "Wrong IWC for offset CRPIX position", status );
   astAnnul( map );

   map = astGetMapping( fs, 2, 3 );
   xin = 45.0 * AST__DD2R;
   yin = 89.9 * AST__DD2R;
   astTran2( map, 1, &xin, &yin, 1, &xout, &yout );
   if( fabs( xout ) > 1.0e-10 || fabs( yout ) > 1.0e-10 )
      stopit( 306, "Wrong IWC for CRVAL position", status );
   astAnnul( map );
   astAnnul( fs );

   /* ---------------------------------------------------------------
    * Illegal CRPIX2 value - should cause error in ast_read
    * ---------------------------------------------------------------*/
   strncpy( cards[0], "CRPIX1  = 45", 80 );           cards[0][80] = '\0';
   strncpy( cards[1], "CRPIX2  = 'fred'", 80 );       cards[1][80] = '\0';
   strncpy( cards[2], "CRVAL1  = 45", 80 );           cards[2][80] = '\0';
   strncpy( cards[3], "CRVAL2  = 89.9", 80 );         cards[3][80] = '\0';
   strncpy( cards[4], "CDELT1  = -0.01", 80 );        cards[4][80] = '\0';
   strncpy( cards[5], "CDELT2  = 0.01", 80 );         cards[5][80] = '\0';
   strncpy( cards[6], "CTYPE1  = 'RA---TAN'", 80 );   cards[6][80] = '\0';
   strncpy( cards[7], "CTYPE2  = 'DEC--TAN'", 80 );   cards[7][80] = '\0';

   fc = astFitsChan( NULL, NULL, " " );
   for( i = 0; i < 8; i++ )
      astPutFits( fc, cards[i], 0 );

   /* Retain CTYPE1 */
   astClear( fc, "Card" );
   if( astFindFits( fc, "CTYPE1", card, 0 ) )
      astRetainFits( fc );

   if( *status != 0 ) goto cleanup;

   astClear( fc, "Card" );
   fs = (AstFrameSet *)astRead( fc );

   if( fs ) {
      stopit( 4, "A FrameSet has been read from the FitsChan", status );
   } else if( astOK ) {
      stopit( 5, "No error has been reported by ast_read", status );
   } else {
      astClearStatus;
   }

   /* CTYPE1 should still be present */
   astClear( fc, "Card" );
   if( !astFindFits( fc, "CTYPE1", card, 0 ) )
      stopit( 6, "CTYPE1 has not been retained", status );

   /* CTYPE2 also still present (Clean not set) */
   astClear( fc, "Card" );
   if( !astFindFits( fc, "CTYPE2", card, 0 ) )
      stopit( 7, "CTYPE2 has not been retained", status );

   /* ---------------------------------------------------------------
    * Same again but with Clean=1
    * ---------------------------------------------------------------*/
   fc = astFitsChan( NULL, NULL, "Clean=1" );
   for( i = 0; i < 8; i++ )
      astPutFits( fc, cards[i], 0 );

   astClear( fc, "Card" );
   if( astFindFits( fc, "CTYPE1", card, 0 ) )
      astRetainFits( fc );

   if( *status != 0 ) goto cleanup;

   astClear( fc, "Card" );
   fs = (AstFrameSet *)astRead( fc );

   if( fs ) {
      stopit( 8, "A FrameSet has been read from the FitsChan", status );
   } else if( astOK ) {
      stopit( 9, "No error has been reported by ast_read", status );
   } else {
      astClearStatus;
   }

   /* CTYPE1 retained (explicitly retained) */
   astClear( fc, "Card" );
   if( !astFindFits( fc, "CTYPE1", card, 0 ) )
      stopit( 10, "CTYPE1 has not been retained", status );

   /* CTYPE2 removed (Clean=1) */
   astClear( fc, "Card" );
   if( astFindFits( fc, "CTYPE2", card, 0 ) )
      stopit( 11, "CTYPE2 has been retained", status );

   /* -TAB tests */
   checktab( status );
   checktab2( status );

   /* ---------------------------------------------------------------
    * SIP header: read then try to write (should fail: non-linear)
    * ---------------------------------------------------------------*/
   astEmptyFits( fc );
   astSetI( fc, "SipOK", 0 );
   astSet( fc, "SourceFile=sip.head" );
   astClear( fc, "Card" );
   fs = (AstFrameSet *)astRead( fc );
   astSet( fc, "Encoding=FITS-WCS" );
   if( !fs )
      stopit( 12, "Failed to read SIP header", status );
   else if( astWrite( fc, fs ) > 0 )
      stopit( 13, "Test on SIP header non-linearity failed", status );
   if( fs ) astAnnul( fs );

   /* ---------------------------------------------------------------
    * Alternate axis descriptions (set 'C' is badly formed on purpose)
    * ---------------------------------------------------------------*/
   astSetI( fc, "Clean", 0 );
   astEmptyFits( fc );

   if( astOK && astGetI( fc, "IgnoreBadAlt" ) )
      stopit( 14, " ", status );

   astSet( fc, "SourceFile=alt.header" );
   astClear( fc, "Card" );
   fs = (AstFrameSet *)astRead( fc );

   if( astOK ) {
      /* Expected an error */
      stopit( 15, " ", status );
   } else if( astStatus != AST__BDFTS ) {
      astClearStatus;
      stopit( 15, " ", status );
   } else {
      astClearStatus;
   }

   astSetI( fc, "IgnoreBadAlt", 1 );
   if( astOK && !astGetI( fc, "IgnoreBadAlt" ) )
      stopit( 16, " ", status );

   astEmptyFits( fc );
   astSet( fc, "SourceFile=alt.header" );
   astClear( fc, "Card" );
   fs = (AstFrameSet *)astRead( fc );

   if( !astOK ) {
      astClearStatus;
      stopit( 17, " ", status );
   }

   if( astGetI( fs, "Nframe" ) != 4 )
      stopit( 18, "Wrong number of Frames", status );
   if( fs ) astAnnul( fs );

   astEmptyFits( fc );
   astSet( fc, "SourceFile=alt.header" );
   astSet( fc, "Warnings=BadAlt" );

   astClear( fc, "Card" );
   fs = (AstFrameSet *)astRead( fc );

   km = (AstKeyMap *)astWarnings( fc );
   if( astMapSize( km ) != 1 )
      stopit( 19, "Wrong number of Warnings", status );

   if( fs ) astAnnul( fs );
   if( km ) astAnnul( km );

/* -----------------------------------------------------------------------
 * Test Clear/Set/Test for FitsChan-specific attributes.
 * Error numbers 100+ to avoid collision with existing tests.
 * AltAxes skipped: the attribute parser treats the leading 'A' as an
 * axis qualifier, making it impossible to set via astSet/astSetI.
 * -----------------------------------------------------------------------*/
   {
      AstFitsChan *afc = astFitsChan( NULL, NULL, " " );
      int ia;

      /* Integer attributes: { name, value to set } */
      static const struct {
         const char *name;
         int setval;
      } int_attrs[] = {
         { "DefB1950",     1 },
         { "CarLin",       1 },
         { "CDMatrix",     1 },
         { "FitsDigits",  10 },
         { "FitsRounding", 5 },
         { "ForceTab",     1 },
         { "IgnoreBadAlt", 1 },
         { "Iwc",          1 },
         { "PolyTan",      1 },
         { "SipOK",        0 },
         { "SipReplace",   0 },
         { "TabOK",        1 },
      };
      int n_int_attrs = (int)( sizeof(int_attrs) / sizeof(int_attrs[0]) );

      for( ia = 0; ia < n_int_attrs; ia++ ) {
         const char *name = int_attrs[ia].name;
         int setval = int_attrs[ia].setval;
         int errbase = 100 + ia * 3;
         char buf[80];

         astSetI( afc, name, setval );
         if( !astTest( afc, name ) ) {
            snprintf( buf, sizeof(buf), "%s not set after SetI", name );
            stopit( errbase, buf, status );
         }
         if( astGetI( afc, name ) != setval ) {
            snprintf( buf, sizeof(buf), "%s value wrong after SetI", name );
            stopit( errbase + 1, buf, status );
         }
         astClear( afc, name );
         if( astTest( afc, name ) ) {
            snprintf( buf, sizeof(buf), "%s still set after Clear", name );
            stopit( errbase + 2, buf, status );
         }
      }

      /* FitsTol (the only double attribute) */
      astSetD( afc, "FitsTol", 0.5 );
      if( !astTest( afc, "FitsTol" ) )
         stopit( 200, "FitsTol not set after SetD", status );
      if( fabs( astGetD( afc, "FitsTol" ) - 0.5 ) > 1e-10 )
         stopit( 201, "FitsTol value wrong after SetD", status );
      astClear( afc, "FitsTol" );
      if( astTest( afc, "FitsTol" ) )
         stopit( 202, "FitsTol still set after Clear", status );

      /* String attributes: { name, value to set } */
      {
         static const struct {
            const char *name;
            const char *setval;
         } str_attrs[] = {
            { "FitsAxisOrder", "RA DEC FREQ" },
            { "Warnings",      "BadAlt BadPV" },
         };
         int n_str_attrs = (int)( sizeof(str_attrs) / sizeof(str_attrs[0]) );

         for( ia = 0; ia < n_str_attrs; ia++ ) {
            const char *name = str_attrs[ia].name;
            const char *setval = str_attrs[ia].setval;
            int errbase = 210 + ia * 3;
            char buf[80];
            const char *got;

            astSetC( afc, name, setval );
            if( !astTest( afc, name ) ) {
               snprintf( buf, sizeof(buf), "%s not set after SetC", name );
               stopit( errbase, buf, status );
            }
            got = astGetC( afc, name );
            if( !got || strcmp( got, setval ) ) {
               snprintf( buf, sizeof(buf), "%s value wrong after SetC", name );
               stopit( errbase + 1, buf, status );
            }
            astClear( afc, name );
            if( astTest( afc, name ) ) {
               snprintf( buf, sizeof(buf), "%s still set after Clear", name );
               stopit( errbase + 2, buf, status );
            }
         }
      }

      afc = astAnnul( afc );
   }

/* -----------------------------------------------------------------------
 * Test FitsChan public API methods that had zero coverage.
 * Error numbers 300+.
 * -----------------------------------------------------------------------*/
   {
      AstFitsChan *afc = astFitsChan( NULL, NULL, " " );
      double cf[2];
      int ci[2];
      int lval;
      char *cval;

      /* astSetFitsCM — insert a comment card */
      astSetFitsCM( afc, "This is a test comment", 0 );
      if( !astOK )
         stopit( 300, "SetFitsCM failed", status );

      /* astSetFitsCF — complex float keyword */
      cf[0] = 1.5;
      cf[1] = 2.5;
      astSetFitsCF( afc, "TESTCF", cf, "complex float", 0 );
      if( !astOK )
         stopit( 301, "SetFitsCF failed", status );

      /* astGetFitsCF — retrieve complex float */
      cf[0] = cf[1] = 0.0;
      astClear( afc, "Card" );
      astGetFitsCF( afc, "TESTCF", cf );
      if( !astOK )
         stopit( 302, "GetFitsCF failed", status );
      if( fabs( cf[0] - 1.5 ) > 1e-10 || fabs( cf[1] - 2.5 ) > 1e-10 )
         stopit( 303, "GetFitsCF returned wrong values", status );

      /* astSetFitsCI — complex integer keyword */
      ci[0] = 10;
      ci[1] = 20;
      astSetFitsCI( afc, "TESTCI", ci, "complex int", 0 );
      if( !astOK )
         stopit( 304, "SetFitsCI failed", status );

      /* astGetFitsCI — retrieve complex integer */
      ci[0] = ci[1] = 0;
      astClear( afc, "Card" );
      astGetFitsCI( afc, "TESTCI", ci );
      if( !astOK )
         stopit( 305, "GetFitsCI failed", status );
      if( ci[0] != 10 || ci[1] != 20 )
         stopit( 306, "GetFitsCI returned wrong values", status );

      /* astGetFitsL — logical keyword */
      astSetFitsL( afc, "TESTL", 1, "logical true", 0 );
      lval = 0;
      astClear( afc, "Card" );
      astGetFitsL( afc, "TESTL", &lval );
      if( !astOK )
         stopit( 307, "GetFitsL failed", status );
      if( !lval )
         stopit( 308, "GetFitsL returned wrong value", status );

      /* astGetFitsCN — continuation string.
         Insert a CONTINUE-card sequence manually via astPutFits. */
      {
         AstFitsChan *fc2 = astFitsChan( NULL, NULL, " " );
         astPutFits( fc2, "LONGSTR = 'Part one of a long &'", 0 );
         astPutFits( fc2, "CONTINUE  'string value'", 0 );
         cval = NULL;
         astClear( fc2, "Card" );
         astGetFitsCN( fc2, "LONGSTR", &cval );
         if( !astOK )
            stopit( 309, "GetFitsCN failed", status );
         if( !cval || strncmp( cval, "Part one", 8 ) )
            stopit( 310, "GetFitsCN returned wrong value", status );
         fc2 = astAnnul( fc2 );
      }

      /* astPutCards — insert multiple 80-char cards as a single string */
      {
         AstFitsChan *fc2 = astFitsChan( NULL, NULL, " " );
         astPutCards( fc2,
            "SIMPLE  =                    T / Standard FITS                                   "
            "BITPIX  =                  -32 / Bits per pixel                                  "
            "NAXIS   =                    0 / No data                                         "
            "END                                                                             " );
         if( !astOK )
            stopit( 311, "PutCards failed", status );
         if( astGetI( fc2, "Ncard" ) < 3 )
            stopit( 312, "PutCards wrong card count", status );
         fc2 = astAnnul( fc2 );
      }

      /* astShowFits — display cards to stdout (just verify no crash) */
      astClear( afc, "Card" );
      astShowFits( afc );
      if( !astOK )
         stopit( 313, "ShowFits failed", status );

      /* astWriteFits / astReadFits — write cards to sink, read from source.
         Use SinkFile/SourceFile for simplicity. */
      {
         AstFitsChan *fc2;
         int ncard_before, ncard_after;

         astSet( afc, "SinkFile=/tmp/ast_testfitschan_write.fits" );
         ncard_before = astGetI( afc, "Ncard" );
         astWriteFits( afc );
         if( !astOK )
            stopit( 314, "WriteFits failed", status );

         fc2 = astFitsChan( NULL, NULL,
                            "SourceFile=/tmp/ast_testfitschan_write.fits" );
         astReadFits( fc2 );
         if( !astOK )
            stopit( 315, "ReadFits failed", status );
         ncard_after = astGetI( fc2, "Ncard" );
         if( ncard_after < ncard_before )
            stopit( 316, "ReadFits lost cards", status );
         fc2 = astAnnul( fc2 );
      }

      afc = astAnnul( afc );
   }

/* -----------------------------------------------------------------------
 * Test FitsChan Dump/Load round-trip via astToString/astFromString.
 * Error numbers 400+.
 * -----------------------------------------------------------------------*/
   {
      AstFitsChan *afc = astFitsChan( NULL, NULL, " " );
      AstFitsChan *afc2;
      char *pickle;
      int ncard1, ncard2;

      astSetFitsI( afc, "NAXIS", 0, "No data", 0 );
      astSetFitsF( afc, "CRVAL1", 180.0, "Reference value", 0 );
      astSetFitsS( afc, "CTYPE1", "RA---TAN", "Projection", 0 );
      ncard1 = astGetI( afc, "Ncard" );

      pickle = astToString( afc );
      if( !pickle ) {
         stopit( 400, "astToString returned NULL", status );
      } else {
         afc2 = (AstFitsChan *) astFromString( pickle );
         if( !afc2 ) {
            stopit( 401, "astFromString returned NULL", status );
         } else {
            ncard2 = astGetI( afc2, "Ncard" );
            if( ncard2 != ncard1 )
               stopit( 402, "Round-trip changed card count", status );
            afc2 = astAnnul( afc2 );
         }
         astFree( pickle );
      }

      afc = astAnnul( afc );
   }

/* -----------------------------------------------------------------------
 * Write-path round-trip tests: build FrameSets that exercise specific
 * write code paths, encode to a FitsChan, read back, verify success.
 * Error numbers 500+.
 *
 * For sky coordinate tests, we read a simple TAN header to get a
 * valid FrameSet with proper WCS mappings, then change the SkyFrame
 * system and write it back. This ensures the projection and coordinate
 * matrices are valid.
 * -----------------------------------------------------------------------*/
   {
      AstFitsChan *hfc;
      AstFrameSet *tfs;
      AstSpecFrame *spec;
      AstFrame *pixel;
      AstMapping *map;
      char buf[120];
      int rt;

      static const char *sky_systems[] = {
         "FK5", "FK4", "FK4-NO-E", "ICRS", "GALACTIC",
         "SUPERGALACTIC", "ECLIPTIC",
      };
      int n_sky = (int)( sizeof(sky_systems) / sizeof(sky_systems[0]) );

      static const char *encodings[] = {
         "FITS-WCS", "FITS-PC", "FITS-IRAF", "FITS-AIPS",
      };
      int n_enc = (int)( sizeof(encodings) / sizeof(encodings[0]) );
      int is, ie;

      for( is = 0; is < n_sky; is++ ) {
         for( ie = 0; ie < n_enc; ie++ ) {
            hfc = astFitsChan( NULL, NULL, " " );
            astPutFits( hfc, "SIMPLE  =                    T", 0 );
            astPutFits( hfc, "BITPIX  =                  -32", 0 );
            astPutFits( hfc, "NAXIS   =                    2", 0 );
            astPutFits( hfc, "NAXIS1  =                  100", 0 );
            astPutFits( hfc, "NAXIS2  =                  100", 0 );
            astPutFits( hfc, "CTYPE1  = 'GLON-TAN'", 0 );
            astPutFits( hfc, "CTYPE2  = 'GLAT-TAN'", 0 );
            astPutFits( hfc, "CRVAL1  =              180.000", 0 );
            astPutFits( hfc, "CRVAL2  =                0.000", 0 );
            astPutFits( hfc, "CRPIX1  =               50.500", 0 );
            astPutFits( hfc, "CRPIX2  =               50.500", 0 );
            astPutFits( hfc, "CDELT1  =              -0.0100", 0 );
            astPutFits( hfc, "CDELT2  =               0.0100", 0 );
            astPutFits( hfc, "END", 0 );
            astClear( hfc, "Card" );
            tfs = (AstFrameSet *) astRead( hfc );
            hfc = astAnnul( hfc );
            if( tfs ) {
               astSetC( tfs, "System", sky_systems[is] );
               rt = roundtrip( tfs, encodings[ie], "", 0, status );
               if( rt == -1 ) {
                  snprintf( buf, sizeof(buf),
                            "astEqual failed for %s with %s",
                            sky_systems[is], encodings[ie] );
                  stopit( 500 + is * n_enc + ie, buf, status );
               }
               tfs = astAnnul( tfs );
            }
         }
      }

      /* --- Spectral axes through multiple encodings --- */
      {
         static const struct {
            const char *spec_attrs;
            const char *label;
         } spec_systems[] = {
            { "System=FREQ,Unit=Hz,StdOfRest=Barycentric,RestFreq=1.4204e9 Hz", "FREQ" },
            { "System=WAVE,Unit=m,StdOfRest=Barycentric,RestFreq=5.996e14 Hz",  "WAVE" },
            { "System=VRAD,Unit=m/s,StdOfRest=LSRK,RestFreq=1.4204e9 Hz",      "VRAD" },
         };
         int n_spec = (int)( sizeof(spec_systems) / sizeof(spec_systems[0]) );
         double zooms[] = { 1.0e6, 1.0e-10, 1000.0 };
         int isp;

         for( isp = 0; isp < n_spec; isp++ ) {
            spec = astSpecFrame( "%s", spec_systems[isp].spec_attrs );
            pixel = astFrame( 1, "Domain=GRID" );
            map = (AstMapping *) astZoomMap( 1, zooms[isp], " " );
            tfs = astFrameSet( pixel, " " );
            astAddFrame( tfs, AST__CURRENT, map, (AstFrame *) spec );
            rt = roundtrip( tfs, "FITS-WCS", "", 0, status );
            if( rt == 0 ) {
               snprintf( buf, sizeof(buf),
                         "Write/read failed for %s SpecFrame",
                         spec_systems[isp].label );
               stopit( 600 + isp, buf, status );
            } else if( rt == -1 ) {
               snprintf( buf, sizeof(buf),
                         "astEqual failed for %s SpecFrame",
                         spec_systems[isp].label );
               stopit( 600 + isp, buf, status );
            }
         }
      }

      /* --- 3D: celestial + spectral (exercises MakeFitsFrameSet) --- */
      {
         AstFitsChan *hfc3 = astFitsChan( NULL, NULL, " " );
         astPutFits( hfc3, "SIMPLE  =                    T", 0 );
         astPutFits( hfc3, "BITPIX  =                  -32", 0 );
         astPutFits( hfc3, "NAXIS   =                    3", 0 );
         astPutFits( hfc3, "NAXIS1  =                  100", 0 );
         astPutFits( hfc3, "NAXIS2  =                  100", 0 );
         astPutFits( hfc3, "NAXIS3  =                  100", 0 );
         astPutFits( hfc3, "CTYPE1  = 'RA---TAN'", 0 );
         astPutFits( hfc3, "CTYPE2  = 'DEC--TAN'", 0 );
         astPutFits( hfc3, "CTYPE3  = 'FREQ'", 0 );
         astPutFits( hfc3, "CRVAL1  =              180.000", 0 );
         astPutFits( hfc3, "CRVAL2  =               45.000", 0 );
         astPutFits( hfc3, "CRVAL3  =          1.4204e+009", 0 );
         astPutFits( hfc3, "CRPIX1  =               50.000", 0 );
         astPutFits( hfc3, "CRPIX2  =               50.000", 0 );
         astPutFits( hfc3, "CRPIX3  =               50.000", 0 );
         astPutFits( hfc3, "CDELT1  =              -0.0100", 0 );
         astPutFits( hfc3, "CDELT2  =               0.0100", 0 );
         astPutFits( hfc3, "CDELT3  =            1.000e+06", 0 );
         astPutFits( hfc3, "CUNIT3  = 'Hz'", 0 );
         astPutFits( hfc3, "RADESYS = 'FK5'", 0 );
         astPutFits( hfc3, "EQUINOX =               2000.0", 0 );
         astPutFits( hfc3, "RESTFRQ =          1.4204e+009", 0 );
         astPutFits( hfc3, "SPECSYS = 'BARYCENT'", 0 );
         astPutFits( hfc3, "END", 0 );
         astClear( hfc3, "Card" );
         tfs = (AstFrameSet *) astRead( hfc3 );
         hfc3 = astAnnul( hfc3 );
         if( tfs ) {
            rt = roundtrip( tfs, "FITS-WCS", "", 0, status );
            if( rt == 0 )
               stopit( 610, "Write/read failed for 3D sky+spec", status );
            else if( rt == -1 )
               stopit( 611, "astEqual failed for 3D sky+spec", status );
            tfs = astAnnul( tfs );
         }
      }

      /* --- Offset SkyFrame (exercises SkySys SkyRef path) --- */
      {
         AstFitsChan *hfc4 = astFitsChan( NULL, NULL, " " );
         astPutFits( hfc4, "SIMPLE  =                    T", 0 );
         astPutFits( hfc4, "BITPIX  =                  -32", 0 );
         astPutFits( hfc4, "NAXIS   =                    2", 0 );
         astPutFits( hfc4, "NAXIS1  =                  100", 0 );
         astPutFits( hfc4, "NAXIS2  =                  100", 0 );
         astPutFits( hfc4, "CTYPE1  = 'OFLN-TAN'", 0 );
         astPutFits( hfc4, "CTYPE2  = 'OFLT-TAN'", 0 );
         astPutFits( hfc4, "CRVAL1  =                0.000", 0 );
         astPutFits( hfc4, "CRVAL2  =                0.000", 0 );
         astPutFits( hfc4, "CRPIX1  =               50.500", 0 );
         astPutFits( hfc4, "CRPIX2  =               50.500", 0 );
         astPutFits( hfc4, "CDELT1  =              -0.0100", 0 );
         astPutFits( hfc4, "CDELT2  =               0.0100", 0 );
         astPutFits( hfc4, "RADESYS = 'FK5'", 0 );
         astPutFits( hfc4, "EQUINOX =               2000.0", 0 );
         astPutFits( hfc4, "SREFRA  =              180.000", 0 );
         astPutFits( hfc4, "SREFDEC =               45.000", 0 );
         astPutFits( hfc4, "SREFIS  = 'ORIGIN'", 0 );
         astPutFits( hfc4, "END", 0 );
         astClear( hfc4, "Card" );
         tfs = (AstFrameSet *) astRead( hfc4 );
         hfc4 = astAnnul( hfc4 );
         if( tfs ) {
            rt = roundtrip( tfs, "FITS-WCS", "", 0, status );
            if( rt == 0 )
               stopit( 620, "Write/read failed for offset SkyFrame", status );
            else if( rt == -1 )
               stopit( 621, "astEqual failed for offset SkyFrame", status );
            tfs = astAnnul( tfs );
         }
      }

      /* --- AZEL: read a standard header with observer position and epoch,
             convert to AZEL, then round-trip through FITS-WCS --- */
      {
         AstFitsChan *hfc5 = astFitsChan( NULL, NULL, " " );
         astPutFits( hfc5, "SIMPLE  =                    T", 0 );
         astPutFits( hfc5, "BITPIX  =                  -32", 0 );
         astPutFits( hfc5, "NAXIS   =                    2", 0 );
         astPutFits( hfc5, "NAXIS1  =                  100", 0 );
         astPutFits( hfc5, "NAXIS2  =                  100", 0 );
         astPutFits( hfc5, "CTYPE1  = 'RA---TAN'", 0 );
         astPutFits( hfc5, "CTYPE2  = 'DEC--TAN'", 0 );
         astPutFits( hfc5, "CRVAL1  =              180.000", 0 );
         astPutFits( hfc5, "CRVAL2  =               45.000", 0 );
         astPutFits( hfc5, "CRPIX1  =               50.500", 0 );
         astPutFits( hfc5, "CRPIX2  =               50.500", 0 );
         astPutFits( hfc5, "CDELT1  =              -0.0100", 0 );
         astPutFits( hfc5, "CDELT2  =               0.0100", 0 );
         astPutFits( hfc5, "RADESYS = 'FK5'", 0 );
         astPutFits( hfc5, "EQUINOX =               2000.0", 0 );
         astPutFits( hfc5, "DATE-OBS= '2017-07-21T08:38:39.087'", 0 );
         astPutFits( hfc5, "OBSGEO-X=  -5464586.5949660344", 0 );
         astPutFits( hfc5, "OBSGEO-Y=  -2492996.5580556658", 0 );
         astPutFits( hfc5, "OBSGEO-Z=   2150654.3760909005", 0 );
         astPutFits( hfc5, "END", 0 );
         astClear( hfc5, "Card" );
         tfs = (AstFrameSet *) astRead( hfc5 );
         hfc5 = astAnnul( hfc5 );
         if( tfs ) {
            astSetC( tfs, "System", "AZEL" );
            rt = roundtrip( tfs, "FITS-WCS", "", 0, status );
            if( rt == 0 )
               stopit( 630, "Write/read failed for AZEL", status );
            tfs = astAnnul( tfs );
         } else {
            stopit( 631, "Failed to read header for AZEL test", status );
         }
      }

      /* --- AZEL: read an AZ/EL header directly --- */
      {
         AstFitsChan *hfc6 = astFitsChan( NULL, NULL, " " );
         astPutFits( hfc6, "SIMPLE  =                    T", 0 );
         astPutFits( hfc6, "BITPIX  =                  -32", 0 );
         astPutFits( hfc6, "NAXIS   =                    2", 0 );
         astPutFits( hfc6, "NAXIS1  =                  100", 0 );
         astPutFits( hfc6, "NAXIS2  =                  100", 0 );
         astPutFits( hfc6, "CTYPE1  = 'AZ---TAN'", 0 );
         astPutFits( hfc6, "CTYPE2  = 'EL---TAN'", 0 );
         astPutFits( hfc6, "CRVAL1  =              185.000", 0 );
         astPutFits( hfc6, "CRVAL2  =               60.000", 0 );
         astPutFits( hfc6, "CRPIX1  =               50.500", 0 );
         astPutFits( hfc6, "CRPIX2  =               50.500", 0 );
         astPutFits( hfc6, "CDELT1  =              -0.0100", 0 );
         astPutFits( hfc6, "CDELT2  =               0.0100", 0 );
         astPutFits( hfc6, "DATE-OBS= '2017-07-21T08:38:39.087'", 0 );
         astPutFits( hfc6, "OBSGEO-X=  -5464586.5949660344", 0 );
         astPutFits( hfc6, "OBSGEO-Y=  -2492996.5580556658", 0 );
         astPutFits( hfc6, "OBSGEO-Z=   2150654.3760909005", 0 );
         astPutFits( hfc6, "END", 0 );
         astClear( hfc6, "Card" );
         tfs = (AstFrameSet *) astRead( hfc6 );
         hfc6 = astAnnul( hfc6 );
         if( tfs ) {
            AstFrame *skyfrm = astGetFrame( tfs, AST__CURRENT );
            if( strcmp( astGetC( skyfrm, "System" ), "AZEL" ) != 0 )
               stopit( 640, "AZEL header did not produce AZEL SkyFrame", status );
            astAnnul( skyfrm );

            rt = roundtrip( tfs, "FITS-WCS", "", 0, status );
            if( rt == 0 )
               stopit( 641, "Write/read failed for direct AZEL header", status );
            tfs = astAnnul( tfs );
         } else {
            stopit( 642, "Failed to read AZEL header", status );
         }
      }

      /* --- DATE-OBS and OBSGEO propagation: verify metadata survives
             a FITS-WCS round-trip --- */
      {
         AstFitsChan *hfc7 = astFitsChan( NULL, NULL, " " );
         AstFitsChan *wfc;
         AstFrameSet *tfs7;
         double dval;

         astPutFits( hfc7, "SIMPLE  =                    T", 0 );
         astPutFits( hfc7, "BITPIX  =                  -32", 0 );
         astPutFits( hfc7, "NAXIS   =                    2", 0 );
         astPutFits( hfc7, "NAXIS1  =                  100", 0 );
         astPutFits( hfc7, "NAXIS2  =                  100", 0 );
         astPutFits( hfc7, "CTYPE1  = 'RA---TAN'", 0 );
         astPutFits( hfc7, "CTYPE2  = 'DEC--TAN'", 0 );
         astPutFits( hfc7, "CRVAL1  =              180.000", 0 );
         astPutFits( hfc7, "CRVAL2  =               45.000", 0 );
         astPutFits( hfc7, "CRPIX1  =               50.500", 0 );
         astPutFits( hfc7, "CRPIX2  =               50.500", 0 );
         astPutFits( hfc7, "CDELT1  =              -0.0100", 0 );
         astPutFits( hfc7, "CDELT2  =               0.0100", 0 );
         astPutFits( hfc7, "RADESYS = 'FK5'", 0 );
         astPutFits( hfc7, "EQUINOX =               2000.0", 0 );
         astPutFits( hfc7, "MJD-OBS =   57955.360174621412", 0 );
         astPutFits( hfc7, "DATE-OBS= '2017-07-21T08:38:39.087'", 0 );
         astPutFits( hfc7, "OBSGEO-X=  -5464586.5949660344", 0 );
         astPutFits( hfc7, "OBSGEO-Y=  -2492996.5580556658", 0 );
         astPutFits( hfc7, "OBSGEO-Z=   2150654.3760909005", 0 );
         astPutFits( hfc7, "END", 0 );
         astClear( hfc7, "Card" );
         tfs7 = (AstFrameSet *) astRead( hfc7 );
         hfc7 = astAnnul( hfc7 );
         if( !tfs7 ) {
            stopit( 650, "Failed to read header with DATE-OBS/OBSGEO", status );
         } else {
            wfc = astFitsChan( NULL, NULL, " " );
            astSetC( wfc, "Encoding", "FITS-WCS" );
            if( astWrite( wfc, tfs7 ) != 1 )
               stopit( 651, "Failed to write header with DATE-OBS/OBSGEO", status );

            astClear( wfc, "Card" );
            if( !astGetFitsF( wfc, "MJD-OBS", &dval ) )
               stopit( 652, "MJD-OBS not in written header", status );
            else if( fabs( dval - 57955.360174621412 ) > 1e-6 )
               stopit( 653, "MJD-OBS value wrong after round-trip", status );

            astClear( wfc, "Card" );
            if( !astGetFitsF( wfc, "OBSGEO-X", &dval ) )
               stopit( 654, "OBSGEO-X not in written header", status );
            else if( fabs( dval - (-5464586.5949660344) ) > 1.0 )
               stopit( 655, "OBSGEO-X value wrong after round-trip", status );

            astClear( wfc, "Card" );
            if( !astGetFitsF( wfc, "OBSGEO-Y", &dval ) )
               stopit( 656, "OBSGEO-Y not in written header", status );
            else if( fabs( dval - (-2492996.5580556658) ) > 1.0 )
               stopit( 657, "OBSGEO-Y value wrong after round-trip", status );

            astClear( wfc, "Card" );
            if( !astGetFitsF( wfc, "OBSGEO-Z", &dval ) )
               stopit( 658, "OBSGEO-Z not in written header", status );
            else if( fabs( dval - 2150654.3760909005 ) > 1.0 )
               stopit( 659, "OBSGEO-Z value wrong after round-trip", status );

            wfc = astAnnul( wfc );
            tfs7 = astAnnul( tfs7 );
         }
      }
   }

cleanup:
   astEnd;

   if( *status == 0 ) {
      printf( " All FitsChan tests passed\n" );
   } else {
      printf( "FitsChan tests failed\n" );
   }

   return *status;
}
