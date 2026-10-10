/*
 *  Test the MocChan class.
 *  Converted from the Fortran test testmocchan.f.
 *
 *  The Fortran version used COMMON-block source/sink callbacks with an
 *  in-memory file array. This C version uses static arrays with C-style
 *  source (returns const char *) and sink (takes const char *) callbacks.
 */
#include "ast.h"
#include <stdio.h>
#include <string.h>

#define MXLINE 500
#define LINELEN 200

static char files[MXLINE][LINELEN];
static int filelen;
static int iline;

static void stopit( const char *text, int *status ) {
   if( *status != 0 ) return;
   *status = 1;
   printf( "%s\n", text );
}

static const char *source( void ) {
   if( iline < filelen ) {
      return files[iline++];
   }
   return NULL;
}

static void sink( const char *line );
static char dump[MXLINE][LINELEN];
static int ndump;

static void dumpsink( const char *line ) {
   if( ndump < MXLINE && line ) {
      strncpy( dump[ndump], line, LINELEN - 1 );
      dump[ndump][LINELEN - 1] = '\0';
      ndump++;
   }
}

static const char *dumpsource( void ) {
   if( iline < ndump ) {
      return dump[iline++];
   }
   return NULL;
}

/* A MocChan with a STRING or JSON MocFormat is dumped with that MocEnc
   name and loaded with that MocFormat. */
static void checkDumpFormat( const char *format, int *status ) {
   AstMocChan *mch;
   AstChannel *dch;
   AstObject *obj;
   char want[ 40 ];
   int i, found;

   if( *status != 0 ) return;
   mch = astMocChan( NULL, NULL, "MocFormat=%s", format );
   dch = astChannel( dumpsource, dumpsink, " " );
   ndump = 0;
   astWrite( dch, mch );
   sprintf( want, "MocEnc = \"%s\"", format );
   found = 0;
   for( i = 0; i < ndump; i++ ) {
      if( strstr( dump[i], want ) ) found = 1;
   }
   if( !found ) stopit( "Error dump 1", status );
   iline = 0;
   obj = astRead( dch );
   if( !obj || strcmp( astGetC( obj, "MocFormat" ), format ) )
      stopit( "Error dump 2", status );
}

/* A MocLineLen of zero, which any negative value becomes, is too short
   for any value and is reported rather than crashing. */
static void checkZeroLineLen( int *status ) {
   AstMocChan *mch;
   AstMoc *moc;
   int ret;

   if( *status != 0 ) return;
   moc = astMoc( " " );
   astAddCell( moc, AST__OR, 3, (int64_t)5 );
   mch = astMocChan( NULL, sink, "MocLineLen=-5" );
   if( astGetI( mch, "MocLineLen" ) != 0 ) stopit( "Error zero 1", status );
   filelen = 0;
   ret = astWrite( mch, moc );
   if( astStatus != AST__SMBUF || ret != 0 ) stopit( "Error zero 2", status );
   astClearStatus;
}

static void sink( const char *line ) {
   if( filelen < MXLINE && line && line[0] ) {
      strncpy( files[filelen], line, LINELEN - 1 );
      files[filelen][LINELEN - 1] = '\0';
      filelen++;
   }
}

int main( void ) {
   int status_value = 0;
   int *status = &status_value;
   AstMoc *moc;
   AstSkyFrame *sf;
   AstRegion *reg1;
   AstMocChan *ch;
   AstObject *obj;
   double point[1], centre[2];

   astWatch( status );
   astBegin;

   /* Create a simple MOC */
   moc = astMoc( "maxorder=9" );
   astAddCell( moc, AST__OR, 7, (int64_t)1000 );
   astAddCell( moc, AST__OR, 7, (int64_t)1001 );
   astAddCell( moc, AST__OR, 7, (int64_t)1002 );
   astAddCell( moc, AST__OR, 7, (int64_t)997 );
   astAddCell( moc, AST__OR, 6, (int64_t)500 );

   /* Write MOC to channel (default STRING format) */
   ch = astMocChan( source, sink, " " );
   filelen = 0;
   if( astWrite( ch, moc ) != 1 )
      stopit( "Error 1", status );

   if( strcmp( files[0], "6/500 7/997,1000-1002 9/" ) != 0 )
      stopit( "Error 2", status );

   if( astTest( ch, "MocFormat" ) )
      stopit( "Error 3", status );

   /* Read back */
   iline = 0;
   obj = astRead( ch );
   if( !obj )
      stopit( "Error 4", status );

   if( obj && !astEqual( obj, moc ) )
      stopit( "Error 5", status );

   if( !astTest( ch, "MocFormat" ) )
      stopit( "Error 6", status );

   if( strcmp( astGetC( ch, "MocFormat" ), "STRING" ) != 0 )
      stopit( "Error 7", status );

   /* Write in JSON format */
   astSetC( ch, "MocFormat", "json" );
   filelen = 0;
   if( astWrite( ch, moc ) != 1 )
      stopit( "Error 8", status );

   if( strcmp( files[0],
              "{\"6\":[500],\"7\":[997,1000,1001,1002],\"9\":[]}" ) != 0 )
      stopit( "Error 9", status );

   /* Read back from JSON */
   astClear( ch, "MocFormat" );
   iline = 0;
   obj = astRead( ch );
   if( !obj )
      stopit( "Error 10", status );

   if( obj && !astEqual( obj, moc ) )
      stopit( "Error 11", status );

   if( strcmp( astGetC( ch, "MocFormat" ), "JSON" ) != 0 )
      stopit( "Error 12", status );

   astClear( ch, "MocFormat" );

   /* Create a MOC from a sky region (circle) */
   moc = astMoc( "maxorder=15,minorder=12" );
   sf = astSkyFrame( "system=icrs" );
   centre[0] = 1.0;
   centre[1] = 1.0;
   point[0] = 0.01;
   reg1 = (AstRegion *)astCircle( sf, 1, centre, point, NULL, " " );
   astAddRegion( moc, AST__OR, reg1 );

   /* Write and read back (STRING format) */
   filelen = 0;
   if( astWrite( ch, moc ) != 1 )
      stopit( "Error 13", status );

   iline = 0;
   obj = astRead( ch );
   if( !obj )
      stopit( "Error 14", status );

   if( obj && !astEqual( obj, moc ) )
      stopit( "Error 15", status );

   /* Write and read back (JSON format) */
   astSetC( ch, "MocFormat", "json" );

   filelen = 0;
   if( astWrite( ch, moc ) != 1 )
      stopit( "Error 16", status );

   iline = 0;
   obj = astRead( ch );
   if( !obj )
      stopit( "Error 17", status );

   if( obj && !astEqual( obj, moc ) )
      stopit( "Error 18", status );

   checkDumpFormat( "STRING", status );
   checkDumpFormat( "JSON", status );
   checkZeroLineLen( status );

   astEnd;

   if( *status == 0 ) {
      printf( " All MocChan tests passed\n" );
   } else {
      printf( "MocChan tests failed\n" );
   }
   return *status;
}
