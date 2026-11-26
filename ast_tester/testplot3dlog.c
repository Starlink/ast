/*
*  Name:
*     testplot3dlog.c
*
*  Purpose:
*     Exercise the Plot3D class through the logging grf3d module, with no
*     external graphics library.
*
*  Usage:
*     testplot3dlog --in FILE [--attr STR] --box "x1 y1 z1 x2 y2 z2"
*                   --bbox "x1 y1 z1 x2 y2 z2" [--mm N] [--log FILE]
*                   [--dump FILE] [--svg FILE [W H]]
*                   [--camera "ex ey ez tx ty tz ux uy uz screen"]
*                   --do ACTION [--do ACTION ...]
*
*  Description:
*     Reads a Frame, FrameSet or Region from FILE (FITS header cards
*     through a FitsChan when FILE ends in ".head", otherwise a Native
*     dump), builds a Plot3D over the given graphics box (--box) and
*     physical box (--bbox) with the attribute settings in --attr,
*     checks that the Plot3D survives a serialisation round trip, writes
*     its Native dump to --dump, and performs the actions in order while
*     the logging grf3d module (grf3d_log.c) records every graphics call
*     in --log.  With --svg the calls are also drawn through a camera:
*     the auto camera for the graphics box, replaced by --camera when
*     given.
*
*     Each action writes an "ACTION <name>" line to the log first.
*     Actions:
*        grid                      astGrid
*        border                    astBorder
*        outline                   astRegionOutline with the input Region
*        refused                   each Plot method a Plot3D refuses, with
*                                  the error written as "ERR <code> <message>"
*        text:STR:x,y,z:ux,uy,uz:JUST
*                                  astText; "~" in STR stands for a space
*        mark:TYPE:x,y,z[/x,y,z...]
*                                  astMark at each position
*        set:ATTRS                 astSet
*        clear:ATTR                astClear
*        setcurrent:N              astSetI Current
*        clearcurrent              astClear Current
*        removeframe:N             astRemoveFrame
*
*     Any AST error outside the "refused" action is a failure.
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ast.h"
#include "grf_log.h"
#include "grf3d_log.h"

static int status = 0;
static FILE *logfp = NULL;

/* The most recent error messages, collected by an error handler so that
   they can be written to the log. */
static char errbuf[4096];

static void collect_error( int status_value, const char *message ) {
   (void)status_value;
   if( strlen( errbuf ) + strlen( message ) + 2 < sizeof(errbuf) ) {
      if( errbuf[0] ) strcat( errbuf, " / " );
      strcat( errbuf, message );
   }
}

static int fail( const char *what ) {
   printf( "testplot3dlog: %s\n", what );
   if( errbuf[0] ) printf( "   %s\n", errbuf );
   return 1;
}

/* Parse "n" whitespace-separated doubles. */
static int parse_doubles( const char *text, int n, double *out ) {
   const char *p = text;
   int i;
   for( i = 0; i < n; i++ ) {
      char *end;
      out[i] = strtod( p, &end );
      if( end == p ) return 0;
      p = end;
   }
   return 1;
}

/* Parse "x,y,z" into three doubles; returns a pointer past the triple. */
static const char *parse_triple( const char *p, double *out ) {
   int i;
   for( i = 0; i < 3; i++ ) {
      char *end;
      out[i] = strtod( p, &end );
      if( end == p ) return NULL;
      p = end;
      if( i < 2 ) {
         if( *p != ',' ) return NULL;
         p++;
      }
   }
   return p;
}

/* Read a whole file into a malloc'd string. */
static char *read_file( const char *name ) {
   FILE *fp = fopen( name, "r" );
   char *buf = NULL;
   size_t total = 0, cap = 0;
   char chunk[1024];
   if( !fp ) return NULL;
   while( fgets( chunk, sizeof(chunk), fp ) ) {
      size_t len = strlen( chunk );
      if( total + len + 1 > cap ) {
         cap = ( total + len + 1 ) * 2;
         buf = realloc( buf, cap );
      }
      memcpy( buf + total, chunk, len );
      total += len;
   }
   fclose( fp );
   if( buf ) buf[total] = '\0';
   return buf;
}

/* Read the input object: FITS cards through a FitsChan, or a Native dump. */
static AstObject *read_input( const char *name ) {
   AstObject *result = NULL;
   const char *ext = strrchr( name, '.' );
   if( ext && strcmp( ext, ".head" ) == 0 ) {
      AstFitsChan *fc = astFitsChan( NULL, NULL, " " );
      FILE *fp = fopen( name, "r" );
      char card[256];
      if( !fp ) return NULL;
      while( fgets( card, sizeof(card), fp ) ) {
         size_t len = strlen( card );
         while( len > 0 && ( card[len-1] == '\n' || card[len-1] == '\r' ) ) {
            card[--len] = '\0';
         }
         astPutFits( fc, card, 0 );
      }
      fclose( fp );
      astClear( fc, "Card" );
      result = astRead( fc );
      fc = astAnnul( fc );
   } else {
      char *text = read_file( name );
      if( text ) {
         result = astFromString( text );
         free( text );
      }
   }
   return result;
}

static void log_action( const char *name ) {
   if( logfp ) fprintf( logfp, "ACTION %s\n", name );
}

/* Call a refused Plot method and log the error it reports. */
static void log_refused( const char *name ) {
   if( logfp ) {
      fprintf( logfp, "ERR %d %s\n", status, errbuf[0] ? errbuf : "(no message)" );
   }
   if( astOK ) printf( "testplot3dlog: %s was not refused\n", name );
   astClearStatus;
   errbuf[0] = '\0';
}

static int do_refused( AstPlot3D *plot ) {
   float lbnd[2], ubnd[2];
   double dlbnd[2] = { 0.0, 0.0 }, dubnd[2] = { 1.0, 1.0 };
   double start[3] = { 0.0, 0.0, 0.0 }, finish[3] = { 1.0, 1.0, 1.0 };
   double poly[6] = { 0.0, 1.0, 0.0, 1.0, 0.0, 1.0 };
   AstUnitMap *um;
   int ok = 1;

   astBoundingBox( plot, lbnd, ubnd );
   if( astOK ) ok = 0;
   log_refused( "astBoundingBox" );

   astClip( plot, AST__CURRENT, dlbnd, dubnd );
   if( astOK ) ok = 0;
   log_refused( "astClip" );

   astCurve( plot, start, finish );
   if( astOK ) ok = 0;
   log_refused( "astCurve" );

   um = astUnitMap( 3, " " );
   astGenCurve( plot, um );
   if( astOK ) ok = 0;
   log_refused( "astGenCurve" );
   um = astAnnul( um );

   astGridLine( plot, 1, start, 1.0 );
   if( astOK ) ok = 0;
   log_refused( "astGridLine" );

   astPolyCurve( plot, 2, 3, 2, poly );
   if( astOK ) ok = 0;
   log_refused( "astPolyCurve" );

   return ok;
}

/* Perform one action; returns 0 on a failure. */
static int do_action( AstPlot3D *plot, AstObject *input, const char *action ) {
   if( strcmp( action, "grid" ) == 0 ) {
      log_action( "grid" );
      astGrid( plot );

   } else if( strcmp( action, "border" ) == 0 ) {
      log_action( "border" );
      astBorder( plot );

   } else if( strcmp( action, "outline" ) == 0 ) {
      log_action( "outline" );
      if( !astIsARegion( input ) ) return fail( "outline needs a Region input" ) == 0;
      astRegionOutline( plot, (AstRegion *) input );

   } else if( strcmp( action, "refused" ) == 0 ) {
      log_action( "refused" );
      if( !do_refused( plot ) ) return fail( "a refused method did not report an error" ) == 0;

   } else if( strncmp( action, "text:", 5 ) == 0 ) {
      char text[256];
      double pos[3], upd[3];
      float up[3];
      const char *p = action + 5;
      const char *colon = strchr( p, ':' );
      size_t i, len;
      if( !colon ) return fail( "bad text action" ) == 0;
      len = (size_t)( colon - p );
      if( len >= sizeof(text) ) return fail( "text too long" ) == 0;
      memcpy( text, p, len );
      text[len] = '\0';
      for( i = 0; i < len; i++ ) if( text[i] == '~' ) text[i] = ' ';
      p = parse_triple( colon + 1, pos );
      if( !p || *p != ':' ) return fail( "bad text position" ) == 0;
      p = parse_triple( p + 1, upd );
      if( !p || *p != ':' ) return fail( "bad text up vector" ) == 0;
      up[0] = (float) upd[0];
      up[1] = (float) upd[1];
      up[2] = (float) upd[2];
      log_action( action );
      astText( plot, text, pos, up, p + 1 );

   } else if( strncmp( action, "mark:", 5 ) == 0 ) {
      double in[3*64];
      int type, n = 0;
      char *end;
      const char *p = action + 5;
      type = (int) strtol( p, &end, 10 );
      if( end == p || *end != ':' ) return fail( "bad mark type" ) == 0;
      p = end + 1;
      while( *p && n < 64 ) {
         double pos[3];
         p = parse_triple( p, pos );
         if( !p ) return fail( "bad mark position" ) == 0;
         in[n] = pos[0];
         in[64 + n] = pos[1];
         in[128 + n] = pos[2];
         n++;
         if( *p == '/' ) p++;
      }
      log_action( action );
      astMark( plot, n, 3, 64, in, type );

   } else if( strncmp( action, "set:", 4 ) == 0 ) {
      log_action( action );
      astSet( plot, "%s", action + 4 );

   } else if( strncmp( action, "clear:", 6 ) == 0 ) {
      log_action( action );
      astClear( plot, action + 6 );

   } else if( strncmp( action, "setcurrent:", 11 ) == 0 ) {
      log_action( action );
      astSetI( plot, "Current", atoi( action + 11 ) );

   } else if( strcmp( action, "clearcurrent" ) == 0 ) {
      log_action( action );
      astClear( plot, "Current" );

   } else if( strncmp( action, "removeframe:", 12 ) == 0 ) {
      log_action( action );
      astRemoveFrame( plot, atoi( action + 12 ) );

   } else {
      printf( "testplot3dlog: unknown action %s\n", action );
      return 0;
   }

   if( !astOK ) {
      char what[300];
      snprintf( what, sizeof(what), "error in action %s", action );
      return fail( what ) == 0;
   }
   return 1;
}

static void usage( void ) {
   printf( "Usage: testplot3dlog --in FILE [--attr STR] --box \"x1 y1 z1 x2 y2 z2\"\n"
           "       --bbox \"x1 y1 z1 x2 y2 z2\" [--mm N] [--log FILE] [--dump FILE]\n"
           "       [--svg FILE [W H]] [--camera \"ex ey ez tx ty tz ux uy uz screen\"]\n"
           "       --do ACTION [--do ACTION ...]\n" );
}

int main( int argc, char **argv ) {
   const char *infile = NULL, *attr = " ", *logfile = NULL, *dumpfile = NULL;
   const char *svgfile = NULL, *camera = NULL;
   const char *actions[64];
   int nact = 0;
   double box[6], bbox[6], mm = 0.0;
   int have_box = 0, have_bbox = 0;
   int svg_w = 720, svg_h = 540;
   FILE *svgfp = NULL;
   AstObject *input, *copy;
   AstPlot3D *plot;
   float gbox[6];
   char *dump1, *dump2;
   int i, rc = 0;

   astWatch( &status );
   astSetPutErr( collect_error );
   errbuf[0] = '\0';

   for( i = 1; i < argc; i++ ) {
      const char *a = argv[i];
      if( strcmp( a, "--in" ) == 0 && i + 1 < argc ) {
         infile = argv[++i];
      } else if( strcmp( a, "--attr" ) == 0 && i + 1 < argc ) {
         attr = argv[++i];
      } else if( strcmp( a, "--box" ) == 0 && i + 1 < argc ) {
         have_box = parse_doubles( argv[++i], 6, box );
      } else if( strcmp( a, "--bbox" ) == 0 && i + 1 < argc ) {
         have_bbox = parse_doubles( argv[++i], 6, bbox );
      } else if( strcmp( a, "--mm" ) == 0 && i + 1 < argc ) {
         mm = atof( argv[++i] );
      } else if( strcmp( a, "--log" ) == 0 && i + 1 < argc ) {
         logfile = argv[++i];
      } else if( strcmp( a, "--dump" ) == 0 && i + 1 < argc ) {
         dumpfile = argv[++i];
      } else if( strcmp( a, "--svg" ) == 0 && i + 1 < argc ) {
         svgfile = argv[++i];
         if( i + 2 < argc && argv[i+1][0] != '-' && argv[i+2][0] != '-' ) {
            svg_w = atoi( argv[++i] );
            svg_h = atoi( argv[++i] );
         }
      } else if( strcmp( a, "--camera" ) == 0 && i + 1 < argc ) {
         camera = argv[++i];
      } else if( strcmp( a, "--do" ) == 0 && i + 1 < argc && nact < 64 ) {
         actions[nact++] = argv[++i];
      } else {
         usage();
         return 1;
      }
   }
   if( !infile || !have_box || !have_bbox || nact == 0 ) {
      usage();
      return 1;
   }

   if( logfile ) {
      logfp = fopen( logfile, "w" );
      if( !logfp ) return fail( "cannot open the log file" );
   }
   if( svgfile ) {
      svgfp = fopen( svgfile, "w" );
      if( !svgfp ) return fail( "cannot open the SVG file" );
   }

   astGrfLogInit( NULL );
   astGrf3DLogInit( logfp );
   if( mm > 0.0 ) astGrf3DLogSetScale( mm );

   for( i = 0; i < 6; i++ ) gbox[i] = (float) box[i];

   if( svgfp ) {
      float lbnd[3] = { gbox[0], gbox[1], gbox[2] };
      float ubnd[3] = { gbox[3], gbox[4], gbox[5] };
      astGrf3DLogSetSvg( svgfp, svg_w, svg_h );
      if( !astGrf3DLogAutoCamera( lbnd, ubnd ) ) return fail( "auto camera failed" );
      if( camera ) {
         double c[10];
         float eye[3], target[3], up[3];
         if( !parse_doubles( camera, 10, c ) ) return fail( "bad --camera" );
         for( i = 0; i < 3; i++ ) {
            eye[i] = (float) c[i];
            target[i] = (float) c[3+i];
            up[i] = (float) c[6+i];
         }
         if( !astGrf3DLogSetCamera( eye, target, up, (float) c[9] ) ) {
            return fail( "camera does not define a view" );
         }
      }
   }

   astBegin;

   input = read_input( infile );
   if( !input || !astOK ) return fail( "cannot read the input object" );

   plot = astPlot3D( input, gbox, bbox, "%s", attr );
   if( !astOK ) return fail( "cannot create the Plot3D" );

   /* The Plot3D must survive a serialisation round trip unchanged. */
   dump1 = astToString( plot );
   copy = astFromString( dump1 );
   if( !copy || !astOK ) return fail( "cannot read back the Plot3D dump" );
   if( !astEqual( plot, copy ) ) return fail( "round-tripped Plot3D is not astEqual" );
   dump2 = astToString( copy );
   if( !dump2 || strcmp( dump1, dump2 ) != 0 ) {
      return fail( "round-tripped Plot3D dumps differently" );
   }
   if( dumpfile ) {
      FILE *fp = fopen( dumpfile, "w" );
      if( !fp ) return fail( "cannot open the dump file" );
      fputs( dump1, fp );
      fclose( fp );
   }
   dump1 = astFree( dump1 );
   dump2 = astFree( dump2 );
   copy = astAnnul( copy );

   for( i = 0; i < nact && rc == 0; i++ ) {
      if( !do_action( plot, input, actions[i] ) ) rc = 1;
   }

   astGrf3DLogClose();
   if( logfp ) fclose( logfp );
   if( svgfp ) fclose( svgfp );

   astEnd;

   if( rc == 0 && !astOK ) rc = fail( "AST error at exit" );
   if( rc == 0 ) printf( "testplot3dlog: all actions completed\n" );
   return rc;
}
