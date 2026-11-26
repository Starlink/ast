/*
*  Name:
*     testgrf3dlog.c
*
*  Purpose:
*     Unit test of the logging grf3d module (grf3d_log.c), independent
*     of the Plot3D class.
*
*  Description:
*     Calls each grf3d function directly and checks the line it writes
*     to the log, its return value, and the values it computes.  The
*     module is the oracle that other implementations of grf3d
*     (in particular the Rust port) are checked against, so the
*     expected strings and numbers here pin its behaviour exactly.
*
*     The 2D logging module (grf_log.c) is linked as well, to show that
*     the two modules keep separate attribute state.
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ast.h"
#include "ast_err.h"
#include "grf.h"
#include "grf3d.h"
#include "grf_log.h"
#include "grf3d_log.h"

static int status = 0;
static int nfail = 0;

static FILE *logfp = NULL;
static long rdpos = 0;

static void fail( const char *what, const char *detail ) {
   nfail++;
   printf( "FAIL: %s: %s\n", what, detail );
}

/* Return the next unread line of the log, without its newline, or ""
   when nothing has been written since the last read. */
static const char *next_line( void ) {
   static char buf[1024];
   buf[0] = '\0';
   fflush( logfp );
   fseek( logfp, rdpos, SEEK_SET );
   if( fgets( buf, sizeof(buf), logfp ) ) {
      size_t len = strlen( buf );
      if( len > 0 && buf[len-1] == '\n' ) buf[len-1] = '\0';
   }
   rdpos = ftell( logfp );
   fseek( logfp, 0, SEEK_END );
   return buf;
}

static void expect_line( const char *what, const char *expected ) {
   const char *got = next_line();
   if( strcmp( got, expected ) != 0 ) {
      char detail[2200];
      snprintf( detail, sizeof(detail), "\n   expected: %s\n   got:      %s",
                expected, got );
      fail( what, detail );
   }
}

static void expect_int( const char *what, int got, int expected ) {
   if( got != expected ) {
      char detail[200];
      snprintf( detail, sizeof(detail), "expected %d, got %d", expected, got );
      fail( what, detail );
   }
}

static void expect_near( const char *what, double got, double expected ) {
   if( fabs( got - expected ) > 1e-6 ) {
      char detail[200];
      snprintf( detail, sizeof(detail), "expected %.9g, got %.9g",
                expected, got );
      fail( what, detail );
   }
}

static void expect_error( const char *what, int expected ) {
   if( status != expected ) {
      char detail[200];
      snprintf( detail, sizeof(detail), "expected status %d, got %d",
                expected, status );
      fail( what, detail );
   }
   astClearStatus;
}

/* Attributes, capabilities, flush and character height. */
static void test_attr_cap_qch( void ) {
   double old = -1.0;
   float ch = 0.0f;

   expect_int( "G3DAttr returns 1",
               astG3DAttr( GRF__COLOUR, 3.0, &old, GRF__LINE ), 1 );
   expect_near( "G3DAttr old value is the default", old, 1.0 );
   expect_line( "G3DAttr log", "G3DAttr COLOUR=3.0000 prim=LINE" );

   expect_int( "G3DAttr query returns 1",
               astG3DAttr( GRF__COLOUR, AST__BAD, &old, GRF__LINE ), 1 );
   expect_near( "G3DAttr query reads the new value", old, 3.0 );
   expect_line( "G3DAttr query logs nothing", "" );

   expect_int( "G3DAttr size", astG3DAttr( GRF__SIZE, 2.0, NULL, GRF__TEXT ), 1 );
   expect_line( "G3DAttr size log", "G3DAttr SIZE=2.0000 prim=TEXT" );

   expect_int( "G3DQch returns 1", astG3DQch( &ch ), 1 );
   expect_near( "G3DQch with size 2 at 200 mm per unit", ch, 8.0 / 200.0 );
   expect_line( "G3DQch log", "G3DQch ch=0.040000" );

   /* The 2D module has its own attribute arrays. */
   expect_int( "GAttr query", astGAttr( GRF__SIZE, AST__BAD, &old, GRF__TEXT ), 1 );
   expect_near( "2D size is untouched by the 3D module", old, 1.0 );

   astGrf3DLogSetScale( 100.0 );
   expect_int( "G3DQch returns 1", astG3DQch( &ch ), 1 );
   expect_near( "G3DQch with size 2 at 100 mm per unit", ch, 8.0 / 100.0 );
   expect_line( "G3DQch log", "G3DQch ch=0.080000" );
   astGrf3DLogSetScale( 200.0 );

   expect_int( "G3DAttr unknown attribute returns 0",
               astG3DAttr( 99, 1.0, NULL, GRF__TEXT ), 0 );
   expect_error( "G3DAttr unknown attribute status", AST__GRFER );
   expect_line( "G3DAttr unknown attribute logs nothing", "" );

   expect_int( "G3DCap ESC", astG3DCap( GRF__ESC, 1 ), 0 );
   expect_int( "G3DCap MJUST", astG3DCap( GRF__MJUST, 1 ), 0 );
   expect_int( "G3DCap SCALES", astG3DCap( GRF__SCALES, 1 ), 0 );
   expect_line( "G3DCap logs nothing", "" );

   expect_int( "G3DFlush returns 1", astG3DFlush(), 1 );
   expect_line( "G3DFlush log", "G3DFlush" );

   /* Init resets the attributes to their defaults. */
   astGrf3DLogInit( logfp );
   expect_int( "G3DQch after init", astG3DQch( &ch ), 1 );
   expect_near( "G3DQch default", ch, 4.0 / 200.0 );
   expect_line( "G3DQch default log", "G3DQch ch=0.020000" );
   expect_int( "G3DAttr colour after init",
               astG3DAttr( GRF__COLOUR, AST__BAD, &old, GRF__LINE ), 1 );
   expect_near( "colour reset by init", old, 1.0 );
}

/* Lines, markers and text. */
static void test_line_mark_text( void ) {
   float x[3] = { 0.0f, 1.0f, 0.5f };
   float y[3] = { 0.0f, 0.5f, 0.5f };
   float z[3] = { 0.0f, 0.25f, 0.125f };
   float norm[3] = { 0.0f, 0.0f, 1.0f };
   float ref[3] = { 0.1f, 0.2f, 0.3f };
   float up[3] = { 0.0f, 1.0f, 0.0f };

   expect_int( "G3DLine returns 1", astG3DLine( 3, x, y, z ), 1 );
   expect_line( "G3DLine log",
      "G3DLine n=3 (0.000,0.000,0.000) (1.000,0.500,0.250) (0.500,0.500,0.125)" );

   expect_int( "G3DMark returns 1", astG3DMark( 2, x, y, z, 2, norm ), 1 );
   expect_line( "G3DMark log",
      "G3DMark type=2 norm=(0.000,0.000,1.000) n=2 (0.000,0.000,0.000) (1.000,0.500,0.250)" );

   expect_int( "G3DText returns 1", astG3DText( "hi", ref, "CC", up, norm ), 1 );
   expect_line( "G3DText log",
      "G3DText \"hi\" ref=(0.100,0.200,0.300) just=CC up=(0.000,1.000,0.000) norm=(0.000,0.000,1.000)" );

   expect_int( "G3DText NULL just returns 1",
               astG3DText( "hi", ref, NULL, up, norm ), 1 );
   expect_line( "G3DText NULL just log",
      "G3DText \"hi\" ref=(0.100,0.200,0.300) just=?? up=(0.000,1.000,0.000) norm=(0.000,0.000,1.000)" );
}

/* Close detaches the log. */
static void test_close( void ) {
   astGrf3DLogClose();
   expect_int( "G3DFlush after close", astG3DFlush(), 1 );
   expect_line( "nothing logged after close", "" );
   astGrf3DLogInit( logfp );
}

int main( void ) {
   astWatch( &status );

   logfp = tmpfile();
   if( !logfp ) {
      printf( "cannot create a temporary log file\n" );
      return 1;
   }

   astGrfLogInit( NULL );
   astGrf3DLogInit( logfp );

   test_attr_cap_qch();
   test_line_mark_text();
   test_close();

   astGrf3DLogClose();
   fclose( logfp );

   if( nfail == 0 ) {
      printf( "All grf3d_log tests passed\n" );
      return 0;
   }
   printf( "%d grf3d_log test(s) failed\n", nfail );
   return 1;
}
