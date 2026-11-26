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

/* Check the four corners and the baseline start of a text box. */
static void expect_box( const char *what, const float *xb, const float *yb,
                        const float *zb, const float *bl,
                        const float expected[5][3] ) {
   char name[200];
   int i;
   for( i = 0; i < 4; i++ ) {
      snprintf( name, sizeof(name), "%s corner %d x", what, i );
      expect_near( name, xb[i], expected[i][0] );
      snprintf( name, sizeof(name), "%s corner %d y", what, i );
      expect_near( name, yb[i], expected[i][1] );
      snprintf( name, sizeof(name), "%s corner %d z", what, i );
      expect_near( name, zb[i], expected[i][2] );
   }
   snprintf( name, sizeof(name), "%s bl x", what );
   expect_near( name, bl[0], expected[4][0] );
   snprintf( name, sizeof(name), "%s bl y", what );
   expect_near( name, bl[1], expected[4][1] );
   snprintf( name, sizeof(name), "%s bl z", what );
   expect_near( name, bl[2], expected[4][2] );
}

/* Text extent in the text's own plane.  With the default scale of 200 mm
   per world unit, a 4 mm character is 0.02 world units high, and the text
   "ab" has phys_hu = 4, phys_hd = -0.8, phys_w = 4.8, padded_w = 5.6 mm,
   so the box is 0.028 wide and 0.024 high in world units. */
static void test_txext( void ) {
   float ref[3] = { 0.0f, 0.0f, 0.0f };
   float up[3] = { 0.0f, 1.0f, 0.0f };
   float norm[3] = { 0.0f, 0.0f, 1.0f };
   float xb[4], yb[4], zb[4], bl[3];
   int i;

   /* BL: the reference point is the left end of the baseline. */
   {
      static const float e[5][3] = {
         { 0.0f,   -0.004f, 0.0f }, { 0.028f, -0.004f, 0.0f },
         { 0.028f,  0.020f, 0.0f }, { 0.0f,    0.020f, 0.0f },
         { 0.0f,    0.0f,   0.0f } };
      expect_int( "TxExt BL returns 1",
                  astG3DTxExt( "ab", ref, "BL", up, norm, xb, yb, zb, bl ), 1 );
      expect_box( "BL", xb, yb, zb, bl, e );
      expect_line( "TxExt BL log",
         "G3DTxExt \"ab\" ref=(0.000,0.000,0.000) just=BL "
         "box=(0.000,-0.004,0.000)-(0.028,-0.004,0.000)-(0.028,0.020,0.000)-"
         "(0.000,0.020,0.000) bl=(0.000,0.000,0.000)" );
   }

   /* CC: the box centre is at cy_mm = 0.5*phys_hd = -0.4 mm = -0.002, so
      the baseline is at -0.002 + 0.0012 - 0.0 ... i.e. at
      cy - 0.5*(hu+hd) = -0.4 - 1.6 = -2.0 mm = -0.010; and bl x is
      cx - half_w = -2.8 mm = -0.014. */
   {
      static const float e[5][3] = {
         { -0.014f, -0.014f, 0.0f }, { 0.014f, -0.014f, 0.0f },
         {  0.014f,  0.010f, 0.0f }, { -0.014f, 0.010f, 0.0f },
         { -0.014f, -0.010f, 0.0f } };
      expect_int( "TxExt CC returns 1",
                  astG3DTxExt( "ab", ref, "CC", up, norm, xb, yb, zb, bl ), 1 );
      expect_box( "CC", xb, yb, zb, bl, e );
      next_line();
   }

   /* NULL justification is CC. */
   {
      static const float e[5][3] = {
         { -0.014f, -0.014f, 0.0f }, { 0.014f, -0.014f, 0.0f },
         {  0.014f,  0.010f, 0.0f }, { -0.014f, 0.010f, 0.0f },
         { -0.014f, -0.010f, 0.0f } };
      expect_int( "TxExt NULL just returns 1",
                  astG3DTxExt( "ab", ref, NULL, up, norm, xb, yb, zb, bl ), 1 );
      expect_box( "NULL just", xb, yb, zb, bl, e );
      expect_line( "TxExt NULL just log",
         "G3DTxExt \"ab\" ref=(0.000,0.000,0.000) just=?? "
         "box=(-0.014,-0.014,0.000)-(0.014,-0.014,0.000)-(0.014,0.010,0.000)-"
         "(-0.014,0.010,0.000) bl=(-0.014,-0.010,0.000)" );
   }

   /* M is not a vertical justification this module knows, so it is C. */
   {
      static const float e[5][3] = {
         { -0.014f, -0.014f, 0.0f }, { 0.014f, -0.014f, 0.0f },
         {  0.014f,  0.010f, 0.0f }, { -0.014f, 0.010f, 0.0f },
         { -0.014f, -0.010f, 0.0f } };
      expect_int( "TxExt MC returns 1",
                  astG3DTxExt( "ab", ref, "MC", up, norm, xb, yb, zb, bl ), 1 );
      expect_box( "MC", xb, yb, zb, bl, e );
      next_line();
   }

   /* TR: the reference point is the top right corner. */
   {
      static const float e[5][3] = {
         { -0.028f, -0.024f, 0.0f }, { 0.0f, -0.024f, 0.0f },
         {  0.0f,    0.0f,   0.0f }, { -0.028f, 0.0f, 0.0f },
         { -0.028f, -0.020f, 0.0f } };
      expect_int( "TxExt TR returns 1",
                  astG3DTxExt( "ab", ref, "TR", up, norm, xb, yb, zb, bl ), 1 );
      expect_box( "TR", xb, yb, zb, bl, e );
      next_line();
   }

   /* TL and BR, mixing the two. */
   {
      static const float e[5][3] = {
         { 0.0f,   -0.024f, 0.0f }, { 0.028f, -0.024f, 0.0f },
         { 0.028f,  0.0f,   0.0f }, { 0.0f,    0.0f,   0.0f },
         { 0.0f,   -0.020f, 0.0f } };
      expect_int( "TxExt TL returns 1",
                  astG3DTxExt( "ab", ref, "TL", up, norm, xb, yb, zb, bl ), 1 );
      expect_box( "TL", xb, yb, zb, bl, e );
      next_line();
   }
   {
      static const float e[5][3] = {
         { -0.028f, -0.004f, 0.0f }, { 0.0f, -0.004f, 0.0f },
         {  0.0f,    0.020f, 0.0f }, { -0.028f, 0.020f, 0.0f },
         { -0.028f,  0.0f,   0.0f } };
      expect_int( "TxExt BR returns 1",
                  astG3DTxExt( "ab", ref, "BR", up, norm, xb, yb, zb, bl ), 1 );
      expect_box( "BR", xb, yb, zb, bl, e );
      next_line();
   }

   /* The reference point offsets every corner. */
   {
      float ref2[3] = { 1.0f, 2.0f, 3.0f };
      static const float e[5][3] = {
         { 1.0f,   1.996f, 3.0f }, { 1.028f, 1.996f, 3.0f },
         { 1.028f, 2.020f, 3.0f }, { 1.0f,   2.020f, 3.0f },
         { 1.0f,   2.0f,   3.0f } };
      expect_int( "TxExt offset ref returns 1",
                  astG3DTxExt( "ab", ref2, "BL", up, norm, xb, yb, zb, bl ), 1 );
      expect_box( "offset ref", xb, yb, zb, bl, e );
      next_line();
   }

   /* An up vector with a component along the normal: only its projection
      onto the text plane counts, so (0,1,1) with norm (0,0,1) is (0,1,0). */
   {
      float up2[3] = { 0.0f, 1.0f, 1.0f };
      static const float e[5][3] = {
         { 0.0f,   -0.004f, 0.0f }, { 0.028f, -0.004f, 0.0f },
         { 0.028f,  0.020f, 0.0f }, { 0.0f,    0.020f, 0.0f },
         { 0.0f,    0.0f,   0.0f } };
      expect_int( "TxExt tilted up returns 1",
                  astG3DTxExt( "ab", ref, "BL", up2, norm, xb, yb, zb, bl ), 1 );
      expect_box( "tilted up", xb, yb, zb, bl, e );
      next_line();
   }

   /* Reversing the normal reverses the baseline direction: tx = up x tz
      is (-1,0,0), so the box runs to negative x. */
   {
      float norm2[3] = { 0.0f, 0.0f, -1.0f };
      static const float e[5][3] = {
         { 0.0f,   -0.004f, 0.0f }, { -0.028f, -0.004f, 0.0f },
         { -0.028f, 0.020f, 0.0f }, { 0.0f,     0.020f, 0.0f },
         { 0.0f,    0.0f,   0.0f } };
      expect_int( "TxExt reversed normal returns 1",
                  astG3DTxExt( "ab", ref, "BL", up, norm2, xb, yb, zb, bl ), 1 );
      expect_box( "reversed normal", xb, yb, zb, bl, e );
      next_line();
   }

   /* A normal that is not a unit vector, along x, with up along z:
      tx = up x tz = (0,1,0), ty = tz x tx = (0,0,1). */
   {
      float norm3[3] = { 5.0f, 0.0f, 0.0f };
      float up3[3] = { 0.0f, 0.0f, 2.0f };
      static const float e[5][3] = {
         { 0.0f, 0.0f,   -0.004f }, { 0.0f, 0.028f, -0.004f },
         { 0.0f, 0.028f,  0.020f }, { 0.0f, 0.0f,    0.020f },
         { 0.0f, 0.0f,    0.0f } };
      expect_int( "TxExt x normal returns 1",
                  astG3DTxExt( "ab", ref, "BL", up3, norm3, xb, yb, zb, bl ), 1 );
      expect_box( "x normal", xb, yb, zb, bl, e );
      next_line();
   }

   /* The text size attribute scales the box. */
   {
      static const float e[5][3] = {
         { 0.0f,   -0.008f, 0.0f }, { 0.056f, -0.008f, 0.0f },
         { 0.056f,  0.040f, 0.0f }, { 0.0f,    0.040f, 0.0f },
         { 0.0f,    0.0f,   0.0f } };
      astG3DAttr( GRF__SIZE, 2.0, NULL, GRF__TEXT );
      next_line();
      expect_int( "TxExt size 2 returns 1",
                  astG3DTxExt( "ab", ref, "BL", up, norm, xb, yb, zb, bl ), 1 );
      expect_box( "size 2", xb, yb, zb, bl, e );
      next_line();
      astG3DAttr( GRF__SIZE, 1.0, NULL, GRF__TEXT );
      next_line();
   }

   /* A zero normal has no meaning without a camera. */
   expect_int( "TxExt zero normal returns 0",
               astG3DTxExt( "ab", ref, "BL", up, (float[3]){ 0.0f, 0.0f, 0.0f },
                            xb, yb, zb, bl ), 0 );
   expect_error( "TxExt zero normal status", AST__GRFER );

   /* An up vector parallel to the normal gives no baseline. */
   expect_int( "TxExt parallel up returns 0",
               astG3DTxExt( "ab", ref, "BL", norm, norm, xb, yb, zb, bl ), 0 );
   expect_error( "TxExt parallel up status", AST__GRFER );
   next_line();

   /* Empty and NULL text give an all-zero box without error. */
   for( i = 0; i < 4; i++ ) { xb[i] = yb[i] = zb[i] = 1.0f; }
   bl[0] = bl[1] = bl[2] = 1.0f;
   {
      static const float e[5][3] = { {0,0,0},{0,0,0},{0,0,0},{0,0,0},{0,0,0} };
      expect_int( "TxExt empty text returns 1",
                  astG3DTxExt( "", ref, "BL", up, norm, xb, yb, zb, bl ), 1 );
      expect_box( "empty text", xb, yb, zb, bl, e );
      expect_line( "TxExt empty text log",
         "G3DTxExt \"\" ref=(0.000,0.000,0.000) just=BL "
         "box=(0.000,0.000,0.000)-(0.000,0.000,0.000)-(0.000,0.000,0.000)-"
         "(0.000,0.000,0.000) bl=(0.000,0.000,0.000)" );
      for( i = 0; i < 4; i++ ) { xb[i] = yb[i] = zb[i] = 1.0f; }
      bl[0] = bl[1] = bl[2] = 1.0f;
      expect_int( "TxExt NULL text returns 1",
                  astG3DTxExt( NULL, ref, "BL", up, norm, xb, yb, zb, bl ), 1 );
      expect_box( "NULL text", xb, yb, zb, bl, e );
      next_line();
   }
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
   test_txext();
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
