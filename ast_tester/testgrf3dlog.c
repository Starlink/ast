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

/* Return the next unread line of "fp", without its newline, or "" when
   nothing has been written since the last read.  "*pos" remembers the
   read position between calls. */
static const char *next_line_from( FILE *fp, long *pos ) {
   static char buf[1024];
   buf[0] = '\0';
   fflush( fp );
   fseek( fp, *pos, SEEK_SET );
   if( fgets( buf, sizeof(buf), fp ) ) {
      size_t len = strlen( buf );
      if( len > 0 && buf[len-1] == '\n' ) buf[len-1] = '\0';
   }
   *pos = ftell( fp );
   fseek( fp, 0, SEEK_END );
   return buf;
}

static const char *next_line( void ) {
   return next_line_from( logfp, &rdpos );
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

/* SVG output through a camera. */
static FILE *svgfp = NULL;
static long svgpos = 0;

static const char *next_svg( void ) {
   return next_line_from( svgfp, &svgpos );
}

static void expect_svg( const char *what, const char *expected ) {
   const char *got = next_svg();
   if( strcmp( got, expected ) != 0 ) {
      char detail[2200];
      snprintf( detail, sizeof(detail), "\n   expected: %s\n   got:      %s",
                expected, got );
      fail( what, detail );
   }
}

static void expect_contains( const char *what, const char *line,
                             const char *needle, int want ) {
   int has = strstr( line, needle ) != NULL;
   if( has != want ) {
      char detail[2200];
      snprintf( detail, sizeof(detail), "%s \"%s\" in: %s",
                want ? "missing" : "unexpected", needle, line );
      fail( what, detail );
   }
}

/* Read the integer attribute "name" from an SVG element line. */
static int svg_attr( const char *line, const char *name ) {
   char key[32];
   const char *p;
   snprintf( key, sizeof(key), " %s=\"", name );
   p = strstr( line, key );
   if( !p ) return -99999;
   return atoi( p + strlen( key ) );
}

/* Count the "x,y" pairs in a polyline points attribute and check that each
   lies on the page. */
static int polyline_points_on_page( const char *line ) {
   const char *p = strstr( line, "points=\"" );
   int n = 0;
   if( !p ) return -1;
   p += strlen( "points=\"" );
   while( *p && *p != '"' ) {
      int x, y, used;
      if( sscanf( p, "%d,%d%n", &x, &y, &used ) != 2 ) return -1;
      if( x < 0 || x > 720 || y < 0 || y > 540 ) return -2;
      n++;
      p += used;
      while( *p == ' ' ) p++;
   }
   return n;
}

static void test_svg( void ) {
   float lbnd[3] = { -1.0f, -1.0f, -1.0f };
   float ubnd[3] = { 1.0f, 1.0f, 1.0f };
   float eye[3] = { 0.0f, 0.0f, 10.0f };
   float target[3] = { 0.0f, 0.0f, 0.0f };
   float camup[3] = { 0.0f, 1.0f, 0.0f };
   float norm[3] = { 0.0f, 0.0f, 1.0f };
   float up[3] = { 0.0f, 1.0f, 0.0f };
   float ref[3] = { 0.0f, 0.0f, 0.0f };
   const char *line;
   long before;
   int cx0, cy0, cx1, cy1, cx2, cxnear, cxfar;

   svgfp = tmpfile();
   if( !svgfp ) {
      fail( "svg", "cannot create a temporary SVG file" );
      return;
   }

   astGrf3DLogSetSvg( svgfp, 720, 540 );
   expect_svg( "svg xml header", "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" );
   expect_svg( "svg root",
      "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"720\" height=\"540\" "
      "viewBox=\"0 0 720 540\">" );
   expect_svg( "svg background",
      "<rect width=\"100%\" height=\"100%\" fill=\"black\"/>" );

   /* The auto camera fits the cube on the page: all eight corners of the
      cube drawn as one polyline land inside it. */
   expect_int( "AutoCamera returns 1", astGrf3DLogAutoCamera( lbnd, ubnd ), 1 );
   {
      float x[8] = { -1, 1, 1, -1, -1, 1, 1, -1 };
      float y[8] = { -1, -1, 1, 1, -1, -1, 1, 1 };
      float z[8] = { -1, -1, -1, -1, 1, 1, 1, 1 };
      expect_int( "G3DLine with svg returns 1", astG3DLine( 8, x, y, z ), 1 );
      expect_line( "G3DLine log is unchanged by the svg",
         "G3DLine n=8 (-1.000,-1.000,-1.000) (1.000,-1.000,-1.000) "
         "(1.000,1.000,-1.000) (-1.000,1.000,-1.000) (-1.000,-1.000,1.000) "
         "(1.000,-1.000,1.000) (1.000,1.000,1.000) (-1.000,1.000,1.000)" );
      line = next_svg();
      expect_contains( "polyline element", line, "<polyline ", 1 );
      expect_contains( "polyline colour", line, "stroke=\"#ff0000\"", 1 );
      expect_int( "polyline has eight points on the page",
                  polyline_points_on_page( line ), 8 );
   }

   /* An orthographic camera on the z axis: x maps to the page x axis
      linearly, and y is unchanged along it. */
   expect_int( "SetCamera orthographic returns 1",
               astGrf3DLogSetCamera( eye, target, camup, 0.0f ), 1 );
   {
      float x[3] = { 0.0f, 0.5f, 1.0f };
      float y[3] = { 0.0f, 0.0f, 0.0f };
      float z[3] = { 0.0f, 0.0f, 0.0f };
      expect_int( "G3DMark with svg returns 1", astG3DMark( 3, x, y, z, 1, norm ), 1 );
      next_line();
      line = next_svg();
      expect_contains( "circle element", line, "<circle ", 1 );
      cx0 = svg_attr( line, "cx" );
      cy0 = svg_attr( line, "cy" );
      line = next_svg();
      cx1 = svg_attr( line, "cx" );
      cy1 = svg_attr( line, "cy" );
      line = next_svg();
      cx2 = svg_attr( line, "cx" );
      if( cx1 <= cx0 ) fail( "orthographic x", "+x does not move right" );
      expect_int( "orthographic y unchanged", cy1, cy0 );
      if( abs( ( cx2 - cx0 ) - 2 * ( cx1 - cx0 ) ) > 1 ) {
         fail( "orthographic x is linear", "1.0 is not twice 0.5" );
      }
   }

   /* A perspective camera: the same offset nearer the eye projects
      further from the centre. */
   expect_int( "SetCamera perspective returns 1",
               astGrf3DLogSetCamera( eye, target, camup, 10.0f ), 1 );
   {
      float x[2] = { 0.5f, 0.5f };
      float y[2] = { 0.0f, 0.0f };
      float z[2] = { 5.0f, -5.0f };
      expect_int( "G3DMark perspective returns 1", astG3DMark( 2, x, y, z, 1, norm ), 1 );
      next_line();
      cxnear = svg_attr( next_svg(), "cx" );
      cxfar = svg_attr( next_svg(), "cx" );
      if( !( cxnear - cx0 > cxfar - cx0 && cxfar > cx0 ) ) {
         fail( "perspective", "nearer point does not project further out" );
      }
   }

   /* Text facing the camera, baseline along +x: no rotation, and the
      justification maps to the SVG anchor and baseline. */
   expect_int( "SetCamera orthographic again",
               astGrf3DLogSetCamera( eye, target, camup, 0.0f ), 1 );
   expect_int( "G3DText with svg returns 1",
               astG3DText( "ab", ref, "BL", up, norm ), 1 );
   next_line();
   line = next_svg();
   expect_contains( "text element", line, "<text ", 1 );
   expect_contains( "text anchor L", line, "text-anchor=\"start\"", 1 );
   expect_contains( "text baseline B", line, "dominant-baseline=\"alphabetic\"", 1 );
   expect_contains( "text not rotated", line, "rotate(", 0 );
   expect_contains( "text content", line, ">ab</text>", 1 );
   expect_contains( "text font size", line, "font-size=\"", 1 );
   if( svg_attr( line, "font-size" ) <= 0 ) fail( "font size", "not positive" );
   /* The text is placed where a marker at the same point was drawn. */
   expect_int( "text x at the reference point", svg_attr( line, "x" ), cx0 );
   expect_int( "text y at the reference point", svg_attr( line, "y" ), cy0 );

   /* Text with its up vector along +x: the baseline runs down the page. */
   {
      float upx[3] = { 1.0f, 0.0f, 0.0f };
      expect_int( "G3DText rotated returns 1",
                  astG3DText( "ab", ref, "BL", upx, norm ), 1 );
      next_line();
      line = next_svg();
      expect_contains( "text rotated 90", line, "rotate(90,", 1 );
   }

   /* Text seen from behind: the normal points away from the eye, so the
      text is flipped about its up vector and L becomes R. */
   {
      float back[3] = { 0.0f, 0.0f, -1.0f };
      expect_int( "G3DText from behind returns 1",
                  astG3DText( "ab", ref, "BL", up, back ), 1 );
      next_line();
      line = next_svg();
      expect_contains( "behind: anchor swapped", line, "text-anchor=\"end\"", 1 );
      expect_contains( "behind: baseline kept", line, "dominant-baseline=\"alphabetic\"", 1 );
      expect_contains( "behind: not rotated", line, "rotate(", 0 );
   }

   /* Text upside down: its up vector points down the page, so it is
      flipped about the normal and both justifications swap. */
   {
      float down[3] = { 0.0f, -1.0f, 0.0f };
      expect_int( "G3DText upside down returns 1",
                  astG3DText( "ab", ref, "BL", down, norm ), 1 );
      next_line();
      line = next_svg();
      expect_contains( "upside down: anchor swapped", line, "text-anchor=\"end\"", 1 );
      expect_contains( "upside down: baseline swapped", line, "dominant-baseline=\"hanging\"", 1 );
      expect_contains( "upside down: not rotated", line, "rotate(", 0 );
   }

   /* Empty text draws nothing. */
   fflush( svgfp );
   fseek( svgfp, 0, SEEK_END );
   before = ftell( svgfp );
   expect_int( "G3DText empty with svg returns 1", astG3DText( "", ref, "BL", up, norm ), 1 );
   next_line();
   fflush( svgfp );
   fseek( svgfp, 0, SEEK_END );
   if( ftell( svgfp ) != before ) fail( "empty text", "wrote to the svg" );

   /* Close finishes the picture and detaches it. */
   astGrf3DLogClose();
   expect_svg( "svg closed", "</svg>" );
   astGrf3DLogInit( logfp );
   fseek( svgfp, 0, SEEK_END );
   before = ftell( svgfp );
   {
      float x[2] = { 0.0f, 1.0f }, y[2] = { 0.0f, 0.0f }, z[2] = { 0.0f, 0.0f };
      expect_int( "G3DLine without svg returns 1", astG3DLine( 2, x, y, z ), 1 );
      expect_line( "G3DLine log without svg",
                   "G3DLine n=2 (0.000,0.000,0.000) (1.000,0.000,0.000)" );
   }
   fflush( svgfp );
   fseek( svgfp, 0, SEEK_END );
   if( ftell( svgfp ) != before ) fail( "no svg set", "still drawing after close" );

   fclose( svgfp );
   svgfp = NULL;
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
   test_svg();
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
