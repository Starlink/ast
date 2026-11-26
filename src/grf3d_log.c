/*
*  Name:
*     grf3d_log.c

*  Purpose:
*     Implement the grf3d module with call logging and an optional SVG
*     picture.

*  Description:
*     This file implements the 3D graphics functions declared in grf3d.h
*     against a virtual 3D space that needs no external graphics library.
*     Every call can be logged to a FILE for analysis or regression
*     testing, in the style of the 2D logging module grf_log.c.  The log
*     is the oracle that other implementations of the interface are
*     checked against, so every number is printed with a fixed format
*     and every arithmetic step is done in the order and precision
*     written here.
*
*     Text is measured in its own plane: the character height, the 0.6
*     width factor, the descender and the padding of grf_log.c are
*     applied along the text's baseline and up directions in 3D.  The
*     scale between the 3D world and physical millimetres is set with
*     astGrf3DLogSetScale and defaults to the 200 mm page of grf_log.c.
*
*     Lines, markers and text can also be drawn into an SVG file through
*     a camera, ported from the PLplot implementation grf3d_plplot.c.
*     The camera is read only when an SVG file is set; the log and the
*     text measurements never depend on it.

*  History:
*     26-NOV-2026 (TIMJ):
*        Original version.

*  Licence:
*     This program is free software: you can redistribute it and/or
*     modify it under the terms of the GNU Lesser General Public
*     License as published by the Free Software Foundation, either
*     version 3 of the License, or (at your option) any later
*     version.
*
*     This program is distributed in the hope that it will be useful,
*     but WITHOUT ANY WARRANTY; without even the implied warranty of
*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*     GNU Lesser General Public License for more details.
*
*     You should have received a copy of the GNU Lesser General
*     License along with this program.  If not, see
*     <http://www.gnu.org/licenses/>.
*/

/* Header files. */
/* ============= */
#include "grf3d.h"
#include "grf3d_log.h"
#include "error.h"
#include "ast_err.h"
#include "pointset.h"

#include <float.h>
#include <math.h>
#include <string.h>
#include <stdio.h>

/* Constants. */
/* ========== */
#define DEFAULT_MM_PER_WORLD 200.0  /* The 200 mm page of grf_log.c */
#define DEFAULT_CHR_HT  4.0         /* mm, matches PLplot default */

#define PAGE_HEIGHT_MM  200.0       /* The page the 14 pixel font is sized for */
#define SVG_FONT_PX     14.0        /* grf_log.c font size, in pixels */
#define CAMERA_OK 123456789         /* Flags that the camera has been set */
#define R2D 57.29577951308232       /* radians to degrees */

/* PLplot-compatible colour map 0 (indices 0-15), as in grf_log.c. */
static const char *cmap0[] = {
   "#000000", "#ff0000", "#ffff00", "#00ff00",
   "#7fffd4", "#ffc0cb", "#f5deb3", "#808080",
   "#a52a2a", "#0000ff", "#8a2be2", "#00ffff",
   "#f0e68c", "#fa8072", "#c0c0c0", "#ffffff"
};
#define NCMAP0 16

/* Type definitions. */
/* ================= */
/* The position and orientation of the camera in 3D world coordinates,
   as in grf3d_plplot.c. */
typedef struct Camera {
   float eye_vector[3];
   float target_vector[3];
   float up_vector[3];
   float w2c_matrix[9];
   float screen_distance;
   int ok_flag;
} Camera;

/* Static state. */
/* ============= */
static FILE *log_fp = NULL;
static double mm_per_world = DEFAULT_MM_PER_WORLD;

/* SVG output: the file, the page size in pixels, the camera, and the
   window in projected screen coordinates that maps onto the page. */
static FILE *svg_fp = NULL;
static int svg_w = 720;
static int svg_h = 540;
static Camera camera = { {0,0,0}, {0,0,0}, {0,0,0},
                         {0,0,0,0,0,0,0,0,0}, 0.0f, 0 };
static float win_x1 = 0.0f;
static float win_x2 = 1.0f;
static float win_y1 = 0.0f;
static float win_y2 = 1.0f;

/* Attribute state: [GRF__TEXT=0, GRF__LINE=1, GRF__MARK=2].  Separate
   from the arrays in grf_log.c, since a program may link both modules. */
static double attr_style[3]  = { 1.0, 1.0, 1.0 };
static double attr_width[3]  = { 1.0, 1.0, 1.0 };
static double attr_size[3]   = { 1.0, 1.0, 1.0 };
static double attr_font[3]   = { 1.0, 1.0, 1.0 };
static double attr_colour[3] = { 1.0, 1.0, 1.0 };

/* Initialisation / shutdown. */
/* ========================== */
void astGrf3DLogInit( FILE *logfile ) {
   int i;
   log_fp = logfile;
   for( i = 0; i < 3; i++ ) {
      attr_style[i]  = 1.0;
      attr_width[i]  = 1.0;
      attr_size[i]   = 1.0;
      attr_font[i]   = 1.0;
      attr_colour[i] = 1.0;
   }
}

void astGrf3DLogSetScale( double scale ) {
   mm_per_world = scale;
}

void astGrf3DLogSetSvg( FILE *svgfile, int width, int height ) {
   svg_fp = svgfile;
   svg_w = width;
   svg_h = height;
   if( svg_fp ) {
      fprintf( svg_fp,
         "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
         "<svg xmlns=\"http://www.w3.org/2000/svg\" "
         "width=\"%d\" height=\"%d\" "
         "viewBox=\"0 0 %d %d\">\n"
         "<rect width=\"100%%\" height=\"100%%\" fill=\"black\"/>\n",
         svg_w, svg_h, svg_w, svg_h );
   }
}

void astGrf3DLogClose( void ) {
   if( log_fp ) fflush( log_fp );
   log_fp = NULL;
   if( svg_fp ) {
      fprintf( svg_fp, "</svg>\n" );
      fflush( svg_fp );
   }
   svg_fp = NULL;
}

/* Helpers. */
/* ======== */

/* Current character height in mm, accounting for size scaling. */
static double ChrHtMM( void ) {
   double scale = attr_size[GRF__TEXT];
   if( scale <= 0.0 ) scale = 1.0;
   return DEFAULT_CHR_HT * scale;
}

/* Vector utilities, as in grf3d_plplot.c. */
/* Vector product of a pair of 3-vectors "a" and "b". */
static void vectorProduct( const float *a, const float *b, float *c ) {
   c[ 0 ] = a[ 1 ]*b[ 2 ] - a[ 2 ]*b[ 1 ];
   c[ 1 ] = a[ 2 ]*b[ 0 ] - a[ 0 ]*b[ 2 ];
   c[ 2 ] = a[ 0 ]*b[ 1 ] - a[ 1 ]*b[ 0 ];
}

/* Dot product of a pair of 3-vectors "a" and "b". */
static float dotProduct( const float *a, const float *b ) {
   return a[ 0 ]*b[ 0 ] + a[ 1 ]*b[ 1 ] + a[ 2 ]*b[ 2 ];
}

/* Vector from "b" to "a" (i.e. a minus b). */
static void vectorSub( const float *a, const float *b, float *c ) {
   c[ 0 ] = a[ 0 ] - b[ 0 ];
   c[ 1 ] = a[ 1 ] - b[ 1 ];
   c[ 2 ] = a[ 2 ] - b[ 2 ];
}

/* Return the length of a vector. */
static float vectorModulus( const float *a ) {
   return sqrtf( a[ 0 ]*a[ 0 ] + a[ 1 ]*a[ 1 ] + a[ 2 ]*a[ 2 ] );
}

/* Normalise a vector to unit length.  Returns 0 if the vector has zero
   length, and 1 otherwise. */
static int vectorNorm( float *a ) {
   float d;
   d = vectorModulus( a );
   if( d > 0.0 ) {
      a[ 0 ] /= d;
      a[ 1 ] /= d;
      a[ 2 ] /= d;
      return 1;
   } else {
      return 0;
   }
}

/* Unit vectors along the three axes of the text plane, expressed in 3D
   world coordinates: "tz" along the supplied normal, "tx" along the
   baseline (perpendicular to both the up vector and the normal) and "ty"
   the projection of the up vector onto the text plane.  Unlike the PLplot
   module, nothing here depends on a camera: a zero normal is an error,
   and the axes are never flipped to face an eye.  Returns 0 after
   reporting an error if the vectors do not define a plane. */
static int TextAxes( const char *caller, const float up[3], const float norm[3],
                     float tx[3], float ty[3], float tz[3] ) {
   tz[ 0 ] = norm[ 0 ];
   tz[ 1 ] = norm[ 1 ];
   tz[ 2 ] = norm[ 2 ];
   if( !vectorNorm( tz ) ) {
      astError( AST__GRFER, "%s: Zero length normal vector supplied.",
                caller );
      return 0;
   }

   vectorProduct( up, tz, tx );
   vectorProduct( tz, tx, ty );

   if( !vectorNorm( tx ) || !vectorNorm( ty ) ) {
      astError( AST__GRFER, "%s: Zero length up-vector supplied (or the "
                "up-vector is parallel to the normal vector).", caller );
      return 0;
   }
   return 1;
}

/* Camera and projection, ported from grf3d_plplot.c. */
/* ================================================== */

int astGrf3DLogSetCamera( float eye[3], float target[3], float up[3],
                          float screen ) {
   Camera *cam = &camera;
   float *u, *v, *w;
   int result = 1;

   memcpy( cam->target_vector, target, 3*sizeof( float ) );
   memcpy( cam->eye_vector, eye, 3*sizeof( float ) );
   cam->screen_distance = screen;

/* The rows of the w2c_matrix are unit vectors along the u, v and w axes
   of the camera, expressed in world coordinates. */
   u = cam->w2c_matrix;
   v = u + 3;
   w = v + 3;

/* The "w" axis points from the eye to the target. */
   vectorSub( target, eye, w );
   if( ! vectorNorm( w ) ) result = 0;

/* The "v" axis is in the plane of "w" and the up vector: take the normal
   to that plane (temporarily in "u"), then the vector normal to both it
   and "w". */
   vectorProduct( w, up, u );
   vectorProduct( u, w, v );
   if( ! vectorNorm( v ) ) result = 0;

/* The "u" axis is perpendicular to both. */
   vectorProduct( v, w, u );
   if( ! vectorNorm( u ) ) result = 0;

/* Store "v" as the camera's up vector. */
   memcpy( cam->up_vector, v, 3*sizeof( float ) );

   cam->ok_flag = result ? CAMERA_OK : CAMERA_OK/2;
   return result;
}

/* Project "n" world positions to screen coordinates (h, r).  Returns 0
   if no usable camera has been set. */
static int transform( int n, const float *x, const float *y, const float *z,
                      float *h, float *r ) {
   Camera *cam = &camera;
   float dx, dy, dz, u, v, w, f;
   int i;

   if( cam->ok_flag != CAMERA_OK ) return 0;

   for( i = 0; i < n; i++ ) {

/* Offset from the position to the eye. */
      dx = x[ i ] - (cam->eye_vector)[ 0 ];
      dy = y[ i ] - (cam->eye_vector)[ 1 ];
      dz = z[ i ] - (cam->eye_vector)[ 2 ];

/* The same vector in the (u,v,w) system. */
      u = (cam->w2c_matrix)[ 0 ]*dx +
          (cam->w2c_matrix)[ 1 ]*dy +
          (cam->w2c_matrix)[ 2 ]*dz;

      v = (cam->w2c_matrix)[ 3 ]*dx +
          (cam->w2c_matrix)[ 4 ]*dy +
          (cam->w2c_matrix)[ 5 ]*dz;

      w = (cam->w2c_matrix)[ 6 ]*dx +
          (cam->w2c_matrix)[ 7 ]*dy +
          (cam->w2c_matrix)[ 8 ]*dz;

/* Screen coordinates by a tangent plane or an orthographic projection.
   A position in the plane of the eye has no perspective projection and
   is flagged with FLT_MAX. */
      if( cam->screen_distance != 0.0 ) {
         if( w != 0.0 ) {
            f = cam->screen_distance/w;
            h[ i ] = -f*u;
            r[ i ] = f*v;
         } else {
            h[ i ] = FLT_MAX;
            r[ i ] = FLT_MAX;
         }
      } else {
         h[ i ] = -u;
         r[ i ] = v;
      }
   }
   return 1;
}

int astGrf3DLogAutoCamera( float lbnd[3], float ubnd[3] ) {
   float target[3], eye[3], up[3], screen, dx, dy, dz, hlo, hhi, rlo, rhi;
   float x[8], y[8], z[8], h[8], r[8];
   int result = 0;
   int i;

/* Look at the middle of the volume. */
   target[ 0 ] = 0.5*( lbnd[ 0 ] + ubnd[ 0 ] );
   target[ 1 ] = 0.5*( lbnd[ 1 ] + ubnd[ 1 ] );
   target[ 2 ] = 0.5*( lbnd[ 2 ] + ubnd[ 2 ] );

/* From slightly off a corner view... */
   eye[ 0 ] = 0.85*ubnd[ 0 ] + 0.15*lbnd[ 0 ];
   eye[ 1 ] = 0.75*ubnd[ 1 ] + 0.25*lbnd[ 1 ];
   eye[ 2 ] = 0.75*ubnd[ 2 ] + 0.25*lbnd[ 2 ];

/* ...seven times the size of the box away from its centre. */
   eye[ 0 ] = 7*(eye[ 0 ] - target[ 0 ] ) + target[ 0 ];
   eye[ 1 ] = 7*(eye[ 1 ] - target[ 1 ] ) + target[ 1 ];
   eye[ 2 ] = 7*(eye[ 2 ] - target[ 2 ] ) + target[ 2 ];

/* With the Z axis up. */
   up[ 0 ] = 0.0;
   up[ 1 ] = 0.0;
   up[ 2 ] = 1.0;

/* The screen passes through the centre of the box. */
   dx = eye[ 0 ] - target[ 0 ];
   dy = eye[ 1 ] - target[ 1 ];
   dz = eye[ 2 ] - target[ 2 ];
   screen = sqrtf( dx*dx + dy*dy + dz*dz );

   if( astGrf3DLogSetCamera( eye, target, up, screen ) ) {

/* The corners of the volume. */
      x[ 0 ] = ubnd[ 0 ];
      x[ 1 ] = ubnd[ 0 ];
      x[ 2 ] = lbnd[ 0 ];
      x[ 3 ] = lbnd[ 0 ];
      x[ 4 ] = ubnd[ 0 ];
      x[ 5 ] = ubnd[ 0 ];
      x[ 6 ] = lbnd[ 0 ];
      x[ 7 ] = lbnd[ 0 ];

      y[ 0 ] = lbnd[ 1 ];
      y[ 1 ] = ubnd[ 1 ];
      y[ 2 ] = ubnd[ 1 ];
      y[ 3 ] = lbnd[ 1 ];
      y[ 4 ] = lbnd[ 1 ];
      y[ 5 ] = ubnd[ 1 ];
      y[ 6 ] = ubnd[ 1 ];
      y[ 7 ] = lbnd[ 1 ];

      z[ 0 ] = lbnd[ 2 ];
      z[ 1 ] = lbnd[ 2 ];
      z[ 2 ] = lbnd[ 2 ];
      z[ 3 ] = lbnd[ 2 ];
      z[ 4 ] = ubnd[ 2 ];
      z[ 5 ] = ubnd[ 2 ];
      z[ 6 ] = ubnd[ 2 ];
      z[ 7 ] = ubnd[ 2 ];

      if( transform( 8, x, y, z, h, r ) ) {

/* The bounds of the projected volume. */
         hlo = FLT_MAX;
         hhi = -FLT_MAX;
         rlo = FLT_MAX;
         rhi = -FLT_MAX;

         for( i = 0; i < 8; i++ ) {
            if( h[ i ] < hlo ) hlo = h[ i ];
            if( h[ i ] > hhi ) hhi = h[ i ];
            if( r[ i ] < rlo ) rlo = r[ i ];
            if( r[ i ] > rhi ) rhi = r[ i ];
         }

/* Extended by 5% at each end. */
         dx = 0.05*( hhi - hlo );
         hhi += dx;
         hlo -= dx;

         dy = 0.05*( rhi - rlo );
         rhi += dy;
         rlo -= dy;

/* A box with area becomes the window mapped onto the page. */
         if( rlo < rhi && hlo < hhi ) {
            win_x1 = hlo;
            win_x2 = hhi;
            win_y1 = rlo;
            win_y2 = rhi;
            result = 1;
         }
      }
   }
   return result;
}

/* SVG helpers. */
/* ============ */

/* Page pixels from projected screen coordinates: the window maps onto
   the page as the unit square does in grf_log.c. */
static int SvgX( float h ) {
   return (int) round( ( h - win_x1 ) / ( win_x2 - win_x1 ) * (float)svg_w );
}

static int SvgY( float r ) {
   return (int) round( (float)svg_h -
                       ( r - win_y1 ) / ( win_y2 - win_y1 ) * (float)svg_h );
}

/* A projected position with no perspective image cannot be drawn. */
static int Projectable( float h, float r ) {
   return h != FLT_MAX && r != FLT_MAX;
}

static const char *SvgColour( int prim ) {
   int idx = (int)round( attr_colour[prim] );
   if( idx < 0 || idx >= NCMAP0 ) idx = 1;
   return cmap0[idx];
}

/* The text axes flipped so that the text reads normally from the camera,
   as getTextAxes does in grf3d_plplot.c, with the justification adjusted
   so that the text keeps its place on the screen.  "newjust" receives
   the two justification characters, defaulting to "CC". */
static void FlipTextAxes( const float ref[3], const char *just,
                          float tx[3], float ty[3], float tz[3],
                          char newjust[3] ) {
   float eye[3];

   if( just ) {
      strncpy( newjust, just, 2 );
      if( !newjust[ 0 ] ) newjust[ 0 ] = 'C';
      if( !newjust[ 1 ] ) newjust[ 1 ] = 'C';
      newjust[ 2 ] = 0;
   } else {
      strcpy( newjust, "CC" );
   }

/* The vector from the reference position to the eye. */
   vectorSub( camera.eye_vector, ref, eye );

/* Text seen from the back is spun about its up vector (ty), so the
   horizontal justification swaps. */
   if( dotProduct( tz, eye ) < 0.0 ) {
      tz[ 0 ] = -tz[ 0 ];
      tz[ 1 ] = -tz[ 1 ];
      tz[ 2 ] = -tz[ 2 ];
      tx[ 0 ] = -tx[ 0 ];
      tx[ 1 ] = -tx[ 1 ];
      tx[ 2 ] = -tx[ 2 ];

      if( newjust[ 1 ] == 'L' ) {
         newjust[ 1 ] = 'R';
      } else if( newjust[ 1 ] == 'R' ) {
         newjust[ 1 ] = 'L';
      }
   }

/* Text seen upside down is spun about the viewing vector (tz), so both
   justifications swap. */
   if( dotProduct( ty, camera.up_vector ) < 0.0 ) {
      ty[ 0 ] = -ty[ 0 ];
      ty[ 1 ] = -ty[ 1 ];
      ty[ 2 ] = -ty[ 2 ];
      tx[ 0 ] = -tx[ 0 ];
      tx[ 1 ] = -tx[ 1 ];
      tx[ 2 ] = -tx[ 2 ];

      if( newjust[ 0 ] == 'B' || newjust[ 0 ] == 'M' ) {
         newjust[ 0 ] = 'T';
      } else if( newjust[ 0 ] == 'T' ) {
         newjust[ 0 ] = 'M';
      }

      if( newjust[ 1 ] == 'L' ) {
         newjust[ 1 ] = 'R';
      } else if( newjust[ 1 ] == 'R' ) {
         newjust[ 1 ] = 'L';
      }
   }
}

/* Draw a text string into the SVG. */
static void SvgText( const char *text, const float ref[3], const char *just,
                     const float up[3], const float norm[3] ) {
   char newjust[3];
   float tx[3], ty[3], tz[3];
   float px[3], py[3], pz[3], h[3], r[3];
   float ch, dxpx, dypx, proj_px, nominal_px;
   double angle;
   const char *anchor;
   const char *baseline;
   int sx, sy, fs, iangle;

   if( !TextAxes( "astG3DText", up, norm, tx, ty, tz ) ) return;
   FlipTextAxes( ref, just, tx, ty, tz, newjust );

/* Project the reference point, a point one unit along the baseline, and
   a point one character height up the text. */
   ch = (float)( ChrHtMM() / mm_per_world );
   px[0] = ref[0];
   py[0] = ref[1];
   pz[0] = ref[2];
   px[1] = ref[0] + tx[0];
   py[1] = ref[1] + tx[1];
   pz[1] = ref[2] + tx[2];
   px[2] = ref[0] + ty[0] * ch;
   py[2] = ref[1] + ty[1] * ch;
   pz[2] = ref[2] + ty[2] * ch;
   if( !transform( 3, px, py, pz, h, r ) ) return;
   if( !Projectable( h[0], r[0] ) || !Projectable( h[1], r[1] ) ||
       !Projectable( h[2], r[2] ) ) return;

   sx = SvgX( h[0] );
   sy = SvgY( r[0] );

/* The baseline angle on the page, in degrees clockwise as SVG rotates. */
   dxpx = ( h[1] - h[0] ) / ( win_x2 - win_x1 ) * (float)svg_w;
   dypx = ( r[1] - r[0] ) / ( win_y2 - win_y1 ) * (float)svg_h;
   angle = -atan2( (double)dypx, (double)dxpx ) * R2D;

/* The font size: grf_log.c's 14 pixels scaled by the projected character
   height relative to the height a default character has on the 200 mm
   page at this page height. */
   dxpx = ( h[2] - h[0] ) / ( win_x2 - win_x1 ) * (float)svg_w;
   dypx = ( r[2] - r[0] ) / ( win_y2 - win_y1 ) * (float)svg_h;
   proj_px = sqrtf( dxpx*dxpx + dypx*dypx );
   nominal_px = (float)( DEFAULT_CHR_HT / PAGE_HEIGHT_MM ) * (float)svg_h;
   fs = (int)round( SVG_FONT_PX * proj_px / nominal_px );
   if( fs < 1 ) fs = 1;

   switch( newjust[1] ) {
      case 'L': anchor = "start";  break;
      case 'R': anchor = "end";    break;
      default:  anchor = "middle"; break;
   }
   switch( newjust[0] ) {
      case 'T': baseline = "hanging";    break;
      case 'B':
      case 'M': baseline = "alphabetic"; break;
      default:  baseline = "central";    break;
   }

   fprintf( svg_fp, "<text x=\"%d\" y=\"%d\" "
            "font-family=\"sans-serif\" font-size=\"%d\" "
            "fill=\"%s\" text-anchor=\"%s\" dominant-baseline=\"%s\"",
            sx, sy, fs, SvgColour(GRF__TEXT), anchor, baseline );
   iangle = (int)round( angle );
   if( iangle != 0 ) {
      fprintf( svg_fp, " transform=\"rotate(%d,%d,%d)\"", iangle, sx, sy );
   }
   fprintf( svg_fp, ">%s</text>\n", text );
}

/* The "i"th character (0 or 1) of a justification string for the log:
   "?" when the string is NULL or too short, as grf_log.c prints a NULL. */
static char JustChar( const char *just, int i ) {
   if( !just ) return '?';
   if( i == 1 && just[0] == '\0' ) return '?';
   return just[i] ? just[i] : '?';
}

/* Write "n" 3D points to the log as " (x,y,z)" groups. */
static void LogPoints( int n, const float *x, const float *y, const float *z ) {
   int i;
   if( x && y && z ) {
      for( i = 0; i < n; i++ ) {
         fprintf( log_fp, " (%.3f,%.3f,%.3f)", x[i], y[i], z[i] );
      }
   }
}

/* grf3d functions. */
/* ================ */

int astG3DAttr( int attr, double value, double *old_value, int prim ) {
   double *arr;
   int p;

   p = prim;
   if( p < 0 || p > 2 ) p = 0;

   if( attr == GRF__STYLE ) {
      arr = attr_style;
   } else if( attr == GRF__WIDTH ) {
      arr = attr_width;
   } else if( attr == GRF__SIZE ) {
      arr = attr_size;
   } else if( attr == GRF__FONT ) {
      arr = attr_font;
   } else if( attr == GRF__COLOUR ) {
      arr = attr_colour;
   } else {
      astError( AST__GRFER,
                "astG3DAttr: Unknown graphics attribute '%d' requested.",
                attr );
      return 0;
   }

   if( old_value ) *old_value = arr[p];

   if( value != AST__BAD ) {
      arr[p] = value;
      if( log_fp ) {
         static const char *aname[] = {
            "STYLE", "WIDTH", "SIZE", "FONT", "COLOUR" };
         static const char *pname[] = { "TEXT", "LINE", "MARK" };
         fprintf( log_fp, "G3DAttr %s=%.4f prim=%s\n",
                  aname[attr], value, pname[p] );
      }
   }

   return 1;
}

int astG3DCap( int cap, int value ) {
   (void)cap;
   (void)value;
   return 0;
}

int astG3DFlush( void ) {
   if( log_fp ) {
      fprintf( log_fp, "G3DFlush\n" );
      fflush( log_fp );
   }
   return 1;
}

int astG3DLine( int n, float *x, float *y, float *z ) {
   if( log_fp ) {
      fprintf( log_fp, "G3DLine n=%d", n );
      LogPoints( n, x, y, z );
      fprintf( log_fp, "\n" );
   }
   if( svg_fp && n > 1 && x && y && z && camera.ok_flag == CAMERA_OK ) {
      float h, r;
      int i;
      double w = attr_width[GRF__LINE];
      if( w < 1.0 ) w = 1.0;
      fprintf( svg_fp, "<polyline fill=\"none\" stroke=\"%s\" "
               "stroke-width=\"%d\" points=\"", SvgColour(GRF__LINE),
               (int)round( w ) );
      for( i = 0; i < n; i++ ) {
         if( transform( 1, x + i, y + i, z + i, &h, &r ) &&
             Projectable( h, r ) ) {
            fprintf( svg_fp, "%d,%d ", SvgX( h ), SvgY( r ) );
         }
      }
      fprintf( svg_fp, "\"/>\n" );
   }
   return 1;
}

int astG3DMark( int n, float *x, float *y, float *z, int type,
                float norm[3] ) {
   if( log_fp ) {
      fprintf( log_fp, "G3DMark type=%d norm=(%.3f,%.3f,%.3f) n=%d",
               type, norm[0], norm[1], norm[2], n );
      LogPoints( n, x, y, z );
      fprintf( log_fp, "\n" );
   }
   if( svg_fp && n > 0 && x && y && z ) {
      float h, r;
      int i;
      for( i = 0; i < n; i++ ) {
         if( transform( 1, x + i, y + i, z + i, &h, &r ) &&
             Projectable( h, r ) ) {
            fprintf( svg_fp, "<circle cx=\"%d\" cy=\"%d\" r=\"2\" "
                     "fill=\"%s\"/>\n", SvgX( h ), SvgY( r ),
                     SvgColour(GRF__MARK) );
         }
      }
   }
   return 1;
}

int astG3DQch( float *ch ) {
   *ch = (float)( ChrHtMM() / mm_per_world );
   if( log_fp ) fprintf( log_fp, "G3DQch ch=%.6f\n", *ch );
   return 1;
}

int astG3DText( const char *text, float ref[3], const char *just,
                float up[3], float norm[3] ) {
   if( log_fp ) {
      fprintf( log_fp, "G3DText \"%s\" ref=(%.3f,%.3f,%.3f) just=%c%c "
               "up=(%.3f,%.3f,%.3f) norm=(%.3f,%.3f,%.3f)\n",
               text ? text : "", ref[0], ref[1], ref[2],
               JustChar( just, 0 ), JustChar( just, 1 ),
               up[0], up[1], up[2], norm[0], norm[1], norm[2] );
   }
   if( svg_fp && text && text[0] != '\0' ) {
      SvgText( text, ref, just, up, norm );
   }
   return 1;
}

int astG3DTxExt( const char *text, float ref[3], const char *just,
                 float up[3], float norm[3], float *xb, float *yb,
                 float *zb, float bl[3] ) {
   char lj[2];
   int i;
   float tx[3], ty[3], tz[3];
   double chr_ht;

   for( i = 0; i < 4; i++ ) {
      xb[i] = 0.0f;
      yb[i] = 0.0f;
      zb[i] = 0.0f;
   }
   bl[0] = 0.0f;
   bl[1] = 0.0f;
   bl[2] = 0.0f;

   if( text && text[0] != '\0' ) {

/* Justification, as grf_log.c reads it: anything other than T, C or B
   vertically (so M, the bottom of the descenders, included) is C, and
   anything other than L, C or R horizontally is C.  A missing character
   is C as well. */
      {
         char j0 = JustChar( just, 0 );
         char j1 = JustChar( just, 1 );
         lj[0] = ( j0 == 'T' || j0 == 'C' || j0 == 'B' ) ? j0 : 'C';
         lj[1] = ( j1 == 'L' || j1 == 'C' || j1 == 'R' ) ? j1 : 'C';
      }

/* The axes of the text plane. */
      if( !TextAxes( "astG3DTxExt", up, norm, tx, ty, tz ) ) return 0;

/* The box in text-plane millimetres, with the expressions of grf_log.c:
   the text rises one character height above its baseline and descends
   0.2 of a height below it, each character is 0.6 of a height wide, and
   the width is padded by 0.2 of a height.  The plane is isotropic, so
   there is no separate scale along each axis and no sign flip of the up
   vector. */
      chr_ht = ChrHtMM();

      float phys_hu = (float)( 1.0 * chr_ht );
      float phys_hd = (float)( -0.2 * chr_ht );
      float phys_w  = (float)( strlen( text ) * chr_ht * 0.6 );
      float padded_w = phys_w + 0.2f * phys_hu;

/* The offset of the box centre from the reference point, in mm along the
   baseline (cx_mm) and the up direction (cy_mm). */
      float cx_mm = 0.0f, cy_mm = 0.0f;

      if( lj[0] == 'B' ) {
         cy_mm = 0.5f * ( phys_hu + phys_hd );
      } else if( lj[0] == 'T' ) {
         cy_mm = -0.5f * ( phys_hu - phys_hd );
      } else if( lj[0] == 'C' ) {
         cy_mm = 0.5f * phys_hd;
      }

      if( lj[1] == 'L' ) {
         cx_mm = 0.5f * padded_w;
      } else if( lj[1] == 'R' ) {
         cx_mm = -0.5f * padded_w;
      }

      float half_w_mm = 0.5f * padded_w;
      float half_h_mm = 0.5f * ( phys_hu - phys_hd );

/* World units per millimetre. */
      float s = (float)( 1.0 / mm_per_world );

/* The box centre in 3D world coordinates. */
      float cx = ref[0] + s * ( cx_mm * tx[0] + cy_mm * ty[0] );
      float cy = ref[1] + s * ( cx_mm * tx[1] + cy_mm * ty[1] );
      float cz = ref[2] + s * ( cx_mm * tx[2] + cy_mm * ty[2] );

/* The corners, anti-clockwise from the bottom left as seen from the
   front of the text plane. */
      xb[0] = cx - s * half_w_mm * tx[0] - s * half_h_mm * ty[0];
      yb[0] = cy - s * half_w_mm * tx[1] - s * half_h_mm * ty[1];
      zb[0] = cz - s * half_w_mm * tx[2] - s * half_h_mm * ty[2];
      xb[1] = cx + s * half_w_mm * tx[0] - s * half_h_mm * ty[0];
      yb[1] = cy + s * half_w_mm * tx[1] - s * half_h_mm * ty[1];
      zb[1] = cz + s * half_w_mm * tx[2] - s * half_h_mm * ty[2];
      xb[2] = cx + s * half_w_mm * tx[0] + s * half_h_mm * ty[0];
      yb[2] = cy + s * half_w_mm * tx[1] + s * half_h_mm * ty[1];
      zb[2] = cz + s * half_w_mm * tx[2] + s * half_h_mm * ty[2];
      xb[3] = cx - s * half_w_mm * tx[0] + s * half_h_mm * ty[0];
      yb[3] = cy - s * half_w_mm * tx[1] + s * half_h_mm * ty[1];
      zb[3] = cz - s * half_w_mm * tx[2] + s * half_h_mm * ty[2];

/* The left end of the baseline.  In text-plane mm relative to the
   reference point, the box spans cx_mm +/- half_w_mm along the baseline
   and cy_mm +/- half_h_mm up it.  Its bottom edge is the descender line,
   phys_hd below the baseline, so the baseline is at

      cy_mm - half_h_mm - phys_hd
    = cy_mm - 0.5*(phys_hu - phys_hd) - phys_hd
    = cy_mm - 0.5*(phys_hu + phys_hd)

   up the text, and its left end is cx_mm - half_w_mm along it.  For
   justification BL both are zero, so bl is the reference point. */
      float blx_mm = cx_mm - half_w_mm;
      float bly_mm = cy_mm - 0.5f * ( phys_hu + phys_hd );

      bl[0] = ref[0] + s * ( blx_mm * tx[0] + bly_mm * ty[0] );
      bl[1] = ref[1] + s * ( blx_mm * tx[1] + bly_mm * ty[1] );
      bl[2] = ref[2] + s * ( blx_mm * tx[2] + bly_mm * ty[2] );
   }

   if( log_fp ) {
      fprintf( log_fp, "G3DTxExt \"%s\" ref=(%.3f,%.3f,%.3f) just=%c%c "
               "box=(%.3f,%.3f,%.3f)-(%.3f,%.3f,%.3f)-(%.3f,%.3f,%.3f)-"
               "(%.3f,%.3f,%.3f) bl=(%.3f,%.3f,%.3f)\n",
               text ? text : "", ref[0], ref[1], ref[2],
               JustChar( just, 0 ), JustChar( just, 1 ),
               xb[0], yb[0], zb[0], xb[1], yb[1], zb[1],
               xb[2], yb[2], zb[2], xb[3], yb[3], zb[3],
               bl[0], bl[1], bl[2] );
   }
   return 1;
}
