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

#include <math.h>
#include <string.h>
#include <stdio.h>

/* Constants. */
/* ========== */
#define DEFAULT_MM_PER_WORLD 200.0  /* The 200 mm page of grf_log.c */
#define DEFAULT_CHR_HT  4.0         /* mm, matches PLplot default */

/* Static state. */
/* ============= */
static FILE *log_fp = NULL;
static double mm_per_world = DEFAULT_MM_PER_WORLD;

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

void astGrf3DLogClose( void ) {
   if( log_fp ) fflush( log_fp );
   log_fp = NULL;
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
               just ? just[0] : '?', just ? just[1] : '?',
               up[0], up[1], up[2], norm[0], norm[1], norm[2] );
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
   anything other than L, C or R horizontally is C. */
      if( just ) {
         lj[0] = ( just[0] == 'T' || just[0] == 'C' || just[0] == 'B' )
                 ? just[0] : 'C';
         lj[1] = ( just[1] == 'L' || just[1] == 'C' || just[1] == 'R' )
                 ? just[1] : 'C';
      } else {
         lj[0] = 'C';
         lj[1] = 'C';
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
               just ? just[0] : '?', just ? just[1] : '?',
               xb[0], yb[0], zb[0], xb[1], yb[1], zb[1],
               xb[2], yb[2], zb[2], xb[3], yb[3], zb[3],
               bl[0], bl[1], bl[2] );
   }
   return 1;
}
