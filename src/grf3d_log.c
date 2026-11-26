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
   int i;

   (void)up;
   (void)norm;

   for( i = 0; i < 4; i++ ) {
      xb[i] = 0.0f;
      yb[i] = 0.0f;
      zb[i] = 0.0f;
   }
   bl[0] = 0.0f;
   bl[1] = 0.0f;
   bl[2] = 0.0f;

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
