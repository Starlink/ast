/*
*  Name:
*     grf3d_log.h

*  Purpose:
*     Declare the set-up functions of the logging grf3d module.

*  Description:
*     The logging grf3d module (grf3d_log.c) implements the standard
*     3D graphics interface declared in grf3d.h with no external
*     dependencies.  Every call can be logged to a FILE for analysis or
*     regression testing, and the calls can optionally be drawn into an
*     SVG file through a camera.  This header declares the non-standard
*     set-up and shutdown functions; the standard grf3d entry points
*     (astG3DLine, astG3DText, etc.) are declared in grf3d.h.
*
*     The log and the text measurements never depend on the camera.
*     The camera is read only when an SVG file has been set.

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

#ifndef GRF3D_LOG_H_INCLUDED
#define GRF3D_LOG_H_INCLUDED

#include <stdio.h>

/* Start logging to "logfile" (NULL for no log) and reset the five
   graphics attributes of every primitive to their defaults.  The scale
   set by astGrf3DLogSetScale is not an attribute and is left alone. */
void astGrf3DLogInit( FILE *logfile );

/* Physical millimetres per 3D world unit, the same along every axis.
   The default of 200 gives a unit cube the character size of the 200 mm
   page used by the 2D logging module. */
void astGrf3DLogSetScale( double mm_per_world );

/* Draw every later line, marker and text into "svgfile" (NULL for no
   picture) on a page of the given size in pixels. */
void astGrf3DLogSetSvg( FILE *svgfile, int width, int height );

/* Place the camera used by the SVG: "screen" zero gives an orthographic
   projection, otherwise a perspective one with the screen at that
   distance from the eye.  Returns 1 on success, 0 if the eye, target
   and up vector do not define a view. */
int astGrf3DLogSetCamera( float eye[3], float target[3], float up[3],
                          float screen );

/* Place the camera so that the given 3D box fills the SVG page, and set
   the page window to the projected box.  Returns 1 on success. */
int astGrf3DLogAutoCamera( float lbnd[3], float ubnd[3] );

/* Set the page window so that the given 3D box, projected through the
   current camera, fills the SVG page with a 5 per cent margin.  Call it
   after astGrf3DLogSetCamera; astGrf3DLogAutoCamera does so itself.
   Returns 1 on success, 0 if there is no usable camera or the projected
   box has no area. */
int astGrf3DLogFitWindow( float lbnd[3], float ubnd[3] );

/* Flush and detach the log, and finish and detach the SVG. */
void astGrf3DLogClose( void );

#endif
