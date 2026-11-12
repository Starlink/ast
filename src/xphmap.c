/*
*class+
*  Name:
*     XphMap

*  Purpose:
*     Transform between different HEALPix projections.

*  Constructor Function:
*     astXphMap

*  Description:
*     The XphMap class implements a Mapping between grid coordinates
*     within different types of HEALPix projection: The available types
*     are :
*
*     - "HPX0":  An HPX projection centred on RA=0h.
*     - "HPX12": An HPX projection centred on RA=12h.
*     - "XPHN":  An XPH ("Butterfly") projection centred on the
*                 north pole.
*     - "XPHS":  An XPH ("Butterfly") projection centred on the
*                south pole.
*
*     The HPX projection is described in "Mapping on the HEALPix grid"
*     by Calabretta and Roukema, A&A, 2005. The XPH projection is
*     described in "Representing the 'butterfly' projection in FITS -
*     projection code XPH" by Calabretta and Lowe, PASA, 2018.
*
*     The forward transformation of an XphMap transforms from grid
*     coordinates within the nominated projection type into grid
*     coordinates within an HPX12 projection.

*  Inheritance:
*     The XphMap class inherits from the Mapping class.

*  Attributes:
*     The XphMap class does not define any new attributes beyond those
*     which are applicable to all Mappings.

*  Functions:
*     The XphMap class does not define any new functions beyond those
*     which are applicable to all Mappings.

*  Copyright:
*     Copyright (C) 2018 East Asian Observatory

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

*  Authors:
*     DSB: David S. Berry (EAO)

*  History:
*     18-OCT-2018 (DSB):
*        Original version.
*     8-SEP-2026 (TIMJ):
*        Free the projection type string read by astLoadXphMap.
*     11-NOV-2026 (TIMJ):
*        astLoadXphMap now checks that the dump describes a Mapping with
*        two inputs and two outputs. The transformation read or wrote beyond
*        the coordinate arrays otherwise.
*     11-NOV-2026 (TIMJ):
*        - The transformation now holds grid cell and facet indices in
*        doubles. They were held in integers, and for a position far from
*        the grid the conversion of the index to an integer, or the sum or
*        difference of two indices, was undefined. The result for a
*        position within range of an integer is unchanged.
*        - astLoadXphMap now rejects an Order outside the range accepted by
*        the constructor. The transformation uses it as a shift count.
*     11-NOV-2026 (TIMJ):
*        The forward transformation of the XPHN and XPHS types now gives
*        bad values for a position outside the grid. Each gore was tested
*        against its sides but not against the edge of the grid where it
*        ends, so a position beyond the grid but in line with a gore was
*        transformed to a position outside the HPX12 grid.
*     12-NOV-2026 (TIMJ):
*        - A position on the boundary of a used grid cell is now
*        transformed, using the transformation of that cell. A cell owns
*        only its low RA edges, so such a position was given bad values if
*        the cell on its high RA side is not used. That was so for a pole,
*        which is the corner of four cells none of which owned it, and for
*        every position on a meridian that divides two polar facets of an
*        HPX grid or the two halves of a facet of an XPH grid. The Moc
*        class uses those positions as test points, and so left out cells
*        beside such a meridian, and could not handle a Region centred on
*        such a position.
*        - A position within the margin that the HPX projection allows for
*        rounding (AST__HPXTOL degrees) of the boundary of a used cell is
*        treated as being on the boundary, since the grid coordinates of a
*        position on the edge of the sky are rarely exact.
*        - The inverse transformation now accepts a position in the unused
*        copy of a cell in the bottom left or top right facet of the HPX12
*        grid, and transforms the equivalent position in the used copy. It
*        gave bad values, although the sky position RA=0, Dec=0 and others
*        near that meridian are projected into the unused copy.
*        - A position less than one facet below or to the left of an HPX
*        grid now gives bad values. The index of its facet was rounded
*        towards zero, which put it in the first facet.
*class-
*/

/* Module Macros. */
/* ============== */
/* Set the name of the class we are implementing. This indicates to
   the header files that define class interfaces that they should make
   "protected" symbols available. */
#define astCLASS XphMap

/* Include files. */
/* ============== */
/* Interface definitions. */
/* ---------------------- */

#include "globals.h"             /* Thread-safe global data access */
#include "error.h"               /* Error reporting facilities */
#include "memory.h"              /* Memory allocation facilities */
#include "globals.h"             /* Thread-safe global data access */
#include "object.h"              /* Base Object class */
#include "pointset.h"            /* Sets of points/coordinates */
#include "mapping.h"             /* Coordinate mappings (parent class) */
#include "channel.h"             /* I/O channels */
#include "unitmap.h"             /* Unit Mappings */
#include "matrixmap.h"           /* Matrix Mappings */
#include "xphmap.h"              /* Interface definition for this class */

/* Error code definitions. */
/* ----------------------- */
#include "ast_err.h"             /* AST error codes */

/* C header files. */
/* --------------- */
#include <float.h>
#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

/* Module Variables. */
/* ================= */

/* Names and descriptions of each type of XphMap. These should correspond
   to the values of constants AST__HPX0, AST__HPX12, AST__XPHN and AST__XPHS
   in xphmap.h */
static const char *proj_name[] = { "HPX0", "HPX12", "XPHN", "XPHS" };
static const char *proj_comm[] = { "HPX projection centred on RA=0h",
                                   "HPX projection centred on RA=12h",
                                   "XPH projection centred on north pole",
                                   "XPH projection centred on south pole" };

/* Address of this static variable is used as a unique identifier for
   member of this class. */
static int class_check;

/* Pointers to parent class methods which are extended by this class. */
static AstPointSet *(* parent_transform)( AstMapping *, AstPointSet *, int, AstPointSet *, int * );

/* Define macros for accessing each item of thread specific global data. */
#ifdef THREAD_SAFE

/* Define how to initialise thread-specific globals. */
#define GLOBAL_inits \
   globals->Class_Init = 0;

/* Create the function that initialises global data for this module. */
astMAKE_INITGLOBALS(XphMap)

/* Define macros for accessing each item of thread specific global data. */
#define class_init astGLOBAL(XphMap,Class_Init)
#define class_vtab astGLOBAL(XphMap,Class_Vtab)
#define getattrib_buff astGLOBAL(XphMap,GetAttrib_Buff)



/* If thread safety is not needed, declare and initialise globals at static
   variables. */
#else

/* Define the class virtual function table and its initialisation flag
   as static variables. */
static AstXphMapVtab class_vtab;   /* Virtual function table */
static int class_init = 0;         /* Virtual function table initialised? */

#endif

/* External Interface Function Prototypes. */
/* ======================================= */
/* The following functions have public prototypes only (i.e. no
   protected prototypes), so we must provide local prototypes for use
   within this module. */
AstXphMap *astXphMapId_( int, int, const char *, ... );

/* Prototypes for Private Member Functions. */
/* ======================================== */
static AstPointSet *Transform( AstMapping *, AstPointSet *, int, AstPointSet *, int * );
static int EdgeCells( double, double, int, double[ 2 ], double[ 2 ] );
static int Equal( AstObject *, AstObject *, int * );
static int FindHpxCell( double *, double *, int, double, int, double * );
static int FindXphCell( int, double *, double *, int, double );
static int GetIsLinear( AstMapping *, int * );
static int HpxCellUsed( double, double, int, int, double * );
static int XphCellGore( int, double, double, int );
static int MapMerge( AstMapping *, int, int, int *, AstMapping ***, int **, int * );
static void Dump( AstObject *, AstChannel *, int * );

/* Member functions. */
/* ================= */

static int EdgeCells( double v, double tol, int upper, double cell[ 2 ],
                      double pos[ 2 ] ) {
/*
*  Name:
*     EdgeCells

*  Purpose:
*     Find the grid cells that hold a position on one axis.

*  Type:
*     Private function.

*  Synopsis:
*     #include "xphmap.h"
*     int EdgeCells( double v, double tol, int upper, double cell[ 2 ],
*                    double pos[ 2 ] )

*  Class Membership:
*     XphMap member function

*  Description:
*     This function returns the index on one axis of the grid cell that
*     owns a position. If the position is on the edge of the cell, it
*     also returns the index of the cell on the other side of the edge.

*  Parameters:
*     v
*        The grid coordinate of the position on the axis. The centre of
*        the first cell is at 1.0, and its lower edge is at 0.5.
*     tol
*        The position is on the edge of a cell if it is within this
*        distance of the edge. Must be less than 0.5.
*     upper
*        If non-zero, a cell owns its upper edge. If zero, a cell owns its
*        lower edge.
*     cell
*        Returned holding the zero-based indices of the cells. The first
*        element is the cell that owns the position.
*     pos
*        Returned holding the grid coordinate to use with each element of
*        "cell". The first element is "v". The second is the coordinate of
*        the edge, so that the position is within the second cell.

*  Returned Value:
*     The number of cells returned: two if the position is on the edge of
*     a cell and one otherwise.

*  Notes:
*     - The indices are whole numbers but are returned as doubles, since
*     the position can be any distance from the grid.
*     - A single cell is returned for a value that is not finite.
*/

/* Local Variables: */
   double edge;                  /* Coordinate of the nearest cell edge */

/* The cell that owns the position. */
   cell[ 0 ] = upper ? ceil( v - 1.5 ) : floor( v - 0.5 );
   pos[ 0 ] = v;

/* The edges of the cells are at half-integer coordinates. Find the edge
   that is nearest to the position. The cells on either side of it have
   zero-based indices of ( edge - 1.5 ) and ( edge - 0.5 ). */
   edge = floor( v ) + 0.5;
   if( fabs( v - edge ) <= tol ) {
      cell[ 1 ] = ( cell[ 0 ] < edge - 1.0 ) ? edge - 0.5 : edge - 1.5;
      pos[ 1 ] = edge;
      return 2;
   }

   return 1;
}

static int Equal( AstObject *this_object, AstObject *that_object, int *status ) {
/*
*  Name:
*     Equal

*  Purpose:
*     Test if two XphMaps are equivalent.

*  Type:
*     Private function.

*  Synopsis:
*     #include "xphmap.h"
*     int Equal( AstObject *this, AstObject *that, int *status )

*  Class Membership:
*     XphMap member function (over-rides the astEqual protected
*     method inherited from the astMapping class).

*  Description:
*     This function returns a boolean result (0 or 1) to indicate whether
*     two XphMaps are equivalent.

*  Parameters:
*     this
*        Pointer to the first Object (a XphMap).
*     that
*        Pointer to the second Object.
*     status
*        Pointer to the inherited status variable.

*  Returned Value:
*     One if the XphMaps are equivalent, zero otherwise.

*  Notes:
*     - A value of zero will be returned if this function is invoked
*     with the global status set, or if it should fail for any reason.
*/

/* Local Variables: */
   AstXphMap *that;
   AstXphMap *this;
   int result;

/* Initialise. */
   result = 0;

/* Check the global error status. */
   if ( !astOK ) return result;

/* Obtain pointers to the two XphMap structures. */
   this = (AstXphMap *) this_object;
   that = (AstXphMap *) that_object;

/* Check the second object is a XphMap. We know the first is a
   XphMap since we have arrived at this implementation of the virtual
   function. */
   if( astIsAXphMap( that ) ) {

/* Check the properties of the two XphMaps are equal. */
      if( astGetInvert( this ) == astGetInvert( that ) ) {
         result = ( this->type == that->type ) &&
                  ( this->order == that->order );
      }
   }

/* Return the result, */
   return result;
}

static int FindHpxCell( double *x, double *y, int nppf, double tol, int wrap,
                        double *s ) {
/*
*  Name:
*     FindHpxCell

*  Purpose:
*     Find the used cell of an HPX grid that holds a position.

*  Type:
*     Private function.

*  Synopsis:
*     #include "xphmap.h"
*     int FindHpxCell( double *x, double *y, int nppf, double tol, int wrap,
*                      double *s )

*  Class Membership:
*     XphMap member function

*  Description:
*     This function finds the used cell that holds a position in the grid
*     of an HPX projection (either HPX0 or HPX12). A cell owns its low RA
*     edges, which are its upper edges on both axes. If the cell that
*     owns the position is not used, the position may be on an edge of
*     that cell, in which case the cells on the other side of the edge
*     are tried.

*  Parameters:
*     x
*        Pointer to the grid X coordinate of the position.
*     y
*        Pointer to the grid Y coordinate of the position.
*     nppf
*        The number of pixels along one edge of a facet.
*     tol
*        The position is on the edge of a cell if it is within this
*        distance of the edge.
*     wrap
*        If non-zero, a position in the unused copy of a cell is moved to
*        the equivalent position in the used copy (see HpxCellUsed).
*     s
*        Returned holding the sum of the zero-based X and Y indices of the
*        used cell, plus one.

*  Returned Value:
*     Non-zero if the position is on the sky. Zero if it is not, in which
*     case the supplied values are left unchanged.

*  Notes:
*     - If the used cell is not the cell that owns the position, the
*     returned position is moved onto the edge of the used cell. It is
*     otherwise unchanged, unless "wrap" is non-zero and the position was
*     in an unused copy of a cell.
*/

/* Local Variables: */
   double ix[ 2 ];               /* Candidate cell indices on the X axis */
   double iy[ 2 ];               /* Candidate cell indices on the Y axis */
   double shift;                 /* Offset to the used copy of a cell */
   double xs[ 2 ];               /* X coordinate to use with each "ix" */
   double ys[ 2 ];               /* Y coordinate to use with each "iy" */
   int i;                        /* Index of candidate on the X axis */
   int j;                        /* Index of candidate on the Y axis */
   int nx;                       /* Number of candidates on the X axis */
   int ny;                       /* Number of candidates on the Y axis */

/* Get the cell that owns the position, and the cell on the other side of
   any edge that the position is on, for each axis. */
   nx = EdgeCells( *x, tol, 1, ix, xs );
   ny = EdgeCells( *y, tol, 1, iy, ys );

/* Try each combination, starting with the cell that owns the position. */
   for( j = 0; j < ny; j++ ) {
      for( i = 0; i < nx; i++ ) {
         if( HpxCellUsed( ix[ i ], iy[ j ], nppf, wrap, &shift ) ) {
            *x = xs[ i ] + shift;
            *y = ys[ j ] + shift;
            *s = ix[ i ] + iy[ j ] + 1 + 2*shift;
            return 1;
         }
      }
   }

   return 0;
}

static int FindXphCell( int type, double *x, double *y, int nppf, double tol ) {
/*
*  Name:
*     FindXphCell

*  Purpose:
*     Find the used cell of an XPH grid that holds a position.

*  Type:
*     Private function.

*  Synopsis:
*     #include "xphmap.h"
*     int FindXphCell( int type, double *x, double *y, int nppf, double tol )

*  Class Membership:
*     XphMap member function

*  Description:
*     This function finds the used cell that holds a position in the grid
*     of an XPH projection, and returns the index of the gore that holds
*     the cell.
*
*     The grid is divided into four quadrants that meet at the pole in
*     the centre of the grid, and each quadrant holds one gore. A gore
*     contains its low RA boundaries but not its high RA boundaries, and
*     a cell owns its low RA edges. Which edges these are depends on the
*     gore. If the cell that owns the position is not a used cell of the
*     gore, the position may be on an edge of that cell, in which case
*     the cells on the other side of the edge are tried. A cell in the
*     same gore is preferred to a cell in a different gore.

*  Parameters:
*     type
*        The type of projection: AST__XPHN or AST__XPHS.
*     x
*        Pointer to the grid X coordinate of the position.
*     y
*        Pointer to the grid Y coordinate of the position.
*     nppf
*        The number of pixels along one edge of a facet.
*     tol
*        The position is on the edge of a cell if it is within this
*        distance of the edge.

*  Returned Value:
*     The index of the gore that holds the used cell, in the range 0 to
*     3, or -1 if the position is not on the sky. The supplied values are
*     left unchanged if -1 is returned.

*  Notes:
*     - If the used cell is not the cell that owns the position, the
*     returned position is moved onto the edge of the used cell. It is
*     otherwise unchanged.
*/

/* Local Variables: */
   double ix[ 2 ];               /* Candidate cell indices on the X axis */
   double iy[ 2 ];               /* Candidate cell indices on the Y axis */
   double mid;                   /* Grid coordinate of the centre of the grid */
   double xs[ 2 ];               /* X coordinate to use with each "ix" */
   double ys[ 2 ];               /* Y coordinate to use with each "iy" */
   int g;                        /* Gore holding a candidate cell */
   int gore;                     /* Gore holding the position */
   int i;                        /* Index of candidate on the X axis */
   int j;                        /* Index of candidate on the Y axis */
   int nx;                       /* Number of candidates on the X axis */
   int ny;                       /* Number of candidates on the Y axis */
   int pass;                     /* Loop counter for passes */
   int upx;                      /* Does a cell own its upper X edge? */
   int upy;                      /* Does a cell own its upper Y edge? */

/* Find the gore that holds the position, and so which edges of a cell
   are its low RA edges. The pole at the centre of the grid is put in
   gore 1. */
   mid = 2*nppf + 0.5;
   if( type == AST__XPHN ) {
      if( *x < mid && *y >= mid ) {
         gore = 0;
         upx = 0;
         upy = 0;
      } else if( *x >= mid && *y > mid ) {
         gore = 3;
         upx = 0;
         upy = 1;
      } else if( *x > mid && *y <= mid ) {
         gore = 2;
         upx = 1;
         upy = 1;
      } else {
         gore = 1;
         upx = 1;
         upy = 0;
      }

   } else {
      if( *x < mid && *y <= mid ) {
         gore = 0;
         upx = 0;
         upy = 1;
      } else if( *x >= mid && *y < mid ) {
         gore = 3;
         upx = 0;
         upy = 0;
      } else if( *x > mid && *y >= mid ) {
         gore = 2;
         upx = 1;
         upy = 0;
      } else {
         gore = 1;
         upx = 1;
         upy = 1;
      }
   }

/* Get the cell that owns the position, and the cell on the other side of
   any edge that the position is on, for each axis. */
   nx = EdgeCells( *x, tol, upx, ix, xs );
   ny = EdgeCells( *y, tol, upy, iy, ys );

/* Try each combination, starting with the cell that owns the position.
   The first pass accepts only a cell in the gore found above, and the
   second pass accepts a cell in any gore. */
   for( pass = 0; pass < 2; pass++ ) {
      for( j = 0; j < ny; j++ ) {
         for( i = 0; i < nx; i++ ) {
            g = XphCellGore( type, ix[ i ], iy[ j ], nppf );
            if( g >= 0 && ( pass || g == gore ) ) {
               *x = xs[ i ];
               *y = ys[ j ];
               return g;
            }
         }
      }
   }

   return -1;
}

static int GetIsLinear( AstMapping *this_mapping, int *status ){
/*
*  Name:
*     GetIsLinear

*  Purpose:
*     Return the value of the IsLinear attribute for a XphMap.

*  Type:
*     Private function.

*  Synopsis:
*     #include "mapping.h"
*     void GetIsLinear( AstMapping *this, int *status )

*  Class Membership:
*     XphMap member function (over-rides the protected astGetIsLinear
*     method inherited from the Mapping class).

*  Description:
*     This function returns the value of the IsLinear attribute for a
*     Frame, which is always one.

*  Parameters:
*     this
*        Pointer to the XphMap.
*     status
*        Pointer to the inherited status variable.
*/
   return (((AstXphMap *) this_mapping)->type == AST__HPX12 );
}

void astInitXphMapVtab_(  AstXphMapVtab *vtab, const char *name, int *status ) {
/*
*+
*  Name:
*     astInitXphMapVtab

*  Purpose:
*     Initialise a virtual function table for a XphMap.

*  Type:
*     Protected function.

*  Synopsis:
*     #include "xphmap.h"
*     void astInitXphMapVtab( AstXphMapVtab *vtab, const char *name )

*  Class Membership:
*     XphMap vtab initialiser.

*  Description:
*     This function initialises the component of a virtual function
*     table which is used by the XphMap class.

*  Parameters:
*     vtab
*        Pointer to the virtual function table. The components used by
*        all ancestral classes will be initialised if they have not already
*        been initialised.
*     name
*        Pointer to a constant null-terminated character string which contains
*        the name of the class to which the virtual function table belongs (it
*        is this pointer value that will subsequently be returned by the Object
*        astClass function).
*-
*/

/* Local Variables: */
   astDECLARE_GLOBALS            /* Pointer to thread-specific global data */
   AstObjectVtab *object;        /* Pointer to Object component of Vtab */
   AstMappingVtab *mapping;      /* Pointer to Mapping component of Vtab */

/* Check the local error status. */
   if ( !astOK ) return;

/* Get a pointer to the thread specific global data structure. */
   astGET_GLOBALS(NULL);

/* Initialize the component of the virtual function table used by the
   parent class. */
   astInitMappingVtab( (AstMappingVtab *) vtab, name );

/* Store a unique "magic" value in the virtual function table. This
   will be used (by astIsAXphMap) to determine if an object belongs
   to this class.  We can conveniently use the address of the (static)
   class_check variable to generate this unique value. */
   vtab->id.check = &class_check;
   vtab->id.parent = &(((AstMappingVtab *) vtab)->id);

/* Initialise member function pointers. */
/* ------------------------------------ */
/* Store pointers to the member functions (implemented here) that provide
   virtual methods for this class. */

/* Save the inherited pointers to methods that will be extended, and
   replace them with pointers to the new member functions. */
   object = (AstObjectVtab *) vtab;
   mapping = (AstMappingVtab *) vtab;

   parent_transform = mapping->Transform;
   mapping->Transform = Transform;

/* Store replacement pointers for methods which will be over-ridden by
   new member functions implemented here. */
   object->Equal = Equal;
   mapping->MapMerge = MapMerge;
   mapping->GetIsLinear = GetIsLinear;

/* Declare the class dump function. There is no copy constructor or
   destructor. */
   astSetDump( vtab, Dump, "XphMap", "HPX variant mapping" );

/* If we have just initialised the vtab for the current class, indicate
   that the vtab is now initialised, and store a pointer to the class
   identifier in the base "object" level of the vtab. */
   if( vtab == &class_vtab ) {
      class_init = 1;
      astSetVtabClassIdentifier( vtab, &(vtab->id) );
   }
}

static int HpxCellUsed( double ix, double iy, int nppf, int wrap,
                        double *shift ) {
/*
*  Name:
*     HpxCellUsed

*  Purpose:
*     See if a cell of an HPX grid is used.

*  Type:
*     Private function.

*  Synopsis:
*     #include "xphmap.h"
*     int HpxCellUsed( double ix, double iy, int nppf, int wrap,
*                      double *shift )

*  Class Membership:
*     XphMap member function

*  Description:
*     This function tests whether a cell in the grid of an HPX projection
*     (either HPX0 or HPX12) is used.
*
*     The grid holds five facets along each axis, as shown in Fig 3 of
*     Calabretta & Roukema 2005, and the facets on its diagonal and on
*     either side of the diagonal cover the sky. The bottom left facet
*     and the top right facet are the same facet. Each of its cells is
*     in both, and only one copy of each cell is used: the copy in the
*     top right half of the bottom left facet or in the bottom left half
*     of the top right facet. The cells on the other diagonal of the
*     facet, which is the boundary between the two halves, are used in
*     the top right facet.

*  Parameters:
*     ix
*        The zero-based X index of the cell.
*     iy
*        The zero-based Y index of the cell.
*     nppf
*        The number of pixels along one edge of a facet.
*     wrap
*        If non-zero, an unused copy of a cell in the bottom left or top
*        right facet is accepted in place of the used copy.
*     shift
*        Returned holding the value to add to both grid coordinates of a
*        position in the cell to get the equivalent position in the used
*        copy of the cell. Zero unless an unused copy has been accepted.

*  Returned Value:
*     Non-zero if the cell is used, or is an unused copy of a used cell
*     and "wrap" is non-zero.

*  Notes:
*     - Zero is returned if either index is not finite.
*/

/* Local Variables: */
   double fx;                    /* Facet offset along X axis */
   double fy;                    /* Facet offset along Y axis */
   double s;                     /* Sum of x and y cell indices, plus one */

/* Initialise */
   *shift = 0.0;

/* Get the zero-based facet (x,y) indices - offsets from the bottom left
   facet - and check that they are legal. */
   fx = floor( ix/nppf );
   fy = floor( iy/nppf );
   if( fx >= 0 && fx < 5 && fy >= 0 && fy < 5 && fabs( fx - fy ) <= 1 ) {

/* Get the sum of the zero-based x and y cell indices, plus one. */
      s = ix + iy + 1;

/* The bottom left half of the bottom left facet and the top right half
   of the top right facet hold the unused copies. The used copy of each
   is four facets away on both axes. */
      if( s > nppf && s <= 9*nppf ) {
         return 1;

      } else if( wrap ) {
         *shift = ( s <= nppf ) ? 4*nppf : -4*nppf;
         return 1;
      }
   }

   return 0;
}

static int MapMerge( AstMapping *this, int where, int series, int *nmap,
                     AstMapping ***map_list, int **invert_list, int *status ) {
/*
*  Name:
*     MapMerge

*  Purpose:
*     Simplify a sequence of Mappings containing a XphMap.

*  Type:
*     Private function.

*  Synopsis:
*     #include "mapping.h"
*     int MapMerge( AstMapping *this, int where, int series, int *nmap,
*                   AstMapping ***map_list, int **invert_list, int *status )

*  Class Membership:
*     XphMap method (over-rides the protected astMapMerge method
*     inherited from the Mapping class).

*  Description:
*     This function attempts to simplify a sequence of Mappings by
*     merging a nominated XphMap in the sequence with its neighbours,
*     so as to shorten the sequence if possible.
*
*     In many cases, simplification will not be possible and the
*     function will return -1 to indicate this, without further
*     action.
*
*     In most cases of interest, however, this function will either
*     attempt to replace the nominated XphMap with one which it
*     considers simpler, or to merge it with the Mappings which
*     immediately precede it or follow it in the sequence (both will
*     normally be considered). This is sufficient to ensure the
*     eventual simplification of most Mapping sequences by repeated
*     application of this function.
*
*     In some cases, the function may attempt more elaborate
*     simplification, involving any number of other Mappings in the
*     sequence. It is not restricted in the type or scope of
*     simplification it may perform, but will normally only attempt
*     elaborate simplification in cases where a more straightforward
*     approach is not adequate.

*  Parameters:
*     this
*        Pointer to the nominated XphMap which is to be merged with
*        its neighbours. This should be a cloned copy of the XphMap
*        pointer contained in the array element "(*map_list)[where]"
*        (see below). This pointer will not be annulled, and the
*        XphMap it identifies will not be modified by this function.
*     where
*        Index in the "*map_list" array (below) at which the pointer
*        to the nominated XphMap resides.
*     series
*        A non-zero value indicates that the sequence of Mappings to
*        be simplified will be applied in series (i.e. one after the
*        other), whereas a zero value indicates that they will be
*        applied in parallel (i.e. on successive sub-sets of the
*        input/output coordinates).
*     nmap
*        Address of an int which counts the number of Mappings in the
*        sequence. On entry this should be set to the initial number
*        of Mappings. On exit it will be updated to record the number
*        of Mappings remaining after simplification.
*     map_list
*        Address of a pointer to a dynamically allocated array of
*        Mapping pointers (produced, for example, by the astMapList
*        method) which identifies the sequence of Mappings. On entry,
*        the initial sequence of Mappings to be simplified should be
*        supplied.
*
*        On exit, the contents of this array will be modified to
*        reflect any simplification carried out. Any form of
*        simplification may be performed. This may involve any of: (a)
*        removing Mappings by annulling any of the pointers supplied,
*        (b) replacing them with pointers to new Mappings, (c)
*        inserting additional Mappings and (d) changing their order.
*
*        The intention is to reduce the number of Mappings in the
*        sequence, if possible, and any reduction will be reflected in
*        the value of "*nmap" returned. However, simplifications which
*        do not reduce the length of the sequence (but improve its
*        execution time, for example) may also be performed, and the
*        sequence might conceivably increase in length (but normally
*        only in order to split up a Mapping into pieces that can be
*        more easily merged with their neighbours on subsequent
*        invocations of this function).
*
*        If Mappings are removed from the sequence, any gaps that
*        remain will be closed up, by moving subsequent Mapping
*        pointers along in the array, so that vacated elements occur
*        at the end. If the sequence increases in length, the array
*        will be extended (and its pointer updated) if necessary to
*        accommodate any new elements.
*
*        Note that any (or all) of the Mapping pointers supplied in
*        this array may be annulled by this function, but the Mappings
*        to which they refer are not modified in any way (although
*        they may, of course, be deleted if the annulled pointer is
*        the final one).
*     invert_list
*        Address of a pointer to a dynamically allocated array which,
*        on entry, should contain values to be assigned to the Invert
*        attributes of the Mappings identified in the "*map_list"
*        array before they are applied (this array might have been
*        produced, for example, by the astMapList method). These
*        values will be used by this function instead of the actual
*        Invert attributes of the Mappings supplied, which are
*        ignored.
*
*        On exit, the contents of this array will be updated to
*        correspond with the possibly modified contents of the
*        "*map_list" array.  If the Mapping sequence increases in
*        length, the "*invert_list" array will be extended (and its
*        pointer updated) if necessary to accommodate any new
*        elements.
*     status
*        Pointer to the inherited status variable.

*  Returned Value:
*     If simplification was possible, the function returns the index
*     in the "map_list" array of the first element which was
*     modified. Otherwise, it returns -1 (and makes no changes to the
*     arrays supplied).

*  Notes:
*     - A value of -1 will be returned if this function is invoked
*     with the global error status set, or if it should fail for any
*     reason.
*/

/* Local Variables: */
   AstMapping *new;              /* Pointer to replacement Mapping */
   const char *class;            /* Pointer to Mapping class string */
   int cancel;                   /* Do mappings cancel? */
   int imap1;                    /* Index of first XphMap */
   int imap2;                    /* Index of last XphMap */
   int imap;                     /* Loop counter for Mappings */
   int old_invert1;              /* Original Invert flag for map1 */
   int old_invert2;              /* Original Invert flag for map2 */
   int result;                   /* Result value to return */

/* Initialise the returned result. */
   result = -1;

/* Check the global error status. */
   if ( !astOK ) return result;

/* Note, an XphMap with a type of AST__HPX12 is not equivalent to a UnitMap
   because it performs clipping of out-of-bounds input positions, even
   though input positions that are within bounds are copied unchanged. */

/* In series. */
/* ---------- */
/* Handle the case where the Mappings are connected in series. */
   if ( series ) {

/* Look for an adjacent XphMap. */
      imap1 = -1;
      imap2 = -1;

      if( where > 0 ) {
         class = astGetClass( ( *map_list )[ where - 1 ] );
         if ( astOK && !strcmp( class, "XphMap" ) ) {
            imap1 = where - 1;
            imap2 = where;
         }
      }

      if( imap1 == -1 && where < *nmap - 1 ) {
         class = astGetClass( ( *map_list )[ where + 1 ] );
         if ( astOK && !strcmp( class, "XphMap" ) ) {
            imap1 = where;
            imap2 = where + 1;
         }
      }

/* If one was found, check to see if the two adjacent XphMaps cancel out. */
      if( imap1 != -1 ) {

/* Ensure they have the required Invert flags, saving the original
   values first so they can be re-instated later. */
         old_invert1 = astGetInvert( ( *map_list )[ imap1 ] );
         old_invert2 = astGetInvert( ( *map_list )[ imap2 ] );
         astSetInvert( ( *map_list )[ imap1 ], ( *invert_list )[ imap1 ] );
         astSetInvert( ( *map_list )[ imap2 ], ( *invert_list )[ imap2 ] );

/* They cancel out if they are equal and opposite. So invert one and then
   compare them for equality. */
         astInvert( ( *map_list )[ imap1 ] );
         cancel = astEqual( ( *map_list )[ imap1 ], ( *map_list )[ imap2 ] );
         astInvert( ( *map_list )[ imap1 ] );

/* Reinstate the original Invert flags. */
         astSetInvert( ( *map_list )[ imap1 ], old_invert1 );
         astSetInvert( ( *map_list )[ imap2 ], old_invert2 );

/* Replace the two XphMaps with a UnitMap if they cancel. */
         if( cancel ) {

/* Annul the old Mapping pointers */
            ( *map_list )[ imap1 ] = astAnnul( ( *map_list )[ imap1 ] );
            ( *map_list )[ imap2 ] = astAnnul( ( *map_list )[ imap2 ] );

/* Create the UnitMap. */
            new = (AstMapping *) astUnitMap( 2, "", status );

/* Insert the pointer to the replacement Mapping and initialise its
   invert flag. */
            ( *map_list )[ imap1 ] = new;
            ( *invert_list )[ imap1 ] = 0;

/* Loop to close the resulting gap by moving subsequent elements down
   in the arrays. */
            for ( imap = imap2 + 1; imap < *nmap; imap++ ) {
               ( *map_list )[ imap - 1 ] = ( *map_list )[ imap ];
               ( *invert_list )[ imap - 1 ] = ( *invert_list )[ imap ];
            }

/* Clear the vacated elements at the end. */
            ( *map_list )[ *nmap - 1 ] = NULL;
            ( *invert_list )[ *nmap - 1 ] = 0;

/* Decrement the Mapping count and return the index of the first
   modified element. */
            ( *nmap )--;
            result = imap1;
         }
      }
   }

/* If an error occurred, clear the returned result. */
   if ( !astOK ) result = -1;

/* Return the result. */
   return result;
}

static AstPointSet *Transform( AstMapping *this, AstPointSet *in,
                               int forward, AstPointSet *out, int *status ) {
/*
*  Name:
*     Transform

*  Purpose:
*     Apply a XphMap to transform a set of points.

*  Type:
*     Private function.

*  Synopsis:
*     #include "xphmap.h"
*     AstPointSet *Transform( AstMapping *this, AstPointSet *in,
*                             int forward, AstPointSet *out, int *status )

*  Class Membership:
*     XphMap member function (over-rides the astTransform protected
*     method inherited from the Mapping class).

*  Description:
*     This function takes a XphMap and a set of points encapsulated in a
*     PointSet and transforms the points from grid coordinates in one type
*     of HEALPix projection to grid coordinates in another. A point that is
*     not on the sky in the projection it is supplied in is given bad
*     values.

*  Parameters:
*     this
*        Pointer to the XphMap.
*     in
*        Pointer to the PointSet holding the input coordinate data.
*     forward
*        A non-zero value indicates that the forward coordinate transformation
*        should be applied, while a zero value requests the inverse
*        transformation.
*     out
*        Pointer to a PointSet which will hold the transformed (output)
*        coordinate values. A NULL value may also be given, in which case a
*        new PointSet will be created by this function.
*     status
*        Pointer to the inherited status variable.

*  Returned Value:
*     Pointer to the output (possibly new) PointSet.

*  Notes:
*     -  A null pointer will be returned if this function is invoked with the
*     global error status set, or if it should fail for any reason.
*     -  The number of coordinate values per point in the input PointSet must
*     match the number of coordinates for the XphMap being applied.
*     -  If an output PointSet is supplied, it must have space for sufficient
*     number of points and coordinate values per point to accommodate the
*     result. Any excess space will be ignored.
*/

/* Local Variables: */
   AstPointSet *result;          /* Pointer to output PointSet */
   AstXphMap *map;               /* Pointer to XphMap to be applied */
   double **ptr_in;              /* Pointer to input coordinate data */
   double **ptr_out;             /* Pointer to output coordinate data */
   double *px_in;                /* Pointer to input X values */
   double *px_out;               /* Pointer to output X values */
   double *py_in;                /* Pointer to input Y values */
   double *py_out;               /* Pointer to output Y values */
   double s;                     /* Sum of x and y cell indices, plus one */
   double tol;                   /* Width of the edge of a cell, in pixels */
   double x;                     /* X value of the current position */
   double xin;                   /* X value before transformation */
   double y;                     /* Y value of the current position */
   double yin;                   /* Y value before transformation */
   int gore;                     /* Index of gore holding the position */
   int n2;                       /* Two times nppf */
   int n3;                       /* Three times nppf */
   int n4;                       /* Four times nppf */
   int n5;                       /* Five times nppf */
   int n6;                       /* Six times nppf */
   int n7;                       /* Seven times nppf */
   int npoint;                   /* Number of points */
   int nppf;                     /* No. of pixels per facet */
   int ok;                       /* Is the current position on the sky? */
   int point;                    /* Loop counter for points */

/* Check the global error status. */
   if ( !astOK ) return NULL;

/* Obtain a pointer to the XphMap. */
   map = (AstXphMap *) this;

/* Apply the parent mapping using the stored pointer to the Transform member
   function inherited from the parent Mapping class. This function validates
   all arguments and generates an output PointSet if necessary, but does not
   actually transform any coordinate values. */
   result = (*parent_transform)( this, in, forward, out, status );

/* We will now extend the parent astTransform method by performing the
   calculations needed to generate the output coordinate values. */

/* Determine the numbers of points and coordinates per point from the input
   PointSet and obtain pointers for accessing the input and output coordinate
   values. */
   npoint = astGetNpoint( in );
   ptr_in = astGetPoints( in );
   ptr_out = astGetPoints( result );

/* Determine whether to apply the forward or inverse mapping, according to the
   direction specified and whether the mapping has been inverted. */
   if ( astGetInvert( map ) ) forward = !forward;

/* Get the number of pixels per facet. */
   nppf = ( 1 << map->order );
   n2 = 2*nppf;
   n3 = 3*nppf;
   n4 = 4*nppf;
   n5 = 5*nppf;
   n6 = 6*nppf;
   n7 = 7*nppf;

/* Rounding errors in the calculation of its grid coordinates can move a
   position that is on the edge of the sky to a position that is just off
   the sky. So a position that is within a small distance of the edge of
   a grid cell is treated as being on the edge. The distance is the margin
   that the HPX projection allows for rounding, which is AST__HPXTOL
   degrees. A step of one pixel along either axis of a grid changes each
   coordinate in the plane of the projection by 45/nppf degrees. */
   tol = AST__HPXTOL*nppf/45.0;

/* Perform coordinate arithmetic. */
/* ------------------------------ */
/* Each grid is divided into square cells, one pixel across, and only
   some of the cells are used: each cell of the sky is one cell of each
   grid. A position is on the sky if it is within a used cell or on its
   boundary, and the transformation of a position is determined by the
   used cell that holds it. A position that is in no used cell, or on the
   boundary of none, is not on the sky and is transformed to bad values.

   A position on the boundary between two cells belongs to the cell on
   its high RA side. If that cell is not used, the position belongs to a
   used cell on another side of it. There may be more than one of these,
   and they give different positions for the same point on the sky. In
   particular, a pole is a single position in the XPH grid centred on it
   but is four positions in each of the other grids, one at the corner of
   each polar facet, and all four are transformed to the single position.

   The indices of cells and facets are whole numbers but are held in
   doubles. An input position can be any distance from the grid, and an
   integer could not hold the index of every cell. */
   if ( astOK && map->type != AST__HPX12 && map->type != AST__HPX0 &&
        map->type != AST__XPHN && map->type != AST__XPHS ) {
      astError( AST__INTER, "astTransform(%s): Invalid projection type "
                "encountered (%d) (internal programming error).", status,
                astGetClass(this), map->type );
   }

   if ( astOK ) {
      px_in = ptr_in[ 0 ];
      py_in = ptr_in[ 1 ];
      px_out = ptr_out[ 0 ];
      py_out = ptr_out[ 1 ];

      for ( point = 0; point < npoint; point++ ) {
         x = px_in[ point ];
         y = py_in[ point ];
         ok = ( x != AST__BAD && y != AST__BAD );

/* Forward transformation - from "type" to HPX12. */
         if( ok && forward ) {

/* Transforming from HPX12 to HPX12 is a unit transformation, except that
   positions that are not on the sky are set bad. The cells in the bottom
   left half of the bottom left facet and the top right half of the top
   right facet of an HPX grid are not used, as each is a second copy of a
   cell in the other of the two facets. */
            if( map->type == AST__HPX12 || map->type == AST__HPX0 ){
               ok = FindHpxCell( &x, &y, nppf, tol, 0, &s );

/* Transforming from HPX0 to HPX12 - swap the top right and bottom left
   halves of the whole sky map, by adding or subtracting 2*nppf on each
   axis. */
               if( ok && map->type == AST__HPX0 ){
                  if( s <= n5 ) {
                     x += n2;
                     y += n2;
                  } else {
                     x -= n2;
                     y -= n2;
                  }
               }

/* Transforming from XPHN or XPHS to HPX12. The transformation depends on
   the gore containing the used cell that holds the position. */
            } else {
               gore = FindXphCell( map->type, &x, &y, nppf, tol );
               ok = ( gore >= 0 );
               xin = x;
               yin = y;

               if( ok && map->type == AST__XPHN ){
                  if( gore == 0 ) {
                     x = n4 + 1.0 - xin;
                     y = n6 + 1.0 - yin;
                  } else if( gore == 3 ) {
                     x = yin - nppf;
                     y = n5 + 1.0 - xin;
                  } else if( gore == 2 ) {
                     x = xin - n2;
                     y = yin;
                  } else {
                     x = n5 + 1.0 - yin;
                     y = n3 + xin;
                  }

               } else if( ok ) {
                  if( gore == 0 ) {
                     x = n2 + yin;
                     y = n4 + 1.0 - xin;
                  } else if( gore == 3 ) {
                     x = n5 + 1.0 - xin;
                     y = n3 + 1.0 - yin;
                  } else if( gore == 2 ) {
                     x = n4 + 1.0 - yin;
                     y = xin - n2;
                  } else {
                     x = n3 + xin;
                     y = nppf + yin;
                  }
               }
            }

/* Inverse transformation - from HPX12 to "type". The supplied position
   may be in the unused copy of a cell in the bottom left or top right
   facet, in which case it is first moved to the equivalent position in
   the used copy of the cell. The transformation then depends on the
   value of "s" - the sum of the zero-based x and y indices of the used
   cell that holds the position, plus one. */
         } else if( ok ) {
            ok = FindHpxCell( &x, &y, nppf, tol, 1, &s );
            xin = x;
            yin = y;

/* Transforming from HPX12 to HPX12 leaves the position unchanged.
   Transforming from HPX12 to HPX0 - swap the top right and bottom left
   halves of the whole sky map. */
            if( ok && map->type == AST__HPX0 ){
               if( s <= n5 ) {
                  x += n2;
                  y += n2;
               } else {
                  x -= n2;
                  y -= n2;
               }

/* Transforming from HPX12 to XPHN */
            } else if( ok && map->type == AST__XPHN ){

/* Gore 2 */
               if( s <= n3 ) {
                  x = n2 + xin;
                  y = yin;
/* Gore 3 */
               } else if( s <= n5 ) {
                  x = n5 + 1.0 - yin;
                  y = nppf + xin;
/* Gore 0 */
               } else if( s <= n7 ) {
                  x = n4 + 1.0 - xin;
                  y = n6 + 1.0 - yin;
/* Gore 1 */
               } else {
                  x = yin - n3;
                  y = n5 + 1.0 - xin;
               }

/* Transforming from HPX12 to XPHS */
            } else if( ok && map->type == AST__XPHS ){

/* Gore 2 */
               if( s <= n3 ) {
                  x = n2 + yin;
                  y = n4 + 1.0 - xin;
/* Gore 3 */
               } else if( s <= n5 ) {
                  x = n5 + 1.0 - xin;
                  y = n3 + 1.0 - yin;
/* Gore 0 */
               } else if( s <= n7 ) {
                  x = n4 + 1.0 - yin;
                  y = xin - n2;
/* Gore 1 */
               } else {
                  x = xin - n3;
                  y = yin - nppf;
               }
            }
         }

/* Store the transformed position, or bad values if the supplied position
   is bad or is not on the sky. */
         if( ok ) {
            px_out[ point ] = x;
            py_out[ point ] = y;
         } else {
            px_out[ point ] = AST__BAD;
            py_out[ point ] = AST__BAD;
         }
      }
   }

/* Return a pointer to the output PointSet. */
   return result;
}

static int XphCellGore( int type, double ix, double iy, int nppf ) {
/*
*  Name:
*     XphCellGore

*  Purpose:
*     Find the gore that holds a cell of an XPH grid.

*  Type:
*     Private function.

*  Synopsis:
*     #include "xphmap.h"
*     int XphCellGore( int type, double ix, double iy, int nppf )

*  Class Membership:
*     XphMap member function

*  Description:
*     This function returns the index of the gore that holds a cell in
*     the grid of an XPH projection, if the cell is used.
*
*     The grid holds four facets along each axis and is divided into four
*     quadrants, each of which holds one gore. The used cells of a gore
*     form a band, two facets wide, that runs along the diagonal of the
*     quadrant from the centre of the grid to the corner of the grid. The
*     cells on each edge of the band are cut in half by the boundary
*     between two gores, and are used in the gore on their low RA side.

*  Parameters:
*     type
*        The type of projection: AST__XPHN or AST__XPHS.
*     ix
*        The zero-based X index of the cell.
*     iy
*        The zero-based Y index of the cell.
*     nppf
*        The number of pixels along one edge of a facet.

*  Returned Value:
*     The index of the gore, in the range 0 to 3, or -1 if the cell is
*     not used.

*  Notes:
*     - A value of -1 is returned if either index is not finite.
*/

/* Local Variables: */
   double d;                     /* Difference between x and y cell indices */
   double s;                     /* Sum of x and y cell indices, plus one */
   int bottom;                   /* Is the cell in the bottom half of the grid? */
   int gore;                     /* The gore in the quadrant holding the cell */
   int left;                     /* Is the cell in the left half of the grid? */
   int n3;                       /* Three times nppf */
   int n5;                       /* Five times nppf */
   int used;                     /* Is the cell used? */

/* Check the cell is in the grid. */
   if( !( ix >= 0 && ix < 4*nppf && iy >= 0 && iy < 4*nppf ) ) return -1;

/* Find the quadrant that holds the cell. */
   left = ( ix < 2*nppf );
   bottom = ( iy < 2*nppf );

/* The band of used cells runs along the diagonal on which the sum of the
   indices is constant in the top left and bottom right quadrants, and
   along the diagonal on which the difference of the indices is constant
   in the other two. */
   n3 = 3*nppf;
   n5 = 5*nppf;
   s = ix + iy + 1;
   d = ix - iy;

   if( type == AST__XPHN ) {
      if( left && !bottom ) {
         gore = 0;
         used = ( s >= n3 && s < n5 );
      } else if( !left && !bottom ) {
         gore = 3;
         used = ( d >= -nppf && d < nppf );
      } else if( !left && bottom ) {
         gore = 2;
         used = ( s > n3 && s <= n5 );
      } else {
         gore = 1;
         used = ( d > -nppf && d <= nppf );
      }

   } else {
      if( left && bottom ) {
         gore = 0;
         used = ( d >= -nppf && d < nppf );
      } else if( !left && bottom ) {
         gore = 3;
         used = ( s >= n3 && s < n5 );
      } else if( !left && !bottom ) {
         gore = 2;
         used = ( d > -nppf && d <= nppf );
      } else {
         gore = 1;
         used = ( s > n3 && s <= n5 );
      }
   }

   return used ? gore : -1;
}

/* Functions which access class attributes. */
/* ---------------------------------------- */
/* Implement member functions to access the attributes associated with
   this class using the macros defined for this purpose in the
   "object.h" file. For a description of each attribute, see the class
   interface (in the associated .h file). */

/* Copy constructor. */
/* ----------------- */
/* No copy constructor is needed, as a byte-by-byte copy suffices. */

/* Destructor. */
/* ----------- */
/* No destructor is needed as no memory, etc. needs freeing. */

/* Dump function. */
/* -------------- */
static void Dump( AstObject *this_object, AstChannel *channel, int *status ) {
/*
*  Name:
*     Dump

*  Purpose:
*     Dump function for XphMap objects.

*  Type:
*     Private function.

*  Synopsis:
*     void Dump( AstObject *this, AstChannel *channel, int *status )

*  Description:
*     This function implements the Dump function which writes out data
*     for the XphMap class to an output Channel.

*  Parameters:
*     this
*        Pointer to the XphMap whose data are being written.
*     channel
*        Pointer to the Channel to which the data are being written.
*     status
*        Pointer to the inherited status variable.
*/

/* Local Variables: */
   AstXphMap *this;             /* Pointer to the XphMap structure */

/* Check the global error status. */
   if ( !astOK ) return;

/* Obtain a pointer to the XphMap structure. */
   this = (AstXphMap *) this_object;

/* Write out values representing the instance variables for the
   XphMap class.  Accompany these with appropriate comment strings,
   possibly depending on the values being written.*/

   astWriteInt( channel, "Order", 1, 1, this->order, "HEALPix order" );
   astWriteString( channel, "Type", 1, 1, proj_name[this->type],
                   proj_comm[this->type]);
}

/* Standard class functions. */
/* ========================= */
/* Implement the astIsAXphMap and astCheckXphMap functions using the macros
   defined for this purpose in the "object.h" header file. */
astMAKE_ISA(XphMap,Mapping)
astMAKE_CHECK(XphMap)

AstXphMap *astXphMap_( int order, int type, const char *options, int *status, ...) {
/*
*+
*  Name:
*     astXphMap

*  Purpose:
*     Create an XphMap.

*  Type:
*     Public function.

*  Synopsis:
*     #include "xphmap.h"
*     AstXphMap *astXphMap( int order, int type, const char *options, ... )

*  Class Membership:
*     XphMap constructor.

*  Description:
*     This function creates a new XphMap and optionally initialises its
*     attributes.
*
*     The XphMap class implements a Mapping between grid coordinates
*     within different types of HEALPix projection: The available types
*     are :
*
*     - "HPX0":  An HPX projection centred on RA=0h.
*     - "HPX12": An HPX projection centred on RA=12h.
*     - "XPHN":  An XPH ("Butterfly") projection centred on the
*                 north pole.
*     - "XPHS":  An XPH ("Butterfly") projection centred on the
*                south pole.
*
*     The HPX projection is described in "Mapping on the HEALPix grid"
*     by Calabretta and Roukema, A&A, 2005. The XPH projection is
*     described in "Representing the 'butterfly' projection in FITS -
*     projection code XPH" by Calabretta and Lowe, PASA, 2018.
*
*     The forward transformation of an XphMap transforms from grid
*     coordinates within the projection type specified by argument
*     "type" into grid coordinates within an HPX0 projection.

*  Parameters:
*     order
*        The HEALPix order of the HPX transformation.
*     type
*        The type of projection corresponding to the inputs of the
*        XphMap. This can be:
*        - AST__HPX0: A basic HPX projection (centred on RA=0h)
*        - AST__HPX12: An HPX projection centred on RA=12h.
*        - AST__XPHN: An XPH ("Butterfly") projection centred on the
*                     north pole.
*        - AST__XPHS: An XPH ("Butterfly") projection centred on the
*                     south pole.
*     options
*        Pointer to a null-terminated string containing an optional
*        comma-separated list of attribute assignments to be used for
*        initialising the new XphMap. The syntax used is identical to
*        that for the astSet function and may include "printf" format
*        specifiers identified by "%" symbols in the normal way.
*     ...
*        If the "options" string contains "%" format specifiers, then
*        an optional list of additional arguments may follow it in
*        order to supply values to be substituted for these
*        specifiers. The rules for supplying these are identical to
*        those for the astSet function (and for the C "printf"
*        function).

*  Returned Value:
*     astXphMap()
*        A pointer to the new XphMap.

*  Notes:
*     - A null Object pointer (AST__NULL) will be returned if this
*     function is invoked with the AST error status set, or if it
*     should fail for any reason.

*  Status Handling:
*     The protected interface to this function includes an extra
*     parameter at the end of the parameter list descirbed above. This
*     parameter is a pointer to the integer inherited status
*     variable: "int *status".

*-
*/

/* Local Variables: */
   astDECLARE_GLOBALS            /* Pointer to thread-specific global data */
   AstXphMap *new;              /* Pointer to new XphMap */
   va_list args;                 /* Variable argument list */

/* Get a pointer to the thread specific global data structure. */
   astGET_GLOBALS(NULL);

/* Check the global status. */
   if ( !astOK ) return NULL;

/* Initialise the XphMap, allocating memory and initialising the
   virtual function table as well if necessary. */
   new = astInitXphMap( NULL, sizeof( AstXphMap ), !class_init, &class_vtab,
                         "XphMap", order, type );

/* If successful, note that the virtual function table has been
   initialised. */
   if ( astOK ) {
      class_init = 1;

/* Obtain the variable argument list and pass it along with the options string
   to the astVSet method to initialise the new XphMap's attributes. */
      va_start( args, status );
      astVSet( new, options, NULL, args );
      va_end( args );

/* If an error occurred, clean up by deleting the new object. */
      if ( !astOK ) new = astDelete( new );
   }

/* Return a pointer to the new XphMap. */
   return new;
}

AstXphMap *astInitXphMap_( void *mem, size_t size, int init,
                           AstXphMapVtab *vtab, const char *name,
                           int order, int type, int *status ) {
/*
*+
*  Name:
*     astInitXphMap

*  Purpose:
*     Initialise a XphMap.

*  Type:
*     Protected function.

*  Synopsis:
*     #include "xphmap.h"
*     AstXphMap *astInitXphMap( void *mem, size_t size, int init,
*                               AstXphMapVtab *vtab, const char *name,
*                               int order, int type )

*  Class Membership:
*     XphMap initialiser.

*  Description:
*     This function is provided for use by class implementations to initialise
*     a new XphMap object. It allocates memory (if necessary) to accommodate
*     the XphMap plus any additional data associated with the derived class.
*     It then initialises a XphMap structure at the start of this memory. If
*     the "init" flag is set, it also initialises the contents of a virtual
*     function table for a XphMap at the start of the memory passed via the
*     "vtab" parameter.

*  Parameters:
*     mem
*        A pointer to the memory in which the XphMap is to be initialised.
*        This must be of sufficient size to accommodate the XphMap data
*        (sizeof(XphMap)) plus any data used by the derived class. If a value
*        of NULL is given, this function will allocate the memory itself using
*        the "size" parameter to determine its size.
*     size
*        The amount of memory used by the XphMap (plus derived class data).
*        This will be used to allocate memory if a value of NULL is given for
*        the "mem" parameter. This value is also stored in the XphMap
*        structure, so a valid value must be supplied even if not required for
*        allocating memory.
*     init
*        A logical flag indicating if the XphMap's virtual function table is
*        to be initialised. If this value is non-zero, the virtual function
*        table will be initialised by this function.
*     vtab
*        Pointer to the start of the virtual function table to be associated
*        with the new XphMap.
*     name
*        Pointer to a constant null-terminated character string which contains
*        the name of the class to which the new object belongs (it is this
*        pointer value that will subsequently be returned by the astGetClass
*        method).
*     order
*        The HEALPix order.
*     type
*        The projection type.

*  Returned Value:
*     A pointer to the new XphMap.

*  Notes:
*     -  A null pointer will be returned if this function is invoked with the
*     global error status set, or if it should fail for any reason.
*-
*/

/* Local Variables: */
   AstXphMap *new;              /* Pointer to new XphMap */

/* Check the global status. */
   if ( !astOK ) return NULL;

/* If necessary, initialise the virtual function table. */
   if ( init ) astInitXphMapVtab( vtab, name );

/* Initialise. */
   new = NULL;


/* Validate inputs */
   if( order < 0 || order > AST__MXORDHPX ) {
      astError( AST__INTER, "astInitXphMap(%s): Invalid order supplied "
                "(%d) (internal programming error).", status, name, order );
   } else if( type < 0 || type > AST__MXPRJHPX ) {
      astError( AST__INTER, "astInitXphMap(%s): Invalid projection type "
                "supplied (%d) (internal programming error).", status,
                 name, type );
   }

/* Initialise a Mapping structure (the parent class) as the first component
   within the XphMap structure, allocating memory if necessary. Specify that
   the Mapping should be defined in both the forward and inverse directions. */
   new = (AstXphMap *) astInitMapping( mem, size, 0, (AstMappingVtab *) vtab,
                                       name, 2, 2, 1, 1 );
   if ( astOK ) {

/* Initialise the XphMap data. */
/* ---------------------------- */
      new->order = order;
      new->type = type;

/* If an error occurred, clean up by deleting the new XphMap. */
      if ( !astOK ) new = astDelete( new );
   }

/* Return a pointer to the new XphMap. */
   return new;
}

AstXphMap *astLoadXphMap_( void *mem, size_t size, AstXphMapVtab *vtab,
                           const char *name, AstChannel *channel, int *status ) {
/*
*+
*  Name:
*     astLoadXphMap

*  Purpose:
*     Load a XphMap.

*  Type:
*     Protected function.

*  Synopsis:
*     #include "xphmap.h"
*     AstXphMap *astLoadXphMap( void *mem, size_t size, AstXphMapVtab *vtab,
*                               const char *name, AstChannel *channel )

*  Class Membership:
*     XphMap loader.

*  Description:
*     This function is provided to load a new XphMap using data read
*     from a Channel. It first loads the data used by the parent class
*     (which allocates memory if necessary) and then initialises a
*     XphMap structure in this memory, using data read from the input
*     Channel.
*
*     If the "init" flag is set, it also initialises the contents of a
*     virtual function table for a XphMap at the start of the memory
*     passed via the "vtab" parameter.


*  Parameters:
*     mem
*        A pointer to the memory into which the XphMap is to be
*        loaded.  This must be of sufficient size to accommodate the
*        XphMap data (sizeof(XphMap)) plus any data used by derived
*        classes. If a value of NULL is given, this function will
*        allocate the memory itself using the "size" parameter to
*        determine its size.
*     size
*        The amount of memory used by the XphMap (plus derived class
*        data).  This will be used to allocate memory if a value of
*        NULL is given for the "mem" parameter. This value is also
*        stored in the XphMap structure, so a valid value must be
*        supplied even if not required for allocating memory.
*
*        If the "vtab" parameter is NULL, the "size" value is ignored
*        and sizeof(AstXphMap) is used instead.
*     vtab
*        Pointer to the start of the virtual function table to be
*        associated with the new XphMap. If this is NULL, a pointer
*        to the (static) virtual function table for the XphMap class
*        is used instead.
*     name
*        Pointer to a constant null-terminated character string which
*        contains the name of the class to which the new object
*        belongs (it is this pointer value that will subsequently be
*        returned by the astGetClass method).
*
*        If the "vtab" parameter is NULL, the "name" value is ignored
*        and a pointer to the string "XphMap" is used instead.

*  Returned Value:
*     A pointer to the new XphMap.

*  Notes:
*     - A null pointer will be returned if this function is invoked
*     with the global error status set, or if it should fail for any
*     reason.
*-
*/

/* Local Variables: */
   AstXphMap *new;              /* Pointer to the new XphMap */
   char *text;                  /* Text for string-valued attribute */
   astDECLARE_GLOBALS           /* Pointer to thread-specific global data */

/* Initialise. */
   new = NULL;

/* Check the global error status. */
   if ( !astOK ) return new;

/* Get a pointer to the thread specific global data structure. */
   astGET_GLOBALS(channel);

/* If a NULL virtual function table has been supplied, then this is
   the first loader to be invoked for this XphMap. In this case the
   XphMap belongs to this class, so supply appropriate values to be
   passed to the parent class loader (and its parent, etc.). */
   if ( !vtab ) {
      size = sizeof( AstXphMap );
      vtab = &class_vtab;
      name = "XphMap";

/* If required, initialise the virtual function table for this class. */
      if ( !class_init ) {
         astInitXphMapVtab( vtab, name );
         class_init = 1;
      }
   }

/* Invoke the parent class loader to load data for all the ancestral
   classes of the current one, returning a pointer to the resulting
   partly-built XphMap. */
   new = astLoadMapping( mem, size, (AstMappingVtab *) vtab, name,
                         channel );

   if ( astOK ) {

/* Read input data. */
/* ================ */
/* Request the input Channel to read all the input data appropriate to
   this class into the internal "values list". */
      astReadClassData( channel, "XphMap" );

/* Now read each individual data item from this list and use it to
   initialise the appropriate instance variable(s) for this class. */
      new->order = astReadInt( channel, "order", 19 );

/* The transformation uses the order as a shift count and holds multiples
   of the shifted value in integers, so the order must be in the range
   accepted by the constructor. */
      if( ( new->order < 0 || new->order > AST__MXORDHPX ) && astOK ) {
         astError( AST__LDERR, "astLoad(%s): The HEALPix order read for a "
                   "%s (%d) is not in the range 0 to %d.", status, name,
                   name, new->order, AST__MXORDHPX );
      }

      text = astReadString( channel, "type", proj_name[AST__HPX12] );

/* COnvert the projection name into an integer projection identifier. */
      for( new->type = 0; new->type <= AST__MXPRJHPX; new->type++ ) {
         if( astChrMatch( text, proj_name[new->type] ) ) break;
      }

/* Report an error if the projection name is unknown. */
      if( new->type > AST__MXPRJHPX && astOK ) {
         astError( AST__OPT, "astRead(XphMap): Illegal value '%s' supplied "
                   "for the XphMap component 'Type'.", status, text );
      }

/* Free the string returned by astReadString. */
      text = astFree( text );

/* An XphMap always has two inputs and two outputs, and the
   transformation relies on this. The numbers read from the dump by the
   parent Mapping loader are not otherwise constrained, so check them. */
      astCheckLoadedCoords( new, name, 2, 2 );

/* If an error occurred, clean up by deleting the new XphMap. */
      if ( !astOK ) new = astDelete( new );
   }

/* Return the new XphMap pointer. */
   return new;
}

/* Virtual function interfaces. */
/* ============================ */
/* These provide the external interface to the virtual functions defined by
   this class. Each simply checks the global error status and then locates and
   executes the appropriate member function, using the function pointer stored
   in the object's virtual function table (this pointer is located using the
   astMEMBER macro defined in "object.h").

   Note that the member function may not be the one defined here, as it may
   have been over-ridden by a derived class. However, it should still have the
   same interface. */





