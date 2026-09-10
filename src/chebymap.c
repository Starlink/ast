/* To do:

   - what to do about input positions that fall outside the bounding box.
   Have an attribute that can be used to select "set bad" or "extrapolate"?

   - what about overriding astRate ?

*/


/*
*class++
*  Name:
*     ChebyMap

*  Purpose:
*     Map coordinates using Chebyshev polynomial functions.

*  Constructor Function:
c     astChebyMap
f     AST_CHEBYMAP

*  Description:
*     A ChebyMap is a form of Mapping which performs a Chebyshev polynomial
*     transformation.  Each output coordinate is a linear combination of
*     Chebyshev polynomials of the first kind, of order zero up to a
*     specified maximum order, evaluated at the input coordinates. The
*     coefficients to be used in the linear combination are specified
*     separately for each output coordinate.
*
*     For a 1-dimensional ChebyMap, the forward transformation is defined
*     as follows:
*
*        f(x) = c0.T0(x') + c1.T1(x') + c2.T2(x') + ...
*
*     where:
*        - Tn(x') is the nth Chebyshev polynomial of the first kind:
*             - T0(x') = 1
*             - T1(x') = x'
*             - Tn+1(x') = 2.x'.Tn(x') - Tn-1(x')
*        - x' is the input axis value, x, offset and scaled to the range
*          [-1, 1] as x ranges over a specified bounding box, given when the
*          ChebyMap is created. The input positions, x,  supplied to the
*          forward transformation must fall within the bounding box - bad
*          axis values (AST__BAD) are generated for points outside the
*          bounding box.
*
*     For an N-dimensional ChebyMap, the forward transformation is a
*     generalisation of the above form. Each output axis value is the sum
c     of "ncoeff"
f     of NCOEFF
*     terms, where each term is the product of a single coefficient
*     value and N factors of the form Tn(x'_i), where "x'_i" is the
*     normalised value of the i'th input axis value.
*
*     The forward and inverse transformations may be defined independently
*     by separate sets of coefficients supplied when the ChebyMap is
*     created. If forward coefficients are supplied, no inverse coefficients
*     are supplied, and the numbers of inputs and outputs are equal, an
*     iterative inverse is provided by default. It uses the analytic
*     Jacobian of the forward series and confines candidate solutions to
*     the forward bounding box. An unsolved position is returned as
*     AST__BAD. See IterInverse, NiterInverse and TolInverse for details.
*
*     Supplied inverse coefficients are used by default. Setting
*     IterInverse to one selects iteration instead; setting it to zero
*     disables iteration. Clearing it restores the default selection.
*     A local iterative inverse does not guarantee a unique solution or
*     convergence at every position.
*
*     Alternatively, the
c     astPolyTran
f     AST_POLYTRAN
*     method can fit an inverse Chebyshev series, choosing coefficients
*     to minimise the residuals of a forward/inverse round trip.

*  Inheritance:
*     The ChebyMap class inherits from the PolyMap class.

*  Attributes:
*     The ChebyMap class does not define any new attributes beyond those
*     which are applicable to all PolyMaps.

*  Functions:
c     In addition to those functions applicable to all PolyMap, the
c     following functions may also be applied to all ChebyMaps:
f     In addition to those routines applicable to all PolyMap, the
f     following routines may also be applied to all ChebyMaps:
*
c     - astChebyDomain: Get the bounds of the domain of the ChebyMap
f     - AST_CHEBYDOMAIN: Get the bounds of the domain of the ChebyMap

*  Copyright:
*     Copyright (C) 2017 East Asian Observatory.
*     All Rights Reserved.

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
*     DSB: D.S. Berry (EAO)

*  History:
*     1-MAR-2017 (DSB):
*        Original version.
*     30-MAR-2017 (DSB):
*        Over-ride the astFitPoly1DInit and astFitPoly2DInit virtual
*        functions inherited form the PolyMap class.
*     5-MAY-2018 (DSB):
*        Correct usage of "forward" argument in astFitPoly1DInit and
*        astFitPoly2DInit.
*     8-SEP-2026 (TIMJ):
*        Over-ride the astGetJacobian method inherited from the PolyMap
*        class to express the derivatives of a Chebyshev series in the
*        basis of the first kind, and the astLinearGuess method to seed
*        the inversion from the normalised forward domain.
*     8-SEP-2026 (TIMJ):
*        Over-ride the astIterInverse method inherited from the PolyMap
*        class with a Newton iteration restricted to the bounding box the
*        forward series is defined over. The bounded algorithm lives in
*        this file; the parent algorithm is used for a forward series
*        without a usable Chebyshev box.
*     9-SEP-2026 (TIMJ):
*        Construct the Jacobian ChebyMaps on the canonical [-1,1] box
*        instead of a physical box reconstructed from the scales and
*        offsets, which are copied into the new Maps in any case.
*     9-SEP-2026 (TIMJ):
*        Report no finite iteration domain if the forward normalisation is
*        zero or non-finite, rather than dividing by it.
*     9-SEP-2026 (TIMJ):
*        Move a recovered iteration bound into the evaluable domain for as
*        long as doing so improves it, and report no finite domain if it
*        cannot be moved far enough, instead of returning a bound that the
*        forward evaluator rejects.
*     9-SEP-2026 (TIMJ):
*        Report an error if a bounding box is omitted for a direction that
*        has coefficients, rather than dereferencing the null pointer.
*     10-SEP-2026 (TIMJ):
*        Remove the astGetIterInverse and astSetIterInverse overrides and
*        the loader comment about an unvalidated recorded IterInverse
*        value. PolyMap now applies one IterInverse validity rule that
*        already covers a ChebyMap with no forward transformation.
*     10-SEP-2026 (TIMJ):
*        Stop checking TolInverse and NiterInverse in IterInverse; the
*        PolyMap setters now reject an unusable value before it can reach
*        here. Only the bounding box is still checked.
*     10-SEP-2026 (TIMJ):
*        Move a position off a singular Jacobian once, towards the side of
*        the box with more room, instead of declaring it unsolved
*        immediately. A forward series with no linear term is singular at
*        the midpoint seed, which the previous behaviour reported as
*        unsolved for every target.
*     10-SEP-2026 (TIMJ):
*        Override GetNiterInverse so that an unset NiterInverse defaults to
*        10 rather than PolyMap's 4. The bounded algorithm checks only the
*        final candidate, so it needs more headroom than the unbounded
*        algorithm, which returns the last iterate.
*     10-SEP-2026 (TIMJ):
*        Rename IterBounds to AxisBounds and use it, rather than a bare
*        scale-and-offset division, in ChebyDomain and in PolyTran's
*        reconstruction of a missing user bounding box, so that a bound
*        returned by astChebyDomain or used by astPolyTran always
*        evaluates without BAD. Remove the unreachable "lo <= hi" test in
*        AxisBounds.
*     10-SEP-2026 (TIMJ):
*        Create the two trial PointSets used by IterSteps once, before the
*        backtrack loop, and shrink them with astSetNpoint instead of
*        re-creating them at every level. Give IterSteps a "work" array and
*        have it copy an accepted trial's forward values into it, so that
*        IterInverse's whole-batch forward transform at the top of each
*        iteration can be skipped except on the first iteration or after a
*        position has been nudged.
*class--
*/

/* Module Macros. */
/* ============== */
/* Set the name of the class we are implementing. This indicates to
   the header files that define class interfaces that they should make
   "protected" symbols available. */
#define astCLASS ChebyMap

/* Include files. */
/* ============== */
/* Interface definitions. */
/* ---------------------- */

#include "globals.h"             /* Thread-safe global data access */
#include "error.h"               /* Error reporting facilities */
#include "memory.h"              /* Memory allocation facilities */
#include "object.h"              /* Base Object class */
#include "pointset.h"            /* Sets of points/coordinates */
#include "polymap.h"             /* Polynomial mappings (parent class) */
#include "cmpmap.h"              /* Compound mappings */
#include "chebymap.h"            /* Interface definition for this class */
#include "unitmap.h"             /* Unit mappings */
#include "matrixmap.h"           /* Affine initial guesses */
#include "shiftmap.h"
#include "permmap.h"
#include "pal.h"                 /* Linear equation solver for Newton steps */
#include "pal.h"                 /* Linear equation solver for Newton steps */

/* Error code definitions. */
/* ----------------------- */
#include "ast_err.h"             /* AST error codes */

/* C header files. */
/* --------------- */
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <float.h>

/* Module Variables. */
/* ================= */

/* Address of this static variable is used as a unique identifier for
   member of this class. */
static int class_check;

/* Pointers to parent class methods which are extended by this class. */
static size_t (* parent_getobjsize)( AstObject *, int * );
static int (* parent_equal)( AstObject *, AstObject *, int * );
static void (* parent_polypowers)( AstPolyMap *, double **, int, const int *, double **, int, int, int * );
static AstPolyMap *(*parent_polytran)( AstPolyMap *, int, double, double, int, const double *, const double *, int * );
static AstPolyMap **(*parent_getjacobian)( AstPolyMap *, int * );
static AstMapping *(*parent_linearguess)( AstPolyMap *, int * );
static void (*parent_iterinverse)( AstPolyMap *, AstPointSet *, AstPointSet *, int * );
static int (*parent_getniterinverse)( AstPolyMap *, int * );

/* A derivative term retains the original orders except on one axis. */
typedef struct ChebyDerivTerm {
   const int *powers;
   int nin;
   int axis;
   int degree;
   int output;
   double coeff;
} ChebyDerivTerm;


#ifdef THREAD_SAFE
/* Define how to initialise thread-specific globals. */
#define GLOBAL_inits \
   globals->Class_Init = 0; \

/* Create the function that initialises global data for this module. */
astMAKE_INITGLOBALS(ChebyMap)

/* Define macros for accessing each item of thread specific global data. */
#define class_init astGLOBAL(ChebyMap,Class_Init)
#define class_vtab astGLOBAL(ChebyMap,Class_Vtab)

#include <pthread.h>


#else

/* Define the class virtual function table and its initialisation flag
   as static variables. */
static AstChebyMapVtab class_vtab;   /* Virtual function table */
static int class_init = 0;       /* Virtual function table initialised? */

#endif


/* External Interface Function Prototypes. */
/* ======================================= */
/* The following functions have public prototypes only (i.e. no
   protected prototypes), so we must provide local prototypes for use
   within this module. */
AstChebyMap *astChebyMapId_( int, int, int, const double[], int, const double[],
                             const double[], const double[], const double[],
                             const double[], const char *, ... );


/* Prototypes for Private Member Functions. */
/* ======================================== */
static AstPolyMap *PolyTran( AstPolyMap *, int, double, double, int, const double *, const double *, int * );
static int Equal( AstObject *, AstObject *, int * );
static AstPolyMap **GetJacobian( AstPolyMap *, int * );
static int CompareDerivTerms( const void *, const void * );
static AstMapping *LinearGuess( AstPolyMap *, int * );
static int GetIterDomain( AstChebyMap *, double *, double *, int * );
static int AxisBounds( double, double, double *, double * );
static void IterInverse( AstPolyMap *, AstPointSet *, AstPointSet *, int * );
static void IterSteps( AstPolyMap *, int, int, double **, double **,
                       double **, const double *, const double *,
                       const double *, const double *, const double *,
                       const double *, int *, int *, int *, int * );
static void MarkUnsolved( double **, int, int, int *, int * );
static int Usable( double );
static double NudgeIntoDomain( double, double, double, double );
static int GetNiterInverse( AstPolyMap *, int * );
static size_t GetObjSize( AstObject *, int * );
static void ChebyDomain( AstChebyMap *, int, double *, double *, int * );
static void Copy( const AstObject *, AstObject *, int * );
static void Delete( AstObject *obj, int * );
static void Dump( AstObject *, AstChannel *, int * );
static void PolyPowers( AstPolyMap *, double **, int, const int *, double **, int, int, int *);
static void FitPoly1DInit( AstPolyMap *, int, double **, AstMinPackData *, double *, int *);
static void FitPoly2DInit( AstPolyMap *, int, double **, AstMinPackData *, double *, int *);

/* Member functions. */
/* ================= */

static void ChebyDomain( AstChebyMap *this, int forward, double *lbnd,
                         double *ubnd, int *status ){
/*
*++
*  Name:
c     astChebyDomain
f     AST_CHEBYDOMAIN

*  Purpose:
*     Returns the bounding box of the domain of a ChebyMap.

*  Type:
*     Public virtual function.

*  Synopsis:
c     #include "chebymap.h"
c     void astChebyDomain( AstChebyMap *this, int forward, double *lbnd,
c                          double *ubnd )
f     CALL AST_CHEBYDOMAIN( THIS, FORWARD, LBND, UBND, STATUS )

*  Class Membership:
*     ChebyMap method.

*  Description:
c     This function
f     This routine
*     returns the upper and lower limits of the box defining the domain
*     of either the forward or inverse transformation of a ChebyMap. These
*     are the values that were supplied when the ChebyMap was created.

*  Parameters:
c     this
f     THIS = INTEGER (Given)
*        Pointer to the ChebyMap.
c     forward
f     FORWARD = LOGICAL (Given)
c        A non-zero
f        A .TRUE.
*        value indicates that the domain of the ChebyMap's
*        forward transformation is to be returned, while a zero
*        value indicates that the domain of the inverse transformation
*        should be returned.
c     lbnd
f     LBND() = DOUBLE PRECISION (Returned)
c        Pointer to an
f        An
*        array in which to return the lower axis bounds of the ChebyMap
*        domain. The number of elements should be at least equal to the
*        number of ChebyMap inputs (if
c        "forward" is non-zero), or outputs (if "forward" is zero).
f        FORWARD is .TRUE.), or outputs (if FORWARD is .FALSE.).
c     ubnd
f     UBND() = DOUBLE PRECISION (Returned)
c        Pointer to an
f        An
*        array in which to return the upper axis bounds of the ChebyMap
*        domain. The number of elements should be at least equal to the
*        number of ChebyMap inputs (if
c        "forward" is non-zero), or outputs (if "forward" is zero).
f        FORWARD is .TRUE.), or outputs (if FORWARD is .FALSE.).
f     STATUS = INTEGER (Given and Returned)
f        The global status.

*  Notes:
*    - If the requested transformation is undefined (i.e. no
*    transformation coefficients were specified when the ChebyMap was
*    created), this method returns a box determined using the
c    astMapBox
f    AST_MAPBOX
*    method on the opposite transformation, if the opposite
*    transformation is defined.
*    - If the above procedure fails to determine a bounding box, the supplied
*    arrays are filled with AST__BAD values but no error is reported.

*--
*/

/* Local Variables: */
   double *lbnd_o;
   double *offset_o;
   double *offset;
   double *scale_o;
   double *scale;
   double *ubnd_o;
   int fwd_o;
   int iax;
   int nax;
   int nax_o;

/* Check the inherited status. */
   if( !astOK ) return;

/* Get the scales and offsets to use, depending on the value of "forward"
   and whether the ChebyMap has been inverted. */
   if( forward != astGetInvert( this ) ) {
      scale = this->scale_f;
      offset = this->offset_f;
      nax = astGetNin( this );
      scale_o = this->scale_i;
      offset_o = this->offset_i;
      nax_o = astGetNout( this );
      fwd_o = 0;
   } else {
      scale = this->scale_i;
      offset = this->offset_i;
      nax = astGetNout( this );
      scale_o = this->scale_f;
      offset_o = this->offset_f;
      nax_o = astGetNin( this );
      fwd_o = 1;
   }

/* Check the domain is defined. */
   if( scale && offset ) {
      for( iax = 0; iax < nax; iax++ ) {
         if( !AxisBounds( scale[ iax ], offset[ iax ], lbnd + iax, ubnd + iax ) ) {
            lbnd[ iax ] = AST__BAD;
            ubnd[ iax ] = AST__BAD;
         }
      }

/* If the requested domain is not defined, see if it can be determined
   by transforming the domain of the other transformation into the
   requested input ot putput space. */
   } else if( scale_o && offset_o ){

/* Allocate memory to hold the bounding box in the other space (input or
   output), and then store the bounding box values. */
      lbnd_o = astMalloc( nax_o*sizeof( *lbnd_o ) );
      ubnd_o = astMalloc( nax_o*sizeof( *ubnd_o ) );
      if( astOK ) {
         for( iax = 0; iax < nax_o; iax++ ) {
            if( !AxisBounds( scale_o[ iax ], offset_o[ iax ], lbnd_o + iax,
                             ubnd_o + iax ) ) {
               lbnd_o[ iax ] = AST__BAD;
               ubnd_o[ iax ] = AST__BAD;
            }
         }

/* Loop round finding the bounds on each input axis of the requested
   transformation. */
         for( iax = 0; iax < nax; iax++ ) {
            astMapBox( this, lbnd_o, ubnd_o, fwd_o, iax, lbnd + iax,
                       ubnd + iax, NULL, NULL );
         }

/* Free resources */
         lbnd_o = astFree( lbnd_o );
         ubnd_o = astFree( ubnd_o );
      }


/* If the domain of the other transformation is not defined, return bad values. */
   } else {
      for( iax = 0; iax < nax; iax++ ) {
         lbnd[ iax ] = AST__BAD;
         ubnd[ iax ] = AST__BAD;
      }
   }
}

static int Equal( AstObject *this_object, AstObject *that_object, int *status ) {
/*
*  Name:
*     Equal

*  Purpose:
*     Test if two ChebyMaps are equivalent.

*  Type:
*     Private function.

*  Synopsis:
*     #include "chebymap.h"
*     int Equal( AstObject *this, AstObject *that, int *status )

*  Class Membership:
*     ChebyMap member function (over-rides the astEqual protected
*     method inherited from the astPolyMap class).

*  Description:
*     This function returns a boolean result (0 or 1) to indicate whether
*     two ChebyMaps are equivalent.

*  Parameters:
*     this
*        Pointer to the first Object (a ChebyMap).
*     that
*        Pointer to the second Object.
*     status
*        Pointer to the inherited status variable.

*  Returned Value:
*     One if the ChebyMaps are equivalent, zero otherwise.

*  Notes:
*     - A value of zero will be returned if this function is invoked
*     with the global status set, or if it should fail for any reason.
*/

/* Local Variables: */
   AstChebyMap *that;
   AstChebyMap *this;
   int i;
   int nin;
   int nout;
   int result;

/* Initialise. */
   result = 0;

/* Check the global error status. */
   if ( !astOK ) return result;

/* Invoke the Equal method inherited from the parent PolyMap class. This
   checks that the PolyMaps are equal. */
   result = (*parent_equal)( this_object, that_object, status );
   if( result ) {

/* Obtain pointers to the two ChebyMap structures. */
      this = (AstChebyMap *) this_object;
      that = (AstChebyMap *) that_object;

/* Check the second object is a ChebyMap. We know the first is a
   ChebyMap since we have arrived at this implementation of the virtual
   function. */
      if( astIsAChebyMap( that ) ) {

/* Get the number of axes in the input bounding box (the original input space). */
         nin = astGetInvert( this ) ? astGetNout( this ) : astGetNin( this );

/* Check the bounding box is the same for both ChebyMaps. */
         if( this->scale_f && that->scale_f &&
             this->offset_f && that->offset_f ) {
            for( i = 0; i < nin && result; i++ ) {
               if( !astEQUAL( this->scale_f[ i ], that->scale_f[ i ] ) ||
                   !astEQUAL( this->offset_f[ i ], that->offset_f[ i ] )){
                  result = 0;
               }
            }
         } else if( this->scale_f || that->scale_f ||
                    this->offset_f || that->offset_f ) {
            result = 0;
         }

/* Get the number of axes in the output bounding box (the original output space). */
         nout = astGetInvert( this ) ? astGetNin( this ) : astGetNout( this );

/* Check the bounding box is the same for both ChebyMaps. */
         if( this->scale_i && that->scale_i &&
             this->offset_i && that->offset_i ) {
            for( i = 0; i < nout && result; i++ ) {
               if( !astEQUAL( this->scale_i[ i ], that->scale_i[ i ] ) ||
                   !astEQUAL( this->offset_i[ i ], that->offset_i[ i ] )){
                  result = 0;
               }
            }
         } else if( this->scale_i || that->scale_i ||
                    this->offset_i || that->offset_i ) {
            result = 0;
         }
      }
   }

/* If an error occurred, clear the result value. */
   if ( !astOK ) result = 0;

/* Return the result, */
   return result;
}

static void FitPoly1DInit( AstPolyMap *this_polymap, int forward, double **table,
                           AstMinPackData *data, double *scales, int *status ){
/*
*  Name:
*     FitPoly1DInit

*  Purpose:
*     Perform initialisation needed for FitPoly1D

*  Type:
*     Private function.

*  Synopsis:
*     #include "chebymap.h"
*     void FitPoly1DInit( AstPolyMap *this, int forward, double **table,
*                         AstMinPackData *data, double *scales, int *status )

*  Class Membership:
*     ChebyMap member function (over-rides the astFitPoly1DInit protected
*     method inherited from the PolyMap class).

*  Description:
*     This function performs initialisation needed for FitPoly1D in the
*     PolyMap class.

*  Parameters:
*     this
*        Pointer to the PolyMap.
*     forward
*        Non-zero if the forward transformation of "this" is being
*        replaced. Zero if the inverse transformation of "this" is being
*        replaced.
*     table
*        Pointer to an array of 2 pointers. Each of these pointers points
*        to an array of "nsamp" doubles, being the scaled and sampled values
*        for x1 and y1 in that order.
*     data
*        Pointer to a structure holding information to pass the the
*        service function invoked by the minimisation function.
*     scales
*        Array holding the scaling factors for the two columns of the table.
*        Multiplying the table values by the scale factor produces PolyMap
*        input or output axis values. The scales are modified on exit to
*        take account of the scaling performed by the ChebyMap Transform
*        method.
*/

/* Local Variables; */
   AstChebyMap *this;
   double *px1;
   double *pxp1;
   double maxx;
   double minx;
   double off;
   double pmax;
   double pmin;
   double scl;
   double x1;
   int k;
   int w1;

/* Check the local error status. */
   if ( !astOK ) return;

/* Get a pointer to the ChebyMap structure. */
   this = (AstChebyMap *) this_polymap;

/* Find the bounds of the supplied x1 values. */
   px1 = table[ 0 ];
   minx = *px1;
   maxx = *px1;
   px1++;
   for( k = 1; k < data->nsamp; k++,px1++ ) {
      if( *px1 > maxx ) {
         maxx = *px1;
      } else if( *px1 < minx ) {
         minx = *px1;
      }
   }

/* Transform the above limits from table values into PolyMap axis values. */
   pmax = maxx*scales[ 0 ];
   pmin = minx*scales[ 0 ];

/* Calculate the scale and offset that map the above bounds onto the range
   [-1,+1], and store them in the ChebyMap. */
   if( pmax != pmin ) {
      scl = 2.0/( pmax - pmin );
      off = -( pmax + pmin )/( pmax - pmin );
   } else if( astOK ){
      astError( AST__BADBX, "astPolyTran(%s): New bounding box has zero width "
                "on axis 1.", status, astGetClass(this));
   }

   if( forward != astGetInvert( this ) ) {
      this->scale_f = (double *) astFree( this->scale_f );
      this->offset_f = (double *) astFree( this->offset_f );

      this->scale_f = (double *) astMalloc( sizeof( double ) );
      this->offset_f = (double *) astMalloc( sizeof( double ) );
      if( astOK ) {
         this->scale_f[ 0 ] = scl;
         this->offset_f[ 0 ] = off;
      }
   } else {
      this->scale_i = (double *) astFree( this->scale_i );
      this->offset_i = (double *) astFree( this->offset_i );

      this->scale_i = (double *) astMalloc( sizeof( double ) );
      this->offset_i = (double *) astMalloc( sizeof( double ) );
      if( astOK ) {
         this->scale_i[ 0 ] = scl;
         this->offset_i[ 0 ] = off;
      }
   }

/* Get pointers to the supplied x1 values. */
   px1 = table[ 0 ];

/* Get pointers to the location for the next "power" of x1. Here "X to
   the power N" is a metaphor for Tn(x). */
   pxp1 = data->xp1;

/* Loop round all samples. */
   for( k = 0; k < data->nsamp; k++ ) {

/* Get the current x1 value, and scale it into the range [-1,+1]. */
      x1 = *(px1++)*scl*scales[0] + off;

/* Find all the required "powers" of x1 and store them in the "xp1"
   component of the data structure. */
      *(pxp1++) = 1.0;
      *(pxp1++) = x1;
      for( w1 = 2; w1 < data->order; w1++,pxp1++ ) {
         pxp1[ 0 ] = 2.0*x1*pxp1[ -1 ] - pxp1[ -2 ];
      }
   }

/* The scaling representing by the scales[0] value will be performed by
   the astTransform method of the ChebyMap class, so reset teh scales[0]
   value to unity, to avoid the scaling being applied twice. */
   scales[ 0 ] = 1.0;

}

static void FitPoly2DInit( AstPolyMap *this_polymap, int forward, double **table,
                           AstMinPackData *data, double *scales, int *status ){
/*
*  Name:
*     FitPoly2DInit

*  Purpose:
*     Perform initialisation needed for FitPoly2D

*  Type:
*     Private function.

*  Synopsis:
*     #include "chebymap.h"
*     void FitPoly2DInit( AstPolyMap *this, int forward, double **table,
*                         AstMinPackData *data, double *scales, int *status )

*  Class Membership:
*     ChebyMap member function (over-rides the astFitPoly2DInit protected
*     method inherited from the PolyMap class).

*  Description:
*     This function performs initialisation needed for FitPoly2D in the
*     PolyMap class..

*  Parameters:
*     this
*        Pointer to the PolyMap.
*     forward
*        Non-zero if the forward transformation of "this" is being
*        replaced. Zero if the inverse transformation of "this" is being
*        replaced.
*     table
*        Pointer to an array of 4 pointers. Each of these pointers points
*        to an array of "nsamp" doubles, being the scaled and sampled values
*        for x1, x2, y1 or y2 in that order.
*     data
*        Pointer to a structure holding information to pass the the
*        service function invoked by the minimisation function.
*     scales
*        Array holding the scaling factors for the four columns of the table.
*        Multiplying the table values by the scale factor produces PolyMap
*        input or output axis values.
*/

/* Local Variables; */
   AstChebyMap *this;
   double *px1;
   double *px2;
   double *pxp1;
   double *pxp2;
   double maxx;
   double maxy;
   double minx;
   double miny;
   double off[ 2 ];
   double pxmax;
   double pxmin;
   double pymax;
   double pymin;
   double scl[ 2 ];
   double x1;
   double x2;
   int k;
   int w1;
   int w2;

/* Check the local error status. */
   if ( !astOK ) return;

/* Get a pointer to the ChebyMap structure. */
   this = (AstChebyMap *) this_polymap;

/* Find the bounds of the supplied x1 and x2 values. */
   px1 = table[ 0 ];
   px2 = table[ 1 ];
   minx = *px1;
   maxx = *px1;
   miny = *px2;
   maxy = *px2;
   px1++;
   px2++;
   for( k = 1; k < data->nsamp; k++,px1++,px2++ ) {
      if( *px1 > maxx ) {
         maxx = *px1;
      } else if( *px1 < minx ) {
         minx = *px1;
      }
      if( *px2 > maxy ) {
         maxy = *px2;
      } else if( *px2 < miny ) {
         miny = *px2;
      }
   }

/* Transform the above limits from table values into PolyMap axis values. */
   pxmax = maxx*scales[ 0 ];
   pxmin = minx*scales[ 0 ];
   pymax = maxy*scales[ 1 ];
   pymin = miny*scales[ 1 ];

/* Calculate the scale and offset that map the above bounds onto the range
   [-1,+1], and store them in the ChebyMap. */
   if( pxmax != pxmin && pymax != pymin ) {
      scl[ 0 ] = 2.0/( pxmax - pxmin );
      off[ 0 ] = -( pxmax + pxmin )/( pxmax - pxmin );
      scl[ 1 ] = 2.0/( pymax - pymin );
      off[ 1 ] = -( pymax + pymin )/( pymax - pymin );
   } else if( astOK ){
      astError( AST__BADBX, "astPolyTran(%s): New bounding box has zero width "
                "on or both axes.", status, astGetClass(this));
   }

   if( forward != astGetInvert( this ) ) {
      this->scale_f = (double *) astFree( this->scale_f );
      this->offset_f = (double *) astFree( this->offset_f );

      this->scale_f = (double *) astMalloc( 2*sizeof( double ) );
      this->offset_f = (double *) astMalloc( 2*sizeof( double ) );
      if( astOK ) {
         this->scale_f[ 0 ] = scl[ 0 ];
         this->offset_f[ 0 ] = off[ 0 ];
         this->scale_f[ 1 ] = scl[ 1 ];
         this->offset_f[ 1 ] = off[ 1 ];
      }
   } else {
      this->scale_i = (double *) astFree( this->scale_i );
      this->offset_i = (double *) astFree( this->offset_i );

      this->scale_i = (double *) astMalloc( 2*sizeof( double ) );
      this->offset_i = (double *) astMalloc( 2*sizeof( double ) );
      if( astOK ) {
         this->scale_i[ 0 ] = scl[ 0 ];
         this->offset_i[ 0 ] = off[ 0 ];
         this->scale_i[ 1 ] = scl[ 1 ];
         this->offset_i[ 1 ] = off[ 1 ];
      }
   }

/* Get pointers to the supplied x1 and x2 values. */
   px1 = table[ 0 ];
   px2 = table[ 1 ];

/* Get pointers to the location for the next "power" of x1 anmd x2. Here "X to
   the power N" is a metaphor for Tn(x). */
   pxp1 = data->xp1;
   pxp2 = data->xp2;

/* Loop round all samples. */
   for( k = 0; k < data->nsamp; k++ ) {

/* Get the current x1 and x2 values, and scale them into the range [-1,+1]. */
      x1 = *(px1++)*scl[0]*scales[0] + off[0];
      x2 = *(px2++)*scl[1]*scales[1] + off[1];

/* Find all the required "powers" of x1 and store them in the "xp1"
   component of the data structure. */
      *(pxp1++) = 1.0;
      *(pxp1++) = x1;
      for( w1 = 2; w1 < data->order; w1++,pxp1++ ) {
         pxp1[ 0 ] = 2.0*x1*pxp1[ -1 ] - pxp1[ -2 ];
      }

/* Find all the required "powers" of x2 and store them in the "xp2"
   component of the data structure. */
      *(pxp2++) = 1.0;
      *(pxp2++) = x2;
      for( w2 = 2; w2 < data->order; w2++,pxp2++ ) {
         pxp2[ 0 ] = 2.0*x2*pxp2[ -1 ] - pxp2[ -2 ];
      }
   }

/* The scaling representing by the scales[0] and scales[1] values will be
   performed by the astTransform method of the ChebyMap class, so reset the
   scales[0] and scales[1] values to unity, to avoid the scaling being
   applied twice. */
   scales[ 0 ] = 1.0;
   scales[ 1 ] = 1.0;

}

static int CompareDerivTerms( const void *a, const void *b ) {
/*
*  Name:
*     CompareDerivTerms

*  Purpose:
*     Compare the output axes and orders of two derivative terms.

*  Type:
*     Private function.

*  Synopsis:
*     int CompareDerivTerms( const void *a, const void *b )

*  Description:
*     This function compares two ChebyDerivTerm structures for use by
*     qsort. Terms are ordered first by output axis and then
*     lexicographically by their polynomial orders. On the
*     differentiated axis the derivative order is used; on each other
*     axis the original order is used.
*
*     Coefficient values are excluded from the comparison, so terms that
*     can be combined into one coefficient compare equal. The comparator
*     uses only information in the supplied structures and no global
*     state.

*  Parameters:
*     a
*        Pointer to the first ChebyDerivTerm structure.
*     b
*        Pointer to the second ChebyDerivTerm structure. Both terms must
*        have the same number of input axes.

*  Returned Value:
*     Minus one if "a" precedes "b", plus one if "a" follows "b", or
*     zero if the terms have the same output axis and polynomial orders.

*  Notes:
*     - This function does not use the inherited status, since its
*     calling sequence is fixed by qsort.
*/
   const ChebyDerivTerm *ta = a;
   const ChebyDerivTerm *tb = b;
   int i, pa, pb;
   if( ta->output != tb->output ) return ta->output < tb->output ? -1 : 1;
   for( i = 0; i < ta->nin; i++ ) {
      pa = i == ta->axis ? ta->degree : ta->powers[ i ];
      pb = i == tb->axis ? tb->degree : tb->powers[ i ];
      if( pa != pb ) return pa < pb ? -1 : 1;
   }
   return 0;
}

static AstPolyMap **GetJacobian( AstPolyMap *map, int *status ) {
/*
*  Name:
*     GetJacobian

*  Purpose:
*     Get the Jacobian of the original forward transformation of a
*     ChebyMap.

*  Type:
*     Private function.

*  Synopsis:
*     #include "polymap.h"
*     AstPolyMap **GetJacobian( AstPolyMap *map, int *status )

*  Class Membership:
*     ChebyMap member function (over-rides the astGetJacobian protected
*     method inherited from the parent PolyMap class).

*  Description:
*     This function returns one derivative Mapping for each original
*     input axis. Each Mapping takes all original inputs and returns the
*     derivatives of all original outputs with respect to that axis,
*     thereby evaluating one column of the Jacobian. The Invert
*     attribute is ignored.
*
*     For a Chebyshev forward series, each derivative is represented as
*     another first-kind ChebyMap. The derivative of T_n contains orders
*     n-1, n-3, ... with coefficients 2*n, except that the coefficient
*     of order zero is n. Each coefficient also includes the physical
*     input normalisation scale. Orders on other axes are unchanged, and
*     terms with equal output axes and orders are combined.
*
*     The derivative Maps retain the exact forward normalisation and
*     have iterative inversion disabled. A zero column is represented by
*     a defined zero transformation. The Maps are cached for subsequent
*     calls. If the original forward series is an ordinary polynomial,
*     the parent PolyMap implementation is used instead.

*  Parameters:
*     map
*        Pointer to the ChebyMap, supplied as a PolyMap pointer. The
*        original forward transformation must be defined.
*     status
*        Pointer to the inherited status variable.

*  Returned Value:
*     Pointer to an array of PolyMap pointers, with one element per
*     original input axis, or NULL if the original forward
*     transformation is undefined or an error occurs. The array and its
*     Maps belong to "map" and must not be modified, freed or annulled
*     by the caller.

*  Notes:
*     - The returned pointer remains valid only while "map" retains its
*     cached Jacobian. Coefficient replacement can invalidate it.
*     - A NULL pointer is returned if the inherited status is set, or if
*     an error occurs.
*/
   AstChebyMap *this = (AstChebyMap *) map;
   AstChebyMap *deriv;
   ChebyDerivTerm *terms = NULL;
   double *coeffs = NULL;
   double *lbnd = NULL;
   double *ubnd = NULL;
   double *pc;
   size_t count, iterm;
   int nin, nout, axis, out, ico, degree, n, i, nco;

   if( !astOK ) return NULL;
   if( !this->scale_f ) return (*parent_getjacobian)( map, status );
   if( map->jacobian ) return map->jacobian;

   nin = ((AstMapping *) map)->nin;
   nout = ((AstMapping *) map)->nout;
   if( !map->ncoeff_f ) return NULL;
   map->jacobian = astCalloc( nin, sizeof( *map->jacobian ) );
   lbnd = astMalloc( nin*sizeof( *lbnd ) );
   ubnd = astMalloc( nin*sizeof( *ubnd ) );
/* The scale and offset derived from this box are overwritten below with
   the values held by "map", so use the canonical box that those values
   normalise to. Reconstructing the physical box would add rounding error,
   and for a box narrower than the resolution of its own centre the two
   reconstructed bounds coincide, which the constructor rejects. */
   if( astOK ) {
      for( i = 0; i < nin; i++ ) {
         lbnd[ i ] = -1.0;
         ubnd[ i ] = 1.0;
      }
   }

   for( axis = 0; axis < nin && astOK; axis++ ) {
      count = 0;
      for( out = 0; out < nout; out++ ) {
         for( ico = 0; ico < map->ncoeff_f[ out ]; ico++ ) {
            n = map->power_f[ out ][ ico ][ axis ];
            count += (size_t) n/2 + n%2;
         }
      }
      if( count > INT_MAX ) {
         astError( AST__INTER, "GetJacobian(%s): Too many derivative terms.",
                   status, astGetClass( this ) );
         break;
      }
      terms = astMalloc( astMAX( count, 1 )*sizeof( *terms ) );
      coeffs = astCalloc( astMAX( count, 1 ),
                         (nin + 2)*sizeof( *coeffs ) );
      if( astOK ) {
         iterm = 0;
         for( out = 0; out < nout; out++ ) {
            for( ico = 0; ico < map->ncoeff_f[ out ]; ico++ ) {
               n = map->power_f[ out ][ ico ][ axis ];
               for( degree = n - 1; degree >= 0; degree -= 2 ) {
                  ChebyDerivTerm *term = terms + iterm++;
                  term->powers = map->power_f[ out ][ ico ];
                  term->nin = nin;
                  term->axis = axis;
                  term->degree = degree;
                  term->output = out + 1;
                  term->coeff = map->coeff_f[ out ][ ico ];
                  if( term->coeff != AST__BAD ) {
                     term->coeff *= this->scale_f[ axis ]*n*
                                    (degree ? 2.0 : 1.0);
                  }
               }
            }
         }

/* Combine equal terms before constructing the derivative Mapping. */
         qsort( terms, count, sizeof( *terms ), CompareDerivTerms );
         nco = 0;
         pc = coeffs;
         for( iterm = 0; iterm < count; iterm++ ) {
            ChebyDerivTerm *term = terms + iterm;
            if( iterm && !CompareDerivTerms( term, term - 1 ) ) {
               pc[ 0 ] = pc[ 0 ] == AST__BAD || term->coeff == AST__BAD ?
                         AST__BAD : pc[ 0 ] + term->coeff;
            } else {
               pc = coeffs + (size_t) nco++*(nin + 2);
               pc[ 0 ] = term->coeff;
               pc[ 1 ] = term->output;
               for( i = 0; i < nin; i++ ) {
                  pc[ i + 2 ] = i == axis ? term->degree : term->powers[ i ];
               }
            }
         }

/* No terms means a defined zero column, not an undefined transformation. */
         if( !nco ) {
            nco = 1;
            coeffs[ 1 ] = 1.0;
         }
         deriv = astChebyMap( nin, nout, nco, coeffs, 0, NULL,
                              lbnd, ubnd, NULL, NULL, "IterInverse=0", status );
         if( astOK ) {
            memcpy( deriv->scale_f, this->scale_f, nin*sizeof( double ) );
            memcpy( deriv->offset_f, this->offset_f, nin*sizeof( double ) );
         }
         map->jacobian[ axis ] = (AstPolyMap *) deriv;
      }
      coeffs = astFree( coeffs );
      terms = astFree( terms );
   }
   lbnd = astFree( lbnd );
   ubnd = astFree( ubnd );
   if( !astOK && map->jacobian ) {
      for( i = 0; i < nin; i++ ) {
         if( map->jacobian[i] ) map->jacobian[i] = astAnnul( map->jacobian[i] );
      }
      map->jacobian = astFree( map->jacobian );
   }
   return astOK ? map->jacobian : NULL;
}

static int GetIterDomain( AstChebyMap *this, double *lbnd, double *ubnd,
                          int *status ) {
/*
*  Name:
*     GetIterDomain

*  Purpose:
*     Get the finite domain for ChebyMap inverse iteration.

*  Type:
*     Private function.

*  Synopsis:
*     #include "chebymap.h"
*     int GetIterDomain( AstChebyMap *this, double *lbnd, double *ubnd,
*                        int *status )

*  Class Membership:
*     ChebyMap member function.

*  Description:
*     This function returns the physical input bounds used to restrict
*     iterative inversion of the original forward Chebyshev series. The
*     bounds are recovered from the stored forward scales and offsets,
*     independently of the Invert attribute.
*
*     A ChebyMap can also contain an ordinary polynomial forward
*     transformation. Such a transformation has no finite iteration
*     domain; in that case both supplied arrays are left unchanged and
*     zero is returned.
*
*     A reconstructed bound can round to a position whose normalised
*     coordinate lies just outside [-1,1]. Each bound is moved towards
*     the other by up to eight adjacent representable values to bring it
*     within the range accepted by the forward evaluator. This does not
*     enable extrapolation of the forward series.

*  Parameters:
*     this
*        Pointer to the ChebyMap.
*     lbnd
*        Pointer to an array in which to return the lower bound on each
*        original input axis. Its length must equal the number of inputs
*        of the uninverted ChebyMap.
*     ubnd
*        Pointer to an array in which to return the upper bound on each
*        original input axis. Its length and ordering must match "lbnd".
*     status
*        Pointer to the inherited status variable.

*  Returned Value:
*     One if bounds have been supplied, or zero if the original forward
*     transformation has no Chebyshev normalisation.

*  Notes:
*     - A value of zero is returned and the bound arrays are left
*     unchanged if the inherited status is set.
*/
   int i, nin = ((AstMapping *) this)->nin;
   double a, b;
   if( !astOK || !this->scale_f ) return 0;

/* Every axis must yield a usable interval before anything is stored, so
   that the supplied arrays are left unchanged when zero is returned. */
   for( i = 0; i < nin; i++ ) {
      if( !AxisBounds( this->scale_f[i], this->offset_f[i], &a, &b ) ) return 0;
   }

   for( i = 0; i < nin; i++ ) {
      (void) AxisBounds( this->scale_f[i], this->offset_f[i], lbnd + i,
                         ubnd + i );
   }
   return 1;
}

static double NudgeIntoDomain( double bound, double towards, double scale,
                               double offset ) {
/*
*  Name:
*     NudgeIntoDomain

*  Purpose:
*     Move a bound to a position the forward evaluator accepts.

*  Type:
*     Private function.

*  Synopsis:
*     #include "chebymap.h"
*     double NudgeIntoDomain( double bound, double towards, double scale,
*                             double offset )

*  Description:
*     The forward Chebyshev series is evaluated only where the normalised
*     coordinate "bound*scale + offset" lies in [-1,+1]. Recovering a bound
*     from the scale and offset can round it to a position just outside
*     that range. This function moves the bound towards the other end of
*     the interval, one representable value at a time, for as long as that
*     reduces the magnitude of the normalised coordinate.
*
*     Stopping as soon as a step makes no improvement leaves the bound
*     unchanged where no representable position is evaluable, which the
*     caller detects. Such normalisations describe an interval narrower
*     than the spacing of the values around it.

*  Parameters:
*     bound
*        The bound to be moved.
*     towards
*        The other end of the interval, giving the direction to move in.
*     scale
*        The scale factor applied to an axis value before evaluation.
*     offset
*        The offset added to an axis value after scaling.

*  Returned Value:
*     The moved bound.
*/
   double next;
   double nresid;
   double resid = fabs( bound*scale + offset );

   while( resid > 1.0 ) {
      next = nextafter( bound, towards );
      if( next == bound ) break;
      nresid = fabs( next*scale + offset );
      if( !( nresid < resid ) ) break;
      bound = next;
      resid = nresid;
   }

   return bound;
}

static int AxisBounds( double scale, double offset, double *lbnd,
                       double *ubnd ) {
/*
*  Name:
*     AxisBounds

*  Purpose:
*     Recover the evaluable interval on one axis from a Chebyshev normalization.

*  Type:
*     Private function.

*  Synopsis:
*     #include "chebymap.h"
*     int AxisBounds( double scale, double offset, double *lbnd,
*                     double *ubnd )

*  Description:
*     This function returns the range of axis values over which a
*     Chebyshev series with the given normalisation can be evaluated. The
*     bounds are the positions whose normalised coordinates are -1 and +1,
*     moved inwards if necessary so that the forward evaluator accepts
*     them.
*
*     Zero is returned, and the supplied bounds are left unchanged, if the
*     normalisation describes no interval at all. That covers a zero or
*     non-finite scale or offset, bounds that overflow, and bounds that
*     cannot be moved to positions the evaluator accepts without crossing.
*     An interval containing a single representable value is returned as
*     such; it is for the caller to decide what to do with it.

*  Parameters:
*     scale
*        The scale factor applied to an axis value before evaluation.
*     offset
*        The offset added to an axis value after scaling.
*     lbnd
*        Pointer to a double in which to return the lower bound.
*     ubnd
*        Pointer to a double in which to return the upper bound.

*  Returned Value:
*     One if bounds have been returned, zero otherwise.
*/
   double a;
   double b;
   double hi;
   double lo;

   if( scale == 0.0 || !isfinite( scale ) || !isfinite( offset ) ) return 0;

   a = ( -1.0 - offset )/scale;
   b = ( 1.0 - offset )/scale;
   if( !isfinite( a ) || !isfinite( b ) ) return 0;

   lo = NudgeIntoDomain( astMIN( a, b ), astMAX( a, b ), scale, offset );
   hi = NudgeIntoDomain( astMAX( a, b ), lo, scale, offset );

   if( fabs( lo*scale + offset ) > 1.0 ) return 0;
   if( fabs( hi*scale + offset ) > 1.0 ) return 0;

   *lbnd = lo;
   *ubnd = hi;
   return 1;
}

static int Usable( double value ) {
/*
*  Name:
*     Usable

*  Purpose:
*     Test whether a coordinate value can take part in the iteration.

*  Type:
*     Private function.

*  Synopsis:
*     #include "chebymap.h"
*     int Usable( double value )

*  Description:
*     This function returns non-zero if "value" is neither AST__BAD nor
*     a non-finite floating point value.

*  Parameters:
*     value
*        The value to test.

*  Returned Value:
*     Non-zero if the value can be used in arithmetic.
*/
   return value != AST__BAD && isfinite( value );
}

static void MarkUnsolved( double **inputs, int ncoord, int ipoint,
                          int *flags, int *nconv ) {
/*
*  Name:
*     MarkUnsolved

*  Purpose:
*     Record that a batch position has no solution.

*  Type:
*     Private function.

*  Synopsis:
*     #include "chebymap.h"
*     void MarkUnsolved( double **inputs, int ncoord, int ipoint,
*                        int *flags, int *nconv )

*  Description:
*     This function sets every coordinate of position "ipoint" in
*     "inputs" to AST__BAD, sets its resolved flag and increments the
*     count of resolved positions, so that the iteration does not touch
*     the position again.

*  Parameters:
*     inputs
*        Array of pointers to input coordinate arrays, indexed as
*        inputs[axis][point].
*     ncoord
*        The number of axes.
*     ipoint
*        The index of the position within the batch.
*     flags
*        Array of resolved flags, one per position.
*     nconv
*        Pointer to the count of resolved positions.
*/
   int i;
   for( i = 0; i < ncoord; i++ ) inputs[ i ][ ipoint ] = AST__BAD;
   flags[ ipoint ] = 1;
   (*nconv)++;
}

static void IterSteps( AstPolyMap *this, int npoint, int ncoord,
                       double **inputs, double **outputs, double **work,
                       const double *lbnd, const double *ubnd,
                       const double *width, const double *scales,
                       const double *steps, const double *norms,
                       int *stepping, int *flags, int *nconv, int *status ) {
/*
*  Name:
*     IterSteps

*  Purpose:
*     Find a bounded Newton step for each unresolved batch position.

*  Type:
*     Private function.

*  Synopsis:
*     void IterSteps( AstPolyMap *this, int npoint, int ncoord,
*                     double **inputs, double **outputs, double **work,
*                     const double *lbnd, const double *ubnd,
*                     const double *width, const double *scales,
*                     const double *steps, const double *norms,
*                     int *stepping, int *flags, int *nconv, int *status )

*  Description:
*     This function applies one bounded Newton update to each position
*     flagged in "stepping", during a batch of inverse transformations.
*     Each trial position is formed by adding the scaled correction to
*     the current input position and projecting it onto the finite
*     forward domain. The original forward transformation is then
*     evaluated at that position.
*
*     The full correction is tried first, followed by successive
*     halvings, with at most 32 trials. A trial is accepted only if its
*     forward values and scaled residuals are finite and its maximum
*     absolute scaled residual is strictly smaller than the position's
*     entry in "norms". The output scaling is held fixed throughout this
*     search. A correction made small by projection is not by itself
*     evidence of convergence.
*
*     An accepted trial replaces the input position, and its forward
*     value is copied into "work" so the caller does not need to
*     re-evaluate the forward transformation for that position before
*     the next iteration. A position for which no acceptable trial is
*     found is set bad and marked as converged. Either way its
*     "stepping" flag is cleared, so the array is left ready for the
*     next iteration.
*
*     All the positions still searching at a given trial size are
*     evaluated by a single call to the forward transformation. The
*     positions that have finished drop out of the following trials.
*
*     The two working PointSets used to hold a trial and its forward
*     value are created once, sized for the largest trial this call will
*     make, and reused at every backtrack level.

*  Parameters:
*     this
*        Pointer to the ChebyMap, supplied as a PolyMap pointer. Its
*        original forward transformation
*        must be defined, with equal numbers of inputs and outputs. The
*        Invert attribute does not change the transformation being
*        evaluated.
*     npoint
*        The number of positions in the batch.
*     ncoord
*        The number of original input axes, which equals the number of
*        original output axes.
*     inputs
*        Array of pointers to input coordinate arrays, indexed as
*        inputs[axis][point]. There must be one array per original input
*        axis. Only flagged positions are modified.
*     outputs
*        Array of pointers to target output coordinate arrays, indexed
*        as outputs[axis][point]. There must be one array per original
*        output axis. These values are not modified.
*     work
*        Array of pointers to the whole batch's original forward values,
*        indexed as work[axis][point]. There must be one array per
*        original output axis. The forward value of an accepted trial is
*        copied here for its position; other positions are not touched.
*     lbnd
*        Array holding the finite lower domain bound on each original
*        input axis.
*     ubnd
*        Array holding the finite upper domain bound on each original
*        input axis, in the same order as "lbnd".
*     width
*        Array holding the positive half-width of the domain on each
*        original input axis, used to convert normalised corrections
*        into physical coordinate offsets.
*     scales
*        Array holding the non-negative residual scale for each original
*        output axis of each position, indexed as
*        scales[point*ncoord+axis]. A positive scale divides that output
*        residual; a zero scale requires an exactly zero residual.
*     steps
*        Array holding the Newton correction on each original input axis
*        of each position, divided by the corresponding value in
*        "width", indexed as steps[point*ncoord+axis].
*     norms
*        Array holding, for each position, the maximum absolute scaled
*        forward residual at its current input position.
*     stepping
*        Array holding a non-zero value for each position that needs an
*        update. Every element is zero on exit.
*     flags
*        Array of convergence flags, set for each position that is
*        resolved here.
*     nconv
*        Pointer to the count of resolved positions, incremented for
*        each position set bad.
*     status
*        Pointer to the inherited status variable.

*  Returned Value:
*     void

*  Notes:
*     - Failure to find an acceptable step does not set the inherited
*     status.
*     - This function returns without action if the inherited status is
*     set.
*/

/* Local Variables: */
   AstPointSet *trial;
   AstPointSet *trial_out;
   double **x;
   double **y;
   double alpha;
   double newnorm;
   double residual;
   double value;
   int *index;
   int backtrack;
   int changed;
   int fwd;
   int i;
   int ipoint;
   int jpoint;
   int nsearch;
   int nstart;
   int ntrial;
   int valid;

/* Check inherited status */
   if( !astOK ) return;

/* Count the positions needing an update, and allocate an array to hold
   the batch index of each position still searching. */
   nsearch = 0;
   for( ipoint = 0; ipoint < npoint; ipoint++ ) {
      if( stepping[ ipoint ] ) nsearch++;
   }
   if( nsearch == 0 ) return;

   index = astMalloc( sizeof( *index )*nsearch );
   if( !astOK ) return;

   fwd = !astGetInvert( this );

/* Create the two trial PointSets once, sized for the largest batch any
   backtrack level will need, and reuse them throughout. */
   trial = astPointSet( nsearch, ncoord, "", status );
   trial_out = astPointSet( nsearch, ncoord, "", status );
   x = astGetPoints( trial );
   if( !astOK ) {
      trial = astAnnul( trial );
      trial_out = astAnnul( trial_out );
      index = astFree( index );
      return;
   }

/* Try the full correction first, then successive halvings. */
   alpha = 1.0;
   for( backtrack = 0; backtrack < 32 && nsearch > 0 && astOK;
        backtrack++, alpha *= 0.5 ) {

/* Gather the trial positions. A position whose trial position is the one
   it already holds cannot be improved by any smaller correction, so it
   stops searching. The number of positions still searching only ever
   shrinks between levels, so the PointSets created above are shrunk to
   fit rather than re-created. */
      nstart = nsearch;
      astSetNpoint( trial, nstart );
      astSetNpoint( trial_out, nstart );

      ntrial = 0;
      for( ipoint = 0; ipoint < npoint; ipoint++ ) {
         if( !stepping[ ipoint ] ) continue;
         changed = 0;
         for( i = 0; i < ncoord; i++ ) {
            value = inputs[ i ][ ipoint ] +
                    alpha*steps[ ipoint*ncoord + i ]*width[ i ];
            value = astMAX( lbnd[ i ], astMIN( ubnd[ i ], value ) );
            x[ i ][ ntrial ] = value;
            if( value != inputs[ i ][ ipoint ] ) changed = 1;
         }
         if( changed ) {
            index[ ntrial++ ] = ipoint;
         } else {
            MarkUnsolved( inputs, ncoord, ipoint, flags, nconv );
            stepping[ ipoint ] = 0;
            nsearch--;
         }
      }

/* Evaluate them all at once. */
      if( ntrial > 0 ) {
         if( ntrial < nstart ) {
            astSetNpoint( trial, ntrial );
            astSetNpoint( trial_out, ntrial );
         }
         (void) astTransform( this, trial, fwd, trial_out );
         y = astGetPoints( trial_out );
      }

/* Accept the first trial that reduces the scaled residual. */
      for( jpoint = 0; jpoint < ntrial && astOK; jpoint++ ) {
         ipoint = index[ jpoint ];
         valid = 1;
         newnorm = 0.0;
         for( i = 0; i < ncoord; i++ ) {
            value = y[ i ][ jpoint ];
            if( !Usable( value ) ) {
               valid = 0;
               break;
            }
            residual = fabs( outputs[ i ][ ipoint ] - value );
            if( scales[ ipoint*ncoord + i ] > 0.0 ) {
               residual /= scales[ ipoint*ncoord + i ];
            }
            if( !isfinite( residual ) ||
                ( scales[ ipoint*ncoord + i ] == 0.0 && residual != 0.0 ) ) {
               valid = 0;
               break;
            }
            newnorm = astMAX( newnorm, residual );
         }
         if( valid && newnorm < norms[ ipoint ] ) {
            for( i = 0; i < ncoord; i++ ) {
               inputs[ i ][ ipoint ] = x[ i ][ jpoint ];
            }
            for( i = 0; i < ncoord; i++ ) {
               work[ i ][ ipoint ] = y[ i ][ jpoint ];
            }
            stepping[ ipoint ] = 0;
            nsearch--;
         }
      }
   }

   trial = astAnnul( trial );
   trial_out = astAnnul( trial_out );

/* Any position that exhausted the trials has no acceptable step. */
   for( ipoint = 0; ipoint < npoint; ipoint++ ) {
      if( stepping[ ipoint ] ) {
         MarkUnsolved( inputs, ncoord, ipoint, flags, nconv );
         stepping[ ipoint ] = 0;
      }
   }

   index = astFree( index );
}

static void IterInverse( AstPolyMap *map, AstPointSet *out,
                         AstPointSet *result, int *status ){
/*
*  Name:
*     IterInverse

*  Purpose:
*     Evaluate the original inverse transformation of a ChebyMap by
*     bounded Newton iteration.

*  Type:
*     Private function.

*  Synopsis:
*     #include "polymap.h"
*     void IterInverse( AstPolyMap *map, AstPointSet *out,
*                       AstPointSet *result, int *status )

*  Class Membership:
*     ChebyMap member function (over-rides the astIterInverse protected
*     method inherited from the PolyMap class).

*  Description:
*     This function transforms a set of original output positions into
*     original input positions using Newton-Raphson iteration restricted
*     to the forward bounding box of the ChebyMap. Initial guesses come
*     from astLinearGuess and are clipped into the box. Each Newton
*     correction is normalised by the box half-widths, each residual by
*     its Jacobian row norm, and IterSteps backtracks corrections that
*     do not reduce the residual. A position converges when both the
*     normalised correction and the scaled residual are within
*     TolInverse; an exactly zero residual converges immediately.
*
*     After NiterInverse updates the final candidate is checked once
*     more. Any position that is still unsolved, was supplied with a bad
*     or non-finite coordinate, or met a singular Jacobian is returned
*     as AST__BAD without setting the inherited status.
*
*     A position whose Jacobian is singular is moved once by a quarter of
*     each half-width before being declared unsolved.
*
*     The whole-batch forward transformation used to form the residual is
*     only re-evaluated on the first iteration and after a position has
*     been nudged off a singular seed; otherwise the forward values IterSteps
*     already obtained for the accepted trial are reused.
*
*     If the forward series has no usable Chebyshev box on every axis,
*     the parent PolyMap algorithm is used instead.

*  Parameters:
*     map
*        Pointer to the ChebyMap, supplied as a PolyMap pointer.
*     out
*        PointSet holding the original output positions to invert.
*     result
*        PointSet to receive the original input positions.
*     status
*        Pointer to the inherited status variable.

*  Notes:
*     - TolInverse measures each input correction as a fraction of the
*     corresponding box half-width. Each output residual is divided by
*     the sum of the absolute Jacobian elements in its row, multiplied
*     by the respective input half-widths.
*     - This function returns without action if the inherited status is
*     set.
*/

/* Local Variables: */
   AstChebyMap *this;
   AstMapping *lintrunc;
   AstPointSet *work;
   AstPointSet **ps_jac;
   AstPolyMap **jacob;
   double ***ptr_jac;
   double **ptr_in;
   double **ptr_out;
   double **ptr_work;
   double *lbnd;
   double *mat;
   double *norms;
   double *pa;
   double *pb;
   double *scale;
   double *scales;
   double *steps;
   double *ubnd;
   double *vec;
   double *width;
   double det;
   double norm;
   double stepnorm;
   double tol;
   double xx;
   int *flags;
   int *iw;
   int *nudged;
   int *stepping;
   int exact;
   int fwd;
   int icol;
   int icoord;
   int ipoint;
   int irow;
   int iter;
   int maxiter;
   int nconv;
   int ncoord;
   int npoint;
   int sing;
   int stale;
   int valid;

/* Check inherited status */
   if( !astOK ) return;

/* Get a pointer to the ChebyMap structure. */
   this = (AstChebyMap *) map;

/* Check the ChebyMap has equal numbers of inputs and outputs. */
   ncoord = astGetNin( map );
   if( ncoord != astGetNout( map ) ) {
      astError( AST__INTER, "astTransform(%s): Supplied %s has unequal numbers"
                " of inputs and outputs and therefore an iterative inverse "
                "cannot be used (internal AST Programming error).", status,
                astGetClass(map), astGetClass(map) );
      return;
   }

/* Without a Chebyshev normalisation or a usable box there is nothing to
   restrict the iteration to, so use the parent algorithm. */
   lbnd = astMalloc( ncoord*sizeof( *lbnd ) );
   ubnd = astMalloc( ncoord*sizeof( *ubnd ) );
   if( !astOK || !this->scale_f ||
       !GetIterDomain( this, lbnd, ubnd, status ) ) {
      lbnd = astFree( lbnd );
      ubnd = astFree( ubnd );
      (*parent_iterinverse)( map, out, result, status );
      return;
   }

/* Get the Jacobian of the forward transformation as a vector of "ncoord"
   Mappings, each giving one column of the matrix. */
   jacob = astGetJacobian( map );

/* Get the number of points to be transformed. */
   npoint = astGetNpoint( out );

/* Allocate the box half-widths and the per-row residual scales. */
   width = astMalloc( ncoord*sizeof( *width ) );
   scale = astMalloc( ncoord*sizeof( *scale ) );

/* Get another PointSet to hold intermediate results. */
   work = astPointSet( npoint, ncoord, " ", status );

/* See if the ChebyMap has been inverted.*/
   fwd = !astGetInvert( map );

/* Get pointers to the data arrays for all PointSets. Note, here "in" and
   "out" refer to inputs and outputs of the forward transformation. These
   are respectively *outputs* and *inputs* of the inverse transformation. */
   ptr_in = astGetPoints( result );  /* Returned input positions */
   ptr_out = astGetPoints( out );    /* Supplied output positions */
   ptr_work = astGetPoints( work );  /* Work space */

/* Allocate an array of PointSets to hold the elements of the Jacobian
   matrix. */
   ptr_jac = astMalloc( sizeof( double ** )*ncoord );
   ps_jac = astCalloc( ncoord, sizeof( AstPointSet * ) );
   if( astOK ) {
      for( icoord = 0; icoord < ncoord; icoord++ ) {
         ps_jac[ icoord ] = astPointSet( npoint, ncoord, " ", status );
         ptr_jac[ icoord ] = astGetPoints( ps_jac[ icoord ] );
      }
   }

/* Allocate an array to hold flags indicating if each position has
   been resolved. Initialise it to hold zero at every element. */
   flags = astCalloc( npoint, sizeof( int ) );

/* Allocate an array to record whether each position has already been
   nudged off a singular seed. Initialise it to hold zero at every
   element. */
   nudged = astCalloc( npoint, sizeof( int ) );

/* Allocate memory to hold the Jacobian matrix at a single point. */
   mat = astMalloc( sizeof( double )*ncoord*ncoord );

/* Allocate memory to hold the offset vector. */
   vec = astMalloc( sizeof( double )*ncoord );

/* Allocate memory to hold work space for palDmat. */
   iw = astMalloc( sizeof( int )*ncoord );

/* Allocate memory to hold the deferred Newton update for each position of
   the batch. */
   scales = astMalloc( sizeof( double )*npoint*ncoord );
   steps = astMalloc( sizeof( double )*npoint*ncoord );
   norms = astMalloc( sizeof( double )*npoint );
   stepping = astCalloc( npoint, sizeof( int ) );

/* Check pointers can be used safely. */
   if( astOK ) {

/* Store the initial guess at the required input positions. These are
   determined by transforming the supplied output positions using the
   inverse of an affine approximation to the forward transformation. */
      lintrunc = astLinearGuess( map );
      (void) astTransform( lintrunc, out, 0, result );
      lintrunc = astAnnul( lintrunc );

/* Get the maximum number of iterations to perform. */
      maxiter = astGetNiterInverse( map );

/* Get the target normalised error for the returned input axis values. */
      tol = astGetTolInverse( map );

/* Initialise the number of positions which have been resolved. */
      nconv = 0;

/* Initialise the count of positions whose forward value at the top of
   the next iteration cannot be taken from the previous iteration's
   trial (none have been nudged off a singular seed yet). */
      stale = 0;

/* The iteration controls are validated by the NiterInverse and TolInverse
   setters, so only the box need be checked here, then every initial guess
   is clipped into it. A position with a bad or non-finite target, or an
   unusable box, is resolved immediately as bad. */
      valid = 1;
      for( icoord = 0; icoord < ncoord; icoord++ ) {
         width[icoord] = 0.5*ubnd[icoord] - 0.5*lbnd[icoord];
         if( !isfinite(lbnd[icoord]) || !isfinite(ubnd[icoord]) ||
             !isfinite(width[icoord]) || width[icoord] <= 0.0 ) valid = 0;
      }
      for( ipoint = 0; ipoint < npoint; ipoint++ ) {
         int good = valid;
         for( icoord = 0; icoord < ncoord; icoord++ ) {
            if( !Usable( ptr_out[icoord][ipoint] ) ) good = 0;
            xx = ptr_in[icoord][ipoint];
            if( !Usable( xx ) ) xx = 0.5*lbnd[icoord] + 0.5*ubnd[icoord];
            ptr_in[icoord][ipoint] = astMAX( lbnd[icoord],
                                            astMIN(ubnd[icoord],xx) );
         }
         if( !good ) MarkUnsolved( ptr_in, ncoord, ipoint, flags, &nconv );
      }

/* Loop round doing iterations of a Newton-Raphson algorithm, until all
   points have been resolved or the maximum number of updates has been
   performed. The final pass only checks the last candidate. */
      for( iter = 0; iter <= maxiter && nconv < npoint && astOK; iter++ ) {

/* Use the original forward transformation to transform the current
   guesses at the required input positions into the corresponding output
   positions. Store the results in the "work" PointSet. Every position
   still iterating was either just stepped by IterSteps, which recorded
   its forward value, or nudged, which did not; only the first iteration
   and a nudge need a fresh evaluation of the whole batch. */
         if( iter == 0 || stale > 0 ) {
            (void) astTransform( map, result, fwd, work );
            stale = 0;
         }

/* Modify the work PointSet so that it holds the offsets from the output
   positions produced by the current input position guesses, and the
   required output positions. */
         for( icoord = 0; icoord < ncoord; icoord++ ) {
            pa = ptr_out[ icoord ];
            pb = ptr_work[ icoord ];
            for( ipoint = 0; ipoint< npoint; ipoint++,pa++,pb++ ) {
               if( *pa != AST__BAD && *pb != AST__BAD ){
                  *pb = *pa - *pb;
               } else {
                  *pb = AST__BAD;
               }
            }
         }

/* Evaluate the elements of the Jacobian matrix at the current input
   position guesses. */
         for( icoord = 0; icoord < ncoord; icoord++ ) {
            (void) astTransform( jacob[ icoord ], result, 1, ps_jac[ icoord ] );
         }

/* For each position, we now invert the matrix equation

    Dy = Jacobian.Dx

   to find a guess at the vector (dx) holding the offsets from the
   current input positions guesses to their required values. Loop over all
   points. */
         for( ipoint = 0; ipoint < npoint; ipoint++ ) {

/* Do not change positions that have already been resolved. */
            if( !flags[ ipoint ] ) {

/* Get the numerical values for the elements of the Jacobian matrix at
   the current point. */
               pa = mat;
               for( irow = 0; irow < ncoord; irow++ ) {
                  for( icol = 0; icol < ncoord; icol++ ) {
                     *(pa++) = ptr_jac[ icol ][ irow ][ ipoint ];
                  }

/* Store the offset from the current output position to the required
   output position. */
                  vec[ irow ] = ptr_work[ irow ][ ipoint ];
               }

/* Check the residual before the Jacobian: an exact solution is usable
   even at a singular point. BAD values must never reach palDmat. */
               sing = 0;
               norm = 0.0;
               valid = 1;
               exact = 1;
               for( irow = 0; irow < ncoord; irow++ ) {
                  if( !Usable( vec[irow] ) ) valid = 0;
                  if( vec[irow] != 0.0 ) exact = 0;
               }
               if( valid && exact ) {
                  flags[ipoint] = 1;
                  nconv++;
                  continue;
               }

/* Solve for dimensionless steps (dx divided by the input half-width).
   Scale each equation by its Jacobian row norm to avoid dependence on
   arbitrary input or output units in the singularity and residual tests. */
               for( irow = 0; irow < ncoord; irow++ ) {
                  scale[irow] = 0.0;
                  for( icol = 0; icol < ncoord; icol++ ) {
                     pa = mat + irow*ncoord + icol;
                     if( !Usable( *pa ) ) valid = 0;
                     *pa *= width[icol];
                     scale[irow] += fabs(*pa);
                  }
                  if( !isfinite(scale[irow]) ) valid = 0;
                  if( scale[irow] > 0.0 ) {
                     for( icol = 0; icol < ncoord; icol++ ) {
                        mat[irow*ncoord+icol] /= scale[irow];
                     }
                     vec[irow] /= scale[irow];
                  }
                  if( !isfinite(vec[irow]) ) valid = 0;
                  norm = astMAX( norm, fabs(vec[irow]) );
               }
               if( !valid ) sing = 1;
               if( !sing ) palDmat( ncoord, mat, vec, &det, &sing, iw );

/* If the matrix was singular, nudge the position once off the seed that
   produced it and try again next iteration. A position singular a second
   time cannot be evaluated, so store a bad value for it and indicate it
   has been resolved. */
               if( sing ) {
                  if( !nudged[ ipoint ] ) {

/* Move once off a stationary point, towards the side of the box with
   more room, and evaluate the forward transformation there next time. */
                     nudged[ ipoint ] = 1;
                     stale++;
                     for( icoord = 0; icoord < ncoord; icoord++ ) {
                        xx = ptr_in[ icoord ][ ipoint ];
                        if( xx - lbnd[ icoord ] >= ubnd[ icoord ] - xx ) {
                           xx -= 0.25*width[ icoord ];
                        } else {
                           xx += 0.25*width[ icoord ];
                        }
                        ptr_in[ icoord ][ ipoint ] = xx;
                     }
                  } else {
                     MarkUnsolved( ptr_in, ncoord, ipoint, flags, &nconv );
                  }

/* Otherwise, see whether the position has converged, has run out of
   updates, or needs a step. */
               } else {
                  valid = 1;
                  stepnorm = 0.0;
                  for( icoord = 0; icoord < ncoord; icoord++ ) {
                     if( !isfinite(vec[icoord]) ) valid = 0;
                     stepnorm = astMAX( stepnorm, fabs(vec[icoord]) );
                  }
                  if( valid && stepnorm <= tol && norm <= tol ) {
                     flags[ipoint] = 1;
                     nconv++;
                  } else if( !valid || iter == maxiter ) {
                     MarkUnsolved( ptr_in, ncoord, ipoint, flags, &nconv );

/* Defer the backtracking search, so that all the positions looking for an
   acceptable step share one evaluation of the forward transformation at
   each trial size. */
                  } else {
                     for( icoord = 0; icoord < ncoord; icoord++ ) {
                        steps[ipoint*ncoord + icoord] = vec[icoord];
                        scales[ipoint*ncoord + icoord] = scale[icoord];
                     }
                     norms[ipoint] = norm;
                     stepping[ipoint] = 1;
                  }
               }
            }
         }

/* Apply the deferred Newton updates for the whole batch. */
         IterSteps( map, npoint, ncoord, ptr_in, ptr_out, ptr_work, lbnd,
                    ubnd, width, scales, steps, norms, stepping, flags,
                    &nconv, status );
      }
   }

/* Free resources. */
   scales = astFree( scales );
   steps = astFree( steps );
   norms = astFree( norms );
   stepping = astFree( stepping );
   vec = astFree( vec );
   iw = astFree( iw );
   mat = astFree( mat );
   flags = astFree( flags );
   nudged = astFree( nudged );
   work = astAnnul( work );

   if( ps_jac ) {
      for( icoord = 0; icoord < ncoord; icoord++ ) {
         ps_jac[ icoord ] = astAnnul( ps_jac[ icoord ] );
      }
      ps_jac = astFree( ps_jac );
   }

   ptr_jac = astFree( ptr_jac );
   lbnd = astFree( lbnd );
   ubnd = astFree( ubnd );
   width = astFree( width );
   scale = astFree( scale );
}

static AstMapping *LinearGuess( AstPolyMap *map, int *status ) {
/*
*  Name:
*     LinearGuess

*  Purpose:
*     Get a Mapping supplying initial guesses for ChebyMap inversion.

*  Type:
*     Private function.

*  Synopsis:
*     #include "polymap.h"
*     AstMapping *LinearGuess( AstPolyMap *map, int *status )

*  Class Membership:
*     ChebyMap member function (over-rides the astLinearGuess protected
*     method inherited from the parent PolyMap class).

*  Description:
*     This function returns a Mapping whose inverse supplies an initial
*     input position for iterative inversion of the original forward
*     transformation, independently of the Invert attribute. It first
*     tries an affine approximation using the complete forward value and
*     Jacobian at the midpoint of the forward domain.
*
*     If that approximation cannot provide a finite inverse, the
*     constant and linear Chebyshev terms are tried, including their
*     physical input normalisation. If neither approximation is usable,
*     a Mapping whose inverse always returns the domain midpoint is
*     supplied. Clipping initial guesses to the forward domain is the
*     responsibility of the caller.
*
*     The Mapping is cached for subsequent calls. If the original
*     forward transformation is an ordinary polynomial, the parent
*     PolyMap implementation is used instead.

*  Parameters:
*     map
*        Pointer to the ChebyMap, supplied as a PolyMap pointer. For a
*        Chebyshev forward series, the original forward transformation
*        must be defined and the numbers of inputs and outputs must be
*        equal.
*     status
*        Pointer to the inherited status variable.

*  Returned Value:
*     A new reference to the cached Mapping, or NULL if an error occurs
*     or the Chebyshev forward transformation does not meet the stated
*     requirements. The caller must annul the reference when it is no
*     longer required.

*  Notes:
*     - A NULL pointer is returned if the inherited status is set, or if
*     an error occurs.
*/
   AstChebyMap *this = (AstChebyMap *) map;
   AstPolyMap **jac;
   AstMapping *mm = NULL, *sm = NULL, *tmp = NULL, *result = NULL;
   double *center, *value, *column, *matrix;
   int *inperm;
   int nin, i, j, attempt, valid, out, ico, axis, order;
   double c;

   if( !astOK ) return NULL;
   if( !this->scale_f ) return (*parent_linearguess)( map, status );
   if( map->lintrunc ) return astClone( map->lintrunc );
   nin = ((AstMapping *) map)->nin;
   if( nin != ((AstMapping *) map)->nout || !map->ncoeff_f ) return NULL;

   center = astMalloc( nin*sizeof( *center ) );
   value = astMalloc( nin*sizeof( *value ) );
   column = astMalloc( nin*sizeof( *column ) );
   matrix = astMalloc( (size_t) nin*nin*sizeof( *matrix ) );
   inperm = astMalloc( nin*sizeof( *inperm ) );
   jac = astGetJacobian( map );
   if( astOK ) {
      for( i = 0; i < nin; i++ ) center[i] = -this->offset_f[i]/this->scale_f[i];
      astTranN( map, 1, nin, 1, center, !astGetInvert(map), nin, 1, value );
      for( j = 0; j < nin; j++ ) {
         astTranN( jac[j], 1, nin, 1, center, 1, nin, 1, column );
         for( i = 0; i < nin; i++ ) matrix[i*nin+j] = column[i];
      }

      for( attempt = 0; attempt < 2 && astOK; attempt++ ) {
         if( attempt ) {
            memset( matrix, 0, (size_t) nin*nin*sizeof( *matrix ) );
            memset( value, 0, nin*sizeof( *value ) );
            for( out = 0; out < nin; out++ ) {
               for( ico = 0; ico < map->ncoeff_f[out]; ico++ ) {
                  axis = -1;
                  order = 0;
                  for( j = 0; j < nin; j++ ) {
                     if( map->power_f[out][ico][j] ) {
                        if( axis >= 0 || map->power_f[out][ico][j] != 1 ) {
                           order = 2;
                           break;
                        }
                        axis = j;
                        order = 1;
                     }
                  }
                  c = map->coeff_f[out][ico];
                  if( order == 0 ) {
                     value[out] += c;
                  } else if( order == 1 ) {
                     matrix[out*nin+axis] += c*this->scale_f[axis];
                     value[out] += c*(this->scale_f[axis]*center[axis] +
                                     this->offset_f[axis]);
                  }
               }
            }
         }
         valid = 1;
         for( i = 0; i < nin; i++ ) {
            if( value[i] == AST__BAD || !isfinite(value[i]) ) valid = 0;
            for( j = 0; j < nin; j++ ) {
               c = matrix[i*nin+j];
               if( c == AST__BAD || !isfinite(c) ) valid = 0;
            }
         }
         if( valid ) {
            mm = (AstMapping *) astMatrixMap( nin, nin, 0, matrix, "", status );
            if( astGetTranInverse( mm ) ) {
/* A diagonal MatrixMap can advertise an inverse but return BAD for a
   zero diagonal element. Check that its inverse actually yields a seed. */
               astTranN( mm, 1, nin, 1, value, 0, nin, 1, column );
               for( i = 0; i < nin; i++ ) {
                  if( column[i] == AST__BAD || !isfinite(column[i]) ) valid = 0;
               }
               if( valid ) break;
            }
            mm = astAnnul( mm );
         }
      }

      if( mm ) {
/* Keep the shifts separate to avoid subtracting a large J*center from
   the function value when the domain is far from the origin. */
         for( i = 0; i < nin; i++ ) column[i] = -center[i];
         sm = (AstMapping *) astShiftMap( nin, column, "", status );
         tmp = (AstMapping *) astCmpMap( sm, mm, 1, "", status );
         sm = astAnnul( sm );
         sm = (AstMapping *) astShiftMap( nin, value, "", status );
         result = (AstMapping *) astCmpMap( tmp, sm, 1, "", status );
      } else {
         for( i = 0; i < nin; i++ ) inperm[i] = -i - 1;
         result = (AstMapping *) astPermMap( nin, inperm, nin, NULL,
                                            center, "", status );
      }
      if( astOK ) map->lintrunc = astClone( result );
   }
   if( mm ) mm = astAnnul( mm );
   if( sm ) sm = astAnnul( sm );
   if( tmp ) tmp = astAnnul( tmp );
   center = astFree( center );
   value = astFree( value );
   column = astFree( column );
   matrix = astFree( matrix );
   inperm = astFree( inperm );
   return result;
}

static size_t GetObjSize( AstObject *this_object, int *status ) {
/*
*  Name:
*     GetObjSize

*  Purpose:
*     Return the in-memory size of an Object.

*  Type:
*     Private function.

*  Synopsis:
*     #include "chebymap.h"
*     size_t GetObjSize( AstObject *this, int *status )

*  Class Membership:
*     ChebyMap member function (over-rides the astGetObjSize protected
*     method inherited from the parent class).

*  Description:
*     This function returns the in-memory size of the supplied ChebyMap,
*     in bytes.

*  Parameters:
*     this
*        Pointer to the ChebyMap.
*     status
*        Pointer to the inherited status variable.

*  Returned Value:
*     The Object size, in bytes.

*  Notes:
*     - A value of zero will be returned if this function is invoked
*     with the global status set, or if it should fail for any reason.
*/

/* Local Variables: */
   AstChebyMap *this;
   int nin;
   int nout;
   size_t result;

/* Initialise. */
   result = 0;

/* Check the global error status. */
   if ( !astOK ) return result;

/* Obtain a pointers to the ChebyMap structure. */
   this = (AstChebyMap *) this_object;

/* Get the number of input and output axes. */
   nin = astGetInvert( this ) ? astGetNout( this ) : astGetNin( this );
   nout = astGetInvert( this ) ? astGetNin( this ) : astGetNout( this );

/* Invoke the GetObjSize method inherited from the parent class, and then
   add on any components of the class structure defined by this class
   which are stored in dynamically allocated memory. */
   result = (*parent_getobjsize)( this_object, status );

   if( this->scale_f ) result += nin*sizeof( *this->scale_f );
   if( this->offset_f ) result += nin*sizeof( *this->offset_f );
   if( this->scale_i ) result += nout*sizeof( *this->scale_i );
   if( this->offset_i ) result += nout*sizeof( *this->offset_i );

/* If an error occurred, clear the result value. */
   if ( !astOK ) result = 0;

/* Return the result, */
   return result;
}

static int GetNiterInverse( AstPolyMap *this, int *status ) {
/*
*  Name:
*     GetNiterInverse

*  Purpose:
*     Return the value of the NiterInverse attribute.

*  Type:
*     Private function.

*  Synopsis:
*     #include "polymap.h"
*     int GetNiterInverse( AstPolyMap *this, int *status )

*  Class Membership:
*     ChebyMap member function (over-rides the astGetNiterInverse
*     protected method inherited from the PolyMap class).

*  Description:
*     This function returns the NiterInverse value. An explicitly set
*     value is returned unchanged. The default for a ChebyMap is ten,
*     because the bounded algorithm checks the final candidate after the
*     last update and returns AST__BAD on exhaustion, so it needs more
*     headroom than the unbounded PolyMap algorithm, whose default is
*     four.

*  Parameters:
*     this
*        Pointer to the PolyMap (in practice always a ChebyMap).
*     status
*        Pointer to the inherited status variable.

*  Returned Value:
*     The NiterInverse value to use.
*/
   if( !astOK ) return 0;
   if( astTestNiterInverse( this ) ) return (*parent_getniterinverse)( this, status );
   return 10;
}

void astInitChebyMapVtab_(  AstChebyMapVtab *vtab, const char *name, int *status ) {
/*
*+
*  Name:
*     astInitChebyMapVtab

*  Purpose:
*     Initialise a virtual function table for a ChebyMap.

*  Type:
*     Protected function.

*  Synopsis:
*     #include "chebymap.h"
*     void astInitChebyMapVtab( AstChebyMapVtab *vtab, const char *name )

*  Class Membership:
*     ChebyMap vtab initialiser.

*  Description:
*     This function initialises the component of a virtual function
*     table which is used by the ChebyMap class.

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
   AstPolyMapVtab *polymap;      /* Pointer to PolyMap component of Vtab */

/* Check the local error status. */
   if ( !astOK ) return;

/* Get a pointer to the thread specific global data structure. */
   astGET_GLOBALS(NULL);

/* Initialize the component of the virtual function table used by the
   parent class. */
   astInitPolyMapVtab( (AstPolyMapVtab *) vtab, name );

/* Store a unique "magic" value in the virtual function table. This
   will be used (by astIsAChebyMap) to determine if an object belongs
   to this class.  We can conveniently use the address of the (static)
   class_check variable to generate this unique value. */
   vtab->id.check = &class_check;
   vtab->id.parent = &(((AstPolyMapVtab *) vtab)->id);

/* Initialise member function pointers. */
/* ------------------------------------ */
/* Store pointers to the member functions (implemented here) that provide
   virtual methods for this class. */
/* none */

/* Save the inherited pointers to methods that will be extended, and
   replace them with pointers to the new member functions. */
   object = (AstObjectVtab *) vtab;
   polymap = (AstPolyMapVtab *) vtab;

   parent_getjacobian = polymap->GetJacobian;
   polymap->GetJacobian = GetJacobian;
   parent_linearguess = polymap->LinearGuess;
   polymap->LinearGuess = LinearGuess;
   parent_iterinverse = polymap->IterInverse;
   polymap->IterInverse = IterInverse;
   parent_getniterinverse = polymap->GetNiterInverse;
   polymap->GetNiterInverse = GetNiterInverse;

   parent_getobjsize = object->GetObjSize;
   object->GetObjSize = GetObjSize;

   parent_polypowers = polymap->PolyPowers;
   polymap->PolyPowers = PolyPowers;

   parent_polytran = polymap->PolyTran;
   polymap->PolyTran = PolyTran;

   parent_equal = object->Equal;
   object->Equal = Equal;

   polymap->FitPoly1DInit = FitPoly1DInit;
   polymap->FitPoly2DInit = FitPoly2DInit;

/* Store pointers to the member functions (implemented here) that
   provide virtual methods for this class. */
   vtab->ChebyDomain = ChebyDomain;

/* Declare the destructor and copy constructor. */
   astSetDelete( (AstObjectVtab *) vtab, Delete );
   astSetCopy( (AstObjectVtab *) vtab, Copy );

/* Declare the class dump function. */
   astSetDump( vtab, Dump, "ChebyMap", "Chebyshev polynomial transformation" );

/* If we have just initialised the vtab for the current class, indicate
   that the vtab is now initialised, and store a pointer to the class
   identifier in the base "object" level of the vtab. */
   if( vtab == &class_vtab ) {
      class_init = 1;
      astSetVtabClassIdentifier( vtab, &(vtab->id) );
   }
}

static void PolyPowers( AstPolyMap *this_polymap, double **work, int ncoord,
                        const int *mxpow, double **ptr, int point, int fwd,
                        int *status ){
/*
*  Name:
*     PolyPowers

*  Purpose:
*     Find the required powers of the input axis values.

*  Type:
*     Private function.

*  Synopsis:
*     #include "chebymap.h"
*     void PolyPowers( AstPolyMap *this, double **work, int ncoord,
*                      const int *mxpow, double **ptr, int point,
*                      int fwd, int *status )

*  Class Membership:
*     ChebyMap member function (over-rides the astPolyPowers protected
*     method inherited from the PolyMap class).

*  Description:
*     This function is used by astTransform to calculate the powers of
*     the axis values for a single input position. In the case of
*     sub-classes, the powers may not be simply powers of the supplied
*     axis values but may be more complex quantities such as a Chebyshev
*     polynomial of the required degree evaluated at the input axis values.

*  Parameters:
*     this
*        Pointer to the PolyMap.
*     work
*        An array of "ncoord" pointers, each pointing to an array of
*        length "max(2,mxpow)". The required values are placed in this
*        array on exit.
*     ncoord
*        The number of axes.
*     mxpow
*        Pointer to an array holding the maximum power required of each
*        axis value. Should have "ncoord" elements.
*     ptr
*        An array of "ncoord" pointers, each pointing to an array holding
*        the axis values. Each of these arrays of axis values must have
*        at least "point+1" elements.
*     point
*        The zero based index of the point within "ptr" that holds the
*        axis values to be exponentiated.
*     fwd
*        Do the supplied coefficients define the foward transformation of
*        the PolyMap?
*/

/* Local Variables; */
   AstChebyMap *this;
   double *scales;
   double *offsets;
   double *pwork;
   double *t;
   double x;
   int coord;
   int ip;

/* Check the local error status. */
   if ( !astOK ) return;

/* Get a pointer to the ChebyMap structure. */
   this = (AstChebyMap *) this_polymap;

/* Either transformation of a ChebyMap (forward or inverse) can be
   defined either as Chebyshev polynomial or as a standard polynomial.
   Chebyshev polynomials always have non-NULL scale array pointers.
   If the scale array pointer is NULL, then the transformation is a
   standard polynomial. If the coefficients relate to a standard
   polynomial, then invoke the astPolyPowers implementation of the parent
   class (PolyMap). */
   if( (fwd && !this->scale_f) || (!fwd && !this->scale_i) ) {
      (*parent_polypowers)( this_polymap, work, ncoord, mxpow, ptr, point,
                            fwd, status );

/* If the coefficients relate to a Chebyshev polynomial... */
   } else {
      scales = fwd ? this->scale_f : this->scale_i;
      offsets = fwd ? this->offset_f : this->offset_i;

/* This method uses a Chebyshev polynomial of the first kind of degree "i"
   evaluated at "x'" instead of "x raised to the power i". Here, "x'" is
   the input axis value scaled and shifted into the range [-1,+1] on each
   axis. Loop over all input axes. */
      for( coord = 0; coord < ncoord; coord++ ) {

/* Get a pointer to the array in which the powers of the current axis
   value are to be returned. */
         pwork = work[ coord ];

/* The Chebyshev function (type 1) of degree zero is always 1.0, regardless
   of the value of x. */
         pwork[ 0 ] = 1.0;

/* Get the input axis value. If it is bad, store bad values for all
   remaining powers. */
         x = ptr[ coord ][ point ];
         if( x == AST__BAD ) {
            for( ip = 1; ip <= mxpow[ coord ]; ip++ ) pwork[ ip ] = AST__BAD;

/* Otherwise, apply the required scaling to the input */
         } else {
            x = x*scales[ coord ] + offsets[ coord ];

/* Return bad values for input positions outside the bounding box
   associated with the transformation. */
            if( fabs( x ) <= 1.0 ) {

/* The Chebyshev function of degree one is equal to x. */
               t = pwork + 1;
               *t = x;

/* Form and store the remaining Chebyshev polynomial values at the input axis value.
   Use the standard recurrence relation: Tn+1(x') = 2.x'.Tn(x') - Tn-1(x'). */
               for( ip = 2; ip <= mxpow[ coord ]; ip++,t++ ) {
                  t[ 1 ] = 2.0*x*t[ 0 ] - t[ -1 ];
               }
            } else {
               for( ip = 1; ip <= mxpow[ coord ]; ip++ ) pwork[ ip ] = AST__BAD;
            }
         }
      }
   }
}

static AstPolyMap *PolyTran( AstPolyMap *this_polymap, int forward, double acc,
                             double maxacc, int maxorder, const double *lbnd,
                             const double *ubnd, int *status ){
/*
*  Name:
*     PolyTran

*  Purpose:
*     Fit a PolyMap inverse or forward transformation.

*  Type:
*     Private function.

*  Synopsis:
*     #include "polymap.h"
*     AstPolyMap *PolyTran( AstPolyMap *this, int forward, double acc,
*                           double maxacc, int maxorder, const double *lbnd,
*                           const double *ubnd )

*  Class Membership:
*     ChebyMap member function (over-rides the astPolyTran method inherited
*     from the PolyMap class).

*  Description:
*     This function creates a new PolyMap which is a copy of the supplied
*     PolyMap, in which a specified transformation (forward or inverse)
*     has been replaced by a new polynomial transformation. The
*     coefficients of the new transformation are estimated by sampling
*     the other transformation and performing a least squares polynomial
*     fit in the opposite direction to the sampled positions and values.
*
*     This method can only be used on (1-input,1-output) or (2-input,2-output)
*     PolyMaps.
*
*     The transformation to create is specified by the "forward" parameter.
*     In what follows "X" refers to the inputs of the PolyMap, and "Y" to
*     the outputs of the PolyMap. The forward transformation transforms
*     input values (X) into output values (Y), and the inverse transformation
*     transforms output values (Y) into input values (X). Within a PolyMap,
*     each transformation is represented by an independent set of
*     polynomials, P_f or P_i: Y=P_f(X) for the forward transformation and
*     X=P_i(Y) for the inverse transformation.
*
*     The "forward" parameter specifies the transformation to be replaced. If
*     it is non-zero, a new forward transformation is created by first finding
*     the input values (X) using the inverse transformation
*     (which must be available) at a regular grid of points (Y) covering a
*     rectangular region of the PolyMap's output space. The coefficients of
*     the required forward polynomial, Y=P_f(X), are chosen in order to
*     minimise the sum of the squared residuals between the sampled values
*     of Y and P_f(X).
*
*     If "forward" is zero (probably the most likely case),
*     a new inverse transformation is created by
*     first finding the output values (Y) using the forward transformation
*     (which must be available) at a regular grid of points (X) covering a
*     rectangular region of the PolyMap's input space. The coefficients of
*     the required inverse polynomial, X=P_i(Y), are chosen in order to
*     minimise the sum of the squared residuals between the sampled values
*     of X and P_i(Y).
*
*     This fitting process is performed repeatedly with increasing
*     polynomial orders (starting with linear) until the target
*     accuracy is achieved, or a specified maximum order is reached. If
*     the target accuracy cannot be achieved even with this maximum-order
*     polynomial, the best fitting maximum-order polynomial is returned so
*     long as its accuracy is better than "maxacc".
*     If it is not, an error is reported.

*  Parameters:
*     this
*        Pointer to the original Mapping.
*     forward
*        If non-zero, the forward PolyMap transformation is replaced.
*        Otherwise the inverse transformation is replaced.
*     acc
*        The target accuracy, expressed as a geodesic distance within
*        the PolyMap's input space (if "forward" is zero)  or output
*        space (if "forward" is non-zero).
*     maxacc
*        The maximum allowed accuracy for an acceptable polynomial,
*        expressed as a geodesic distance within the PolyMap's input
*        space (if "forward" is zero)  or output space (if "forward" is
*        non-zero).
*     maxorder
*        The maximum allowed polynomial order. This is one more than the
*        maximum power of either input axis. So for instance, a value of
*        3 refers to a quadratic polynomial. Note, cross terms with total
*        powers greater than or equal to maxorder are not inlcuded in the
*        fit. So the maximum number of terms in each of the fitted
*        polynomials is maxorder*(maxorder+1)/2.
*     lbnd
*        Pointer to an array holding the lower bounds of a rectangular
*        region within the PolyMap's input space (if "forward" is zero)
*        or output space (if "forward" is non-zero). The new polynomial
*        will be evaluated over this rectangle. The length of this array
*        should equal the value of the PolyMap's Nin or Nout attribute,
*        depending on "forward". If a NULL pointer is supplied, the lower
*        bounds of the box supplied when the ChebyMap was constructed is
*        used.
*     ubnd
*        Pointer to an
*        array holding the upper bounds of a rectangular region within
*        the PolyMap's input space (if "forward" is zero)  or output space
*        (if "forward" is non-zero). The new polynomial will be evaluated
*        over this rectangle.  The length of this array should equal the
*        value of the PolyMap's Nin or Nout attribute, depending on "forward".
*        If a NULL pointer is supplied, the upper bounds of the box supplied
*        when the ChebyMap was constructed is used.

*  Returned Value:
*     astPolyTran()
*        A pointer to the new PolyMap. A NULL pointer will be returned if
*        the fit fails to achieve the accuracy specified by "maxacc", but
*        no error will be reported.

*  Notes:
*     - This function can only be used on 1D or 2D PolyMaps which have
*     the same number of inputs and outputs.
*     - A null Object pointer (AST__NULL) will be returned if this
*     function is invoked with the AST error status set, or if it
*     should fail for any reason.
*/

/* Local Variables: */
   AstChebyMap *this;
   AstPolyMap *result;
   const char *word;
   double *offset;
   double *scale;
   double this_lbnd[ 2 ];
   double this_ubnd[ 2 ];
   int inverted;
   int k;
   int nax;

/* Initialise. */
   result = NULL;

/* Check the inherited status. */
   if ( !astOK ) return result;

/* Get a pointer to the CHebyMap structure. */
   this = (AstChebyMap *) this_polymap;

/* Select the ChebyMap scales and offsets to be used. */
   inverted = astGetInvert( this );
   if( ( inverted && !forward ) || ( !inverted && forward ) ) {
      word = "inverse";
      scale = this->scale_i;
      offset = this->offset_i;
      nax = ((AstMapping *)this)->nout;
   } else {
      word = "forward";
      scale = this->scale_f;
      offset = this->offset_f;
      nax = ((AstMapping *)this)->nin;
   }

/* The scaled box for a Chebyshev polynomial spans [-1,+1] on each axis.
   Create the corresponding unscaled box. If the user supplies both bounds
   arrays, use them in preference to the bounds in the ChebyMap. Otherwise
   reconstruct the missing bound(s) from the ChebyMap's own normalization,
   through the same AxisBounds helper used elsewhere in this file, so that
   a bound used here always evaluates without BAD. */
   if( lbnd && ubnd ) {
      for( k = 0; k < nax; k++ ) {
         this_lbnd[ k ] = lbnd[ k ];
         this_ubnd[ k ] = ubnd[ k ];
      }

   } else if( scale && offset ) {
      for( k = 0; k < nax; k++ ) {
         if( !AxisBounds( scale[ k ], offset[ k ], this_lbnd + k,
                          this_ubnd + k ) && astOK ) {
            astError( AST__NOBOX, "astPolyTran(%s): The %s transformation "
                      "has no usable bounding box on axis %d.", status,
                      astGetClass( this ), word, k + 1 );
         }
      }
      if( lbnd ) {
         for( k = 0; k < nax; k++ ) this_lbnd[ k ] = lbnd[ k ];
      }
      if( ubnd ) {
         for( k = 0; k < nax; k++ ) this_ubnd[ k ] = ubnd[ k ];
      }

   } else {
      if( !lbnd && astOK ) {
         astError( AST__NOBOX, "astPolyTran(%s): The %s transformation is "
                   "not a Chebyshev polynomial and therefore requires a "
                   "user-supplied bounding box. But no lower bounds were "
                   "supplied. ", status, astGetClass( this ), word );
      } else if( lbnd ) {
         for( k = 0; k < nax; k++ ) this_lbnd[ k ] = lbnd[ k ];
      }
      if( !ubnd && astOK ) {
         astError( AST__NOBOX, "astPolyTran(%s): The %s transformation is "
                   "not a Chebyshev polynomial and therefore requires a "
                   "user-supplied bounding box. But no upper bounds were "
                   "supplied. ", status, astGetClass( this ), word );
      } else if( ubnd ) {
         for( k = 0; k < nax; k++ ) this_ubnd[ k ] = ubnd[ k ];
      }
   }

/* Invoke the parent astPolyMap method, using the bounding box selected
   above. */
   result = (*parent_polytran)( this_polymap, forward, acc, maxacc, maxorder,
                                this_lbnd, this_ubnd, status );

/* Return the new ChebyMap. */
   return result;
}


/* Functions which access class attributes. */
/* ---------------------------------------- */
/* Implement member functions to access the attributes associated with
   this class using the macros defined for this purpose in the
   "object.h" file. For a description of each attribute, see the class
   interface (in the associated .h file). */

/* Copy constructor. */
/* ----------------- */
static void Copy( const AstObject *objin, AstObject *objout, int *status ) {
/*
*  Name:
*     Copy

*  Purpose:
*     Copy constructor for ChebyMap objects.

*  Type:
*     Private function.

*  Synopsis:
*     void Copy( const AstObject *objin, AstObject *objout, int *status )

*  Description:
*     This function implements the copy constructor for ChebyMap objects.

*  Parameters:
*     objin
*        Pointer to the object to be copied.
*     objout
*        Pointer to the object being constructed.
*     status
*        Pointer to the inherited status variable.

*  Returned Value:
*     void

*  Notes:
*     -  This constructor makes a deep copy, including a copy of the
*     coefficients associated with the input ChebyMap.
*/


/* Local Variables: */
   AstChebyMap *in;               /* Pointer to input ChebyMap */
   AstChebyMap *out;              /* Pointer to output ChebyMap */
   int nin;                       /* No. of input coordinates */
   int nout;                      /* No. of output coordinates */

/* Check the global error status. */
   if ( !astOK ) return;

/* Obtain pointers to the input and output ChebyMaps. */
   in = (AstChebyMap *) objin;
   out = (AstChebyMap *) objout;

/* Nullify the pointers stored in the output object since these will
   currently be pointing at the input data (since the output is a simple
   byte-for-byte copy of the input). Otherwise, the input data could be
   freed by accidient if the output object is deleted due to an error
   occuring in this function. */
   out->scale_f = NULL;
   out->offset_f = NULL;
   out->scale_i = NULL;
   out->offset_i = NULL;

/* Get the number of inputs and outputs of the uninverted Mapping. */
   nin = ( (AstMapping *) in )->nin;
   nout = ( (AstMapping *) in )->nout;

/* Copy the bounding box arrays. */
   if( in->scale_f ) out->scale_f = (double *) astStore( NULL,
                                       (void *) in->scale_f,
                                       sizeof( double )*nin );
   if( in->offset_f ) out->offset_f = (double *) astStore( NULL,
                                       (void *) in->offset_f,
                                       sizeof( double )*nin );
   if( in->scale_i ) out->scale_i = (double *) astStore( NULL,
                                       (void *) in->scale_i,
                                       sizeof( double )*nout );
   if( in->offset_i ) out->offset_i = (double *) astStore( NULL,
                                       (void *) in->offset_i,
                                       sizeof( double )*nout );
}

/* Destructor. */
/* ----------- */
static void Delete( AstObject *obj, int *status ) {
/*
*  Name:
*     Delete

*  Purpose:
*     Destructor for ChebyMap objects.

*  Type:
*     Private function.

*  Synopsis:
*     void Delete( AstObject *obj, int *status )

*  Description:
*     This function implements the destructor for ChebyMap objects.

*  Parameters:
*     obj
*        Pointer to the object to be deleted.
*     status
*        Pointer to the inherited status variable.

*  Returned Value:
*     void

*  Notes:
*     This function attempts to execute even if the global error status is
*     set.
*/

/* Local Variables: */
   AstChebyMap *this;

/* Obtain a pointer to the ChebyMap structure. */
   this = (AstChebyMap *) obj;

/* Free the boundib box arrays. */
   this->scale_f = astFree( this->scale_f );
   this->offset_f = astFree( this->offset_f );
   this->scale_i = astFree( this->scale_i );
   this->offset_i = astFree( this->offset_i );
}

/* Dump function. */
/* -------------- */
static void Dump( AstObject *this_object, AstChannel *channel, int *status ) {
/*
*  Name:
*     Dump

*  Purpose:
*     Dump function for ChebyMap objects.

*  Type:
*     Private function.

*  Synopsis:
*     void Dump( AstObject *this, AstChannel *channel, int *status )

*  Description:
*     This function implements the Dump function which writes out data
*     for the ChebyMap class to an output Channel.

*  Parameters:
*     this
*        Pointer to the ChebyMap whose data are being written.
*     channel
*        Pointer to the Channel to which the data are being written.
*     status
*        Pointer to the inherited status variable.
*/

#define KEY_LEN 50               /* Maximum length of a keyword */

/* Local Variables: */
   AstChebyMap *this;             /* Pointer to the ChebyMap structure */
   char buff[ KEY_LEN + 1 ];     /* Buffer for keyword string */
   char comm[ 100 ];             /* Buffer for comment string */
   int i;                        /* Loop index */
   int nin;                      /* No. of input coords */
   int nout;                     /* No. of output coords */

/* Check the global error status. */
   if ( !astOK ) return;

/* Obtain a pointer to the ChebyMap structure. */
   this = (AstChebyMap *) this_object;

/* Find the number of inputs and outputs of the uninverted Mapping. */
   nin = ( (AstMapping *) this )->nin;
   nout = ( (AstMapping *) this )->nout;

/* Write out values representing the instance variables for the
   ChebyMap class.  */

/* The input axis scale factors. */
   if( this->scale_f ){
      for( i = 0; i < nin; i++ ){
         (void) sprintf( buff, "FSCL%d", i + 1 );
         (void) sprintf( comm, "Scale factor on input %d", i + 1 );
         astWriteDouble( channel, buff, 1, 1, (this->scale_f)[ i ], comm );
      }
   }

/* The input axis offsets. */
   if( this->offset_f ){
      for( i = 0; i < nin; i++ ){
         (void) sprintf( buff, "FOFF%d", i + 1 );
         (void) sprintf( comm, "Offset on input %d", i + 1 );
         astWriteDouble( channel, buff, 1, 1, (this->offset_f)[ i ], comm );
      }
   }

/* The output axis scale factors. */
   if( this->scale_i ){
      for( i = 0; i < nout; i++ ){
         (void) sprintf( buff, "ISCL%d", i + 1 );
         (void) sprintf( comm, "Scale factor on output %d", i + 1 );
         astWriteDouble( channel, buff, 1, 1, (this->scale_i)[ i ], comm );
      }
   }

/* The output axis offsets. */
   if( this->offset_i ){
      for( i = 0; i < nout; i++ ){
         (void) sprintf( buff, "IOFF%d", i + 1 );
         (void) sprintf( comm, "Offset on output %d", i + 1 );
         astWriteDouble( channel, buff, 1, 1, (this->offset_i)[ i ], comm );
      }
   }

/* Undefine macros local to this function. */
#undef KEY_LEN
}

/* Standard class functions. */
/* ========================= */
/* Implement the astIsAChebyMap and astCheckChebyMap functions using the macros
   defined for this purpose in the "object.h" header file. */
astMAKE_ISA(ChebyMap,Mapping)
astMAKE_CHECK(ChebyMap)

AstChebyMap *astChebyMap_( int nin, int nout, int ncoeff_f, const double coeff_f[],
                           int ncoeff_i, const double coeff_i[],
                           const double lbnd_f[], const double ubnd_f[],
                           const double lbnd_i[], const double ubnd_i[],
                           const char *options, int *status, ...){
/*
*++
*  Name:
c     astChebyMap
f     AST_CHEBYMAP

*  Purpose:
*     Create a ChebyMap.

*  Type:
*     Public function.

*  Synopsis:
c     #include "chebymap.h"
c     AstChebyMap *astChebyMap( int nin, int nout, int ncoeff_f, const double coeff_f[],
c                               int ncoeff_i, const double coeff_i[],
c                               const double lbnd_f[], const double ubnd_f[],
c                               const double lbnd_i[], const double ubnd_i[],
c                               const char *options, ... )
f     RESULT = AST_CHEBYMAP( NIN, NOUT, NCOEFF_F, COEFF_F, NCOEFF_I, COEFF_I,
f                            LBND_F, UBND_F, LBND_I, UBND_I, OPTIONS, STATUS )

*  Class Membership:
*     ChebyMap constructor.

*  Description:
*     This function creates a new ChebyMap and optionally initialises
*     its attributes.
*
*     A ChebyMap is a form of Mapping which performs a Chebyshev polynomial
*     transformation.  Each output coordinate is a linear combination of
*     Chebyshev polynomials of the first kind, of order zero up to a
*     specified maximum order, evaluated at the input coordinates. The
*     coefficients to be used in the linear combination are specified
*     separately for each output coordinate.
*
*     For a 1-dimensional ChebyMap, the forward transformation is defined
*     as follows:
*
*        f(x) = c0.T0(x') + c1.T1(x') + c2.T2(x') + ...
*
*     where:
*        - Tn(x') is the nth Chebyshev polynomial of the first kind:
*             - T0(x') = 1
*             - T1(x') = x'
*             - Tn+1(x') = 2.x'.Tn(x') - Tn-1(x')
*        - x' is the input axis value, x, offset and scaled to the range
*          [-1, 1] as x ranges over a specified bounding box, given when the
*          ChebyMap is created. The input positions, x,  supplied to the
*          forward transformation must fall within the bounding box - bad
*          axis values (AST__BAD) are generated for points outside the
*          bounding box.
*
*     For an N-dimensional ChebyMap, the forward transformation is a
*     generalisation of the above form. Each output axis value is the sum
c     of "ncoeff"
f     of NCOEFF
*     terms, where each term is the product of a single coefficient
*     value and N factors of the form Tn(x'_i), where "x'_i" is the
*     normalised value of the i'th input axis value.
*
*     The forward and inverse transformations may be defined independently
*     by separate sets of coefficients supplied when the ChebyMap is
*     created. If forward coefficients are supplied, no inverse coefficients
*     are supplied, and the numbers of inputs and outputs are equal, an
*     iterative inverse is provided by default. It uses the analytic
*     Jacobian of the forward series and confines candidate solutions to
*     the forward bounding box. An unsolved position is returned as
*     AST__BAD. See IterInverse, NiterInverse and TolInverse for details.
*
*     Supplied inverse coefficients are used by default. Setting
*     IterInverse to one selects iteration instead; setting it to zero
*     disables iteration. Clearing it restores the default selection.
*     A local iterative inverse does not guarantee a unique solution or
*     convergence at every position.
*
*     Alternatively, the
c     astPolyTran
f     AST_POLYTRAN
*     method can fit an inverse Chebyshev series, choosing coefficients
*     to minimise the residuals of a forward/inverse round trip.

*  Parameters:
c     nin
f     NIN = INTEGER (Given)
*        The number of input coordinates.
c     nout
f     NOUT = INTEGER (Given)
*        The number of output coordinates.
c     ncoeff_f
f     NCOEFF_F = INTEGER (Given)
*        The number of non-zero coefficients necessary to define the
*        forward transformation of the ChebyMap. If zero is supplied, the
*        forward transformation will be undefined.
c     coeff_f
f     COEFF_F( * ) = DOUBLE PRECISION (Given)
*        An array containing
c        "ncoeff_f*( 2 + nin )" elements. Each group of "2 + nin"
f        "NCOEFF_F*( 2 + NIN )" elements. Each group of "2 + NIN"
*        adjacent elements describe a single coefficient of the forward
*        transformation. Within each such group, the first element is the
*        coefficient value; the next element is the integer index of the
*        ChebyMap output which uses the coefficient within its defining
*        expression (the first output has index 1); the remaining elements
*        of the group give the integer powers to use with each input
*        coordinate value (powers must not be negative, and floating
*        point values are rounded to the nearest integer).
c        If "ncoeff_f" is zero, a NULL pointer may be supplied for "coeff_f".
*
*        For instance, if the ChebyMap has 3 inputs and 2 outputs, each group
*        consisting of 5 elements, A groups such as "(1.2, 2.0, 1.0, 3.0, 0.0)"
*        describes a coefficient with value 1.2 which is used within the
*        definition of output 2. The output value is incremented by the
*        product of the coefficient value, the value of the Chebyshev
*        polynomial of power 1 evaluated at input coordinate 1, and the
*        value of the Chebyshev polynomial of power 3 evaluated at input
*        coordinate 2. Input coordinate 3 is not used since its power is
*        specified as zero. As another example, the group "(-1.0, 1.0,
*        0.0, 0.0, 0.0 )" adds a constant value -1.0 onto output 1 (it is
*        a constant value since the power for every input axis is given as
*        zero).
*
c        Each final output coordinate value is the sum of the "ncoeff_f" terms
c        described by the "ncoeff_f" groups within the supplied array.
f        Each final output coordinate value is the sum of the "NCOEFF_F" terms
f        described by the "NCOEFF_F" groups within the supplied array.
c     ncoeff_i
f     NCOEFF_I = INTEGER (Given)
*        The number of non-zero coefficients necessary to define the
*        inverse transformation of the ChebyMap. If zero is supplied,
*        an iterative inverse is provided when forward coefficients exist
*        and the numbers of inputs and outputs are equal (see IterInverse).
c     coeff_i
f     COEFF_I( * ) = DOUBLE PRECISION (Given)
*        An array containing
c        "ncoeff_i*( 2 + nout )" elements. Each group of "2 + nout"
f        "NCOEFF_I*( 2 + NOUT )" elements. Each group of "2 + NOUT"
*        adjacent elements describe a single coefficient of the inverse
c        transformation, using the same schame as "coeff_f",
f        transformation, using the same schame as "COEFF_F",
*        except that "inputs" and "outputs" are transposed.
c        If "ncoeff_i" is zero, a NULL pointer may be supplied for "coeff_i".
c     lbnd_f
f     LBND_F( * ) = DOUBLE PRECISION (Given)
*        An array containing the lower bounds of the input bounding box within
*        which the ChebyMap is defined. This argument is not used or
*        accessed if
c        ncoeff_f is zero, and so a NULL pointer may be supplied.
f        NCOEFF_F is zero.
*        If supplied, the array should contain
c        "nin" elements.
f        "NIN" elements.
c     ubnd_f
f     UBND_F( * ) = DOUBLE PRECISION (Given)
*        An array containing the upper bounds of the input bounding box within
*        which the ChebyMap is defined. This argument is not used or
*        accessed if
c        ncoeff_f is zero, and so a NULL pointer may be supplied.
f        NCOEFF_F is zero.
*        If supplied, the array should contain
c        "nin" elements.
f        "NIN" elements.
c     lbnd_i
f     LBND_I( * ) = DOUBLE PRECISION (Given)
*        An array containing the lower bounds of the output bounding box within
*        which the ChebyMap is defined. This argument is not used or
*        accessed if
c        ncoeff_i is zero, and so a NULL pointer may be supplied.
f        NCOEFF_I is zero.
*        If supplied, the array should contain
c        "nout" elements.
f        "NOUT" elements.
c     ubnd_i
f     UBND_I( * ) = DOUBLE PRECISION (Given)
*        An array containing the upper bounds of the output bounding box within
*        which the ChebyMap is defined. This argument is not used or
*        accessed if
c        ncoeff_i is zero, and so a NULL pointer may be supplied.
f        NCOEFF_I is zero.
*        If supplied, the array should contain
c        "nout" elements.
f        "NOUT" elements.
c     options
f     OPTIONS = CHARACTER * ( * ) (Given)
c        Pointer to a null-terminated string containing an optional
c        comma-separated list of attribute assignments to be used for
c        initialising the new ChebyMap. The syntax used is identical to
c        that for the astSet function and may include "printf" format
c        specifiers identified by "%" symbols in the normal way.
f        A character string containing an optional comma-separated
f        list of attribute assignments to be used for initialising the
f        new ChebyMap. The syntax used is identical to that for the
f        AST_SET routine.
c     ...
c        If the "options" string contains "%" format specifiers, then
c        an optional list of additional arguments may follow it in
c        order to supply values to be substituted for these
c        specifiers. The rules for supplying these are identical to
c        those for the astSet function (and for the C "printf"
c        function).
f     STATUS = INTEGER (Given and Returned)
f        The global status.

*  Returned Value:
c     astChebyMap()
f     AST_CHEBYMAP = INTEGER
*        A pointer to the new ChebyMap.

*  Notes:
*     - A null Object pointer (AST__NULL) will be returned if this
c     function is invoked with the AST error status set, or if it
f     function is invoked with STATUS set to an error value, or if it
*     should fail for any reason.
*--
*/

/* Local Variables: */
   astDECLARE_GLOBALS          /* Pointer to thread-specific global data */
   AstChebyMap *new;            /* Pointer to new ChebyMap */
   va_list args;               /* Variable argument list */

/* Check the global status. */
   if ( !astOK ) return NULL;

/* Get a pointer to the thread specific global data structure. */
   astGET_GLOBALS(NULL);

/* Initialise the ChebyMap, allocating memory and initialising the
   virtual function table as well if necessary. */
   new = astInitChebyMap( NULL, sizeof( AstChebyMap ), !class_init,
                          &class_vtab, "ChebyMap", nin, nout,
                          ncoeff_f, coeff_f, ncoeff_i, coeff_i,
                          lbnd_f, ubnd_f, lbnd_i, ubnd_i );

/* If successful, note that the virtual function table has been
   initialised. */
   if ( astOK ) {
      class_init = 1;

/* Obtain the variable argument list and pass it along with the options string
   to the astVSet method to initialise the new ChebyMap's attributes. */
      va_start( args, status );
      astVSet( new, options, NULL, args );
      va_end( args );

/* If an error occurred, clean up by deleting the new object. */
      if ( !astOK ) new = astDelete( new );
   }

/* Return a pointer to the new ChebyMap. */
   return new;
}

AstChebyMap *astChebyMapId_( int nin, int nout, int ncoeff_f, const double coeff_f[],
                             int ncoeff_i, const double coeff_i[],
                             const double lbnd_f[], const double ubnd_f[],
                             const double lbnd_i[], const double ubnd_i[],
                             const char *options, ... ){
/*
*  Name:
*     astChebyMapId_

*  Purpose:
*     Create a ChebyMap.

*  Type:
*     Private function.

*  Synopsis:
*     #include "chebymap.h"
*     AstChebyMap *astChebyMap( int nin, int nout, int ncoeff_f, const double coeff_f[],
*                             int ncoeff_i, const double coeff_i[], const
*                             double lbnd_f[], const double ubnd_f[],
*                             double lbnd_i[], const double ubnd_i[],
*                             const char *options, ... )

*  Class Membership:
*     ChebyMap constructor.

*  Description:
*     This function implements the external (public) interface to the
*     astChebyMap constructor function. It returns an ID value (instead
*     of a true C pointer) to external users, and must be provided
*     because astChebyMap_ has a variable argument list which cannot be
*     encapsulated in a macro (where this conversion would otherwise
*     occur).
*
*     The variable argument list also prevents this function from
*     invoking astChebyMap_ directly, so it must be a re-implementation
*     of it in all respects, except for the final conversion of the
*     result to an ID value.

*  Parameters:
*     As for astChebyMap_.

*  Returned Value:
*     The ID value associated with the new ChebyMap.
*/

/* Local Variables: */
   astDECLARE_GLOBALS            /* Pointer to thread-specific global data */
   AstChebyMap *new;              /* Pointer to new ChebyMap */
   va_list args;                 /* Variable argument list */
   int *status;                  /* Pointer to inherited status value */

/* Get a pointer to the inherited status value. */
   status = astGetStatusPtr;

/* Get a pointer to the thread specific global data structure. */
   astGET_GLOBALS(NULL);

/* Check the global status. */
   if ( !astOK ) return NULL;

/* Initialise the ChebyMap, allocating memory and initialising the
   virtual function table as well if necessary. */
   new = astInitChebyMap( NULL, sizeof( AstChebyMap ), !class_init,
                         &class_vtab, "ChebyMap", nin, nout,
                         ncoeff_f, coeff_f, ncoeff_i, coeff_i,
                         lbnd_f, ubnd_f, lbnd_i, ubnd_i );

/* If successful, note that the virtual function table has been
   initialised. */
   if ( astOK ) {
      class_init = 1;

/* Obtain the variable argument list and pass it along with the options string
   to the astVSet method to initialise the new ChebyMap's attributes. */
      va_start( args, options );
      astVSet( new, options, NULL, args );
      va_end( args );

/* If an error occurred, clean up by deleting the new object. */
      if ( !astOK ) new = astDelete( new );
   }

/* Return an ID value for the new ChebyMap. */
   return astMakeId( new );
}

AstChebyMap *astInitChebyMap_( void *mem, size_t size, int init,
                             AstChebyMapVtab *vtab, const char *name,
                             int nin, int nout, int ncoeff_f, const double coeff_f[],
                             int ncoeff_i, const double coeff_i[],
                             const double lbnd_f[], const double ubnd_f[],
                             const double lbnd_i[], const double ubnd_i[],
                             int *status ){
/*
*+
*  Name:
*     astInitChebyMap

*  Purpose:
*     Initialise a ChebyMap.

*  Type:
*     Protected function.

*  Synopsis:
*     #include "chebymap.h"
*     AstChebyMap *astInitChebyMap( void *mem, size_t size, int init,
*                                 AstChebyMapVtab *vtab, const char *name,
*                                 int nin, int nout, int ncoeff_f,
*                                 const double coeff_f[], int ncoeff_i,
*                                 const double coeff_i[]
*                                 const double lbnd_f[], const double ubnd_f[],
*                                 const double lbnd_i[], const double ubnd_i[] )

*  Class Membership:
*     ChebyMap initialiser.

*  Description:
*     This function is provided for use by class implementations to initialise
*     a new ChebyMap object. It allocates memory (if necessary) to accommodate
*     the ChebyMap plus any additional data associated with the derived class.
*     It then initialises a ChebyMap structure at the start of this memory. If
*     the "init" flag is set, it also initialises the contents of a virtual
*     function table for a ChebyMap at the start of the memory passed via the
*     "vtab" parameter.

*  Parameters:
*     mem
*        A pointer to the memory in which the ChebyMap is to be initialised.
*        This must be of sufficient size to accommodate the ChebyMap data
*        (sizeof(ChebyMap)) plus any data used by the derived class. If a value
*        of NULL is given, this function will allocate the memory itself using
*        the "size" parameter to determine its size.
*     size
*        The amount of memory used by the ChebyMap (plus derived class data).
*        This will be used to allocate memory if a value of NULL is given for
*        the "mem" parameter. This value is also stored in the ChebyMap
*        structure, so a valid value must be supplied even if not required for
*        allocating memory.
*     init
*        A logical flag indicating if the ChebyMap's virtual function table is
*        to be initialised. If this value is non-zero, the virtual function
*        table will be initialised by this function.
*     vtab
*        Pointer to the start of the virtual function table to be associated
*        with the new ChebyMap.
*     name
*        Pointer to a constant null-terminated character string which contains
*        the name of the class to which the new object belongs (it is this
*        pointer value that will subsequently be returned by the astGetClass
*        method).
*     nin
*        The number of input coordinate values per point. This is the
*        same as the number of columns in the matrix.
*     nout
*        The number of output coordinate values per point. This is the
*        same as the number of rows in the matrix.
*     ncoeff_f
*        The number of non-zero coefficients necessary to define the
*        forward transformation of the ChebyMap. If zero is supplied, the
*        forward transformation will be undefined.
*     coeff_f
*        An array containing "ncoeff_f*( 2 + nin )" elements. Each group
*	 of "2 + nin" adjacent elements describe a single coefficient of
*	 the forward transformation. Within each such group, the first
*	 element is the coefficient value; the next element is the
*	 integer index of the ChebyMap output which uses the coefficient
*	 within its defining polynomial (the first output has index 1);
*	 the remaining elements of the group give the integer powers to
*	 use with each input coordinate value (powers must not be
*	 negative)
*
*        For instance, if the ChebyMap has 3 inputs and 2 outputs, each group
*        consisting of 5 elements, A groups such as "(1.2, 2.0, 1.0, 3.0, 0.0)"
*        describes a coefficient with value 1.2 which is used within the
*        definition of output 2. The output value is incremented by the
*        product of the coefficient value, the value of input coordinate
*        1 raised to the power 1, and the value of input coordinate 2 raised
*        to the power 3. Input coordinate 3 is not used since its power is
*        specified as zero. As another example, the group "(-1.0, 1.0,
*        0.0, 0.0, 0.0 )" describes adds a constant value -1.0 onto
*        output 1 (it is a constant value since the power for every input
*        axis is given as zero).
*
*        Each final output coordinate value is the sum of the "ncoeff_f" terms
*        described by the "ncoeff_f" groups within the supplied array.
*     ncoeff_i
*        The number of non-zero coefficients necessary to define the
*        inverse transformation of the ChebyMap. If zero is supplied,
*        an iterative inverse is provided when forward coefficients exist
*        and the numbers of inputs and outputs are equal (see IterInverse).
*     coeff_i
*        An array containing
*        "ncoeff_i*( 2 + nout )" elements. Each group of "2 + nout"
*        adjacent elements describe a single coefficient of the inverse
*        transformation, using the same schame as "coeff_f", except that
*        "inputs" and "outputs" are transposed.
*     lbnd_f
*        An array containing the lower bounds of the input bounding box within
*        which the ChebyMap is defined. The array should contain "nin" elements.
*     ubnd_f
*        An array containing the upper bounds of the input bounding box within
*        which the ChebyMap is defined. The array should contain "nin" elements.
*     lbnd_i
*        An array containing the lower bounds of the output bounding box within
*        which the ChebyMap is defined. The array should contain "nout" elements.
*     ubnd_i
*        An array containing the upper bounds of the output bounding box within
*        which the ChebyMap is defined. The array should contain "nout" elements.

*  Returned Value:
*     A pointer to the new ChebyMap.

*  Notes:
*     -  A null pointer will be returned if this function is invoked with the
*     global error status set, or if it should fail for any reason.
*-
*/

/* Local Variables: */
   AstChebyMap *new;
   int i;

/* Check the global status. */
   if ( !astOK ) return NULL;

/* If necessary, initialise the virtual function table. */
   if ( init ) astInitChebyMapVtab( vtab, name );

/* Initialise a PolyMap structure (the parent class) as the first component
   within the ChebyMap structure, allocating memory if necessary. */
   new = (AstChebyMap *) astInitPolyMap( mem, size, 0,
                                        (AstPolyMapVtab *) vtab, name,
                                        nin, nout, ncoeff_f, coeff_f,
                                        ncoeff_i, coeff_i );
   if ( astOK ) {

/* Initialise the ChebyMap data. */
/* ---------------------------- */

/* First initialise the pointers in case of errors. */
      new->scale_f = NULL;
      new->offset_f = NULL;
      new->scale_i = NULL;
      new->offset_i = NULL;

/* A bounding box is needed for each direction that has coefficients. The
   checks below leave the status set, so the boxes are not read. */
      if( ncoeff_f > 0 && ( !lbnd_f || !ubnd_f ) ) {
         astError( AST__NOBOX, "astInitChebyMap(%s): No input bounding box "
                   "supplied, but the forward transformation is defined.",
                   status, name );

      } else if( ncoeff_i > 0 && ( !lbnd_i || !ubnd_i ) ) {
         astError( AST__NOBOX, "astInitChebyMap(%s): No output bounding box "
                   "supplied, but the inverse transformation is defined.",
                   status, name );
      }

/* Calculate the scales and offsets that map the supplied input bounding box
   onto the range [-1,+1] on each input axis, and store them. */
      if( ncoeff_f > 0 ) {
         new->scale_f = (double *) astMalloc( sizeof( double )*nin );
         new->offset_f = (double *) astMalloc( sizeof( double )*nin );
         if( astOK ) {
            for( i = 0; i < nin; i++ ) {
               if( ubnd_f[ i ] != lbnd_f[ i ] ) {
                  new->scale_f[ i ] = 2.0/( ubnd_f[ i ] - lbnd_f[ i ] );
                  new->offset_f[ i ] = -( ubnd_f[ i ] + lbnd_f[ i ] )/( ubnd_f[ i ] - lbnd_f[ i ] );
               } else if( astOK ){
                  astError( AST__BADBX, "astInitChebyMap(%s): Input bounding box "
                            "has zero width on input axis %d.", status, name, i + 1 );
                  break;
               }
            }
         }
      }

/* Calculate the scales and offsets that map the supplied output bounding box
   onto the range [-1,+1] on each output axis, and store them. */
      if( ncoeff_i > 0 ) {
         new->scale_i = (double *) astMalloc( sizeof( double )*nout );
         new->offset_i = (double *) astMalloc( sizeof( double )*nout );
         if( astOK ) {
            for( i = 0; i < nout; i++ ) {
               if( ubnd_i[ i ] != lbnd_i[ i ] ) {
                  new->scale_i[ i ] = 2.0/( ubnd_i[ i ] - lbnd_i[ i ] );
                  new->offset_i[ i ] = -( ubnd_i[ i ] + lbnd_i[ i ] )/( ubnd_i[ i ] - lbnd_i[ i ] );
               } else if( astOK ){
                  astError( AST__BADBX, "astInitChebyMap(%s): Output bounding box "
                            "has zero width on output axis %d.", status, name, i + 1 );
                  break;
               }
            }
         }
      }

/* If an error occurred, clean up by deleting the new ChebyMap. */
      if ( !astOK ) new = astDelete( new );
   }

/* Return a pointer to the new ChebyMap. */
   return new;
}

AstChebyMap *astLoadChebyMap_( void *mem, size_t size,
                               AstChebyMapVtab *vtab, const char *name,
                               AstChannel *channel, int *status ) {
/*
*+
*  Name:
*     astLoadChebyMap

*  Purpose:
*     Load a ChebyMap.

*  Type:
*     Protected function.

*  Synopsis:
*     #include "chebymap.h"
*     AstChebyMap *astLoadChebyMap( void *mem, size_t size,
*                                   AstChebyMapVtab *vtab, const char *name,
*                                   AstChannel *channel )

*  Class Membership:
*     ChebyMap loader.

*  Description:
*     This function is provided to load a new ChebyMap using data read
*     from a Channel. It first loads the data used by the parent class
*     (which allocates memory if necessary) and then initialises a
*     ChebyMap structure in this memory, using data read from the input
*     Channel.
*
*     If the "init" flag is set, it also initialises the contents of a
*     virtual function table for a ChebyMap at the start of the memory
*     passed via the "vtab" parameter.


*  Parameters:
*     mem
*        A pointer to the memory into which the ChebyMap is to be
*        loaded.  This must be of sufficient size to accommodate the
*        ChebyMap data (sizeof(ChebyMap)) plus any data used by derived
*        classes. If a value of NULL is given, this function will
*        allocate the memory itself using the "size" parameter to
*        determine its size.
*     size
*        The amount of memory used by the ChebyMap (plus derived class
*        data).  This will be used to allocate memory if a value of
*        NULL is given for the "mem" parameter. This value is also
*        stored in the ChebyMap structure, so a valid value must be
*        supplied even if not required for allocating memory.
*
*        If the "vtab" parameter is NULL, the "size" value is ignored
*        and sizeof(AstChebyMap) is used instead.
*     vtab
*        Pointer to the start of the virtual function table to be
*        associated with the new ChebyMap. If this is NULL, a pointer
*        to the (static) virtual function table for the ChebyMap class
*        is used instead.
*     name
*        Pointer to a constant null-terminated character string which
*        contains the name of the class to which the new object
*        belongs (it is this pointer value that will subsequently be
*        returned by the astGetClass method).
*
*        If the "vtab" parameter is NULL, the "name" value is ignored
*        and a pointer to the string "ChebyMap" is used instead.

*  Returned Value:
*     A pointer to the new ChebyMap.

*  Notes:
*     - A null pointer will be returned if this function is invoked
*     with the global error status set, or if it should fail for any
*     reason.
*-
*/

#define KEY_LEN 50               /* Maximum length of a keyword */

   astDECLARE_GLOBALS            /* Pointer to thread-specific global data */
/* Local Variables: */
   AstChebyMap *new;              /* Pointer to the new ChebyMap */
   char buff[ KEY_LEN + 1 ];     /* Buffer for keyword string */
   int i;                        /* Loop index */
   int ngood;                    /* No. of non-bad values */
   int nin;                      /* No. of input coords */
   int nout;                     /* No. of output coords */

/* Get a pointer to the thread specific global data structure. */
   astGET_GLOBALS(channel);

/* Initialise. */
   new = NULL;

/* Check the global error status. */
   if ( !astOK ) return new;

/* If a NULL virtual function table has been supplied, then this is
   the first loader to be invoked for this ChebyMap. In this case the
   ChebyMap belongs to this class, so supply appropriate values to be
   passed to the parent class loader (and its parent, etc.). */
   if ( !vtab ) {
      size = sizeof( AstChebyMap );
      vtab = &class_vtab;
      name = "ChebyMap";

/* If required, initialise the virtual function table for this class. */
      if ( !class_init ) {
         astInitChebyMapVtab( vtab, name );
         class_init = 1;
      }
   }

/* Invoke the parent class loader to load data for all the ancestral
   classes of the current one, returning a pointer to the resulting
   partly-built ChebyMap. */
   new = astLoadPolyMap( mem, size, (AstPolyMapVtab *) vtab, name,
                         channel );

   if ( astOK ) {

/* Get the number of inputs and outputs for the uninverted Mapping. */
   nin = ( (AstMapping *) new )->nin;
   nout = ( (AstMapping *) new )->nout;

/* Initialise values */
   new->scale_f = NULL;
   new->offset_f = NULL;
   new->scale_i = NULL;
   new->offset_i = NULL;

/* Read input data. */
/* ================ */

/* Request the input Channel to read all the input data appropriate to
   this class into the internal "values list". */
      astReadClassData( channel, "ChebyMap" );

/* Is the forward transformation defined? */
      if( ((AstPolyMap *) new)->ncoeff_f ) {

/* Allocate memory to hold the scales and offsets for the input axes. */
         new->scale_f = astMalloc( sizeof( double )*(size_t) nin );
         new->offset_f = astMalloc( sizeof( double )*(size_t) nin );
         if( astOK ) {

/* Get the scale factors. */
            ngood = 0;
            for( i = 0; i < nin; i++ ){
               (void) sprintf( buff, "fscl%d", i + 1 );
               (new->scale_f)[ i ] = astReadDouble( channel, buff, AST__BAD );
               if( (new->scale_f)[ i ] != AST__BAD ) ngood++;
            }

/* Get the offsets of the bounding box. */
            for( i = 0; i < nin; i++ ){
               (void) sprintf( buff, "foff%d", i + 1 );
               (new->offset_f)[ i ] = astReadDouble( channel, buff, AST__BAD );
               if( (new->offset_f)[ i ] != AST__BAD ) ngood++;
            }

/* The scale and offset values should all be AST__BAD if the transformation
   is a standard polynomial. Anull the scale and offset arrays to
   indicate this. */
            if( ngood == 0 ) {
               new->scale_f = astFree( new->scale_f );
               new->offset_f = astFree( new->offset_f );

/* Otherwise, there should be no bad values. */
            } else if( ngood != 2*nin && astOK ) {
               astError( AST__OBJIN, "astLoadChebyMap: insufficient scale "
                         "and offset values for the forward transformation "
                         "in loaded ChebyMap.", status );
            }
         }
      }

/* Is the inverse transformation defined? */
      if( ((AstPolyMap *) new)->ncoeff_i ) {

/* Allocate memory to hold the scales and offsets for the output axes. */
         new->scale_i = astMalloc( sizeof( double )*(size_t) nout );
         new->offset_i = astMalloc( sizeof( double )*(size_t) nout );
         if( astOK ) {

/* Get the scale factors. */
            ngood = 0;
            for( i = 0; i < nout; i++ ){
               (void) sprintf( buff, "iscl%d", i + 1 );
               (new->scale_i)[ i ] = astReadDouble( channel, buff, AST__BAD );
               if( (new->scale_i)[ i ] != AST__BAD ) ngood++;
            }

/* Get the offsets of the bounding box. */
            for( i = 0; i < nout; i++ ){
               (void) sprintf( buff, "ioff%d", i + 1 );
               (new->offset_i)[ i ] = astReadDouble( channel, buff, AST__BAD );
               if( (new->offset_i)[ i ] != AST__BAD ) ngood++;
            }

/* The scale and offset values should all be AST__BAD if the transformation
   is a standard polynomial. Anull the scale and offset arrays to
   indicate this. */
            if( ngood == 0 ) {
               new->scale_i = astFree( new->scale_i );
               new->offset_i = astFree( new->offset_i );

/* Otherwise, there should be no bad values. */
            } else if( ngood != 2*nout && astOK ) {
               astError( AST__OBJIN, "astLoadChebyMap: insufficient scale "
                         "and offset values for the inverse transformation "
                         "in loaded ChebyMap.", status );
            }
         }
      }

/* If an error occurred, clean up by deleting the new ChebyMap. */
      if ( !astOK ) new = astDelete( new );
   }

/* Return the new ChebyMap pointer. */
   return new;

/* Undefine macros local to this function. */
#undef KEY_LEN
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


void astChebyDomain_( AstChebyMap *this, int forward, double *lbnd, double *ubnd, int *status ){
   if ( !astOK ) return;
   (**astMEMBER(this,ChebyMap,ChebyDomain))( this, forward, lbnd, ubnd, status );
}
