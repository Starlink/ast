/*
*class++
*  Name:
*     TranMap

*  Purpose:
*     Mapping with specified forward and inverse transformations.

*  Constructor Function:
c     astTranMap
f     AST_TRANMAP

*  Description:
*     A TranMap is a Mapping which combines the forward transformation of
*     a supplied Mapping with the inverse transformation of another
*     supplied Mapping, ignoring the un-used transformation in each
*     Mapping (indeed the un-used transformation need not exist).
*
*     When the forward transformation of the TranMap is referred to, the
*     transformation actually used is the forward transformation of the
*     first Mapping supplied when the TranMap was constructed. Likewise, 
*     when the inverse transformation of the TranMap is referred to, the
*     transformation actually used is the inverse transformation of the
*     second Mapping supplied when the TranMap was constructed.

*  Inheritance:
*     The TranMap class inherits from the Mapping class.

*  Attributes:
*     The TranMap class does not define any new attributes beyond those
*     which are applicable to all Mappings.

*  Functions:
c     The TranMap class does not define any new functions beyond those
f     The TranMap class does not define any new routines beyond those
*     which are applicable to all Mappings.

*  Copyright:
*     <COPYRIGHT_STATEMENT>

*  Authors:
*     DSB: David S. Berry (Starlink)

*  History:
*     10-FEB-2004 (DSB):
*        Original version.
*     19-JAN-2005 (DSB):
*        Fix memory leak.
*class--
*/

/* Module Macros. */
/* ============== */
/* Set the name of the class we are implementing. This indicates to
   the header files that define class interfaces that they should make
   "protected" symbols available. */
#define astCLASS TranMap

/* Include files. */
/* ============== */
/* Interface definitions. */
/* ---------------------- */
#include "error.h"               /* Error reporting facilities */
#include "memory.h"              /* Memory allocation facilities */
#include "object.h"              /* Base Object class */
#include "pointset.h"            /* Sets of points/coordinates */
#include "mapping.h"             /* Coordinate Mappings (parent class) */
#include "channel.h"             /* I/O channels */
#include "permmap.h"             /* Coordinate permutation Mappings */
#include "cmpmap.h"              /* Compound Mappings */
#include "unitmap.h"             /* Unit Mappings */
#include "tranmap.h"             /* Interface definition for this class */

/* Error code definitions. */
/* ----------------------- */
#include "ast_err.h"             /* AST error codes */

/* C header files. */
/* --------------- */
#include <stdarg.h>
#include <stddef.h>
#include <string.h>

/* Module Variables. */
/* ================= */
/* Define the class virtual function table and its initialisation flag
   as static variables. */
static AstTranMapVtab class_vtab; /* Virtual function table */
static int class_init = 0;       /* Virtual function table initialised? */

/* Pointers to parent class methods which are extended by this class. */
static AstPointSet *(* parent_transform)( AstMapping *, AstPointSet *, int, AstPointSet * );
static int *(* parent_mapsplit)( AstMapping *, int, int *, AstMapping ** );

/* External Interface Function Prototypes. */
/* ======================================= */
/* The following functions have public prototypes only (i.e. no
   protected prototypes), so we must provide local prototypes for use
   within this module. */
AstTranMap *astTranMapId_( void *, void *, const char *, ... );

/* Prototypes for Private Member Functions. */
/* ======================================== */
static AstPointSet *Transform( AstMapping *, AstPointSet *, int, AstPointSet * );
static double Rate( AstMapping *, double *, int, int );
static int *MapSplit( AstMapping *, int, int *, AstMapping ** );
static int MapMerge( AstMapping *, int, int, int *, AstMapping ***, int ** );
static void Copy( const AstObject *, AstObject * );
static void Delete( AstObject * );
static void Dump( AstObject *, AstChannel * );

/* Member functions. */
/* ================= */
void astInitTranMapVtab_(  AstTranMapVtab *vtab, const char *name ) {
/*
*+
*  Name:
*     astInitTranMapVtab

*  Purpose:
*     Initialise a virtual function table for a TranMap.

*  Type:
*     Protected function.

*  Synopsis:
*     #include "tranmap.h"
*     void astInitTranMapVtab( AstTranMapVtab *vtab, const char *name )

*  Class Membership:
*     TranMap vtab initialiser.

*  Description:
*     This function initialises the component of a virtual function
*     table which is used by the TranMap class.

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
   AstMappingVtab *mapping;      /* Pointer to Mapping component of Vtab */

/* Check the local error status. */
   if ( !astOK ) return;

/* Initialize the component of the virtual function table used by the
   parent class. */
   astInitMappingVtab( (AstMappingVtab *) vtab, name );

/* Store a unique "magic" value in the virtual function table. This
   will be used (by astIsATranMap) to determine if an object belongs to
   this class.  We can conveniently use the address of the (static)
   class_init variable to generate this unique value. */
   vtab->check = &class_init;

/* Initialise member function pointers. */
/* ------------------------------------ */
/* Store pointers to the member functions (implemented here) that
   provide virtual methods for this class. */

/* None. */

/* Save the inherited pointers to methods that will be extended, and
   replace them with pointers to the new member functions. */
   mapping = (AstMappingVtab *) vtab;

   parent_transform = mapping->Transform;
   mapping->Transform = Transform;

   parent_mapsplit = mapping->MapSplit;
   mapping->MapSplit = MapSplit;

/* Store replacement pointers for methods which will be over-ridden by
   new member functions implemented here. */
   mapping->MapMerge = MapMerge;
   mapping->Rate = Rate;

/* Declare the copy constructor, destructor and class dump function. */
   astSetCopy( vtab, Copy );
   astSetDelete( vtab, Delete );
   astSetDump( vtab, Dump, "TranMap", "Compound Transformation Mapping" );
}

static int MapMerge( AstMapping *this, int where, int series, int *nmap,
                     AstMapping ***map_list, int **invert_list ) {
/*
*  Name:
*     MapMerge

*  Purpose:
*     Simplify a sequence of Mappings containing a TranMap.

*  Type:
*     Private function.

*  Synopsis:
*     #include "mapping.h"
*     int MapMerge( AstMapping *this, int where, int series, int *nmap,
*                   AstMapping ***map_list, int **invert_list )

*  Class Membership:
*     TranMap method (over-rides the protected astMapMerge method
*     inherited from the Mapping class).

*  Description:
*     This function attempts to simplify a sequence of Mappings by
*     merging a nominated TranMap in the sequence with its neighbours,
*     so as to shorten the sequence if possible.
*
*     In many cases, simplification will not be possible and the
*     function will return -1 to indicate this, without further
*     action.
*
*     In most cases of interest, however, this function will either
*     attempt to replace the nominated TranMap with one which it
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
*        Pointer to the nominated TranMap which is to be merged with
*        its neighbours. This should be a cloned copy of the TranMap
*        pointer contained in the array element "(*map_list)[where]"
*        (see below). This pointer will not be annulled, and the
*        TranMap it identifies will not be modified by this function.
*     where
*        Index in the "*map_list" array (below) at which the pointer
*        to the nominated TranMap resides.
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
   AstCmpMap *cmap;              /* Pointer to compound Mapping */
   AstMapping *cmap_f;           /* Pointer to compound Mapping */
   AstMapping *cmap_i;           /* Pointer to compound Mapping */
   AstMapping *hmap1;            /* Pointer to 1st comp of higher TranMap */
   AstMapping *hmap2;            /* Pointer to 2nd comp of higher TranMap */
   AstMapping *hmap_f;           /* Pointer to fwd Mapping of higher TranMap */
   AstMapping *hmap_i;           /* Pointer to inv Mapping of higher TranMap */
   AstMapping *map1;             /* Pointer to 1st comp of nominated TranMap */
   AstMapping *map2;             /* Pointer to 2nd comp of nominated TranMap */
   AstMapping *map_f;            /* Pointer to fwd Mapping of nominated TranMap */
   AstMapping *map_i;            /* Pointer to inv Mapping of nominated TranMap */
   AstMapping *smap;             /* Pointer to simplified Mapping */
   AstMapping *smap_f;           /* Pointer to simplified Mapping */
   AstMapping *smap_i;           /* Pointer to simplified Mapping */
   AstTranMap *hmap;             /* Pointer to higher TranMap */
   AstTranMap *map;              /* Pointer to this TranMap */
   AstTranMap *new;              /* Pointer to merged TranMap */
   int i;                        /* Loop count */
   int old_hinv1;                /* Original Invert flag for hmap->map1 */
   int old_hinv2;                /* Original Invert flag for hmap->map2 */
   int old_inv1;                 /* Original Invert flag for this->map1 */
   int old_inv2;                 /* Original Invert flag for this->map2 */
   int result;                   /* The value to return */

/* Initialise.*/
   result = -1;

/* Check the inherited status. */
   if ( !astOK ) return result;

/* Get a pointer to this TranMap. */
   map = (AstTranMap *) this;

/* Get the two component Mappings,and temporarily set their Invert
   attributes back to the values they had when the TranMap was created,
   saving their current Invert values so that they can be re-instated later. */
   map1 = map->map1;
   old_inv1 = astGetInvert( map1 );
   astSetInvert( map1, map->invert1 );

   map2 = map->map2;
   old_inv2 = astGetInvert( map2 );
   astSetInvert( map2, map->invert2 );

/* Simplify the TranMap on its own. */
/* ================================ */

/* If the TranMap is inverted, creat an equal TranMap which is not inverted. 
   To do this, invert and swap the component Mappings. */
   if( ( *invert_list )[ where ] ) {
      astInvert( map1 );
      astInvert( map2 );
      new = astTranMap( map2, map1, "" );
      astInvert( map1 );
      astInvert( map2 );

      (void) astAnnul( ( *map_list )[ where ] );
      ( *map_list )[ where ] = (AstMapping *) new;
      ( *invert_list )[ where ] = 0;
      result = where;

/* Otherwise, try to simplify each of the component Mappings. */
   } else {
      smap_f = astSimplify( map1 );
      smap_i = astSimplify( map2 );

/* Assume some simplification took place if the pointers have changed. */
      if( smap_f != map1 || smap_i != map2 ) { 

/* Construct a new TranMap from these simplifgied Mappings. */
         (void) astAnnul( ( *map_list )[ where ] );
         ( *map_list )[ where ] = (AstMapping *) astTranMap( smap_f, smap_i, "" );
         result = where;

/* Otherwise, if the both component Mappings are defined in both directions... */
      } else if( astGetTranForward( map1 ) && astGetTranInverse( map1 ) &&    
                 astGetTranForward( map2 ) && astGetTranInverse( map2 ) ) {

/* Form a series CmpMap from the two component Mappings, with the second
   Mapping inverted. */
         astInvert( map2 );
         cmap = astCmpMap( map1, map2, 1, "" );
         astInvert( map2 );

/* If this CmpMap simplifies to a UnitMap, then the two components of the 
   TranMap are equal, and so we can replace the entire TranMap with either
   of its components. Note, we leave the supplied invert flag unchanged,
   since the copycreated below refers to the Mapping as it was when the
   TranMap was created. However, we invert the returned Mapping if
   necessary. */
         smap = astSimplify( cmap );
         if( astIsAUnitMap( smap ) ) {
            (void) astAnnul( ( *map_list )[ where ] );
            ( *map_list )[ where ] = astCopy( map1 );
            if( ( *invert_list )[ where ] ) astInvert( ( *map_list )[ where ] );
            result = where;
         }

/* Release resources. */
         smap = astAnnul( smap );
         cmap = astAnnul( cmap );
      }

/* Release resources. */
      smap_f = astAnnul( smap_f );
      smap_i = astAnnul( smap_i );
   }

/* Merge the TranMap with a neighbouring TranMap. */
/* ============================================== */
/* Only do this if no change was made above, and we are combining the
   Mappings in series. */
   if( result == -1 && series ) {

/* Is the higher neighbour a TranMap? */
      if( where < ( *nmap - 1 ) && 
          astIsATranMap( ( *map_list )[ where + 1 ] ) ){

/* Get the two component Mappings of the higher TranMap, and temporarily set 
   their Invert attributes back to the values they had when the TranMap was 
   created, saving their current Invert values so that they can be re-instated 
   later. */
         hmap = (AstTranMap *) ( *map_list )[ where + 1 ];

         hmap1 = hmap->map1;
         old_hinv1 = astGetInvert( hmap1 );
         astSetInvert( hmap1, hmap->invert1 );

         hmap2 = hmap->map2;
         old_hinv2 = astGetInvert( hmap2 );
         astSetInvert( hmap2, hmap->invert2 );

/* Get the Mappings which defines the forward and inverse transformation of 
   the lower TranMap ("this"). Then, map_f and map_i are pointers to
   Mappings which could be used to construct a new TranMap which would be
   equivalent to "this" with the supplied invert setting. */
         if( ( *invert_list )[ where ] ) {
            map_f = map2;
            map_i = map1;
            astInvert( map_f );
            astInvert( map_i );
         } else {
            map_f = map1;
            map_i = map2;
         }

/* Likewise, get the Mappings which defines the forward and inverse 
   transformation of the higher TranMap. */
         if( ( *invert_list )[ where + 1 ] ) {
            hmap_f = hmap2;
            hmap_i = hmap1;
            astInvert( hmap_f );
            astInvert( hmap_i );
         } else {
            hmap_f = hmap1;
            hmap_i = hmap2;
         }

/* Combine the two forward Mappings together into a series CmpMap, and
   simplify it. */
         cmap_f = (AstMapping *) astCmpMap( map_f,  hmap_f, 1, "" );
         smap_f = astSimplify( cmap_f );

/* Do the same for the inverse Mappings */
         cmap_i = (AstMapping *) astCmpMap( map_i,  hmap_i, 1, "" );
         smap_i = astSimplify( cmap_i );

/* Was any simplification performed? We assume this is the case if the
   either of the simplied pointer differs from the original pointer. */
         if( cmap_f != smap_f || cmap_i != smap_i ) {

/* In which case,construct a new TranMap from the simplified Mappings. */
            new = astTranMap( smap_f, smap_i, "" );

         } else {
            new = NULL;
         }

/* Free resources.*/
         cmap_f = astAnnul( cmap_f );
         smap_f = astAnnul( smap_f );
         cmap_i = astAnnul( cmap_i );
         smap_i = astAnnul( smap_i );

/* Re-instate the original Invert values for the component Mappings of
   the higher TranMap. */
         astSetInvert( hmap1, old_hinv1 );
         astSetInvert( hmap2, old_hinv2 );

/* If we have a new TranMap, annul the first of the two Mappings, and replace 
   it with the merged TranMap. Also set the invert flag. */ 
         if( new ) {
            (void) astAnnul( ( *map_list )[ where ] );
            ( *map_list )[ where ] = (AstMapping *) new;
            ( *invert_list )[ where ] = 0;

/* Annul the second of the two Mappings, and shuffle down the rest of the 
   list to fill the gap. */
            (void) astAnnul( ( *map_list )[ where + 1 ] );
            for ( i = where + 2; i < *nmap; i++ ) {
               ( *map_list )[ i - 1 ] = ( *map_list )[ i ];
               ( *invert_list )[ i - 1 ] = ( *invert_list )[ i ];
            }

/* Clear the vacated element at the end. */
            ( *map_list )[ *nmap - 1 ] = NULL;
            ( *invert_list )[ *nmap - 1 ] = 0;

/* Decrement the Mapping count and return the index of the first
   modified element. */
            ( *nmap )--;
            result = where;
         }
      }
   }

/* Re-instate the original Invert values for the component Mappings. */
   astSetInvert( map1, old_inv1 );
   astSetInvert( map2, old_inv2 );

/* If an error occurred, clear the result value. */
   if ( !astOK ) result = -1;

/* Return the result. */
   return result;
}

static int *MapSplit( AstMapping *this_map, int nin, int *in, AstMapping **map ){
/*
*  Name:
*     MapSplit

*  Purpose:
*     Create a Mapping representing a subset of the inputs of an existing
*     TranMap.

*  Type:
*     Private function.

*  Synopsis:
*     #include "tranmap.h"
*     int *MapSplit( AstMapping *this, int nin, int *in, AstMapping **map )

*  Class Membership:
*     TranMap method (over-rides the protected astMapSplit method
*     inherited from the Mapping class).

*  Description:
*     This function creates a new Mapping by picking specified inputs from 
*     an existing TranMap. This is only possible if the specified inputs
*     correspond to some subset of the TranMap outputs. That is, there
*     must exist a subset of the TranMap outputs for which each output
*     depends only on the selected TranMap inputs, and not on any of the
*     inputs which have not been selected. If this condition is not met
*     by the supplied TranMap, then a NULL Mapping is returned.

*  Parameters:
*     this
*        Pointer to the TranMap to be split (the TranMap is not actually 
*        modified by this function).
*     nin
*        The number of inputs to pick from "this".
*     in
*        Pointer to an array of indices (zero based) for the inputs which
*        are to be picked. This array should have "nin" elements. If "Nin"
*        is the number of inputs of the supplied TranMap, then each element 
*        should have a value in the range zero to Nin-1.
*     map
*        Address of a location at which to return a pointer to the new
*        Mapping. This Mapping will have "nin" inputs (the number of
*        outputs may be different to "nin"). A NULL pointer will be
*        returned if the supplied TranMap has no subset of outputs which 
*        depend only on the selected inputs.

*  Returned Value:
*     A pointer to a dynamically allocated array of ints. The number of
*     elements in this array will equal the number of outputs for the 
*     returned Mapping. Each element will hold the index of the
*     corresponding output in the supplied TranMap. The array should be
*     freed using astFree when no longer needed. A NULL pointer will
*     be returned if no output Mapping can be created.

*  Notes:
*     - If this function is invoked with the global error status set,
*     or if it should fail for any reason, then NULL values will be
*     returned as the function value and for the "map" pointer.
*/

/* Local Variables: */
   AstMapping *fmap;          /* Pointer to forward Mapping in supplied TranMap */
   AstMapping *imap;          /* Pointer to inverse Mapping in supplied TranMap */
   AstMapping *rfmap;         /* Pointer to split forward Mapping */
   AstMapping *rimap;         /* Pointer to split inverse Mapping */
   AstTranMap *this;          /* Pointer to TranMap structure */
   int *ires;                 /* I/ps of inv Mapping dependent on selected o/ps */
   int *out;                  /* O/ps of fwd Mapping dependent on selected i/ps */
   int *result;               /* Pointer to returned array */
   int finv;                  /* Invert flag to use with fmap */
   int i;                     /* Loop count */
   int iinv;                  /* Invert flag to use with imap */
   int nout;                  /* No. of outputs dependent on selected inputs */
   int ok;                    /* Can required Mapping be created? */
   int old_finv;              /* Original Invert flag for fmap */
   int old_iinv;              /* Original Invert flag for imap */

/* Initialise */
   result = NULL;
   *map = NULL;

/* Check the global error status. */
   if ( !astOK ) return result;

/* Invoke the parent astMapSplit method to see if it can do the job. */
   result = (*parent_mapsplit)( this_map, nin, in, map );

/* If not, we provide a special implementation here. */
   if( !result ) {

/* Get a pointer to the TranMap structure. */
      this = (AstTranMap *) this_map;

/* Get pointers to the forward and inverse Mappings, taking into account
   whether the TranMap has been inverted. */
      if( !astGetInvert( this ) ) {
         fmap = this->map1;
         finv = this->invert1;
         imap = this->map2;
         iinv = this->invert2;
      } else {
         imap = this->map1;
         iinv = this->invert1;
         fmap = this->map2;
         finv = this->invert2;
      }

/* Temporarily set the Invert flag of both Mappings back to their
   original values. */ 
      old_finv = astGetInvert( fmap );
      astSetInvert( fmap, finv );
      old_iinv = astGetInvert( imap );
      astSetInvert( imap, iinv );

/* Try to split the forward Mapping. */
      out = astMapSplit( fmap, nin, in, &rfmap );

/* Check the split could be done. */
      if( out ) {

/* Get the number of outputs which are fed by the selected inputs. */
         nout = astGetNout( rfmap );

/* See if the inverse Mapping can be split using these outputs as inputs. */
         astInvert( imap );
         ires = astMapSplit( imap, nout, out, &rimap );
         astInvert( imap );
         if( ires ) {
            astInvert( rimap );

/* Check that the resulting inputs are the same as the supplied inputs. */
            if( astGetNin( rimap ) == nin ) {
               ok = 1;
               for( i = 0; i < nin; i++ ) {
                  if( in[ i ] != ires[ i ] ) {
                     ok = 0;
                     break;
                  }
               }

/* If so create the required new TranMap. */
               if( ok ) {
                  *map = (AstMapping *) astTranMap( rfmap, rimap, "" );
                  result = out;
               }
            }

/* Free resources. */
            ires = astFree( ires );
            rimap = astAnnul( rimap );
         }

         if( !result ) out = astFree( out );
         rfmap = astAnnul( rfmap );
      }

/* Re-instate the Invert flags of the component Mappings. */
      astSetInvert( fmap, old_finv );
      astSetInvert( imap, old_iinv );
   }

/* Free returned resources if an error has occurred. */
   if( !astOK ) {
      result = astFree( result );
      *map = astAnnul( *map );
   }

/* Return the list of output indices. */
   return result;
}

static double Rate( AstMapping *this, double *at, int ax1, int ax2 ){
/*
*  Name:
*     Rate

*  Purpose:
*     Calculate the rate of change of a Mapping output.

*  Type:
*     Private function.

*  Synopsis:
*     #include "tranmap.h"
*     result = Rate( AstMapping *this, double *at, int ax1, int ax2 )

*  Class Membership:
*     TranMap member function (overrides the astRate method inherited
*     from the Mapping class ).

*  Description:
*     This function returns the rate of change of a specified output of 
*     the supplied Mapping with respect to a specified input, at a 
*     specified input position. Also evaluates the second derivative.

*  Parameters:
*     this
*        Pointer to the Mapping to be applied.
*     at
*        The address of an array holding the axis values at the position 
*        at which the rate of change is to be evaluated. The number of 
*        elements in this array should equal the number of inputs to the 
*        Mapping.
*     ax1
*        The index of the Mapping output for which the rate of change is to 
*        be found (output numbering starts at 0 for the first output).
*     ax2
*        The index of the Mapping input which is to be varied in order to
*        find the rate of change (input numbering starts at 0 for the first 
*        input).

*  Returned Value:
*     The rate of change of Mapping output "ax1" with respect to input 
*     "ax2", evaluated at "at", or AST__BAD if the value cannot be 
*     calculated.

*/

/* Local Variables: */
   AstTranMap *map;
   AstMapping *cmap;
   double result;
   int cinv;
   int old_inv;

/* Check inherited status */
   if( !astOK ) return AST__BAD;

/* Get a pointer to the TranMap structure. */
   map = (AstTranMap *) this;

/* Choose the component Mapping to use, and get its original Invert
   value. Invert this if the TranMap itself has been inverted (this is
   because the astRate function has no "invert" argument so we need to
   invert the Mapping before calling astRate). */
   if( astGetInvert( this ) ) {
      cmap = map->map2;
      cinv = !(map->invert2);
   } else {
      cmap = map->map1;
      cinv = map->invert1;
   }

/* Temporarily set the Invert flag of the component Mapping back to its
   original value. */ 
   old_inv = astGetInvert( cmap );
   astSetInvert( cmap, cinv );

/* Use the astRate method of the component Mapping. */
   result = astRate( cmap, at, ax1, ax2 );

/* Re-instate the Invert flag of the component Mapping. */
   astSetInvert( cmap, old_inv );

/* Return the result. */
   return result;
}

static AstPointSet *Transform( AstMapping *this, AstPointSet *in,
                               int forward, AstPointSet *out ) {
/*
*  Name:
*     Transform

*  Purpose:
*     Apply a TranMap to transform a set of points.

*  Type:
*     Private function.

*  Synopsis:
*     #include "tranmap.h"
*     AstPointSet *Transform( AstMapping *this, AstPointSet *in,
*                             int forward, AstPointSet *out )

*  Class Membership:
*     TranMap member function (over-rides the astTransform method inherited
*     from the Mapping class).

*  Description:
*     This function takes a TranMap and a set of points encapsulated in a
*     PointSet and transforms the points so as to apply the required Mapping.
*     This implies applying each of the TranMap's component Mappings in turn,
*     either in series or in parallel.

*  Parameters:
*     this
*        Pointer to the TranMap.
*     in
*        Pointer to the PointSet associated with the input coordinate values.
*     forward
*        A non-zero value indicates that the forward coordinate transformation
*        should be applied, while a zero value requests the inverse
*        transformation.
*     out
*        Pointer to a PointSet which will hold the transformed (output)
*        coordinate values. A NULL value may also be given, in which case a
*        new PointSet will be created by this function.

*  Returned Value:
*     Pointer to the output (possibly new) PointSet.

*  Notes:
*     -  A null pointer will be returned if this function is invoked with the
*     global error status set, or if it should fail for any reason.
*     -  The number of coordinate values per point in the input PointSet must
*     match the number of coordinates for the TranMap being applied.
*     -  If an output PointSet is supplied, it must have space for sufficient
*     number of points and coordinate values per point to accommodate the
*     result. Any excess space will be ignored.
*/

/* Local Variables: */
   AstMapping *cmap;             /* Mapping which defines the required transformation */
   AstPointSet *result;          /* Pointer to output PointSet */
   AstTranMap *map;              /* Pointer to TranMap to be applied */
   int cinv;                     /* Invert flag when TranMap was created */
   int old_inv;                  /* Invert flag on entry to this function */

/* Check the global error status. */
   if ( !astOK ) return NULL;

/* Obtain a pointer to the TranMap. */
   map = (AstTranMap *) this;

/* Apply the parent Mapping using the stored pointer to the Transform member
   function inherited from the parent Mapping class. This function validates
   all arguments and generates an output PointSet if necessary, but does not
   actually transform any coordinate values. */
   result = (*parent_transform)( this, in, forward, out );

/* We now extend the parent astTransform method by applying the component
   Mappings of the TranMap to generate the output coordinate values. */

/* Determine whether to apply the forward or inverse Mapping, according to the
   direction specified and whether the Mapping has been inverted. */
   if ( astGetInvert( map ) ) forward = !forward;

/* Choose the component Mapping to use, and get its original Invert value. */
   if( forward ) {
      cmap = map->map1;
      cinv = map->invert1;
   }else {
      cmap = map->map2;
      cinv = map->invert2;
   }

/* Temporarily set the Invert flag of the component Mapping back to its
   original value. */ 
   old_inv = astGetInvert( cmap );
   astSetInvert( cmap, cinv );

/* Use the Transform method of the component Mapping. */
   result = astTransform( cmap, in, forward, out );

/* Re-instate the Invert flag of the component Mapping. */
   astSetInvert( cmap, old_inv );

/* If an error occurred, clean up by deleting the output PointSet (if
   allocated by this function) and setting a NULL result pointer. */
   if ( !astOK ) {
      if ( !out ) result = astDelete( result );
      result = NULL;
   }

/* Return a pointer to the output PointSet. */
   return result;
}

/* Copy constructor. */
/* ----------------- */
static void Copy( const AstObject *objin, AstObject *objout ) {
/*
*  Name:
*     Copy

*  Purpose:
*     Copy constructor for TranMap objects.

*  Type:
*     Private function.

*  Synopsis:
*     void Copy( const AstObject *objin, AstObject *objout )

*  Description:
*     This function implements the copy constructor for TranMap objects.

*  Parameters:
*     objin
*        Pointer to the object to be copied.
*     objout
*        Pointer to the object being constructed.

*  Returned Value:
*     void

*  Notes:
*     -  This constructor makes a deep copy, including a copy of the component
*     Mappings within the TranMap.
*/

/* Local Variables: */
   AstTranMap *in;                /* Pointer to input TranMap */
   AstTranMap *out;               /* Pointer to output TranMap */

/* Check the global error status. */
   if ( !astOK ) return;

/* Obtain pointers to the input and output TranMaps. */
   in = (AstTranMap *) objin;
   out = (AstTranMap *) objout;

/* For safety, start by clearing any references to the input component
   Mappings from the output TranMap. */
   out->map1 = NULL;
   out->map2 = NULL;

/* Make copies of these Mappings and store pointers to them in the output
   TranMap structure. */
   out->map1 = astCopy( in->map1 );
   out->map2 = astCopy( in->map2 );
}

/* Destructor. */
/* ----------- */
static void Delete( AstObject *obj ) {
/*
*  Name:
*     Delete

*  Purpose:
*     Destructor for TranMap objects.

*  Type:
*     Private function.

*  Synopsis:
*     void Delete( AstObject *obj )

*  Description:
*     This function implements the destructor for TranMap objects.

*  Parameters:
*     obj
*        Pointer to the object to be deleted.

*  Returned Value:
*     void

*  Notes:
*     This function attempts to execute even if the global error status is
*     set.
*/

/* Local Variables: */
   AstTranMap *this;              /* Pointer to TranMap */

/* Obtain a pointer to the TranMap structure. */
   this = (AstTranMap *) obj;

/* Annul the pointers to the component Mappings. */
   this->map1 = astAnnul( this->map1 );
   this->map2 = astAnnul( this->map2 );

/* Clear the remaining TranMap variables. */
   this->invert1 = 0;
   this->invert2 = 0;
}

/* Dump function. */
/* -------------- */
static void Dump( AstObject *this_object, AstChannel *channel ) {
/*
*  Name:
*     Dump

*  Purpose:
*     Dump function for TranMap objects.

*  Type:
*     Private function.

*  Synopsis:
*     void Dump( AstObject *this, AstChannel *channel )

*  Description:
*     This function implements the Dump function which writes out data
*     for the TranMap class to an output Channel.

*  Parameters:
*     this
*        Pointer to the TranMap whose data are being written.
*     channel
*        Pointer to the Channel to which the data are being written.
*/

/* Local Variables: */
   AstTranMap *this;              /* Pointer to the TranMap structure */
   int ival;                     /* Integer value */
   int set;                      /* Attribute value set? */

/* Check the global error status. */
   if ( !astOK ) return;

/* Obtain a pointer to the TranMap structure. */
   this = (AstTranMap *) this_object;

/* Write out values representing the instance variables for the TranMap
   class.  Accompany these with appropriate comment strings, possibly
   depending on the values being written.*/

/* In the case of attributes, we first use the appropriate (private)
   Test...  member function to see if they are set. If so, we then use
   the (private) Get... function to obtain the value to be written
   out.

   For attributes which are not set, we use the astGet... method to
   obtain the value instead. This will supply a default value
   (possibly provided by a derived class which over-rides this method)
   which is more useful to a human reader as it corresponds to the
   actual default attribute value.  Since "set" will be zero, these
   values are for information only and will not be read back. */

/* First Invert flag. */
/* ------------------ */
   ival = this->invert1;
   set = ( ival != 0 );
   astWriteInt( channel, "InvA", set, 0, ival,
                ival ? "First Mapping used in inverse direction" :
                       "First Mapping used in forward direction" );

/* Second Invert flag. */
/* ------------------- */
   ival = this->invert2;
   set = ( ival != 0 );
   astWriteInt( channel, "InvB", set, 0, ival,
                ival ? "Second Mapping used in inverse direction" :
                       "Second Mapping used in forward direction" );

/* First Mapping. */ 
/* -------------- */
   astWriteObject( channel, "MapA", 1, 1, this->map1, 
                   "Mapping for forward transformation" );

/* Second Mapping. */
/* --------------- */
   astWriteObject( channel, "MapB", 1, 1, this->map2,
                   "Mapping for inverse transformation" );
}

/* Standard class functions. */
/* ========================= */
/* Implement the astIsATranMap and astCheckTranMap functions using the
   macros defined for this purpose in the "object.h" header file. */
astMAKE_ISA(TranMap,Mapping,check,&class_init)
astMAKE_CHECK(TranMap)

AstTranMap *astTranMap_( void *map1_void, void *map2_void, const char *options, ... ) {
/*
*+
*  Name:
*     astTranMap

*  Purpose:
*     Create a TranMap.

*  Type:
*     Protected function.

*  Synopsis:
*     #include "tranmap.h"
*     AstTranMap *astTranMap( AstMapping *map1, AstMapping *map2, const char *options, ... )

*  Class Membership:
*     TranMap constructor.

*  Description:
*     This function creates a new TranMap and optionally initialises its
*     attributes.

*  Parameters:
*     map1
*        Pointer to the first Mapping (which deinfes the forward
*        transformation).
*     map2
*        Pointer to the second Mapping (which deinfes the inverse
*        transformation).
*     options
*        Pointer to a null terminated string containing an optional
*        comma-separated list of attribute assignments to be used for
*        initialising the new TranMap. The syntax used is the same as for the
*        astSet method and may include "printf" format specifiers identified
*        by "%" symbols in the normal way.
*     ...
*        If the "options" string contains "%" format specifiers, then an
*        optional list of arguments may follow it in order to supply values to
*        be substituted for these specifiers. The rules for supplying these
*        are identical to those for the astSet method (and for the C "printf"
*        function).

*  Returned Value:
*     A pointer to the new TranMap.

*  Notes:
*     - A null pointer will be returned if this function is invoked
*     with the global error status set, or if it should fail for any
*     reason.
*-

*  Implementation Notes:
*     - This function implements the basic TranMap constructor which is
*     available via the protected interface to the TranMap class.  A
*     public interface is provided by the astTranMapId_ function.
*     - Because this function has a variable argument list, it is
*     invoked by a macro that evaluates to a function pointer (not a
*     function invocation) and no checking or casting of arguments is
*     performed before the function is invoked. Because of this, the
*     "map1" and "map2" parameters are of type (void *) and are
*     converted and validated within the function itself.
*/

/* Local Variables: */
   AstTranMap *new;              /* Pointer to new TranMap */
   AstMapping *map1;             /* Pointer to first Mapping structure */
   AstMapping *map2;             /* Pointer to second Mapping structure */
   va_list args;                 /* Variable argument list */

/* Initialise. */
   new = NULL;

/* Check the global status. */
   if ( !astOK ) return new;

/* Obtain and validate pointers to the Mapping structures provided. */
   map1 = astCheckMapping( map1_void );
   map2 = astCheckMapping( map2_void );
   if ( astOK ) {

/* Initialise the TranMap, allocating memory and initialising the
   virtual function table as well if necessary. */
      new = astInitTranMap( NULL, sizeof( AstTranMap ), !class_init, &class_vtab,
                           "TranMap", map1, map2 );

/* If successful, note that the virtual function table has been
   initialised. */
      if ( astOK ) {
         class_init = 1;

/* Obtain the variable argument list and pass it along with the
   options string to the astVSet method to initialise the new TranMap's
   attributes. */
         va_start( args, options );
         astVSet( new, options, args );
         va_end( args );

/* If an error occurred, clean up by deleting the new object. */
         if ( !astOK ) new = astDelete( new );
      }
   }

/* Return a pointer to the new TranMap. */
   return new;
}

AstTranMap *astTranMapId_( void *map1_void, void *map2_void, 
                           const char *options, ... ) {
/*
*++
*  Name:
c     astTranMap
f     AST_TRANMAP

*  Purpose:
*     Create a TranMap.

*  Type:
*     Public function.

*  Synopsis:
c     #include "tranmap.h"
c     AstTranMap *astTranMap( AstMapping *map1, AstMapping *map2,
c                           const char *options, ... )
f     RESULT = AST_TRANMAP( MAP1, MAP2, OPTIONS, STATUS )

*  Class Membership:
*     TranMap constructor.

*  Description:
*     This function creates a new TranMap and optionally initialises
*     its attributes.
*
*     A TranMap is a Mapping which combines the forward transformation of
*     a supplied Mapping with the inverse transformation of another
*     supplied Mapping, ignoring the un-used transformation in each
*     Mapping (indeed the un-used transformation need not exist).
*
*     When the forward transformation of the TranMap is referred to, the
*     transformation actually used is the forward transformation of the
*     first Mapping supplied when the TranMap was constructed. Likewise, 
*     when the inverse transformation of the TranMap is referred to, the
*     transformation actually used is the inverse transformation of the
*     second Mapping supplied when the TranMap was constructed.

*  Parameters:
c     map1
f     MAP1 = INTEGER (Given)
*        Pointer to the first component Mapping, which defines the
*        forward transformation.
c     map2
f     MAP2 = INTEGER (Given)
*        Pointer to the second component Mapping, which defines the
*        inverse transformation.
c     options
f     OPTIONS = CHARACTER * ( * ) (Given)
c        Pointer to a null-terminated string containing an optional
c        comma-separated list of attribute assignments to be used for
c        initialising the new TranMap. The syntax used is identical to
c        that for the astSet function and may include "printf" format
c        specifiers identified by "%" symbols in the normal way.
f        A character string containing an optional comma-separated
f        list of attribute assignments to be used for initialising the
f        new TranMap. The syntax used is identical to that for the
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
c     astTranMap()
f     AST_TRANMAP = INTEGER
*        A pointer to the new TranMap.

*  Notes:
*     - The number of output coordinates generated by the two Mappings
*     (their Nout attribute) must be equal, as must the number of input
*     coordinates accepted by each Mapping (their Nin attribute).
*     - The forward transformation of the first Mapping must exist.
*     - The inverse transformation of the second Mapping must exist.
c     - Note that the component Mappings supplied are not copied by
c     astTranMap (the new TranMap simply retains a reference to
c     them). They may continue to be used for other purposes, but
c     should not be deleted. If a TranMap containing a copy of its
c     component Mappings is required, then a copy of the TranMap should
c     be made using astCopy.
f     - Note that the component Mappings supplied are not copied by
f     AST_TRANMAP (the new TranMap simply retains a reference to
f     them). They may continue to be used for other purposes, but
f     should not be deleted. If a TranMap containing a copy of its
f     component Mappings is required, then a copy of the TranMap should
f     be made using AST_COPY.
*     - A null Object pointer (AST__NULL) will be returned if this
c     function is invoked with the AST error status set, or if it
f     function is invoked with STATUS set to an error value, or if it
*     should fail for any reason.
*--

*  Implementation Notes:
*     - This function implements the external (public) interface to
*     the astTranMap constructor function. It returns an ID value
*     (instead of a true C pointer) to external users, and must be
*     provided because astTranMap_ has a variable argument list which
*     cannot be encapsulated in a macro (where this conversion would
*     otherwise occur).
*     - Because no checking or casting of arguments is performed
*     before the function is invoked, the "map1" and "map2" parameters
*     are of type (void *) and are converted from an ID value to a
*     pointer and validated within the function itself.
*     - The variable argument list also prevents this function from
*     invoking astTranMap_ directly, so it must be a re-implementation
*     of it in all respects, except for the conversions between IDs
*     and pointers on input/output of Objects.
*/

/* Local Variables: */
   AstTranMap *new;               /* Pointer to new TranMap */
   AstMapping *map1;             /* Pointer to first Mapping structure */
   AstMapping *map2;             /* Pointer to second Mapping structure */
   va_list args;                 /* Variable argument list */

/* Initialise. */
   new = NULL;

/* Check the global status. */
   if ( !astOK ) return new;

/* Obtain the Mapping pointers from the ID's supplied and validate the
   pointers to ensure they identify valid Mappings. */
   map1 = astCheckMapping( astMakePointer( map1_void ) );
   map2 = astCheckMapping( astMakePointer( map2_void ) );
   if ( astOK ) {

/* Initialise the TranMap, allocating memory and initialising the
   virtual function table as well if necessary. */
      new = astInitTranMap( NULL, sizeof( AstTranMap ), !class_init, &class_vtab,
                           "TranMap", map1, map2 );

/* If successful, note that the virtual function table has been initialised. */
      if ( astOK ) {
         class_init = 1;

/* Obtain the variable argument list and pass it along with the
   options string to the astVSet method to initialise the new TranMap's
   attributes. */
         va_start( args, options );
         astVSet( new, options, args );
         va_end( args );

/* If an error occurred, clean up by deleting the new object. */
         if ( !astOK ) new = astDelete( new );
      }
   }

/* Return an ID value for the new TranMap. */
   return astMakeId( new );
}

AstTranMap *astInitTranMap_( void *mem, size_t size, int init,
                           AstTranMapVtab *vtab, const char *name,
                           AstMapping *map1, AstMapping *map2 ) {
/*
*+
*  Name:
*     astInitTranMap

*  Purpose:
*     Initialise a TranMap.

*  Type:
*     Protected function.

*  Synopsis:
*     #include "tranmap.h"
*     AstTranMap *astInitTranMap( void *mem, size_t size, int init,
*                               AstTranMapVtab *vtab, const char *name,
*                               AstMapping *map1, AstMapping *map2 )

*  Class Membership:
*     TranMap initialiser.

*  Description:
*     This function is provided for use by class implementations to initialise
*     a new TranMap object. It allocates memory (if necessary) to
*     accommodate the TranMap plus any additional data associated with the
*     derived class. It then initialises a TranMap structure at the start
*     of this memory. If the "init" flag is set, it also initialises the
*     contents of a virtual function table for a TranMap at the start of
*     the memory passed via the "vtab" parameter.

*  Parameters:
*     mem
*        A pointer to the memory in which the TranMap is to be initialised.
*        This must be of sufficient size to accommodate the TranMap data
*        (sizeof(TranMap)) plus any data used by the derived class. If a
*        value of NULL is given, this function will allocate the memory itself
*        using the "size" parameter to determine its size.
*     size
*        The amount of memory used by the TranMap (plus derived class
*        data). This will be used to allocate memory if a value of NULL is
*        given for the "mem" parameter. This value is also stored in the
*        TranMap structure, so a valid value must be supplied even if not
*        required for allocating memory.
*     init
*        A logical flag indicating if the TranMap's virtual function table
*        is to be initialised. If this value is non-zero, the virtual function
*        table will be initialised by this function.
*     vtab
*        Pointer to the start of the virtual function table to be associated
*        with the new TranMap.
*     name
*        Pointer to a constant null-terminated character string which contains
*        the name of the class to which the new object belongs (it is this
*        pointer value that will subsequently be returned by the Object
*        astClass function).
*     map1
*        Pointer to the first Mapping.
*     map2
*        Pointer to the second Mapping.

*  Returned Value:
*     A pointer to the new TranMap.

*  Notes:
*     -  A null pointer will be returned if this function is invoked with the
*     global error status set, or if it should fail for any reason.
*-
*/

/* Local Variables: */
   AstTranMap *new;               /* Pointer to new TranMap */
   int nin;                      /* No. input coordinates for TranMap */
   int nout;                     /* No. output coordinates for TranMap */

/* Check the global status. */
   if ( !astOK ) return NULL;

/* If necessary, initialise the virtual function table. */
   if ( init ) astInitTranMapVtab( vtab, name );

/* Initialise. */
   new = NULL;

/* Report an error if map1 has no forward transformation. */
   if( !astGetTranForward( map1 ) && astOK ) {
      astError( AST__INTRD, "astInitTranMap(%s): The first supplied Mapping "
              "is not able to transform coordinates in the forward direction.",
              name );
   }

/* Report an error if map2 has no inverse transformation. */
   if( !astGetTranInverse( map2 ) && astOK ) {
      astError( AST__INTRD, "astInitTranMap(%s): The second supplied Mapping "
              "is not able to transform coordinates in the inverse direction.",
              name );
   }

/* Check that the number of coordinates are compatible and report an error if 
   they are not. */
   nout = astGetNout( map1 );
   if ( astGetNout( map2 ) != nout && astOK ) {
      astError( AST__INNCO, "astInitTranMap(%s): The number of output "
                      "coordinates per point (%d) for the first Mapping "
                      "supplied does not match the number of output "
                      "coordinates (%d) for the second Mapping.", name, nout,
                      astGetNout( map2 ) );
   }

   nin = astGetNin( map1 );
   if ( astGetNin( map2 ) != nin && astOK ) {
      astError( AST__INNCO, "astInitTranMap(%s): The number of input "
                      "coordinates per point (%d) for the first Mapping "
                      "supplied does not match the number of input "
                      "coordinates (%d) for the second Mapping.", name, nin,
                      astGetNin( map2 ) );
   }

/* Initialise a Mapping structure (the parent class) as the first component
   within the TranMap structure, allocating memory if necessary. Specify
   the number of input and output coordinates and in which directions the
   Mapping should be defined. */
   if ( astOK ) {
      new = (AstTranMap *) astInitMapping( mem, size, 0,
                                          (AstMappingVtab *) vtab, name,
                                          nin, nout, 1, 1 );

      if ( astOK ) {

/* Initialise the TranMap data. */
/* --------------------------- */
/* Store pointers to the component Mappings. */
         new->map1 = astClone( map1 );
         new->map2 = astClone( map2 );

/* Save the initial values of the inversion flags for these Mappings. */
         new->invert1 = astGetInvert( map1 );
         new->invert2 = astGetInvert( map2 );

/* If an error occurred, clean up by annulling the Mapping pointers and
   deleting the new object. */
         if ( !astOK ) {
            new->map1 = astAnnul( new->map1 );
            new->map2 = astAnnul( new->map2 );
            new = astDelete( new );
         }
      }
   }

/* Return a pointer to the new object. */
   return new;
}

AstTranMap *astLoadTranMap_( void *mem, size_t size,
                           AstTranMapVtab *vtab, const char *name,
                           AstChannel *channel ) {
/*
*+
*  Name:
*     astLoadTranMap

*  Purpose:
*     Load a TranMap.

*  Type:
*     Protected function.

*  Synopsis:
*     #include "tranmap.h"
*     AstTranMap *astLoadTranMap( void *mem, size_t size,
*                               AstTranMapVtab *vtab, const char *name,
*                               AstChannel *channel )

*  Class Membership:
*     TranMap loader.

*  Description:
*     This function is provided to load a new TranMap using data read
*     from a Channel. It first loads the data used by the parent class
*     (which allocates memory if necessary) and then initialises a
*     TranMap structure in this memory, using data read from the input
*     Channel.
*
*     If the "init" flag is set, it also initialises the contents of a
*     virtual function table for a TranMap at the start of the memory
*     passed via the "vtab" parameter.


*  Parameters:
*     mem
*        A pointer to the memory into which the TranMap is to be
*        loaded.  This must be of sufficient size to accommodate the
*        TranMap data (sizeof(TranMap)) plus any data used by derived
*        classes. If a value of NULL is given, this function will
*        allocate the memory itself using the "size" parameter to
*        determine its size.
*     size
*        The amount of memory used by the TranMap (plus derived class
*        data).  This will be used to allocate memory if a value of
*        NULL is given for the "mem" parameter. This value is also
*        stored in the TranMap structure, so a valid value must be
*        supplied even if not required for allocating memory.
*
*        If the "vtab" parameter is NULL, the "size" value is ignored
*        and sizeof(AstTranMap) is used instead.
*     vtab
*        Pointer to the start of the virtual function table to be
*        associated with the new TranMap. If this is NULL, a pointer to
*        the (static) virtual function table for the TranMap class is
*        used instead.
*     name
*        Pointer to a constant null-terminated character string which
*        contains the name of the class to which the new object
*        belongs (it is this pointer value that will subsequently be
*        returned by the astGetClass method).
*
*        If the "vtab" parameter is NULL, the "name" value is ignored
*        and a pointer to the string "TranMap" is used instead.

*  Returned Value:
*     A pointer to the new TranMap.

*  Notes:
*     - A null pointer will be returned if this function is invoked
*     with the global error status set, or if it should fail for any
*     reason.
*-
*/

/* Local Variables: */
   AstTranMap *new;               /* Pointer to the new TranMap */

/* Initialise. */
   new = NULL;

/* Check the global error status. */
   if ( !astOK ) return new;

/* If a NULL virtual function table has been supplied, then this is
   the first loader to be invoked for this TranMap. In this case the
   TranMap belongs to this class, so supply appropriate values to be
   passed to the parent class loader (and its parent, etc.). */
   if ( !vtab ) {
      size = sizeof( AstTranMap );
      vtab = &class_vtab;
      name = "TranMap";

/* If required, initialise the virtual function table for this class. */
      if ( !class_init ) {
         astInitTranMapVtab( vtab, name );
         class_init = 1;
      }
   }

/* Invoke the parent class loader to load data for all the ancestral
   classes of the current one, returning a pointer to the resulting
   partly-built TranMap. */
   new = astLoadMapping( mem, size, (AstMappingVtab *) vtab, name,
                         channel );

   if ( astOK ) {

/* Read input data. */
/* ================ */
/* Request the input Channel to read all the input data appropriate to
   this class into the internal "values list". */
      astReadClassData( channel, "TranMap" );

/* Now read each individual data item from this list and use it to
   initialise the appropriate instance variable(s) for this class. */

/* In the case of attributes, we first read the "raw" input value,
   supplying the "unset" value as the default. If a "set" value is
   obtained, we then use the appropriate (private) Set... member
   function to validate and set the value properly. */

/* First Invert flag. */
/* ------------------ */
      new->invert1 = astReadInt( channel, "inva", 0 );
      new->invert1 = ( new->invert1 != 0 );

/* Second Invert flag. */
/* ------------------- */
      new->invert2 = astReadInt( channel, "invb", 0 );
      new->invert2 = ( new->invert2 != 0 );

/* First Mapping. */
/* -------------- */
      new->map1 = astReadObject( channel, "mapa", NULL );

/* Second Mapping. */
/* --------------- */
      new->map2 = astReadObject( channel, "mapb", NULL );

/* If an error occurred, clean up by deleting the new TranMap. */
      if ( !astOK ) new = astDelete( new );
   }

/* Return the new TranMap pointer. */
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

/* None. */
