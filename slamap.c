/*
*class++
*  Name:
*     SlaMap

*  Purpose:
*     Sequence of celestial coordinate conversions.

*  Constructor Function:
c     astSlaMap (also see astSlaAdd)
f     AST_SLAMAP (also see AST_SLAADD)

*  Description:
*     An SlaMap is a specialised form of Mapping which can be used to
*     represent a sequence of conversions between standard celestial
*     (longitude, latitude) coordinate systems.
*
*     When an SlaMap is first created, it simply performs a unit
c     (null) Mapping on a pair of coordinates. Using the astSlaAdd
f     (null) Mapping on a pair of coordinates. Using the AST_SLAADD
c     function, a series of coordinate conversion steps may then be
f     routine, a series of coordinate conversion steps may then be
*     added, selected from those provided by the SLALIB Positional
*     Astronomy Library (Starlink User Note SUN/67). This allows
*     multi-step conversions between a variety of celestial coordinate
*     systems to be assembled out of the building blocks provided by
*     SLALIB.
*
*     For details of the individual coordinate conversions available,
c     see the description of the astSlaAdd function.
f     see the description of the AST_SLAADD routine.

*  Inheritance:
*     The SlaMap class inherits from the Mapping class.

*  Attributes:
*     The SlaMap class does not define any new attributes beyond those
*     which are applicable to all Mappings.

*  Functions:
c     In addition to those functions applicable to all Mappings, the
c     following function may also be applied to all SlaMaps:
f     In addition to those routines applicable to all Mappings, the
f     following routine may also be applied to all SlaMaps:
*
c     - astSlaAdd: Add a celestial coordinate conversion to an SlaMap
f     - AST_SLAADD: Add a celestial coordinate conversion to an SlaMap

*  Copyright:
*     <COPYRIGHT_STATEMENT>

*  Authors:
*     RFWS: R.F. Warren-Smith (Starlink)

*  History:
*     25-APR-1996 (RFWS):
*        Original version.
*     28-MAY-1996 (RFWS):
*        Fixed bug in argument order to slaMappa for AST__SLA_AMP case.
*     26-SEP-1996 (RFWS):
*        Added external interface and I/O facilities.
*     23-MAY-1997 (RFWS):
*        Over-ride the astMapMerge method.
*     28-MAY-1997 (RFWS):
*        Use strings to specify conversions for the public interface
*        and convert to macros (from an enumerated type) for the
*        internal representation. Tidy the public prologues.
*class--
*/

/* Module Macros. */
/* ============== */
/* Set the name of the class we are implementing. This indicates to
   the header files that define class interfaces that they should make
   "protected" symbols available. */
#define astCLASS SlaMap

/* Codes to identify SLALIB sky coordinate conversions. */
#define AST__SLA_NULL    0       /* Null value */
#define AST__SLA_ADDET   1       /* Add E-terms of aberration */
#define AST__SLA_SUBET   2       /* Subtract E-terms of aberration */
#define AST__SLA_PREBN   3       /* Bessel-Newcomb (FK4) precession */
#define AST__SLA_PREC    4       /* Apply IAU 1975 (FK5) precession model */
#define AST__SLA_FK45Z   5       /* FK4 to FK5, no proper motion or parallax */
#define AST__SLA_FK54Z   6       /* FK5 to FK4, no proper motion or parallax */
#define AST__SLA_AMP     7       /* Geocentric apparent to mean place */
#define AST__SLA_MAP     8       /* Mean place to geocentric apparent */
#define AST__SLA_ECLEQ   9       /* Ecliptic to J2000.0 equatorial */
#define AST__SLA_EQECL  10       /* Equatorial J2000.0 to ecliptic */
#define AST__SLA_GALEQ  11       /* Galactic to J2000.0 equatorial */
#define AST__SLA_EQGAL  12       /* J2000.0 equatorial to galactic */
#define AST__SLA_GALSUP 13       /* Galactic to supergalactic */
#define AST__SLA_SUPGAL 14       /* Supergalactic to galactic */

/* Maximum number of arguments required by an SLALIB conversion. */
#define MAX_SLA_ARGS 2

/* The alphabet (used for generating keywords for arguments). */
#define ALPHABET "abcdefghijklmnopqrstuvwxyz"

/* Include files. */
/* ============== */
/* Interface definitions. */
/* ---------------------- */
#include "slalib.h"              /* SLALIB interface */
#include "error.h"               /* Error reporting facilities */
#include "memory.h"              /* Memory allocation facilities */
#include "object.h"              /* Base Object class */
#include "pointset.h"            /* Sets of points/coordinates */
#include "mapping.h"             /* Coordinate Mappings (parent class) */
#include "unitmap.h"             /* Unit (null) Mappings */
#include "slamap.h"              /* Interface definition for this class */

/* Error code definitions. */
/* ----------------------- */
#include "ast_err.h"             /* AST error codes */

/* C header files. */
/* --------------- */
#include <ctype.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

/* Module Variables. */
/* ================= */
/* Define the class virtual function table and its initialisation flag
   as static variables. */
static AstSlaMapVtab class_vtab; /* Virtual function table */
static int class_init = 0;       /* Virtual function table initialised? */

/* Pointers to parent class methods which are extended by this class. */
static AstPointSet *(* parent_transform)( AstMapping *, AstPointSet *, int, AstPointSet * );

/* External Interface Function Prototypes. */
/* ======================================= */
/* The following functions have public prototypes only (i.e. no
   protected prototypes), so we must provide local prototypes for use
   within this module. */
AstSlaMap *astSlaMapId_( int, const char *, ... );

/* Prototypes for Private Member Functions. */
/* ======================================== */
static AstPointSet *Transform( AstMapping *, AstPointSet *, int, AstPointSet * );
static const char *CvtString( int, const char **, int *, const char *[ MAX_SLA_ARGS ] );
static int ChrMatch( const char *, const char * );
static int CvtCode( const char * );
static int MapMerge( AstMapping *, int, int, int *, AstMapping ***, int ** );
static void AddSlaCvt( AstSlaMap *, int, const double * );
static void Copy( const AstObject *, AstObject * );
static void Delete( AstObject * );
static void Dump( AstObject *, AstChannel * );
static void InitVtab( AstSlaMapVtab * );
static void SlaAdd( AstSlaMap *, const char *, const double[] );

/* Member functions. */
/* ================= */
static void AddSlaCvt( AstSlaMap *this, int cvttype, const double *args ) {
/*
*  Name:
*     AddSlaCvt

*  Purpose:
*     Add a coordinate conversion step to an SlaMap.

*  Type:
*     Private function.

*  Synopsis:
*     #include "slamap.h"
*     void AddSlaCvt( AstSlaMap *this, int cvttype, const double *args )

*  Class Membership:
*     SlaMap member function.

*  Description:
*     This function allows one of the sky coordinate conversions
*     supported by SLALIB to be appended to an SlaMap. When an SlaMap
*     is first created (using astSlaMap), it simply performs a unit
*     mapping. By using AddSlaCvt repeatedly, a series of sky
*     coordinate conversions may then be specified which the SlaMap
*     will subsequently perform in sequence. This allows a complex
*     coordinate conversion to be assembled out of the basic building
*     blocks provided by SLALIB. The SlaMap will also perform the
*     inverse coordinate conversion (applying the individual
*     conversion steps in reverse) if required.

*  Parameters:
*     this
*        Pointer to the SlaMap.
*     cvttype
*        A code to identify which sky coordinate conversion is to be
*        appended.  See the "SLALIB Coordinate Conversions" section
*        for details of those available.
*     args
*        Pointer to an array of double containing the argument values
*        required to fully specify the required coordinate
*        conversion. The number of arguments depends on the conversion
*        (see the "SLALIB Coordinate Conversions" section for
*        details). This value is ignored and may be NULL if no
*        arguments are required.

*  Returned Value:
*     void.

*  SLALIB Coordinate Conversions:
*     The following values may be supplied for the "cvttype" parameter
*     in order to specify the sky coordinate conversion to be
*     performed. In each case the value is named after the SLALIB
*     routine that performs the conversion, and the relevant SLALIB
*     documentation should be consulted for full details.
*
*     The argument(s) required to fully specify each conversion are
*     indicated in parentheses after each value. Values for these
*     should be given in the array pointed at by "args". The argument
*     names given match the corresponding SLALIB function arguments
*     (in the Fortran 77 documentation - SUN/67) and their values
*     should be given using the same units, time scale, calendar,
*     etc. as in SLALIB.
*
*        AST__SLA_ADDET( EQ )
*           Add E-terms of aberration.
*        AST__SLA_SUBET( EQ )
*           Subtract E-terms of aberration.
*        AST__SLA_PREBN( BEP0, BEP1 )
*           Apply Bessel-Newcomb pre-IAU 1976 (FK4) precession model.
*        AST__SLA_PREC( EP0, EP1 )
*           Apply IAU 1975 (FK5) precession model.
*        AST__SLA_FK45Z( BEPOCH )
*           Convert FK4 to FK5 (no proper motion or parallax).
*        AST__SLA_FK54Z( BEPOCH )
*           Convert FK5 to FK4 (no proper motion or parallax).
*        AST__SLA_AMP( DATE, EQ )
*           Convert geocentric apparent to mean place.
*        AST__SLA_MAP( EQ, DATE )
*           Convert mean place to geocentric apparent.
*        AST__SLA_ECLEQ( DATE )
*           Convert ecliptic coordinates to J2000.0 equatorial.
*        AST__SLA_EQECL( DATE )
*           Convert equatorial J2000.0 to ecliptic coordinates.
*        AST__SLA_GALEQ( )
*           Convert galactic coordinates to J2000.0 equatorial.
*        AST__SLA_EQGAL( )
*           Convert J2000.0 equatorial to galactic coordinates.
*        AST__SLA_GALSUP( )
*           Convert galactic to supergalactic coordinates.
*        AST__SLA_SUPGAL( )
*           Convert supergalactic coordinates to galactic.

*  Notes:
*     - The specified conversion is appended only if the SlaMap's
*     Invert attribute is zero. If it is non-zero, this function
*     effectively prefixes the inverse of the conversion specified
*     instead.
*     - Sky coordinate values are in radians (as for SLALIB) and all
*     conversions are performed using double arithmetic.
*/

/* Local Variables: */
   const char *argdesc[ MAX_SLA_ARGS ]; /* Pointers to argument descriptions */
   const char *comment;          /* Pointer to comment string */
   const char *cvt_string;       /* Pointer to conversion type string */
   int nargs;                    /* Number of arguments */
   int ncvt;                     /* Number of coordinate conversions */

/* Check the global error status. */
   if ( !astOK ) return;

/* Validate the coordinate conversion type and obtain the number of
   required arguments. */
   cvt_string = CvtString( cvttype, &comment, &nargs, argdesc );

/* If the sky coordinate conversion type was not valid, then report an
   error. */
   if ( astOK && !cvt_string ) {
      astError( AST__SLAIN,
                "Invalid SLALIB sky coordinate conversion type (%d).",
                astGetClass( this ), (int) cvttype );
   }

/* Note the number of coordinate conversions already stored in the SlaMap. */
   if ( astOK ) {
      ncvt = this->ncvt;

/* Extend the array of conversion types and the array of pointers to
   their argument lists to accommodate the new one. */
      this->cvttype = (int *) astGrow( this->cvttype, ncvt + 1,
                                       sizeof( int ) );
      this->cvtargs = (double **) astGrow( this->cvtargs, ncvt + 1,
                                           sizeof( double * ) );

/* If OK, allocate memory and store a copy of the argument list,
   putting a pointer to the copy into the SlaMap. */
      if ( astOK ) {
         this->cvtargs[ ncvt ] = astStore( NULL, args,
                                           sizeof( double ) * (size_t) nargs );
      }

/* Store the conversion type and increment the conversion count. */
      if ( astOK ) {
         this->cvttype[ ncvt ] = cvttype;
         this->ncvt++;
      }
   }
}

static int ChrMatch( const char *str1, const char *str2 ) {
/*
*  Name:
*     ChrMatch

*  Purpose:
*     Case insensitive string comparison.

*  Type:
*     Private function.

*  Synopsis:
*     #include "slamap.h"
*     int ChrMatch( const char *str1, const char *str2 )

*  Class Membership:
*     SlaMap member function.

*  Description:
*     This function compares two null terminated strings for equality,
*     discounting differences in case and any trailing white space in
*     either string.

*  Parameters:
*     str1
*        Pointer to the first string.
*     str2
*        Pointer to the second string.

*  Returned Value:
*     Non-zero if the two strings match, otherwise zero.

*  Notes:
*     - A value of zero is returned if this function is invoked with
*     the global error status set or if it should fail for any reason.
*/

/* Local Variables: */
   int match;                    /* Strings match? */

/* Check the global error status. */
   if ( !astOK ) return 0;

/* Initialise. */
   match = 1;

/* Loop to compare characters in the two strings until a mis-match
   occurs or we reach the end of the longer string. */
   while ( match && ( *str1 || *str2 ) ) {

/* Two characters match if (a) we are at the end of one string and the
   other string contains white space or (b) both strings contain the
   same character when converted to lower case. */
      match = ( !*str1 && isspace( *str2 ) ) ||
              ( !*str2 && isspace( *str1 ) ) ||
              ( tolower( *str1 ) == tolower( *str2 ) );

/* Step through each string a character at a time until its end is
   reached. */
      if ( *str1 ) str1++;
      if ( *str2 ) str2++;
   }

/* Return the result. */
   return match;
}
   
static int CvtCode( const char *cvt_string ) {
/*
*  Name:
*     CvtCode

*  Purpose:
*     Convert a conversion type from a string representation to a code value.

*  Type:
*     Private function.

*  Synopsis:
*     #include "slamap.h"
*     int CvtCode( const char *cvt_string )

*  Class Membership:
*     SlaMap member function.

*  Description:
*     This function accepts a string used to repersent one of the
*     SLALIB sky coordinate conversions and converts it into a code
*     value for internal use.

*  Parameters:
*     cvt_string
*        Pointer to a constant null-terminated string representing a
*        sky coordinate conversion. This is case sensitive and should
*        contain no unnecessary white space.

*  Returned Value:
*     The equivalent conversion code. If the string was not
*     recognised, the code AST__SLA_NULL is returned, without error.

*  Notes:
*     - A value of AST__SLA_NULL will be returned if this function is
*     invoked with the global error status set, or if it should fail
*     for any reason.
*/

/* Local Variables: */
   int result;                   /* Result value to return */

/* Initialise. */
   result = AST__SLA_NULL;

/* Check the global error status. */
   if ( !astOK ) return result;

/* Test the string against each recognised value in turn and assign
   the result. */
   if ( ChrMatch( cvt_string, "ADDET" ) ) {
      result = AST__SLA_ADDET;

   } else if ( ChrMatch( cvt_string, "SUBET" ) ) {
      result = AST__SLA_SUBET;

   } else if ( ChrMatch( cvt_string, "PREBN" ) ) {
      result = AST__SLA_PREBN;

   } else if ( ChrMatch( cvt_string, "PREC" ) ) {
      result = AST__SLA_PREC;

   } else if ( ChrMatch( cvt_string, "FK45Z" ) ) {
      result = AST__SLA_FK45Z;

   } else if ( ChrMatch( cvt_string, "FK54Z" ) ) {
      result = AST__SLA_FK54Z;

   } else if ( ChrMatch( cvt_string, "AMP" ) ) {
      result = AST__SLA_AMP;

   } else if ( ChrMatch( cvt_string, "MAP" ) ) {
      result = AST__SLA_MAP;

   } else if ( ChrMatch( cvt_string, "ECLEQ" ) ) {
      result = AST__SLA_ECLEQ;

   } else if ( ChrMatch( cvt_string, "EQECL" ) ) {
      result = AST__SLA_EQECL;

   } else if ( ChrMatch( cvt_string, "GALEQ" ) ) {
      result = AST__SLA_GALEQ;

   } else if ( ChrMatch( cvt_string, "EQGAL" ) ) {
      result = AST__SLA_EQGAL;

   } else if ( ChrMatch( cvt_string, "GALSUP" ) ) {
      result = AST__SLA_GALSUP;

   } else if ( ChrMatch( cvt_string, "SUPGAL" ) ) {
      result = AST__SLA_SUPGAL;
   }

/* Return the result. */
   return result;
}

static const char *CvtString( int cvt_code, const char **comment,
                              int *nargs, const char *arg[ MAX_SLA_ARGS ] ) {
/*
*  Name:
*     CvtString

*  Purpose:
*     Convert a conversion type from a code value to a string representation.

*  Type:
*     Private function.

*  Synopsis:
*     #include "slamap.h"
*     const char *CvtString( int cvt_code, const char **comment,
*                            int *nargs, const char *arg[ MAX_SLA_ARGS ] )

*  Class Membership:
*     SlaMap member function.

*  Description:
*     This function accepts a code value used to represent one of the
*     SLALIB sky coordinate conversions and converts it into an
*     equivalent string representation. It also returns a descriptive
*     comment and information about the arguments required in order to
*     perform the conversion.

*  Parameters:
*     cvt_code
*        The conversion code.
*     comment
*        Address of a location to return a pointer to a constant
*        null-terminated string containing a description of the
*        conversion.
*     nargs
*        Address of an int in which to return the number of arguments
*        required in order to perform the conversion (may be zero).
*     arg
*        An array in which to return a pointer to a constant
*        null-terminated string for each argument (above) containing a
*        description of what each argument represents.

*  Returned Value:
*     Pointer to a constant null-terminated string representation of
*     the conversion code value supplied. If the code supplied is not
*     valid, a NULL pointer will be returned, without error.

*  Notes:
*     - A NULL pointer value will be returned if this function is
*     invoked with the global error status set, or if it should fail
*     for any reason.
*/

/* Local Variables: */
   const char *result;           /* Result pointer to return */

/* Initialise the returned values. */
   *comment = NULL;
   *nargs = 0;
   result = NULL;

/* Check the global error status. */
   if ( !astOK ) return result;
      
/* Test for each valid code value in turn and assign the appropriate
   return values. */
   switch ( cvt_code ) {

   case AST__SLA_ADDET:
      result = "ADDET";
      *comment = "Add E-terms of aberration";
      *nargs = 1;
      arg[ 0 ] = "Besselian epoch of mean equinox (FK4)";
      break;

   case AST__SLA_SUBET:
      result = "SUBET";
      *comment = "Subtract E-terms of aberration";
      *nargs = 1;
      arg[ 0 ] = "Besselian epoch of mean equinox (FK4)";
      break;

   case AST__SLA_PREBN:
      result = "PREBN";
      *comment = "Apply Bessel-Newcomb (FK4) precession";
      *nargs = 2;
      arg[ 0 ] = "From Besselian epoch";
      arg[ 1 ] = "To Besselian epoch";
      break;

   case AST__SLA_PREC:
      result = "PREC";
      *comment = "Apply IAU 1975 (FK5) precession";
      *nargs = 2;
      arg[ 0 ] = "From Julian epoch";
      arg[ 1 ] = "To Julian epoch";
      break;

   case AST__SLA_FK45Z:
      result = "FK45Z";
      *comment = "FK4 to FK5 J2000.0 (no PM or parallax)";
      arg[ 0 ] = "Besselian epoch of FK4 coordinates";
      *nargs = 1;
      break;

   case AST__SLA_FK54Z:
      result = "FK54Z";
      *comment = "FK5 J2000.0 to FK4 (no PM or parallax)";
      *nargs = 1;
      arg[ 0 ] = "Besselian epoch of FK4 system";
      break;

   case AST__SLA_AMP:
      result = "AMP";
      *comment = "Geocentric apparent to mean place (FK5)";
      *nargs = 2;
      arg[ 0 ] = "TDB of apparent place (as MJD)";
      arg[ 1 ] = "Julian epoch of mean equinox (FK5)";
      break;

   case AST__SLA_MAP:
      result = "MAP";
      *comment = "Mean place (FK5) to geocentric apparent";
      *nargs = 2;
      arg[ 0 ] = "Julian epoch of mean equinox (FK5)";
      arg[ 1 ] = "TDB of apparent place (as MJD)";
      break;

   case AST__SLA_ECLEQ:
      result = "ECLEQ";
      *comment = "Ecliptic (IAU 1980) to J2000.0 equatorial (FK5)";
      *nargs = 1;
      arg[ 0 ] = "TDB of mean ecliptic (as MJD)";
      break;

   case AST__SLA_EQECL:
      result = "EQECL";
      *comment = "Equatorial J2000.0 (FK5) to ecliptic (IAU 1980)";
      *nargs = 1;
      arg[ 0 ] = "TDB of mean ecliptic (as MJD)";
      break;

   case AST__SLA_GALEQ:
      result = "GALEQ";
      *comment = "Galactic (IAU 1958) to J2000.0 equatorial (FK5)";
      *nargs = 0;
      break;

   case AST__SLA_EQGAL:
      result = "EQGAL";
      *comment = "J2000.0 equatorial (FK5) to galactic (IAU 1958)";
      *nargs = 0;
      break;

   case AST__SLA_GALSUP:
      result = "GALSUP";
      *comment = "Galactic (IAU 1958) to supergalactic";
      *nargs = 0;
      break;

   case AST__SLA_SUPGAL:
      result = "SUPGAL";
      *comment = "Supergalactic to galactic (IAU 1958)";
      *nargs = 0;
      break;
   }

/* Return the result. */
   return result;
}

static void InitVtab( AstSlaMapVtab *vtab ) {
/*
*  Name:
*     InitVtab

*  Purpose:
*     Initialise a virtual function table for a SlaMap.

*  Type:
*     Private function.

*  Synopsis:
*     #include "slamap.h"
*     void InitVtab( AstSlaMapVtab *vtab )

*  Class Membership:
*     SlaMap member function.

*  Description:
*     This function initialises the component of a virtual function
*     table which is used by the SlaMap class.

*  Parameters:
*     vtab
*        Pointer to the virtual function table. The components used by
*        all ancestral classes should already have been initialised.
*/

/* Local Variables: */
   AstMappingVtab *mapping;      /* Pointer to Mapping component of Vtab */

/* Check the local error status. */
   if ( !astOK ) return;

/* Store a unique "magic" value in the virtual function table. This
   will be used (by astIsASlaMap) to determine if an object belongs to
   this class.  We can conveniently use the address of the (static)
   class_init variable to generate this unique value. */
   vtab->check = &class_init;

/* Initialise member function pointers. */
/* ------------------------------------ */
/* Store pointers to the member functions (implemented here) that
   provide virtual methods for this class. */
   vtab->SlaAdd = SlaAdd;

/* Save the inherited pointers to methods that will be extended, and
   replace them with pointers to the new member functions. */
   mapping = (AstMappingVtab *) vtab;

   parent_transform = mapping->Transform;
   mapping->Transform = Transform;

/* Store replacement pointers for methods which will be over-ridden by
   new member functions implemented here. */
   mapping->MapMerge = MapMerge;

/* Declare the copy constructor, destructor and class dump
   function. */
   astSetCopy( vtab, Copy );
   astSetDelete( vtab, Delete );
   astSetDump( vtab, Dump, "SlaMap",
               "Conversion between sky coordinate systems" );
}

static int MapMerge( AstMapping *this, int where, int series, int *nmap,
                     AstMapping ***map_list, int **invert_list ) {
/*
*  Name:
*     MapMerge

*  Purpose:
*     Simplify a sequence of Mappings containing an SlaMap.

*  Type:
*     Private function.

*  Synopsis:
*     #include "mapping.h"
*     int MapMerge( AstMapping *this, int where, int series, int *nmap,
*                   AstMapping ***map_list, int **invert_list )

*  Class Membership:
*     SlaMap method (over-rides the protected astMapMerge method
*     inherited from the Mapping class).

*  Description:
*     This function attempts to simplify a sequence of Mappings by
*     merging a nominated SlaMap in the sequence with its neighbours,
*     so as to shorten the sequence if possible.
*
*     In many cases, simplification will not be possible and the
*     function will return -1 to indicate this, without further
*     action.
*
*     In most cases of interest, however, this function will either
*     attempt to replace the nominated SlaMap with one which it
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
*        Pointer to the nominated SlaMap which is to be merged with
*        its neighbours. This should be a cloned copy of the SlaMap
*        pointer contained in the array element "(*map_list)[where]"
*        (see below). This pointer will not be annulled, and the
*        SlaMap it identifies will not be modified by this function.
*     where
*        Index in the "*map_list" array (below) at which the pointer
*        to the nominated SlaMap resides.
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
   AstMapping *new;              /* Pointer to replacement Mapping */
   AstSlaMap *slamap;            /* Pointer to SlaMap */
   const char *argdesc[ MAX_SLA_ARGS ]; /* Argument descriptions (junk) */
   const char *class;            /* Pointer to Mapping class string */
   const char *comment;          /* Pointer to comment string (junk) */
   double (*cvtargs)[ MAX_SLA_ARGS ]; /* Pointer to argument arrays */
   int *cvttype;                 /* Pointer to transformation type codes */
   int *narg;                    /* Pointer to argument count array */
   int done;                     /* Finished (no further simplification)? */
   int iarg;                     /* Loop counter for arguments */
   int icvt1;                    /* Loop initial value */
   int icvt2;                    /* Loop final value */
   int icvt;                     /* Loop counter for transformation steps */
   int ikeep;                    /* Index to store step being kept */
   int imap1;                    /* Index of first SlaMap to merge */
   int imap2;                    /* Index of last SlaMap to merge */
   int imap;                     /* Loop counter for Mappings */
   int inc;                      /* Increment for transformation step loop */
   int invert;                   /* SlaMap applied in inverse direction? */
   int istep;                    /* Loop counter for transformation steps */
   int keep;                     /* Keep transformation step? */
   int ngone;                    /* Number of Mappings eliminated */
   int nstep0;                   /* Original number of transformation steps */
   int nstep;                    /* Total number of transformation steps */
   int result;                   /* Result value to return */
   int simpler;                  /* Simplification possible? */
   int unit;                     /* Replacement Mapping is a UnitMap? */

/* Initialise. */
   result = -1;

/* Check the global error status. */
   if ( !astOK ) return result;

/* SlaMaps can only be merged if they are in series (or if there is
   only one Mapping present, in which case it makes no difference), so
   do nothing if they are not. */
   if ( series || ( *nmap == 1 ) ) {

/* Initialise the number of transformation steps to be merged to equal
   the number in the nominated SlaMap. */
      nstep = ( (AstSlaMap *) ( *map_list )[ where ] )->ncvt;

/* Search adjacent lower-numbered Mappings until one is found which is
   not an SlaMap. Accumulate the number of transformation steps
   involved in any SlaMaps found. */
      imap1 = where;
      while ( ( imap1 - 1 >= 0 ) && astOK ) {
         class = astGetClass( ( *map_list )[ imap1 - 1 ] );
         if ( !astOK || strcmp( class, "SlaMap" ) ) break;
         nstep += ( (AstSlaMap *) ( *map_list )[ imap1 - 1 ] )->ncvt;
         imap1--;
      }

/* Similarly search adjacent higher-numbered Mappings. */
      imap2 = where;
      while ( ( imap2 + 1 < *nmap ) && astOK ) {
         class = astGetClass( ( *map_list )[ imap2 + 1 ] );
         if ( !astOK || strcmp( class, "SlaMap" ) ) break;
         nstep += ( (AstSlaMap *) ( *map_list )[ imap2 + 1 ] )->ncvt;
         imap2++;
      }

/* Remember the initial number of transformation steps. */
      nstep0 = nstep;

/* Allocate memory for accumulating a list of all the transformation
   steps involved in all the SlaMaps found. */
      cvttype = astMalloc( sizeof( int ) * (size_t) nstep );
      cvtargs = astMalloc( sizeof( double[ MAX_SLA_ARGS ] ) * (size_t) nstep );
      narg = astMalloc( sizeof( int ) * (size_t) nstep );

/* Loop to obtain the transformation data for each SlaMap being merged. */
      nstep = 0;
      for ( imap = imap1; astOK && ( imap <= imap2 ); imap++ ) {

/* Obtain a pointer to the SlaMap and note if it is being applied in
   its inverse direction. */
         slamap = (AstSlaMap *) ( *map_list )[ imap ];
         invert = ( *invert_list )[ imap ];

/* Set up loop limits and an increment to scan the transformation
   steps in each SlaMap in either the forward or reverse direction, as
   dictated by the associated "invert" value. */
         icvt1 = invert ? slamap->ncvt - 1 : 0;
         icvt2 = invert ? -1 : slamap->ncvt;
         inc = invert ? -1 : 1;

/* Loop through each transformation step in the SlaMap. */
         for ( icvt = icvt1; icvt != icvt2; icvt += inc ) {

/* Store the transformation type code and use "CvtString" to determine
   the associated number of arguments. Then store these arguments. */
            cvttype[ nstep ] = slamap->cvttype[ icvt ];
            (void) CvtString( cvttype[ nstep ], &comment, narg + nstep,
                              argdesc );
            if ( !astOK ) break;
            for ( iarg = 0; iarg < narg[ nstep ]; iarg++ ) {
               cvtargs[ nstep ][ iarg ] = slamap->cvtargs[ icvt ][ iarg ];
            }

/* If the SlaMap is inverted, we must not only accumulate its
   transformation steps in reverse, but also apply them in
   reverse. For some steps this means swapping arguments, for some it
   means changing the transformation type code to a complementary
   value, and for others it means both.  Define macros to perform each
   of these changes. */

/* Macro to swap the values of two nominated arguments if the
   transformation type code matches "code". */
#define SWAP_ARGS( code, arg1, arg2 ) \
            if ( cvttype[ nstep ] == code ) { \
               double tmp = cvtargs[ nstep ][ arg1 ]; \
               cvtargs[ nstep ][ arg1 ] = cvtargs[ nstep ][ arg2 ]; \
               cvtargs[ nstep ][ arg2 ] = tmp; \
            }

/* Macro to exchange a transformation type code for its inverse (and
   vice versa). */
#define SWAP_CODES( code1, code2 ) \
            if ( cvttype[ nstep ] == code1 ) { \
               cvttype[ nstep ] = code2; \
            } else if ( cvttype[ nstep ] == code2 ) { \
               cvttype[ nstep ] = code1; \
            }

/* Use these macros to apply the changes where needed. */
            if ( invert ) {

/* E-terms of aberration. */
/* ---------------------- */
/* Exchange addition and subtraction of E-terms. */
               SWAP_CODES( AST__SLA_ADDET, AST__SLA_SUBET )

/* Bessel-Newcomb pre-IAU 1976 (FK4) precession model. */
/* --------------------------------------------------- */
/* Exchange the starting and ending Besselian epochs. */
               SWAP_ARGS( AST__SLA_PREBN, 0, 1 )

/* IAU 1975 (FK5) precession model. */
/* -------------------------------- */
/* Exchange the starting and ending epochs. */
               SWAP_ARGS( AST__SLA_PREC, 0, 1 )

/* FK4 to FK5 (no proper motion or parallax). */
/* ------------------------------------------ */
/* Exchange FK5 to FK4 conversion for its inverse, and vice versa. */
               SWAP_CODES( AST__SLA_FK54Z, AST__SLA_FK45Z )

/* Geocentric apparent to mean place. */
/* ---------------------------------- */
/* Exchange the transformation code for its inverse and also exchange
   the order of the date and equinox arguments. */
               SWAP_CODES( AST__SLA_AMP, AST__SLA_MAP )
               SWAP_ARGS( AST__SLA_AMP, 0, 1 )
               SWAP_ARGS( AST__SLA_MAP, 0, 1 )

/* Ecliptic coordinates to J2000.0 equatorial. */
/* ------------------------------------------- */
/* Exchange the transformation code for its inverse. */
               SWAP_CODES( AST__SLA_ECLEQ, AST__SLA_EQECL )

/* Galactic coordinates to J2000.0 equatorial. */
/* ------------------------------------------- */
/* Exchange the transformation code for its inverse. */
               SWAP_CODES( AST__SLA_GALEQ, AST__SLA_EQGAL )

/* Galactic to supergalactic coordinates. */
/* -------------------------------------- */
/* Exchange the transformation code for its inverse. */
               SWAP_CODES( AST__SLA_GALSUP, AST__SLA_SUPGAL )
            }

/* Undefine the local macros. */
#undef SWAP_ARGS
#undef SWAP_CODES

/* Count the transformation steps. */
            nstep++;
         }
      }

/* Loop to simplify the sequence of transformation steps until no
   further improvement is possible. */
      done = 0;
      while ( astOK && !done ) {

/* Examine each remaining transformation step in turn.  */
         ikeep = -1;
         for ( istep = 0; istep < nstep; istep++ ) {

/* Initially assume we will retain the current step. */
            keep = 1;

/* Eliminate redundant precession corrections. */
/* ------------------------------------------- */
/* First check if this is a redundant precession transformation
   (i.e. the starting and ending epochs are the same). If so, then
   note that it should not be kept. */
            if ( ( ( cvttype[ istep ] == AST__SLA_PREBN ) ||
                   ( cvttype[ istep ] == AST__SLA_PREC ) ) &&
                 ( cvtargs[ istep ][ 0 ] == cvtargs[ istep ][ 1 ] ) ) {
               keep = 0;

/* The remaining simplifications act to combine adjacent
   transformation steps, so only apply them while there are at least 2
   steps left. */
            } else if ( istep < ( nstep - 1 ) ) {

/* Define a macro to test if two adjacent transformation type codes
   have specified values. */
#define PAIR_CVT( code1, code2 ) \
               ( ( cvttype[ istep ] == code1 ) && \
                 ( cvttype[ istep + 1 ] == code2 ) )

/* Combine adjacent precession corrections. */
/* ---------------------------------------- */
/* If two precession corrections are adjacent, and have an equinox
   value in common, then they may be combined into a single correction
   by eliminating the common equinox. */
               if ( ( PAIR_CVT( AST__SLA_PREBN, AST__SLA_PREBN ) ||
                      PAIR_CVT( AST__SLA_PREC, AST__SLA_PREC ) ) &&
                    ( cvtargs[ istep ][ 1 ] == cvtargs[ istep + 1 ][ 0 ] ) ) {

/* Retain the second correction, changing its first argument, and
   eliminate the first correction. */
                  cvtargs[ istep + 1 ][ 0 ] = cvtargs[ istep ][ 0 ];
                  istep++;

/* Eliminate redundant E-term handling. */
/* ------------------------------------ */
/* Check if adjacent steps implement a matching pair of corrections
   for the E-terms of aberration with the same argument value. If so,
   they will cancel, so eliminate them both. */
               } else if ( ( PAIR_CVT( AST__SLA_SUBET, AST__SLA_ADDET ) ||
                             PAIR_CVT( AST__SLA_ADDET, AST__SLA_SUBET ) ) &&
                           ( cvtargs[ istep ][ 0 ] ==
                             cvtargs[ istep + 1 ][ 0 ] ) ) {
                  istep++;
                  keep = 0;

/* Eliminate redundant FK4/FK5 conversions. */
/* ---------------------------------------- */
/* Similarly, check for a matching pair of FK4/FK5 conversions with
   the same argument value and eliminate them both if possible. */
               } else if ( ( PAIR_CVT( AST__SLA_FK45Z, AST__SLA_FK54Z ) ||
                             PAIR_CVT( AST__SLA_FK54Z, AST__SLA_FK45Z ) ) &&
                           ( cvtargs[ istep ][ 0 ] ==
                             cvtargs[ istep + 1 ][ 0 ] ) ) {
                  istep++;
                  keep = 0;

/* Eliminate redundant geocentric apparent conversions. */
/* ---------------------------------------------------- */
/* As above, check for a matching pair of conversions with matching
   argument values (note the argument order reverses for the two
   directions) and eliminate them if possible. */
               } else if ( ( PAIR_CVT( AST__SLA_AMP, AST__SLA_MAP ) ||
                             PAIR_CVT( AST__SLA_MAP, AST__SLA_AMP ) ) &&
                           ( cvtargs[ istep ][ 0 ] ==
                             cvtargs[ istep + 1 ][ 1 ] ) &&
                           ( cvtargs[ istep ][ 1 ] ==
                             cvtargs[ istep + 1 ][ 0 ] ) ) {
                  istep++;
                  keep = 0;

/* Eliminate redundant ecliptic coordinate conversions. */
/* ---------------------------------------------------- */
/* This is handled in the same way as the FK4/FK5 case. */
               } else if ( ( PAIR_CVT( AST__SLA_ECLEQ, AST__SLA_EQECL ) ||
                             PAIR_CVT( AST__SLA_EQECL, AST__SLA_ECLEQ ) ) &&
                           ( cvtargs[ istep ][ 0 ] ==
                             cvtargs[ istep + 1 ][ 0 ] ) ) {
                  istep++;
                  keep = 0;

/* Eliminate redundant galactic coordinate conversions. */
/* ---------------------------------------------------- */
/* This is handled as above, except that there are no arguments to
   check. */
               } else if ( PAIR_CVT( AST__SLA_GALEQ, AST__SLA_EQGAL ) ||
                           PAIR_CVT( AST__SLA_EQGAL, AST__SLA_GALEQ ) ) {
                  istep++;
                  keep = 0;

/* Eliminate redundant supergalactic coordinate conversions. */
/* --------------------------------------------------------- */
/* This is handled as above. */
               } else if ( PAIR_CVT( AST__SLA_GALSUP, AST__SLA_SUPGAL ) ||
                           PAIR_CVT( AST__SLA_SUPGAL, AST__SLA_GALSUP ) ) {
                  istep++;
                  keep = 0;
               }

/* Undefine the local macro. */
#undef PAIR_CVT
            }

/* If the current transformation (possibly modified above) is being
   kept, then increment the index that identifies its new location in
   the list of transformation steps. */
            if ( keep ) {
               ikeep++;

/* If the new location is different to its current location, copy the
   transformation data into the new location. */
               if ( ikeep != istep ) {
                  cvttype[ ikeep ] = cvttype[ istep ];
                  for ( iarg = 0; iarg < narg[ istep ]; iarg++ ) {
                     cvtargs[ ikeep ][ iarg ] = cvtargs[ istep ][ iarg ];
                  }
                  narg[ ikeep ] = narg[ istep ];
               }
            }
         }

/* Note if no simplification was achieved on this iteration (i.e. the
   number of transformation steps was not reduced). This is the signal
   to quit. */
         done = ( ( ikeep + 1 ) >= nstep );

/* Note how many transformation steps now remain. */
         nstep = ikeep + 1;
      }

/* Determine how many Mappings can be eliminated by condensing all
   those considered above into a single Mapping. */
      if ( astOK ) {
         ngone = imap2 - imap1;

/* Determine if the replacement Mapping can be a UnitMap (a null
   Mapping). This will only be the case if all the transformation
   steps were eliminated above. */
         unit = ( nstep == 0 );

/* Determine if simplification is possible. This will be the case if
   (a) Mappings were eliminated ("ngone" is non-zero), or (b) the
   number of transformation steps was reduced, or (c) the SlaMap(s)
   can be replaced by a UnitMap, or (d) if there was initially only
   one SlaMap present, its invert flag was set (this flag will always
   be cleared in the replacement Mapping). */
         simpler = ngone || ( nstep < nstep0 ) || unit ||
                   ( *invert_list )[ where ];

/* Do nothing more unless simplification is possible. */
         if ( simpler ) {

/* If the replacement Mapping is a UnitMap, then create it. */
            if ( unit ) {
               new = (AstMapping *)
                        astUnitMap( astGetNin( ( *map_list )[ where ] ), "" );

/* Otherwise, create a replacement SlaMap and add each of the
   remaining transformation steps to it. */
            } else {
               new = (AstMapping *) astSlaMap( 0, "" );
               for ( istep = 0; istep < nstep; istep++ ) {
                  AddSlaCvt( (AstSlaMap *) new, cvttype[ istep ],
                             cvtargs[ istep ] );
               }
            }

/* Annul the pointers to the Mappings being eliminated. */
            if ( astOK ) {
               for ( imap = imap1; imap <= imap2; imap++ ) {
                  ( *map_list )[ imap ] = astAnnul( ( *map_list )[ imap ] );
               }

/* Insert the pointer and invert value for the new Mapping. */
               ( *map_list )[ imap1 ] = new;
               ( *invert_list )[ imap1 ] = 0;

/* Move any subsequent Mapping information down to close the gap. */
               for ( imap = imap2 + 1; imap < *nmap; imap++ ) {
                  ( *map_list )[ imap - ngone ] = ( *map_list )[ imap ];
                  ( *invert_list )[ imap - ngone ] = ( *invert_list )[ imap ];
               }

/* Blank out any information remaining at the end of the arrays. */
               for ( imap = ( *nmap - ngone ); imap < *nmap; imap++ ) {
                  ( *map_list )[ imap ] = NULL;
                  ( *invert_list )[ imap ] = 0;
               }

/* Decrement the Mapping count and return the index of the first
   Mapping which was eliminated. */
               ( *nmap ) -= ngone;
               result = imap1;

/* If an error occurred, annul the new Mapping pointer. */
            } else {
               new = astAnnul( new );
            }
         }
      }

/* Free the memory used for the transformation steps. */
      cvttype = astFree( cvttype );
      cvtargs = astFree( cvtargs );
      narg = astFree( narg );
   }

/* If an error occurred, clear the returned value. */
   if ( !astOK ) result = -1;

/* Return the result. */
   return result;
}

static void SlaAdd( AstSlaMap *this, const char *cvt, const double args[] ) {
/*
*++
*  Name:
c     astSlaAdd
f     AST_SLAADD

*  Purpose:
*     Add a celestial coordinate conversion to an SlaMap.

*  Type:
*     Public virtual function.

*  Synopsis:
c     #include "slamap.h"
c     void astSlaAdd( AstSlaMap *this, const char *cvt, const double args[] )
f     CALL AST_SLAADD( THIS, CVT, ARGS, STATUS )

*  Class Membership:
*     SlaMap method.

*  Description:
c     This function adds one of the standard celestial coordinate
f     This routine adds one of the standard celestial coordinate
*     system conversions provided by the SLALIB Positional Astronomy
*     Library (Starlink User Note SUN/67) to an existing SlaMap.
*
c     When an SlaMap is first created (using astSlaMap), it simply
f     When an SlaMap is first created (using AST_SLAMAP), it simply
c     performs a unit (null) Mapping. By using astSlaAdd (repeatedly
f     performs a unit (null) Mapping. By using AST_SLAADD (repeatedly
*     if necessary), one or more coordinate conversion steps may then
*     be added, which the SlaMap will perform in sequence. This allows
*     multi-step conversions between a variety of celestial coordinate
*     systems to be assembled out of the building blocks provided by
*     SLALIB.
*
*     Normally, if an SlaMap's Invert attribute is zero (the default),
*     then its forward transformation is performed by carrying out
*     each of the individual coordinate conversions specified by
c     astSlaAdd in the order given (i.e. with the most recently added
f     AST_SLAADD in the order given (i.e. with the most recently added
*     conversion applied last).
*
*     This order is reversed if the SlaMap's Invert attribute is
*     non-zero (or if the inverse transformation is requested by any
*     other means) and each individual coordinate conversion is also
*     replaced by its own inverse. This process inverts the overall
*     effect of the SlaMap. In this case, the first conversion to be
*     applied would be the inverse of the one most recently added.

*  Parameters:
c     this
f     THIS = INTEGER (Given)
*        Pointer to the SlaMap.
c     cvt
f     CVT = CHARACTER * ( * ) (Given)
c        Pointer to a null-terminated string which identifies the
f        A character string which identifies the
*        celestial coordinate conversion to be added to the
*        SlaMap. See the "SLALIB Conversions" section for details of
*        those available.
c     args
f     ARGS( * ) = DOUBLE PRECISION (Given)
*        An array containing argument values for the celestial
*        coordinate conversion. The number of arguments required, and
*        hence the number of array elements used, depends on the
*        conversion specified (see the "SLALIB Conversions"
*        section). This array is ignored
c        and a NULL pointer may be supplied
*        if no arguments are needed.
f     STATUS = INTEGER (Given and Returned)
f        The global status.

*  Notes:
*     - All coordinate values processed by an SlaMap are in
*     radians. The first coordinate is the celestial longitude and the
*     second coordinate is the celestial latitude.
*     - When assembling a multi-stage conversion, it can sometimes be
*     difficult to determine the most economical conversion path. For
*     example, converting to the standard FK5 coordinate system as an
*     intermediate stage is often sensible in formulating the problem,
*     but may introduce unnecessary extra conversion steps. A solution
*     to this is to include all the steps which are (logically)
c     necessary, but then to use astSimplify to simplify the resulting
f     necessary, but then to use AST_SIMPLIFY to simplify the resulting
*     SlaMap. The simplification process will eliminate any steps
*     which turn out not to be needed.
c     - This function does not check to ensure that the sequence of
f     - This routine does not check to ensure that the sequence of
*     coordinate conversions added to an SlaMap is physically
*     meaningful.

*  SLALIB Conversions:
*     The following strings (which are case-insensitive) may be supplied
c     via the "cvt" parameter to indicate which celestial coordinate
f     via the CVT argument to indicate which celestial coordinate
*     conversion is to be added to the SlaMap. Each string is derived
*     from the name of the SLALIB routine that performs the
*     conversion and the relevant documentation (SUN/67) should be
*     consulted for details.  Where arguments are needed by
*     the conversion, they are listed in parentheses. Values for
c     these arguments should be given, via the "args" array, in the
f     these arguments should be given, via the ARGS array, in the
*     order indicated. The argument names match the corresponding
*     SLALIB routine arguments and their values should be given using
*     exactly the same units, time scale, calendar, etc. as described
*     in SUN/67:
*
*     - "ADDET" (EQ): Add E-terms of aberration.
*     - "SUBET" (EQ): Subtract E-terms of aberration.
*     - "PREBN" (BEP0,BEP1): Apply Bessel-Newcomb pre-IAU 1976 (FK4)
*     precession model.
*     - "PREC" (EP0,EP1): Apply IAU 1975 (FK5) precession model.
*     - "FK45Z" (BEPOCH): Convert FK4 to FK5 (no proper motion or parallax).
*     - "FK54Z" (BEPOCH): Convert FK5 to FK4 (no proper motion or parallax).
*     - "AMP" (DATE,EQ): Convert geocentric apparent to mean place.
*     - "MAP" (EQ,DATE): Convert mean place to geocentric apparent.
*     - "ECLEQ" (DATE): Convert ecliptic coordinates to J2000.0 equatorial.
*     - "EQECL" (DATE): Convert equatorial J2000.0 to ecliptic coordinates.
*     - "GALEQ": Convert galactic coordinates to J2000.0 equatorial.
*     - "EQGAL": Convert J2000.0 equatorial to galactic coordinates.
*     - "GALSUP": Convert galactic to supergalactic coordinates.
*     - "SUPGAL": Convert supergalactic coordinates to galactic.
*
*     For example, to use the "ADDET" conversion, which takes a single
*     argument EQ, you should consult the documentation for the SLALIB
*     routine SLA_ADDET. This describes the conversion in detail and
*     shows that EQ is the Besselian epoch of the mean equator and
*     equinox.
c     This value should then be supplied to astSlaAdd in args[0].
f     This value should then be supplied to AST_SLAADD in ARGS(1).
*--
*/

/* Local Variables: */
   int cvttype;                  /* Conversion type code */

/* Check the inherited status. */
   if ( !astOK ) return;

/* Validate the type string supplied and obtain the equivalent
   conversion type code. */
   cvttype = CvtCode( cvt );

/* If the string was not recognised, then report an error. */
   if ( astOK && ( cvttype == AST__SLA_NULL ) ) {
      astError( AST__SLAIN,
                "astSlaAdd(%s): Invalid SLALIB sky coordinate conversion "
                "type \"%s\".", astGetClass( this ), cvt );
   }

/* Add the new conversion to the SlaMap. */
   AddSlaCvt( this, cvttype, args );
}

static AstPointSet *Transform( AstMapping *this, AstPointSet *in,
                               int forward, AstPointSet *out ) {
/*
*  Name:
*     Transform

*  Purpose:
*     Apply an SlaMap to transform a set of points.

*  Type:
*     Private function.

*  Synopsis:
*     #include "slamap.h"
*     AstPointSet *Transform( AstMapping *this, AstPointSet *in,
*                             int forward, AstPointSet *out )

*  Class Membership:
*     SlaMap member function (over-rides the astTransform method inherited
*     from the Mapping class).

*  Description:
*     This function takes an SlaMap and a set of points encapsulated
*     in a PointSet and transforms the points so as to perform the
*     sequence of SLALIB sky coordinate conversions specified by
*     previous invocations of astSlaAdd.

*  Parameters:
*     this
*        Pointer to the SlaMap.
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

*  Returned Value:
*     Pointer to the output (possibly new) PointSet.

*  Notes:
*     -  A null pointer will be returned if this function is invoked with the
*     global error status set, or if it should fail for any reason.
*     -  The number of coordinate values per point in the input PointSet must
*     match the number of coordinates for the SlaMap being applied.
*     -  If an output PointSet is supplied, it must have space for sufficient
*     number of points and coordinate values per point to accommodate the
*     result. Any excess space will be ignored.
*/

/* Local Variables: */
   AstPointSet *result;          /* Pointer to output PointSet */
   AstSlaMap *map;               /* Pointer to SlaMap to be applied */
   double **ptr_in;              /* Pointer to input coordinate data */
   double **ptr_out;             /* Pointer to output coordinate data */
   double *alpha;                /* Pointer to longitude array */
   double *args;                 /* Pointer to argument list for conversion */
   double *delta;                /* Pointer to latitude array */
   int coord;                    /* Loop counter for coordinates */
   int cvt;                      /* Loop counter for conversions */
   int end;                      /* Termination index for conversion loop */
   int inc;                      /* Increment for conversion loop */
   int ncoord_in;                /* Number of coordinates per input point */
   int npoint;                   /* Number of points */
   int point;                    /* Loop counter for points */
   int start;                    /* Starting index for conversion loop */

/* Check the global error status. */
   if ( !astOK ) return NULL;

/* Obtain a pointer to the SlaMap. */
   map = (AstSlaMap *) this;

/* Apply the parent mapping using the stored pointer to the Transform member
   function inherited from the parent Mapping class. This function validates
   all arguments and generates an output PointSet if necessary, but does not
   actually transform any coordinate values. */
   result = (*parent_transform)( this, in, forward, out );

/* We will now extend the parent astTransform method by performing the
   coordinate conversions needed to generate the output coordinate values. */

/* Determine the numbers of points and coordinates per point from the input
   PointSet and obtain pointers for accessing the input and output coordinate
   values. */
   ncoord_in = astGetNcoord( in );
   npoint = astGetNpoint( in );
   ptr_in = astGetPoints( in );      
   ptr_out = astGetPoints( result );

/* Determine whether to apply the forward or inverse transformation, according
   to the direction specified and whether the mapping has been inverted. */
   if ( astGetInvert( this ) ) forward = !forward;

/* Transform the coordinate values. */
/* -------------------------------- */
/* Use "alpha" and "delta" as synonyms for the arrays of longitude and latitude
   coordinate values stored in the output PointSet. */
   if ( astOK ) {
      alpha = ptr_out[ 0 ];
      delta = ptr_out[ 1 ];

/* Initialise the output coordinate values by copying the input ones. */
      (void) memcpy( alpha, ptr_in[ 0 ], sizeof( double ) * (size_t) npoint );
      (void) memcpy( delta, ptr_in[ 1 ], sizeof( double ) * (size_t) npoint );

/* We will loop to apply each SLALIB sky coordinate conversion in turn to the
   (alpha,delta) arrays. However, if the inverse transformation was requested,
   we must loop through these transformations in reverse order, so set up
   appropriate limits and an increment to control this loop. */
      start = forward ? 0 : map->ncvt - 1;
      end = forward ? map->ncvt : -1;
      inc = forward ? 1 : -1;

/* Loop through the coordinate conversions in the required order and obtain a
   pointer to the argument list for the current conversion. */
      for ( cvt = start; cvt != end; cvt += inc ) {
         args = map->cvtargs[ cvt ];

/* Define a local macro as a shorthand to apply the code given as "function"
   (the macro argument) to each element of the (alpha,delta) arrays in turn.
   Before applying this conversion function, each element is first checked for
   "bad" coordinates (indicated by the value AST__BAD) and appropriate "bad"
   result values are assigned if necessary. */
#define TRAN_ARRAY(function) \
        for ( point = 0; point < npoint; point++ ) { \
           if ( ( alpha[ point ] == AST__BAD ) || \
                ( delta[ point ] == AST__BAD ) ) { \
              alpha[ point ] == AST__BAD; \
              delta[ point ] == AST__BAD; \
	   } else { \
              function \
	   } \
        }

/* Classify the SLALIB sky coordinate conversion to be applied. */
         switch ( map->cvttype[ cvt ] ) {

/* Add E-terms of aberration. */
/* -------------------------- */
/* Add or subtract (for the inverse) the E-terms from each coordinate pair
   in turn, returning the results to the same arrays. */
            case AST__SLA_ADDET:
               if ( forward ) {
                  TRAN_ARRAY(slaAddet( alpha[ point ], delta[ point ],
                                       args[ 0 ],
                                       alpha + point, delta + point );)
	       } else {
                  TRAN_ARRAY(slaSubet( alpha[ point ], delta[ point ],
                                       args[ 0 ],
                                       alpha + point, delta + point );)
               }
               break;

/* Subtract E-terms of aberration. */
/* ------------------------------- */
/* This is the same as above, but with the forward and inverse cases
   transposed. */
            case AST__SLA_SUBET:
               if ( forward ) {
                  TRAN_ARRAY(slaSubet( alpha[ point ], delta[ point ],
                                       args[ 0 ],
                                       alpha + point, delta + point );)
	       } else {
                  TRAN_ARRAY(slaAddet( alpha[ point ], delta[ point ],
                                       args[ 0 ],
                                       alpha + point, delta + point );)
	       }
               break;

/* Apply Bessel-Newcomb pre-IAU 1976 (FK4) precession model. */
/* --------------------------------------------------------- */
/* Since we are transforming a sequence of points, first set up the required
   precession matrix, swapping the argument order to get the inverse matrix
   if required. */
            case AST__SLA_PREBN:
               {
                  double epoch1 = forward ? args[ 0 ] : args[ 1 ];
                  double epoch2 = forward ? args[ 1 ] : args[ 0 ];
                  double precess_matrix[ 3 ][ 3 ];
                  double vec1[ 3 ];
                  double vec2[ 3 ];
                  slaPrebn( epoch1, epoch2, precess_matrix );

/* For each point in the (alpha,delta) arrays, convert to Cartesian
   coordinates, apply the precession matrix, convert back to polar coordinates
   and then constrain the longitude result to lie in the range 0 to 2*pi
   (slaDcc2s doesn't do this itself). */
                  TRAN_ARRAY(slaDcs2c( alpha[ point ], delta[ point ], vec1 );
                             slaDmxv( precess_matrix, vec1, vec2 );
                             slaDcc2s( vec2, alpha + point, delta + point );
                             alpha[ point ] = slaDranrm( alpha[ point ] );)
	       }
               break;

/* Apply IAU 1975 (FK5) precession model. */
/* -------------------------------------- */
/* This is handled in the same way as above, but using the appropriate FK5
   precession matrix. */
            case AST__SLA_PREC:
               {
                  double epoch1 = forward ? args[ 0 ] : args[ 1 ];
                  double epoch2 = forward ? args[ 1 ] : args[ 0 ];
                  double precess_matrix[ 3 ][ 3 ];
                  double vec1[ 3 ];
                  double vec2[ 3 ];
                  slaPrec( epoch1, epoch2, precess_matrix );
                  TRAN_ARRAY(slaDcs2c( alpha[ point ], delta[ point ], vec1 );
                             slaDmxv( precess_matrix, vec1, vec2 );
                             slaDcc2s( vec2, alpha + point, delta + point );
                             alpha[ point ] = slaDranrm( alpha[ point ] );)
	       }
               break;

/* Convert FK4 to FK5 (no proper motion or parallax). */
/* -------------------------------------------------- */
/* Apply the conversion to each point. */
	    case AST__SLA_FK45Z:
               if ( forward ) {
                  TRAN_ARRAY(slaFk45z( alpha[ point ], delta[ point ],
                                       args[ 0 ],
                                       alpha + point, delta + point );)

/* The inverse transformation is also straightforward, except that we need a
   couple of dummy variables as function arguments. */
	       } else {
                  double dr1950;
                  double dd1950;
                  TRAN_ARRAY(slaFk54z( alpha[ point ], delta[ point ],
                                       args[ 0 ],
                                       alpha + point, delta + point,
                                       &dr1950, &dd1950 );)
	       }
               break;

/* Convert FK5 to FK4 (no proper motion or parallax). */
/* -------------------------------------------------- */
/* This is the same as above, but with the forward and inverse cases
   transposed. */
	    case AST__SLA_FK54Z:
               if ( forward ) {
                  double dr1950;
                  double dd1950;
                  TRAN_ARRAY(slaFk54z( alpha[ point ], delta[ point ],
                                       args[ 0 ],
                                       alpha + point, delta + point,
                                       &dr1950, &dd1950 );)
	       } else {
                  TRAN_ARRAY(slaFk45z( alpha[ point ], delta[ point ],
                                       args[ 0 ],
                                       alpha + point, delta + point );)
               }
               break;

/* Convert geocentric apparent to mean place. */
/* ------------------------------------------ */
/* Since we are transforming a sequence of points, first set up the required
   parameter array. Than apply this to each point in turn. */
	    case AST__SLA_AMP:
               {
                  double amprms[ 21 ];
                  slaMappa( args[ 1 ], args[ 0 ], amprms );
                  if ( forward ) {
                     TRAN_ARRAY(slaAmpqk( alpha[ point ], delta[ point ],
                                          amprms,
                                          alpha + point, delta + point );)

/* The inverse uses the same parameter array but converts from mean place
   to geocentric apparent. */
                  } else {
                     TRAN_ARRAY(slaMapqkz( alpha[ point ], delta[ point ],
                                           amprms,
                                           alpha + point, delta + point );)
		  }
               }
               break;

/* Convert mean place to geocentric apparent. */
/* ------------------------------------------ */
/* This is the same as above, but with the forward and inverse cases
   transposed. */
	    case AST__SLA_MAP:
               {
                  double amprms[ 21 ];
                  slaMappa( args[ 0 ], args[ 1 ], amprms );
                  if ( forward ) {
                     TRAN_ARRAY(slaMapqkz( alpha[ point ], delta[ point ],
                                           amprms,
                                           alpha + point, delta + point );)
                  } else {
                     TRAN_ARRAY(slaAmpqk( alpha[ point ], delta[ point ],
                                          amprms,
                                          alpha + point, delta + point );)
		  }
               }
               break;

/* Convert ecliptic coordinates to J2000.0 equatorial. */
/* --------------------------------------------------- */
/* Since we are transforming a sequence of points, first set up the required
   conversion matrix (the conversion is a rotation). */
	    case AST__SLA_ECLEQ:
               {
                  double convert_matrix[ 3 ][ 3 ];
                  double precess_matrix[ 3 ][ 3 ];
                  double rotate_matrix[ 3 ][ 3 ];
                  double vec1[ 3 ];
                  double vec2[ 3 ];

/* Obtain the matrix that precesses equatorial coordinates from J2000.0 to the
   required date. Also obtain the rotation matrix that converts from
   equatorial to ecliptic coordinates.  */
                  slaPrec( 2000.0, slaEpj( args[ 0 ] ), precess_matrix );
                  slaEcmat( args[ 0 ], rotate_matrix );

/* Multiply these matrices to give the overall matrix that converts from
   equatorial J2000.0 coordinates to ecliptic coordinates for the required
   date. */
                  slaDmxm( rotate_matrix, precess_matrix, convert_matrix );

/* Apply the conversion by transforming from polar to Cartesian coordinates,
   multiplying by the inverse conversion matrix and converting back to polar
   coordinates. Then constrain the longitude result to lie in the range
   0 to 2*pi (slaDcc2s doesn't do this itself). */
                  if ( forward ) {
                     TRAN_ARRAY(slaDcs2c( alpha[ point ], delta[ point ],
                                          vec1 );
                                slaDimxv( convert_matrix, vec1, vec2 );
                                slaDcc2s( vec2, alpha + point, delta + point );
                                alpha[ point ] = slaDranrm ( alpha[ point ] );)

/* The inverse conversion is the same except that we multiply by the forward
   conversion matrix (slaDmxv instead of slaDimxv). */
                  } else {
                     TRAN_ARRAY(slaDcs2c( alpha[ point ], delta[ point ],
                                          vec1 );
                                slaDmxv( convert_matrix, vec1, vec2 );
                                slaDcc2s( vec2, alpha + point, delta + point );
                                alpha[ point ] = slaDranrm ( alpha[ point ] );)
                  }
	       }
               break;

/* Convert equatorial J2000.0 to ecliptic coordinates. */
/* --------------------------------------------------- */
/* This is the same as above, but with the forward and inverse cases
   transposed. */
	    case AST__SLA_EQECL:
               {
                  double convert_matrix[ 3 ][ 3 ];
                  double precess_matrix[ 3 ][ 3 ];
                  double rotate_matrix[ 3 ][ 3 ];
                  double vec1[ 3 ];
                  double vec2[ 3 ];

/* Create the conversion matrix. */
                  slaPrec( 2000.0, slaEpj( args[ 0 ] ), precess_matrix );
                  slaEcmat( args[ 0 ], rotate_matrix );
                  slaDmxm( rotate_matrix, precess_matrix, convert_matrix );

/* Apply it. */
                  if ( forward ) {
                     TRAN_ARRAY(slaDcs2c( alpha[ point ], delta[ point ],
                                          vec1 );
                                slaDmxv( convert_matrix, vec1, vec2 );
                                slaDcc2s( vec2, alpha + point, delta + point );
                                alpha[ point ] = slaDranrm ( alpha[ point ] );)
                  } else {
                     TRAN_ARRAY(slaDcs2c( alpha[ point ], delta[ point ],
                                          vec1 );
                                slaDimxv( convert_matrix, vec1, vec2 );
                                slaDcc2s( vec2, alpha + point, delta + point );
                                alpha[ point ] = slaDranrm ( alpha[ point ] );)
                  }
	       }
               break;

/* Convert galactic coordinates to J2000.0 equatorial. */
/* --------------------------------------------------- */
/* Apply the conversion to each point. */
	    case AST__SLA_GALEQ:
               if ( forward ) {
                  TRAN_ARRAY(slaGaleq( alpha[ point ], delta[ point ],
                                       alpha + point, delta + point );)

/* The inverse simply uses the inverse SLALIB function. */
	       } else {
                  TRAN_ARRAY(slaEqgal( alpha[ point ], delta[ point ],
                                       alpha + point, delta + point );)
	       }
               break;

/* Convert J2000.0 equatorial to galactic coordinates. */
/* --------------------------------------------------- */
/* This is the same as above, but with the forward and inverse cases
   transposed. */
	    case AST__SLA_EQGAL:
               if ( forward ) {
                  TRAN_ARRAY(slaEqgal( alpha[ point ], delta[ point ],
                                       alpha + point, delta + point );)
	       } else {
                  TRAN_ARRAY(slaGaleq( alpha[ point ], delta[ point ],
                                       alpha + point, delta + point );)
               }
               break;

/* Convert galactic to supergalactic coordinates. */
/* ---------------------------------------------- */
/* Apply the conversion to each point. */
	    case AST__SLA_GALSUP:
               if ( forward ) {
                  TRAN_ARRAY(slaGalsup( alpha[ point ], delta[ point ],
                                        alpha + point, delta + point );)

/* The inverse simply uses the inverse SLALIB function. */
               } else {
                  TRAN_ARRAY(slaSupgal( alpha[ point ], delta[ point ],
                                        alpha + point, delta + point );)
               }
               break;

/* Convert supergalactic coordinates to galactic. */
/* ---------------------------------------------- */
/* This is the same as above, but with the forward and inverse cases
   transposed. */
	    case AST__SLA_SUPGAL:
               if ( forward ) {
                  TRAN_ARRAY(slaSupgal( alpha[ point ], delta[ point ],
                                        alpha + point, delta + point );)
               } else {
                  TRAN_ARRAY(slaGalsup( alpha[ point ], delta[ point ],
                                        alpha + point, delta + point );)
               }
               break;

/* If the conversion type was not recognised, then report an error
   (this should not happen unless validation in astSlaAdd has failed
   to detect a bad value previously). */
            default:
               astError( AST__SLAIN, "astTransform(%s): Corrupt %s contains "
                         "invalid SLALIB sky coordinate conversion code (%d).",
                         astGetClass( this ), astGetClass( this ),
                         (int) map->cvttype[ cvt ] );
               break;
         }
      }
   }

/* If an error has occurred and a new PointSet may have been created, then
   clean up by annulling it. In any case, ensure that a NULL result is
   returned.*/
   if ( !astOK ) {
      if ( !out ) result = astAnnul( result );
      result = NULL;
   }

/* Return a pointer to the output PointSet. */
   return result;

/* Undefine macros local to this function. */
#undef TRAN_ARRAY
}

/* Copy constructor. */
/* ----------------- */
static void Copy( const AstObject *objin, AstObject *objout ) {
/*
*  Name:
*     Copy

*  Purpose:
*     Copy constructor for SlaMap objects.

*  Type:
*     Private function.

*  Synopsis:
*     void Copy( const AstObject *objin, AstObject *objout )

*  Description:
*     This function implements the copy constructor for SlaMap objects.

*  Parameters:
*     objin
*        Pointer to the object to be copied.
*     objout
*        Pointer to the object being constructed.

*  Returned Value:
*     void

*  Notes:
*     -  This constructor makes a deep copy.
*/

/* Local Variables: */
   AstSlaMap *in;                /* Pointer to input SlaMap */
   AstSlaMap *out;               /* Pointer to output SlaMap */
   int cvt;                      /* Loop counter for coordinate conversions */

/* Check the global error status. */
   if ( !astOK ) return;

/* Obtain pointers to the input and output SlaMap structures. */
   in = (AstSlaMap *) objin;
   out = (AstSlaMap *) objout;

/* For safety, first clear any references to the input memory from the output
   SlaMap. */
   out->cvtargs = NULL;
   out->cvttype = NULL;

/* Allocate memory for the output array of argument list pointers. */
   out->cvtargs = astMalloc( sizeof( double * ) * (size_t) in->ncvt );

/* If necessary, allocate memory and make a copy of the input array of sky
   coordinate conversion codes. */
   if ( in->cvttype ) out->cvttype = astStore( NULL, in->cvttype,
                                               sizeof( int )
                                               * (size_t) in->ncvt );

/* If OK, loop through each conversion in the input SlaMap and make a copy of
   its argument list, storing the new pointer in the output argument list
   array. */
   if ( astOK ) {
      for ( cvt = 0; cvt < in->ncvt; cvt++ ) {
         out->cvtargs[ cvt ] = astStore( NULL, in->cvtargs[ cvt ],
                                         astSizeOf( in->cvtargs[ cvt ] ) );
      }

/* If an error occurred while copying the argument lists, loop through the
   conversions again and clean up by ensuring that the new memory allocated for
   each argument list is freed. */
      if ( !astOK ) {
         for ( cvt = 0; cvt < in->ncvt; cvt++ ) {
            out->cvtargs[ cvt ] = astFree( out->cvtargs[ cvt ] );
	 }
      }
   }

/* If an error occurred, free all other memory allocated above. */
   if ( !astOK ) {
      out->cvtargs = astFree( out->cvtargs );
      out->cvttype = astFree( out->cvttype );
   }
}

/* Destructor. */
/* ----------- */
static void Delete( AstObject *obj ) {
/*
*  Name:
*     Delete

*  Purpose:
*     Destructor for SlaMap objects.

*  Type:
*     Private function.

*  Synopsis:
*     void Delete( AstObject *obj )

*  Description:
*     This function implements the destructor for SlaMap objects.

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
   AstSlaMap *this;              /* Pointer to SlaMap */
   int cvt;                      /* Loop counter for coordinate conversions */

/* Obtain a pointer to the SlaMap structure. */
   this = (AstSlaMap *) obj;

/* Loop to free the memory containing the argument list for each sky coordinate
   conversion. */
   for ( cvt = 0; cvt < this->ncvt; cvt++ ) {
      this->cvtargs[ cvt ] = astFree( this->cvtargs[ cvt ] );
   }

/* Free the memory holding the array of conversion types and the array of
   argument list pointers. */
   this->cvtargs = astFree( this->cvtargs );
   this->cvttype = astFree( this->cvttype );
}

/* Dump function. */
/* -------------- */
static void Dump( AstObject *this_object, AstChannel *channel ) {
/*
*  Name:
*     Dump

*  Purpose:
*     Dump function for SlaMap objects.

*  Type:
*     Private function.

*  Synopsis:
*     #include "slamap.h"
*     void Dump( AstObject *this, AstChannel *channel )

*  Description:
*     This function implements the Dump function which writes out data
*     for the SlaMap class to an output Channel.

*  Parameters:
*     this
*        Pointer to the SlaMap whose data are being written.
*     channel
*        Pointer to the Channel to which the data are being written.
*/

/* Local Constants: */
#define KEY_LEN 50               /* Maximum length of a keyword */

/* Local Variables: */
   AstSlaMap *this;              /* Pointer to the SlaMap structure */
   char key[ KEY_LEN + 1 ];      /* Buffer for keyword string */
   const char *argdesc[ MAX_SLA_ARGS ]; /* Pointers to argument descriptions */
   const char *comment;          /* Pointer to comment string */
   const char *sval;             /* Pointer to string value */
   int iarg;                     /* Loop counter for arguments */
   int icvt;                     /* Loop counter for conversion steps */
   int ival;                     /* Integer value */
   int nargs;                    /* Number of conversion arguments */
   int set;                      /* Attribute value set? */

/* Check the global error status. */
   if ( !astOK ) return;

/* Obtain a pointer to the SlaMap structure. */
   this = (AstSlaMap *) this_object;

/* Write out values representing the instance variables for the SlaMap
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

/* Number of conversion steps. */
/* --------------------------- */
/* Regard this as "set" if it is non-zero. */
   ival = this->ncvt;
   set = ( ival != 0 );
   astWriteInt( channel, "Nsla", set, 0, ival, "Number of conversion steps" );

/* Write out data for each conversion step... */
   for ( icvt = 0; icvt < this->ncvt; icvt++ ) {

/* Conversion type. */
/* ---------------- */
/* Change each conversion type code into an equivalent string and
   obtain associated descriptive information. If the conversion code
   was not recognised, report an error and give up. */
      if ( astOK ) {
         sval = CvtString( this->cvttype[ icvt ], &comment, &nargs, argdesc );
         if ( astOK && !sval ) {
            astError( AST__SLAIN,
                      "astWrite(%s): Corrupt %s contains invalid SLALIB "
                      "sky coordinate conversion code (%d).",
                      astGetClass( channel ), astGetClass( this ),
                      (int) this->cvttype[ icvt ] );
            break;
         }

/* Create an appropriate keyword and write out the conversion code
   information. */
         (void) sprintf( key, "Sla%d", icvt + 1 );
         astWriteString( channel, key, 1, 1, sval, comment );

/* Write out data for each conversion argument... */
         for ( iarg = 0; iarg < nargs; iarg++ ) {

/* Arguments. */
/* ---------- */
/* Create an appropriate keyword and write out the argument value,
   accompanied by the descriptive comment obtained above. */
            (void) sprintf( key, "Sla%d%c", icvt + 1, ALPHABET[ iarg ] );
            astWriteDouble( channel, key, 1, 1, this->cvtargs[ icvt ][ iarg ],
                            argdesc[ iarg ] );
         }

/* Quit looping if an error occurs. */
         if ( !astOK ) break;
      }
   }

/* Undefine macros local to this function. */
#undef KEY_LEN
}

/* Standard class functions. */
/* ========================= */
/* Implement the astIsASlaMap and astCheckSlaMap functions using the macros
   defined for this purpose in the "object.h" header file. */
astMAKE_ISA(SlaMap,Mapping,check,&class_init)
astMAKE_CHECK(SlaMap)

AstSlaMap *astSlaMap_( int flags, const char *options, ... ) {
/*
*++
*  Name:
c     astSlaMap
f     AST_SLAMAP

*  Purpose:
*     Create an SlaMap.

*  Type:
*     Public function.

*  Synopsis:
c     #include "slamap.h"
c     AstSlaMap *astSlaMap( int flags, const char *options, ... )
f     RESULT = AST_SLAMAP( FLAGS, OPTIONS, STATUS )

*  Class Membership:
*     SlaMap constructor.

*  Description:
*     This function creates a new SlaMap and optionally initialises
*     its attributes.
*
*     An SlaMap is a specialised form of Mapping which can be used to
*     represent a sequence of conversions between standard celestial
*     (longitude, latitude) coordinate systems.
*
*     When an SlaMap is first created, it simply performs a unit
c     (null) Mapping on a pair of coordinates. Using the astSlaAdd
f     (null) Mapping on a pair of coordinates. Using the AST_SLAADD
c     function, a series of coordinate conversion steps may then be
f     routine, a series of coordinate conversion steps may then be
*     added, selected from those provided by the SLALIB Positional
*     Astronomy Library (Starlink User Note SUN/67). This allows
*     multi-step conversions between a variety of celestial coordinate
*     systems to be assembled out of the building blocks provided by
*     SLALIB.
*
*     For details of the individual coordinate conversions available,
c     see the description of the astSlaAdd function.
f     see the description of the AST_SLAADD routine.

*  Parameters:
c     flags
f     FLAGS = INTEGER (Given)
c        This parameter is reserved for future use and should currently
f        This argument is reserved for future use and should currently
*        always be set to zero.
c     options
f     OPTIONS = CHARACTER * ( * ) (Given)
c        Pointer to a null-terminated string containing an optional
c        comma-separated list of attribute assignments to be used for
c        initialising the new SlaMap. The syntax used is identical to
c        that for the astSet function and may include "printf" format
c        specifiers identified by "%" symbols in the normal way.
c        If no initialisation is required, a zero-length string may be
c        supplied.
f        A character string containing an optional comma-separated
f        list of attribute assignments to be used for initialising the
f        new SlaMap. The syntax used is identical to that for the
f        AST_SET routine. If no initialisation is required, a blank
f        value may be supplied.
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
c     astSlaMap()
f     AST_SLAMAP = INTEGER
*        A pointer to the new SlaMap.

*  Notes:
*     - The Nin and Nout attributes (number of input and output
*     coordinates) for an SlaMap are both equal to 2. The first
*     coordinate is the celestial longitude and the second coordinate
*     is the celestial latitude. All coordinate values are in radians.
*     - A null Object pointer (AST__NULL) will be returned if this
c     function is invoked with the AST error status set, or if it
f     function is invoked with STATUS set to an error value, or if it
*     should fail for any reason.
*--
*/

/* Local Variables: */
   AstSlaMap *new;               /* Pointer to the new SlaMap */
   va_list args;                 /* Variable argument list */

/* Check the global status. */
   if ( !astOK ) return NULL;

/* Initialise the SlaMap, allocating memory and initialising the virtual
   function table as well if necessary. */
   new = astInitSlaMap( NULL, sizeof( AstSlaMap ), !class_init, &class_vtab,
                        "SlaMap", flags );

/* If successful, note that the virtual function table has been initialised. */
   if ( astOK ) {
      class_init = 1;

/* Obtain the variable argument list and pass it along with the options string
   to the astVSet method to initialise the new SlaMap's attributes. */
      va_start( args, options );
      astVSet( new, options, args );
      va_end( args );

/* If an error occurred, clean up by deleting the new object. */
      if ( !astOK ) new = astDelete( new );
   }

/* Return a pointer to the new SlaMap. */
   return new;
}

AstSlaMap *astSlaMapId_( int flags, const char *options, ... ) {
/*
*  Name:
*     astSlaMapId_

*  Purpose:
*     Create an SlaMap.

*  Type:
*     Private function.

*  Synopsis:
*     #include "slamap.h"
*     AstSlaMap *astSlaMapId_( int flags, const char *options, ... )

*  Class Membership:
*     SlaMap constructor.

*  Description:
*     This function implements the external (public) interface to the
*     astSlaMap constructor function. It returns an ID value (instead
*     of a true C pointer) to external users, and must be provided
*     because astSlaMap_ has a variable argument list which cannot be
*     encapsulated in a macro (where this conversion would otherwise
*     occur).
*
*     The variable argument list also prevents this function from
*     invoking astSlaMap_ directly, so it must be a re-implementation
*     of it in all respects, except for the final conversion of the
*     result to an ID value.

*  Parameters:
*     As for astSlaMap_.

*  Returned Value:
*     The ID value associated with the new SlaMap.
*/

/* Local Variables: */
   AstSlaMap *new;               /* Pointer to the new SlaMap */
   va_list args;                 /* Variable argument list */

/* Check the global status. */
   if ( !astOK ) return NULL;

/* Initialise the SlaMap, allocating memory and initialising the virtual
   function table as well if necessary. */
   new = astInitSlaMap( NULL, sizeof( AstSlaMap ), !class_init, &class_vtab,
                        "SlaMap", flags );

/* If successful, note that the virtual function table has been initialised. */
   if ( astOK ) {
      class_init = 1;

/* Obtain the variable argument list and pass it along with the options string
   to the astVSet method to initialise the new SlaMap's attributes. */
      va_start( args, options );
      astVSet( new, options, args );
      va_end( args );

/* If an error occurred, clean up by deleting the new object. */
      if ( !astOK ) new = astDelete( new );
   }

/* Return an ID value for the new SlaMap. */
   return astMakeId( new );
}

AstSlaMap *astInitSlaMap_( void *mem, size_t size, int init,
                           AstSlaMapVtab *vtab, const char *name,
                           int flags ) {
/*
*+
*  Name:
*     astInitSlaMap

*  Purpose:
*     Initialise an SlaMap.

*  Type:
*     Protected function.

*  Synopsis:
*     #include "slamap.h"
*     AstSlaMap *astInitSlaMap( void *mem, size_t size, int init,
*                               AstSlaMapVtab *vtab, const char *name,
*                               int flags )

*  Class Membership:
*     SlaMap initialiser.

*  Description:
*     This function is provided for use by class implementations to initialise
*     a new SlaMap object. It allocates memory (if necessary) to accommodate
*     the SlaMap plus any additional data associated with the derived class.
*     It then initialises an SlaMap structure at the start of this memory. If
*     the "init" flag is set, it also initialises the contents of a virtual
*     function table for an SlaMap at the start of the memory passed via the
*     "vtab" parameter.

*  Parameters:
*     mem
*        A pointer to the memory in which the SlaMap is to be initialised.
*        This must be of sufficient size to accommodate the SlaMap data
*        (sizeof(SlaMap)) plus any data used by the derived class. If a value
*        of NULL is given, this function will allocate the memory itself using
*        the "size" parameter to determine its size.
*     size
*        The amount of memory used by the SlaMap (plus derived class data).
*        This will be used to allocate memory if a value of NULL is given for
*        the "mem" parameter. This value is also stored in the SlaMap
*        structure, so a valid value must be supplied even if not required for
*        allocating memory.
*     init
*        A logical flag indicating if the SlaMap's virtual function table is
*        to be initialised. If this value is non-zero, the virtual function
*        table will be initialised by this function.
*     vtab
*        Pointer to the start of the virtual function table to be associated
*        with the new SlaMap.
*     name
*        Pointer to a constant null-terminated character string which contains
*        the name of the class to which the new object belongs (it is this
*        pointer value that will subsequently be returned by the astClass
*        method).
*     flags
*        This parameter is reserved for future use. It is currently ignored.

*  Returned Value:
*     A pointer to the new SlaMap.

*  Notes:
*     -  A null pointer will be returned if this function is invoked with the
*     global error status set, or if it should fail for any reason.
*-
*/

/* Local Variables: */
   AstSlaMap *new;               /* Pointer to the new SlaMap */

/* Check the global status. */
   if ( !astOK ) return NULL;

/* Initialise a Mapping structure (the parent class) as the first component
   within the SlaMap structure, allocating memory if necessary. Specify that
   the Mapping should be defined in both the forward and inverse directions. */
   new = (AstSlaMap *) astInitMapping( mem, size, init,
                                       (AstMappingVtab *) vtab, name,
                                       2, 2, 1, 1 );

/* If necessary, initialise the virtual function table. */
/* ---------------------------------------------------- */
   if ( init ) InitVtab( vtab );
   if ( astOK ) {

/* Initialise the SlaMap data. */
/* --------------------------- */
/* The initial state is with no SLALIB conversions set, in which condition the
   SlaMap simply implements a unit mapping. */
      new->ncvt = 0;
      new->cvtargs = NULL;
      new->cvttype = NULL;

/* If an error occurred, clean up by deleting the new object. */
      if ( !astOK ) new = astDelete( new );
   }

/* Return a pointer to the new object. */
   return new;
}

AstSlaMap *astLoadSlaMap_( void *mem, size_t size, int init,
                           AstSlaMapVtab *vtab, const char *name,
                           AstChannel *channel ) {
/*
*+
*  Name:
*     astLoadSlaMap

*  Purpose:
*     Load a SlaMap.

*  Type:
*     Protected function.

*  Synopsis:
*     #include "slamap.h"
*     AstSlaMap *astLoadSlaMap( void *mem, size_t size, int init,
*                               AstSlaMapVtab *vtab, const char *name,
*                               AstChannel *channel )

*  Class Membership:
*     SlaMap loader.

*  Description:
*     This function is provided to load a new SlaMap using data read
*     from a Channel. It first loads the data used by the parent class
*     (which allocates memory if necessary) and then initialises a
*     SlaMap structure in this memory, using data read from the input
*     Channel.
*
*     If the "init" flag is set, it also initialises the contents of a
*     virtual function table for a SlaMap at the start of the memory
*     passed via the "vtab" parameter.

*  Parameters:
*     mem
*        A pointer to the memory into which the SlaMap is to be
*        loaded.  This must be of sufficient size to accommodate the
*        SlaMap data (sizeof(SlaMap)) plus any data used by derived
*        classes. If a value of NULL is given, this function will
*        allocate the memory itself using the "size" parameter to
*        determine its size.
*     size
*        The amount of memory used by the SlaMap (plus derived class
*        data).  This will be used to allocate memory if a value of
*        NULL is given for the "mem" parameter. This value is also
*        stored in the SlaMap structure, so a valid value must be
*        supplied even if not required for allocating memory.
*
*        If the "vtab" parameter is NULL, the "size" value is ignored
*        and sizeof(AstSlaMap) is used instead.
*     init
*        A boolean flag indicating if the SlaMap's virtual function
*        table is to be initialised. If this value is non-zero, the
*        virtual function table will be initialised by this function.
*
*        If the "vtab" parameter is NULL, the "init" value is ignored
*        and the (static) virtual function table initialisation flag
*        for the SlaMap class is used instead.
*     vtab
*        Pointer to the start of the virtual function table to be
*        associated with the new SlaMap. If this is NULL, a pointer to
*        the (static) virtual function table for the SlaMap class is
*        used instead.
*     name
*        Pointer to a constant null-terminated character string which
*        contains the name of the class to which the new object
*        belongs (it is this pointer value that will subsequently be
*        returned by the astGetClass method).
*
*        If the "vtab" parameter is NULL, the "name" value is ignored
*        and a pointer to the string "SlaMap" is used instead.

*  Returned Value:
*     A pointer to the new SlaMap.

*  Notes:
*     - A null pointer will be returned if this function is invoked
*     with the global error status set, or if it should fail for any
*     reason.
*-
*/

/* Local Constants: */
#define KEY_LEN 50               /* Maximum length of a keyword */

/* Local Variables: */
   AstSlaMap *new;               /* Pointer to the new SlaMap */
   char *sval;                   /* Pointer to string value */
   char key[ KEY_LEN + 1 ];      /* Buffer for keyword string */
   const char *argdesc[ MAX_SLA_ARGS ]; /* Pointers to argument descriptions */
   const char *comment;          /* Pointer to comment string */
   int iarg;                     /* Loop counter for arguments */
   int icvt;                     /* Loop counter for conversion steps */
   int nargs;                    /* Number of conversion arguments */

/* Initialise. */
   new = NULL;

/* Check the global error status. */
   if ( !astOK ) return new;

/* If a NULL virtual function table has been supplied, then this is
   the first loader to be invoked for this SlaMap. In this case the
   SlaMap belongs to this class, so supply appropriate values to be
   passed to the parent class loader (and its parent, etc.). */
   if ( !vtab ) {
      size = sizeof( AstSlaMap );
      init = !class_init;
      vtab = &class_vtab;
      name = "SlaMap";
   }

/* Invoke the parent class loader to load data for all the ancestral
   classes of the current one, returning a pointer to the resulting
   partly-built SlaMap. */
   new = astLoadMapping( mem, size, init, (AstMappingVtab *) vtab, name,
                         channel );

/* If required, initialise the part of the virtual function table used
   by this class. */
   if ( init ) InitVtab( vtab );

/* Note if we have successfully initialised the (static) virtual
   function table owned by this class (so that this is done only
   once). */
   if ( astOK ) {
      if ( ( vtab == &class_vtab ) && init ) class_init = 1;

/* Read input data. */
/* ================ */
/* Request the input Channel to read all the input data appropriate to
   this class into the internal "values list". */
      astReadClassData( channel, "SlaMap" );

/* Now read each individual data item from this list and use it to
   initialise the appropriate instance variable(s) for this class. */

/* In the case of attributes, we first read the "raw" input value,
   supplying the "unset" value as the default. If a "set" value is
   obtained, we then use the appropriate (private) Set... member
   function to validate and set the value properly. */

/* Number of conversion steps. */
/* --------------------------- */
/* Read the number of conversion steps and allocate memory to hold
   data for each step. */
      new->ncvt = astReadInt( channel, "nsla", 0 );
      if ( new->ncvt < 0 ) new->ncvt = 0;
      new->cvttype = astMalloc( sizeof( int ) * (size_t) new->ncvt );
      new->cvtargs = astMalloc( sizeof( double * ) * (size_t) new->ncvt );

/* If an error occurred, ensure that all allocated memory is freed. */
      if ( !astOK ) {
         new->cvttype = astFree( new->cvttype );
         new->cvtargs = astFree( new->cvtargs );

/* Otherwise, initialise the argument pointer array. */
      } else {
         for ( icvt = 0; icvt < new->ncvt; icvt++ ) {
            new->cvtargs[ icvt ] = NULL;
         }

/* Read in data for each conversion step... */
         for ( icvt = 0; icvt < new->ncvt; icvt++ ) {

/* Conversion type. */
/* ---------------- */
/* Create an appropriate keyword and read the string representation of
   the conversion type. */
            (void) sprintf( key, "sla%d", icvt + 1 );
            sval = astReadString( channel, key, NULL );

/* If no value was read, report an error. */
            if ( astOK ) {
               if ( !sval ) {
                  astError( AST__BADIN,
                            "astRead(%s): An SLALIB sky coordinate conversion "
                            "type is missing from the input SlaMap data.",
                            astGetClass( channel ) );

/* Otherwise, convert the string representation into the required
   conversion type code. */
               } else {
                  new->cvttype[ icvt ] = CvtCode( sval );

/* If the string was not recognised, report an error. */
                  if ( new->cvttype[ icvt ] == AST__SLA_NULL ) {
                     astError( AST__BADIN,
                              "astRead(%s): Invalid SLALIB sky conversion "
                              "type \"%s\" in SlaMap data.",
                              astGetClass( channel ), sval );
                  }
               }

/* Free the memory holding the string value. */
               sval = astFree( sval );
            }

/* Obtain the number of arguments associated with the conversion and
   allocate memory to hold them. */
            (void) CvtString( new->cvttype[ icvt ], &comment, &nargs,
                              argdesc );
            new->cvtargs[ icvt ] = astMalloc( sizeof( double ) *
                                              (size_t) nargs );

/* Read in data for each argument... */
            if ( astOK ) {
               for ( iarg = 0; iarg < nargs; iarg++ ) {

/* Arguments. */
/* ---------- */
/* Create an appropriate keyword and read each argument value. */
                  (void) sprintf( key, "sla%d%c", icvt + 1, ALPHABET[ iarg ] );
                  new->cvtargs[ icvt ][ iarg ] = astReadDouble( channel, key,
                                                                AST__BAD );
               }
            }

/* Quit looping if an error occurs. */
            if ( !astOK ) break;
         }
      }

/* If an error occurred, clean up by deleting the new SlaMap. */
      if ( !astOK ) new = astDelete( new );
   }

/* Return the new SlaMap pointer. */
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
void astSlaAdd_( AstSlaMap *this, const char *cvt, const double args[] ) {
   if ( !astOK ) return;
   (**astMEMBER(this,SlaMap,SlaAdd))( this, cvt, args );
}
