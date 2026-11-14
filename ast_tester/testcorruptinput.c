/*
*  Name:
*     testcorruptinput

*  Purpose:
*     Check that corrupt or malicious input is rejected without harm.

*  Description:
*     An object dump or a FITS header can come from an untrusted source, so
*     reading one must never crash, corrupt memory, or run for a time or
*     use an amount of memory out of proportion to its size.  Each case
*     here is a small input of the kind produced by a fuzzer, which changes
*     one value in a valid input, and each made the library do one of those
*     things.  The test reads every input and requires either the stated
*     error or, where the input is merely unusual, a usable object.  It is
*     intended to be run with the address and undefined-behavior sanitizers
*     enabled, which is how most of the faults are detected.  The cases
*     cover:
*
*       - Integers: a value too large for an int, which used to wrap round
*         to another value.
*       - Axis limit: a Mapping or Frame cannot have more than AST__MXAXES
*         inputs, outputs or axes, whether it comes from a constructor or
*         from a dump.  The limit exists so that a file that declares an
*         axis count of billions cannot make the library allocate memory
*         and loop in proportion to that count.
*       - Components: an object that is missing from a dump, or belongs
*         to the wrong class, where another object needs it as a component.
*       - Regions: defining points that are missing or the wrong number,
*         components with the wrong number of axes, and a PointSet that
*         claims more values than a dump can hold.
*       - FrameSets: more Frames or nodes claimed than are present, node
*         links that form a loop, and Mappings and Frames that disagree
*         about the number of axes.
*       - Frames: the wrong number of axes for a SkyFrame, SpecFrame,
*         TimeFrame or FluxFrame, a SkyFrame whose axes are not described,
*         and the wrong components for a SpecFluxFrame.
*       - Mappings: a loader that carries on after an error, a PermMap that
*         refers to a constant it does not have, a full MatrixMap with more
*         elements than can be counted or named, a LutMap given an input
*         far beyond its table, and a KeyMap that asks for a huge hash
*         table.
*       - Large values: a sky position or an Epoch that is too large to
*         format or to bring into range one turn at a time.
*       - Text of any length: values, Labels, Domains and Units that are
*         far longer than the fixed buffers into which error messages,
*         warnings, titles and formatted values used to be written without
*         a check, and a Format string that asks "sprintf" for values that
*         are not there.
*       - Splitting: a compound Mapping in which a MatrixMap or WcsMap
*         couples two of many axes, for which astMapSplit used to try
*         every arrangement of every subset of the axes.
*       - FITS headers: a WCSAXES value that is not a usable number of
*         axes, CTYPE values of fewer than three characters, an IRAF WAT
*         card with no string value, a header whose transformation fails
*         part of the way through a compound Mapping, a pixel size so
*         small that its reciprocal is infinite, a rest frequency so small
*         that a fraction of it is zero, and a floating point value that is
*         too large to convert to an integer.

*  Notes:
*     - This is a native C test with no Fortran original.
*     - Public C API only.  Errors expected from a deliberately invalid
*       input are caught through astStatus and cleared with astClearStatus
*       so the remaining cases can run; failures of the test itself are
*       counted and reported at the end.
*/

#include <stdio.h>
#include <string.h>
#include "sae_par.h"
#include "ast.h"
#include "ast_err.h"

/* Check that the error "expected" has been reported, for the operation
   described by "what" that was meant to fail and return "obj" (which is
   annulled if it is not NULL).  Returns 0 if so, and 1 otherwise.  Any
   raised status is cleared before returning. */
static int check_failed( const char *what, void *obj, int expected ) {
   int code = astStatus;
   if ( !astOK ) astClearStatus;
   if ( obj ) {
      obj = astAnnul( obj );
      fprintf( stderr, "testcorruptinput: %s gave an object but should "
               "have failed\n", what );
      return 1;
   }
   if ( code != expected ) {
      fprintf( stderr, "testcorruptinput: %s gave status %d but status %d "
               "was expected\n", what, code, expected );
      return 1;
   }
   return 0;
}

/* Attempt to load a dump that should be rejected with the error code
   "expected".  Returns 0 if no object was returned and that error was
   reported, 1 otherwise. */
static int check_rejected( const char *what, const char *dump, int expected ) {
   return check_failed( what, astFromString( dump ), expected );
}

/* Check that the operation described by "what" gave an object, "obj", and
   no error.  The object is annulled.  Returns 0 if so, and 1 otherwise. */
static int check_made( const char *what, void *obj ) {
   int code = astStatus;
   if ( !astOK ) astClearStatus;
   if ( obj ) obj = astAnnul( obj );
   if ( code != SAI__OK ) {
      fprintf( stderr, "testcorruptinput: %s gave status %d but should "
               "have worked\n", what, code );
      return 1;
   }
   return 0;
}

/* Attempt to load a dump that should be accepted.  Returns 0 if an object
   was returned without error, 1 otherwise. */
static int check_accepted( const char *what, const char *dump ) {
   AstObject *obj = astFromString( dump );
   if ( !obj && astOK ) {
      fprintf( stderr, "testcorruptinput: %s gave no object\n", what );
      return 1;
   }
   return check_made( what, obj );
}

/* Return a copy of "dump" in which the first occurrence of "old" is
   replaced by "new".  The copy is allocated with astMalloc and should be
   freed by the caller.  Returns NULL if "old" is not found. */
static char *edit( const char *dump, const char *old, const char *new ) {
   const char *pos = dump ? strstr( dump, old ) : NULL;
   size_t pre;
   char *out;
   if ( !pos ) return NULL;
   pre = (size_t) ( pos - dump );
   out = astMalloc( strlen( dump ) - strlen( old ) + strlen( new ) + 1 );
   if ( out ) {
      memcpy( out, dump, pre );
      strcpy( out + pre, new );
      strcat( out, pos + strlen( old ) );
   }
   return out;
}

/* As check_rejected, but for the result of editing the dump "dump" with
   the function above.  A NULL result is reported as a failure. */
static int check_edit( const char *what, const char *dump, const char *old,
                       const char *new, int expected ) {
   char *text = edit( dump, old, new );
   int nfail;
   if ( !text ) {
      fprintf( stderr, "testcorruptinput: could not build the dump for "
               "%s\n", what );
      if ( !astOK ) astClearStatus;
      return 1;
   }
   nfail = check_rejected( what, text, expected );
   text = astFree( text );
   return nfail;
}

/* The dump of a UnitMap, for use where an object of another class is
   needed, and the dumps of Frames with one and two axes. */
#define UNITMAP " Begin UnitMap\n Nin = 2\n IsA Mapping\n End UnitMap\n"
#define FRAME1 " Begin Frame\n Naxes = 1\n Ax1 =\n Begin Axis\n End Axis\n" \
               " End Frame\n"
#define FRAME2 " Begin Frame\n Naxes = 2\n Ax1 =\n Begin Axis\n End Axis\n" \
               " Ax2 =\n Begin Axis\n End Axis\n End Frame\n"

/* Integers that cannot be held in an int. */
static int test_integers( void ) {
   int nfail = 0;

   nfail += check_accepted( "Frame with the largest int as MaxAxes",
      " Begin Frame\n Naxes = 1\n MaxAx = 2147483647\n End Frame\n" );
   nfail += check_accepted( "Frame with the smallest int as MatchEnd",
      " Begin Frame\n Naxes = 1\n MchEnd = -2147483648\n End Frame\n" );
   nfail += check_rejected( "UnitMap with Nin one more than the largest int",
      " Begin UnitMap\n Nin = 2147483648\n IsA Mapping\n End UnitMap\n",
      AST__BADIN );
   nfail += check_rejected( "UnitMap with Nin that wraps to a small int",
      " Begin UnitMap\n Nin = 4294967298\n IsA Mapping\n End UnitMap\n",
      AST__BADIN );
   nfail += check_rejected( "UnitMap with Nin beyond a 64 bit int",
      " Begin UnitMap\n Nin = 99999999999999999999\n IsA Mapping\n"
      " End UnitMap\n", AST__BADIN );
   nfail += check_rejected( "UnitMap with Nin one less than the smallest int",
      " Begin UnitMap\n Nin = -2147483649\n IsA Mapping\n End UnitMap\n",
      AST__BADIN );
   nfail += check_rejected( "UnitMap with a floating point Nin",
      " Begin UnitMap\n Nin = 1e308\n IsA Mapping\n End UnitMap\n",
      AST__BADIN );
   nfail += check_rejected( "UnitMap with text after Nin",
      " Begin UnitMap\n Nin = 2 3\n IsA Mapping\n End UnitMap\n",
      AST__BADIN );
   return nfail;
}

/* The limit on the number of axes. */
static int test_axis_limit( void ) {
   AstFrame *big;
   AstFrame *half;
   char dump[ 200 ];
   int nfail = 0;

/* Constructors accept the limit and reject one more. */
   nfail += check_made( "astUnitMap with AST__MXAXES inputs",
                        astUnitMap( AST__MXAXES, " " ) );
   nfail += check_failed( "astUnitMap with too many inputs",
                          astUnitMap( AST__MXAXES + 1, " " ), AST__BADNI );
   nfail += check_failed( "astUnitMap with INT_MAX inputs",
                          astUnitMap( 2147483647, " " ), AST__BADNI );
   nfail += check_failed( "astFrame with too many axes",
                          astFrame( AST__MXAXES + 1, " " ), AST__NAXIN );
   nfail += check_failed( "astFrame with INT_MAX axes",
                          astFrame( 2147483647, " " ), AST__NAXIN );

/* Two Frames that are each within the limit cannot be combined into a
   CmpFrame that is beyond it. */
   big = astFrame( AST__MXAXES, " " );
   half = astFrame( AST__MXAXES/2, " " );
   nfail += check_made( "astCmpFrame with AST__MXAXES axes",
                        astCmpFrame( half, half, " " ) );
   nfail += check_failed( "astCmpFrame with too many axes",
                          astCmpFrame( big, half, " " ), AST__NAXIN );
   if ( big ) big = astAnnul( big );
   if ( half ) half = astAnnul( half );
   if ( !astOK ) astClearStatus;

/* Loaders do the same. */
   snprintf( dump, sizeof( dump ), " Begin UnitMap\n Nin = %d\n"
             " IsA Mapping\n End UnitMap\n", AST__MXAXES );
   nfail += check_accepted( "UnitMap dump with AST__MXAXES inputs", dump );
   snprintf( dump, sizeof( dump ), " Begin UnitMap\n Nin = %d\n"
             " IsA Mapping\n End UnitMap\n", AST__MXAXES + 1 );
   nfail += check_rejected( "UnitMap dump with too many inputs", dump,
                            AST__BADNI );
   snprintf( dump, sizeof( dump ), " Begin PermMap\n Nin = 1\n Nout = %d\n"
             " IsA Mapping\n End PermMap\n", AST__MXAXES + 1 );
   nfail += check_rejected( "PermMap dump with too many outputs", dump,
                            AST__BADNO );
   nfail += check_rejected( "ShiftMap dump with INT_MAX inputs",
      " Begin ShiftMap\n Nin = 2147483647\n IsA Mapping\n End ShiftMap\n",
      AST__BADNI );
   snprintf( dump, sizeof( dump ), " Begin Frame\n Naxes = %d\n End Frame\n",
             AST__MXAXES + 1 );
   nfail += check_rejected( "Frame dump with too many axes", dump,
                            AST__NAXIN );
   nfail += check_rejected( "Frame dump with INT_MAX axes",
      " Begin Frame\n Naxes = 2147483647\n End Frame\n", AST__NAXIN );
   nfail += check_rejected( "TimeFrame dump with INT_MAX axes",
      " Begin TimeFrame\n Naxes = 2147483647\n IsA Frame\n"
      " End TimeFrame\n", AST__NAXIN );
   return nfail;
}

/* Component objects that are missing or of the wrong class. */
static int test_components( void ) {
   int nfail = 0;

   nfail += check_accepted( "CmpFrame of two Frames",
      " Begin CmpFrame\n FrameA =\n" FRAME2 " FrameB =\n" FRAME1
      " End CmpFrame\n" );
   nfail += check_rejected( "CmpFrame with a UnitMap as its first Frame",
      " Begin CmpFrame\n FrameA =\n" UNITMAP " FrameB =\n" FRAME1
      " End CmpFrame\n", AST__LDERR );
   nfail += check_rejected( "CmpFrame with a UnitMap as its second Frame",
      " Begin CmpFrame\n FrameA =\n" FRAME2 " FrameB =\n" UNITMAP
      " End CmpFrame\n", AST__LDERR );
   nfail += check_rejected( "CmpFrame with the smallest int as an axis",
      " Begin CmpFrame\n FrameA =\n" FRAME2 " FrameB =\n" FRAME1
      " Axp1 = -2147483648\n End CmpFrame\n", AST__PRMIN );
   nfail += check_rejected( "Frame with a UnitMap as an Axis",
      " Begin Frame\n Naxes = 1\n Ax1 =\n" UNITMAP " End Frame\n",
      AST__LDERR );
   nfail += check_rejected( "Frame with the smallest int as an axis",
      " Begin Frame\n Naxes = 1\n Prm1 = -2147483648\n End Frame\n",
      AST__PRMIN );
   nfail += check_rejected( "Frame with a UnitMap as its variants",
      " Begin Frame\n Naxes = 1\n Vrnts =\n" UNITMAP " End Frame\n",
      AST__LDERR );
   nfail += check_rejected( "UnitMap with a Frame as its KeyMap",
      " Begin UnitMap\n KeyMap =\n" FRAME1 " IsA Object\n Nin = 2\n"
      " IsA Mapping\n End UnitMap\n", AST__LDERR );
   nfail += check_rejected( "FluxFrame with a Frame as its SpecFrame",
      " Begin FluxFrame\n Naxes = 1\n IsA Frame\n SpcFr =\n" FRAME1
      " End FluxFrame\n", AST__LDERR );
   nfail += check_rejected( "Table with a UnitMap as its columns",
      " Begin Table\n IsA KeyMap\n Columns =\n" UNITMAP " End Table\n",
      AST__LDERR );
   nfail += check_rejected( "Table with a UnitMap as its parameters",
      " Begin Table\n IsA KeyMap\n Params =\n" UNITMAP " End Table\n",
      AST__LDERR );
   nfail += check_rejected( "FitsTable with no header",
      " Begin FitsTable\n IsA KeyMap\n IsA Table\n End FitsTable\n",
      AST__LDERR );
   return nfail;
}

/* Regions. */
static int test_regions( void ) {
   AstFrame *frm1 = astFrame( 1, " " );
   AstFrame *frm2 = astFrame( 2, " " );
   AstRegion *box;
   AstRegion *circle;
   AstRegion *interval;
   char *dump;
   double p1[ 2 ] = { 0.0, 0.0 };
   double p2[ 2 ] = { 1.0, 1.0 };
   double p3[ 2 ] = { 0.5, 0.0 };
   double poly[ 6 ] = { 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 };
   double radius = 1.0;
   int nfail = 0;

/* The dump of a Circle holds two points with two coordinates each. */
   circle = (AstRegion *) astCircle( frm2, 1, p1, &radius, NULL, " " );
   dump = astToString( circle );
   nfail += check_accepted( "Circle", dump );
   nfail += check_edit( "Circle with no points", dump,
      " Points =\n Begin PointSet\n Npoint = 2\n Ncoord = 2\n X1 = 0\n"
      " X2 = 0\n X3 = 1\n X4 = 0\n End PointSet\n", "", AST__LDERR );
   nfail += check_edit( "Circle with one point", dump, " Npoint = 2\n",
                        " Npoint = 1\n", AST__LDERR );
   nfail += check_edit( "Circle with three points", dump, " Npoint = 2\n",
                        " Npoint = 3\n", AST__LDERR );
   nfail += check_edit( "Circle with INT_MAX points", dump, " Npoint = 2\n",
                        " Npoint = 2147483647\n", AST__LDERR );
   nfail += check_edit( "Circle with points of three coordinates", dump,
                        " Ncoord = 2\n", " Ncoord = 3\n", AST__REGIN );
   nfail += check_edit( "Circle with a UnitMap as its points", dump,
      " Begin PointSet\n Npoint = 2\n Ncoord = 2\n X1 = 0\n X2 = 0\n"
      " X3 = 1\n X4 = 0\n End PointSet\n", UNITMAP, AST__REGIN );
   nfail += check_edit( "Circle with a UnitMap as its Frame", dump,
                        " Frm =\n" FRAME2, " Frm =\n" UNITMAP, AST__LDERR );
   nfail += check_edit( "Circle with a Frame as its FrameSet", dump,
                        " Frm =\n", " FrmSet =\n", AST__LDERR );
   nfail += check_edit( "Circle with a number as its uncertainty", dump,
                        " IsA Region\n", " Unc = 0\n IsA Region\n",
                        AST__BADIN );
   nfail += check_edit( "Circle with a UnitMap as its uncertainty", dump,
                        " IsA Region\n", " Unc =\n" UNITMAP " IsA Region\n",
                        AST__LDERR );
   nfail += check_edit( "Circle with an uncertainty on one axis", dump,
      " IsA Region\n",
      " Unc =\n Begin Interval\n Points =\n Begin PointSet\n Npoint = 2\n"
      " Ncoord = 1\n X1 = 0\n X2 = 1\n End PointSet\n IsA Region\n"
      " End Interval\n IsA Region\n", AST__REGIN );
   nfail += check_edit( "Circle with an uncertainty in a Frame of one axis",
      dump, " IsA Region\n",
      " Unc =\n Begin Interval\n Frm =\n" FRAME1 " Points =\n"
      " Begin PointSet\n Npoint = 2\n Ncoord = 1\n X1 = 0\n X2 = 1\n"
      " End PointSet\n IsA Region\n End Interval\n IsA Region\n",
      AST__LDERR );
   nfail += check_edit( "Circle with a System in its Frame data", dump,
                        " Begin Circle\n",
                        " Begin Circle\n System = \"FK5\"\n IsA Frame\n",
                        AST__ATTIN );
   dump = astFree( dump );

/* A Box and an Interval also need two points, an Ellipse needs three, and
   a Polygon and a PointList need at least one. */
   box = (AstRegion *) astBox( frm2, 1, p1, p2, NULL, " " );
   dump = astToString( box );
   nfail += check_accepted( "Box", dump );
   nfail += check_edit( "Box with one point", dump, " Npoint = 2\n",
                        " Npoint = 1\n", AST__LDERR );
   nfail += check_edit( "Box with no Region data", dump, " IsA Region\n", "",
                        AST__LDERR );
   dump = astFree( dump );

   interval = (AstRegion *) astInterval( frm1, p1, p2, NULL, " " );
   dump = astToString( interval );
   nfail += check_accepted( "Interval", dump );
   nfail += check_edit( "Interval with three points", dump, " Npoint = 2\n",
                        " Npoint = 3\n", AST__LDERR );
   dump = astFree( dump );

   dump = astToString( astEllipse( frm2, 1, p1, p2, p3, NULL, " " ) );
   nfail += check_accepted( "Ellipse", dump );
   nfail += check_edit( "Ellipse with two points", dump, " Npoint = 3\n",
                        " Npoint = 2\n", AST__LDERR );
   dump = astFree( dump );

   dump = astToString( astPolygon( frm2, 3, 3, poly, NULL, " " ) );
   nfail += check_accepted( "Polygon", dump );
   nfail += check_edit( "Polygon with a UnitMap as its points", dump,
                        " Points =\n Begin PointSet\n", " Points =\n" UNITMAP
                        " Junk =\n Begin PointSet\n", AST__REGIN );
   nfail += check_rejected( "Polygon with no points",
      " Begin Polygon\n Frm =\n" FRAME2 " IsA Region\n End Polygon\n",
      AST__LDERR );
   nfail += check_rejected( "PointList with no points",
      " Begin PointList\n Frm =\n" FRAME2 " IsA Region\n End PointList\n",
      AST__LDERR );
   nfail += check_rejected( "Polygon on one axis",
      " Begin Polygon\n Frm =\n" FRAME1 " Points =\n Begin PointSet\n"
      " Npoint = 3\n Ncoord = 1\n X1 = 0\n X2 = 1\n X3 = 2\n"
      " End PointSet\n IsA Region\n End Polygon\n", AST__LDERR );

/* A Prism has all the axes of both its components. */
   dump = astToString( astPrism( circle, interval, " " ) );
   nfail += check_accepted( "Prism", dump );
   nfail += check_edit( "Prism with no second Region", dump,
                        " RegionB =\n", " Junk =\n", AST__LDERR );
   nfail += check_edit( "Prism with a UnitMap as its first Region", dump,
                        " RegionA =\n", " RegionA =\n" UNITMAP " Junk =\n",
                        AST__LDERR );
   nfail += check_edit( "Prism with too many axes", dump,
                        " FrameB =\n" FRAME1, " FrameB =\n" FRAME2,
                        AST__LDERR );
   nfail += check_edit( "Prism with a component that has too few "
                        "coordinates", dump, " Ncoord = 2\n X1 = 0\n",
                        " Ncoord = 1\n X1 = 0\n", AST__LDERR );
   nfail += check_rejected( "NullRegion with INT_MAX axes",
      " Begin NullRegion\n RegAxes = 2147483647\n IsA Region\n"
      " End NullRegion\n", AST__NAXIN );
   dump = astFree( dump );

/* Both components of a CmpRegion have the axes of the CmpRegion. */
   dump = astToString( astCmpRegion( circle, box, AST__OR, " " ) );
   nfail += check_accepted( "CmpRegion", dump );
   nfail += check_edit( "CmpRegion with an unknown operator", dump,
                        " Operator = 2\n", " Operator = 7\n", AST__LDERR );
   nfail += check_edit( "CmpRegion with no second Region", dump,
                        " RegionB =\n", " Junk =\n", AST__LDERR );
   nfail += check_edit( "CmpRegion with a UnitMap as its first Region", dump,
                        " RegionA =\n", " RegionA =\n" UNITMAP " Junk =\n",
                        AST__LDERR );
   nfail += check_edit( "CmpRegion with a component on one axis", dump,
      " RegionA =\n",
      " RegionA =\n Begin Interval\n Points =\n Begin PointSet\n"
      " Npoint = 2\n Ncoord = 1\n X1 = 0\n X2 = 1\n End PointSet\n"
      " IsA Region\n End Interval\n Junk =\n", AST__LDERR );
   dump = astFree( dump );

/* An Stc encapsulates a Region with the same axes, and describes each of
   its AstroCoords with a KeyMap. */
   dump = astToString( astStcSearchLocation( circle, 0, NULL, " " ) );
   nfail += check_accepted( "StcSearchLocation", dump );
   nfail += check_edit( "Stc with no Region", dump, " Region =\n",
                        " Junk =\n", AST__LDERR );
   nfail += check_edit( "Stc with a UnitMap as its Region", dump,
                        " Region =\n", " Region =\n" UNITMAP " Junk =\n",
                        AST__LDERR );
   nfail += check_edit( "Stc with INT_MAX AstroCoords", dump, " IsA Stc\n",
                        " Ncoord = 2147483647\n IsA Stc\n", AST__LDERR );
   nfail += check_edit( "Stc with a UnitMap as an AstroCoords", dump,
                        " IsA Stc\n", " Ncoord = 1\n Coord1 =\n" UNITMAP
                        " IsA Stc\n", AST__LDERR );
   dump = astFree( dump );

/* A PointSet cannot hold more values than a dump can name. */
   nfail += check_accepted( "PointSet with two bad points",
      " Begin PointSet\n Npoint = 2\n Ncoord = 1\n End PointSet\n" );
   nfail += check_rejected( "PointSet with INT_MAX points",
      " Begin PointSet\n Npoint = 2147483647\n Ncoord = 1\n"
      " End PointSet\n", AST__LDERR );
   nfail += check_rejected( "PointSet with INT_MAX coordinates",
      " Begin PointSet\n Npoint = 1\n Ncoord = 2147483647\n"
      " End PointSet\n", AST__LDERR );
   nfail += check_rejected( "PointSet with INT_MAX points and coordinates",
      " Begin PointSet\n Npoint = 2147483647\n Ncoord = 2147483647\n"
      " End PointSet\n", AST__LDERR );

   if ( !astOK ) astClearStatus;
   return nfail;
}

/* FrameSets. */
static int test_framesets( void ) {
   AstFrameSet *fs = astFrameSet( astFrame( 2, " " ), " " );
   char *dump;
   int nfail = 0;

/* A FrameSet with two Frames of two axes, joined by a UnitMap. */
   astAddFrame( fs, AST__BASE, astUnitMap( 2, " " ),
                astFrame( 2, "Domain=B" ) );
   dump = astToString( fs );
   nfail += check_accepted( "FrameSet", dump );

/* Frames and Mappings that are not there. */
   nfail += check_edit( "FrameSet with three Frames stated", dump,
                        " Nframe = 2\n", " Nframe = 3\n", AST__LDERR );
   nfail += check_edit( "FrameSet with INT_MAX Frames stated", dump,
                        " Nframe = 2\n", " Nframe = 2147483647\n",
                        AST__LDERR );
   nfail += check_edit( "FrameSet with INT_MAX nodes stated", dump,
                        " Nframe = 2\n",
                        " Nframe = 2\n Nnode = 2147483647\n", AST__LDERR );
   nfail += check_edit( "FrameSet with no Mapping", dump, " Map2 =\n",
                        " Junk =\n", AST__LDERR );
   nfail += check_edit( "FrameSet with a UnitMap as a Frame", dump,
                        " Frm1 =\n", " Frm1 =\n" UNITMAP " Junk =\n",
                        AST__LDERR );

/* Node links that are out of range or form a loop. */
   nfail += check_edit( "FrameSet with the smallest int as a link", dump,
                        " Lnk2 = 1\n", " Lnk2 = -2147483648\n",
                        AST__FRSIN );
   nfail += check_edit( "FrameSet with a node derived from itself", dump,
                        " Lnk2 = 1\n", " Lnk2 = 2\n", AST__FRSIN );
   nfail += check_edit( "FrameSet with two nodes derived from each other",
                        dump, " Lnk2 = 1\n",
                        " Nnode = 3\n Lnk2 = 3\n Lnk3 = 2\n Map3 =\n"
                        UNITMAP, AST__FRSIN );

/* Mappings and Frames that disagree about the number of axes. */
   nfail += check_edit( "FrameSet with a Mapping of three inputs", dump,
                        " Begin UnitMap\n Nin = 2\n",
                        " Begin UnitMap\n Nin = 3\n", AST__FRSIN );
   nfail += check_edit( "FrameSet with a Frame of one axis", dump,
                        " Frm2 =\n", " Frm2 =\n" FRAME1 " Junk =\n",
                        AST__FRSIN );
   nfail += check_edit( "FrameSet with Mappings that disagree", dump,
                        " Lnk2 = 1\n",
                        " Nnode = 3\n Lnk2 = 1\n Lnk3 = 1\n Map3 =\n"
                        " Begin UnitMap\n Nin = 3\n IsA Mapping\n"
                        " End UnitMap\n", AST__FRSIN );
   dump = astFree( dump );

   if ( !astOK ) astClearStatus;
   return nfail;
}

/* Frames that always have a particular number of axes. */
static int test_frames( void ) {
   AstFrame *frm;
   const char *text;
   int nfail = 0;

   nfail += check_rejected( "SkyFrame with one axis",
      " Begin SkyFrame\n Naxes = 1\n IsA Frame\n End SkyFrame\n",
      AST__LDERR );
   nfail += check_rejected( "SkyFrame with a negative number of axes",
      " Begin SkyFrame\n Naxes = -2\n IsA Frame\n End SkyFrame\n",
      AST__LDERR );
   nfail += check_rejected( "SkyFrame with three axes",
      " Begin SkyFrame\n Naxes = 3\n IsA Frame\n End SkyFrame\n",
      AST__LDERR );
   nfail += check_rejected( "SpecFrame with two axes",
      " Begin SpecFrame\n Naxes = 2\n IsA Frame\n End SpecFrame\n",
      AST__LDERR );
   nfail += check_rejected( "TimeFrame with no axes",
      " Begin TimeFrame\n Naxes = 0\n IsA Frame\n End TimeFrame\n",
      AST__LDERR );
   nfail += check_rejected( "FluxFrame with two axes",
      " Begin FluxFrame\n Naxes = 2\n IsA Frame\n End FluxFrame\n",
      AST__LDERR );
   nfail += check_rejected( "SpecFluxFrame made from two simple Frames",
      " Begin SpecFluxFrame\n IsA Frame\n FrameA =\n" FRAME1 " FrameB =\n"
      FRAME1 " IsA CmpFrame\n End SpecFluxFrame\n", AST__LDERR );

/* A SkyFrame read from a dump that does not describe its axes still has
   sky axes, so it formats a value as hours, minutes and seconds. */
   frm = (AstFrame *) astFromString( " Begin SkyFrame\n Naxes = 2\n"
                                     " IsA Frame\n End SkyFrame\n" );
   text = frm ? astFormat( frm, 1, 1.0 ) : NULL;
   if ( !astOK || !text || !strchr( text, ':' ) ) {
      fprintf( stderr, "testcorruptinput: a SkyFrame read without its axes "
               "formatted a longitude as \"%s\"\n", text ? text : "" );
      nfail++;
   }
   if ( !astOK ) astClearStatus;
   if ( frm ) frm = astAnnul( frm );
   return nfail;
}

/* Mappings. */
static int test_mappings( void ) {
   AstLutMap *lutmap;
   double element = 1.0;
   double far[ 2 ] = { 1.0E308, -1.0E308 };
   double farout[ 2 ];
   double lut[ 3 ] = { 1.0, 2.0, 4.0 };
   int nfail = 0;

/* A loader must not use what it failed to read. Here the parent loader
   fails because the dump ends with the wrong class name. */
   nfail += check_rejected( "WcsMap that ends as a WinMap",
      " Begin WcsMap\n Nin = 2\n IsA Mapping\n End WinMap\n",
      AST__BADIN );
   nfail += check_rejected( "UnitNormMap that ends as an Axis",
      " Begin UnitNormMap\n End Axis\n", AST__BADIN );
   nfail += check_rejected( "UnitNormMap with no inputs",
      " Begin UnitNormMap\n Nin = 0\n Nout = 1\n IsA Mapping\n"
      " End UnitNormMap\n", AST__LDERR );

/* A PermMap refers to its constants by number. */
#define PERMMAP( extra ) " Begin PermMap\n Nin = 2\n IsA Mapping\n" \
   " Out1 = -1\n Out2 = 1\n In1 = 2\n In2 = 0\n" extra " End PermMap\n"
   nfail += check_accepted( "PermMap with one constant",
                            PERMMAP( " Nconst = 1\n Con1 = 5\n" ) );
   nfail += check_accepted( "PermMap with no constants", PERMMAP( "" ) );
   nfail += check_rejected( "PermMap that refers to a second constant",
      " Begin PermMap\n Nin = 2\n IsA Mapping\n Out1 = -2\n Out2 = 1\n"
      " In1 = 2\n In2 = 0\n Nconst = 1\n Con1 = 5\n End PermMap\n",
      AST__LDERR );
   nfail += check_rejected( "PermMap that refers to constant INT_MAX + 1",
      " Begin PermMap\n Nin = 2\n IsA Mapping\n Out1 = -1\n Out2 = 1\n"
      " In1 = -2147483648\n In2 = 0\n End PermMap\n", AST__LDERR );
   nfail += check_rejected( "PermMap with INT_MAX constants",
                            PERMMAP( " Nconst = 2147483647\n" ), AST__LDERR );
#undef PERMMAP

/* A count of conversion steps far larger than the number of steps in the
   dump. The arrays that describe the steps were sized from the count. */
   nfail += check_rejected( "SpecMap with INT_MAX steps stated",
      " Begin SpecMap\n Nin = 1\n IsA Mapping\n Nspec = 2147483647\n"
      " End SpecMap\n", AST__BADIN );
   nfail += check_rejected( "SlaMap with INT_MAX steps stated",
      " Begin SlaMap\n Nin = 2\n IsA Mapping\n Nsla = 2147483647\n"
      " End SlaMap\n", AST__BADIN );
   nfail += check_rejected( "TimeMap with INT_MAX steps stated",
      " Begin TimeMap\n Nin = 1\n IsA Mapping\n Ntime = 2147483647\n"
      " End TimeMap\n", AST__BADIN );
   nfail += check_accepted( "SlaMap with one step",
      " Begin SlaMap\n Nin = 2\n IsA Mapping\n Nsla = 1\n"
      " Sla1 = \"FK45Z\"\n Sla1a = 1950\n End SlaMap\n" );

/* The elements of a full matrix are counted with an int, and a dump can
   name fewer than ten million of them. */
   nfail += check_accepted( "MatrixMap with four elements",
      " Begin MatrixMap\n Nin = 2\n IsA Mapping\n M0 = 1\n M1 = 2\n"
      " M2 = 3\n M3 = 4\n Form = \"Full\"\n End MatrixMap\n" );
   nfail += check_rejected( "MatrixMap with 2**32 elements",
      " Begin MatrixMap\n Nin = 65536\n IsA Mapping\n Form = \"Full\"\n"
      " End MatrixMap\n", AST__LDERR );
   nfail += check_rejected( "MatrixMap with sixteen million elements",
      " Begin MatrixMap\n Nin = 4000\n IsA Mapping\n Form = \"Full\"\n"
      " End MatrixMap\n", AST__LDERR );
   nfail += check_failed( "astMatrixMap with 2**32 elements",
                          astMatrixMap( 65536, 65536, 0, &element, " " ),
                          AST__TOOBG );
   nfail += check_made( "astMatrixMap with 65536 diagonal elements",
                        astMatrixMap( 65536, 65536, 2, NULL, " " ) );

/* An axis index that cannot be converted to a zero-based index. */
   nfail += check_rejected( "WcsMap with the smallest int as an axis",
      " Begin WcsMap\n Nin = 2\n IsA Mapping\n Type = \"TAN\"\n"
      " WcsAx1 = -2147483648\n End WcsMap\n", AST__WCSAX );

/* A LutMap finds a table index from its input, which can be far beyond
   the range of an int. */
   lutmap = astLutMap( 3, lut, 0.0, 1.0, " " );
   astTran1( lutmap, 2, far, 1, farout );
   astSet( lutmap, "LutInterp=1" );
   astTran1( lutmap, 2, far, 1, farout );
   nfail += check_made( "LutMap given inputs of 1.0E308", lutmap );

/* The size of the hash table of a KeyMap is not taken on trust. */
   nfail += check_accepted( "KeyMap that asks for a table of INT_MAX",
      " Begin KeyMap\n MapSz = 2147483647\n End KeyMap\n" );
   nfail += check_accepted( "KeyMap with a size guess of INT_MAX",
      " Begin KeyMap\n SzGss = 2147483647\n End KeyMap\n" );
   return nfail;
}

/* Values too large to treat one turn at a time. */
static int test_large_values( void ) {
   AstFrameSet *fs;
   AstSkyFrame *azel;
   AstSkyFrame *sky = astSkyFrame( " " );
   char symbol[ 201 ];
   char *dump;
   const char *text;
   double xin = 1.0;
   double xout;
   double yin = 0.5;
   double yout;
   int nfail = 0;

/* A sky position of 1.0E300 radians has far more hours than the format
   buffer has room for. */
   text = astFormat( sky, 1, 1.0E300 );
   if ( !astOK || !text || strcmp( text, "<bad>" ) ) {
      fprintf( stderr, "testcorruptinput: a longitude of 1.0E300 was "
               "formatted as \"%s\"\n", text ? text : "" );
      nfail++;
   }
   text = astFormat( sky, 2, -1.0E300 );
   if ( !astOK || !text || strcmp( text, "<bad>" ) ) {
      fprintf( stderr, "testcorruptinput: a latitude of -1.0E300 was "
               "formatted as \"%s\"\n", text ? text : "" );
      nfail++;
   }
   text = astFormat( sky, 1, 1.0 );
   if ( !astOK || !text || strcmp( text, "3:49:11.0" ) ) {
      fprintf( stderr, "testcorruptinput: a longitude of 1.0 was formatted "
               "as \"%s\"\n", text ? text : "" );
      nfail++;
   }

/* The comments that a SkyFrame adds to its dump include the Symbol of
   each axis and the formatted reference position. */
   memset( symbol, 'x', sizeof( symbol ) - 1 );
   symbol[ sizeof( symbol ) - 1 ] = '\0';
   astSetC( sky, "Symbol(1)", symbol );
   astSetD( sky, "SkyRef(1)", 1.0 );
   astSetD( sky, "SkyRef(2)", 0.5 );
   astSetD( sky, "SkyRefP(1)", 2.0 );
   astSetD( sky, "SkyRefP(2)", 0.25 );
   dump = astToString( sky );
   nfail += check_accepted( "SkyFrame with a long Symbol", dump );
   dump = astFree( dump );

/* The local apparent sidereal time for an Epoch far in the future is a
   very large number of turns. */
   azel = astSkyFrame( "System=AZEL,Epoch=99999999999" );
   fs = astConvert( azel, sky, " " );
   if ( fs ) astTran2( fs, 1, &xin, &yin, 1, &xout, &yout );
   if ( !astOK ) astClearStatus;

   return nfail;
}

/* Return a dump made by joining "head", "nfill" copies of the character
   "fill", and "tail".  It should be freed using astFree. */
static char *long_dump( const char *head, int nfill, char fill,
                        const char *tail ) {
   size_t nh = strlen( head );
   char *dump = astMalloc( nh + (size_t) nfill + strlen( tail ) + 1 );
   if ( dump ) {
      memcpy( dump, head, nh );
      memset( dump + nh, fill, (size_t) nfill );
      strcpy( dump + nh + nfill, tail );
   }
   return dump;
}

/* Text that is longer than the buffers it used to be copied into. */
static int test_long_text( void ) {
   AstFrame *frm;
   AstFrame *frm2;
   AstFrame *cmp;
   AstObject *obj;
   char *dump;
   const char *text;
   int i;
   int nfail = 0;

/* A value that cannot be read as an integer is quoted in the error
   message, which was formatted on the stack in 1024 characters. */
   dump = long_dump( " Begin Frame\n Naxes = ", 5000, '7', "\n End Frame\n" );
   nfail += check_rejected( "Frame with 5000 digits as its Naxes", dump,
                            AST__BADIN );
   dump = astFree( dump );
   dump = long_dump( " Begin Frame\n Naxes = 1\n System = \"", 5000, 'S',
                     "\"\n End Frame\n" );
   nfail += check_rejected( "Frame with 5000 characters as its System", dump,
                            AST__ATTIN );
   dump = astFree( dump );

/* An item that is not recognised is quoted in a warning, which was
   formatted on the stack in 201 characters. */
   dump = long_dump( " Begin Frame\n Naxes = 1\n Junk = \"", 5000, 'J',
                     "\"\n End Frame\n" );
   obj = astFromString( dump );
   if ( !astOK ) astClearStatus;
   if ( obj ) obj = astAnnul( obj );
   dump = astFree( dump );

/* The default Domain of a CmpFrame joins the Domains of its components,
   and its default Title and Labels are made in the same way, all in
   buffers of a fixed size. */
   dump = long_dump( " Begin Frame\n Naxes = 1\n Domain = \"", 300, 'D',
                     "\"\n Lbl1 = \"A long label\"\n Uni1 = \"km\"\n"
                     " End Frame\n" );
   frm = (AstFrame *) astFromString( dump );
   dump = astFree( dump );
   cmp = frm ? (AstFrame *) astCmpFrame( frm, frm, " " ) : NULL;
   for ( i = 0; i < 4 && cmp; i++ ) {
      frm2 = (AstFrame *) astCmpFrame( cmp, frm, " " );
      cmp = astAnnul( cmp );
      cmp = frm2;
   }
   text = cmp ? astGetC( cmp, "Domain" ) : NULL;
   if ( !astOK || !text || strlen( text ) > 100 || text[ 0 ] != 'D' ) {
      fprintf( stderr, "testcorruptinput: the Domain of a CmpFrame of Frames "
               "with long Domains was \"%s\"\n", text ? text : "" );
      nfail++;
   }
   if ( cmp ) {
      (void) astGetC( cmp, "Title" );
      (void) astGetC( cmp, "Label(1)" );
      dump = astToString( cmp );
      dump = astFree( dump );
      cmp = astAnnul( cmp );
   }
   if ( frm ) frm = astAnnul( frm );
   if ( !astOK ) astClearStatus;

/* The default Title of a SpecFluxFrame joins the Labels of its
   components, and that of a FluxFrame includes a Unit. */
   dump = long_dump( " Begin SpecFluxFrame\n IsA Frame\n FrameA =\n"
                     " Begin SpecFrame\n Naxes = 1\n Lbl1 = \"", 400, 'L',
                     "\"\n IsA Frame\n End SpecFrame\n FrameB =\n"
                     " Begin FluxFrame\n Naxes = 1\n IsA Frame\n"
                     " End FluxFrame\n IsA CmpFrame\n End SpecFluxFrame\n" );
   obj = astFromString( dump );
   dump = astFree( dump );
   if ( obj ) {
      (void) astGetC( obj, "Title" );
      obj = astAnnul( obj );
   }
   if ( !astOK ) astClearStatus;

/* A Format is a "printf" format for the one axis value. A Format that asks
   for anything else is refused when it is used, as is one that needs more
   room than there is. */
   frm = astFrame( 1, " " );
   astSetC( frm, "Format(1)", "%.3lf metres" );
   text = astFormat( frm, 1, 1.5 );
   if ( !astOK || !text || strcmp( text, "1.500 metres" ) ) {
      fprintf( stderr, "testcorruptinput: 1.5 was formatted as \"%s\"\n",
               text ? text : "" );
      nfail++;
   }
   if ( !astOK ) astClearStatus;
   astSetC( frm, "Format(1)", "%5d%%" );
   text = astFormat( frm, 1, 42.0 );
   if ( !astOK || !text || strcmp( text, "   42%" ) ) {
      fprintf( stderr, "testcorruptinput: 42 was formatted as \"%s\"\n",
               text ? text : "" );
      nfail++;
   }
   if ( !astOK ) astClearStatus;
   {
      const char *bad[] = { "%s", "%n", "%s%s%s%s", "%.2f and %s", "%*d",
                            "%Lf", "%ld", "%.2f %.2f", "%p", "%.500f",
                            "%500.1f" };
      for ( i = 0; i < (int) ( sizeof( bad )/sizeof( bad[ 0 ] ) ); i++ ) {
         astSetC( frm, "Format(1)", bad[ i ] );
         if ( !astOK ) astClearStatus;
         text = astFormat( frm, 1, 1.5 );
         if ( astOK || astStatus != AST__FMTER ) {
            fprintf( stderr, "testcorruptinput: the Format \"%s\" gave "
                     "status %d and \"%s\"\n", bad[ i ], astStatus,
                     text ? text : "" );
            nfail++;
         }
         if ( !astOK ) astClearStatus;
      }
   }
   frm = astAnnul( frm );

/* The same Format read from a dump. */
   frm = (AstFrame *) astFromString( " Begin Frame\n Naxes = 1\n"
                                     " Fmt1 = \"%s%n\"\n End Frame\n" );
   if ( frm ) {
      (void) astFormat( frm, 1, 1.5 );
      nfail += check_failed( "astFormat with a Format of %s%n from a dump",
                             NULL, AST__FMTER );
      frm = astAnnul( frm );
   }
   if ( !astOK ) astClearStatus;
   return nfail;
}

/* Splitting a Mapping in which two of many axes are coupled. */
static int test_split( void ) {
#define NAX 12
   AstCmpMap *cmp;
   AstMapping *map;
   AstMapping *part;
   double matrix[ NAX*NAX ];
   int inperm[ NAX ];
   int outperm[ 2 ];
   int in[ 2 ];
   int out[ NAX ];
   int i;
   int nfail = 0;
   int which;

/* The first of two Mappings in series has NAX axes, of which the second
   and third are coupled, and the rest are independent. It is a MatrixMap
   that rotates them, and then a WcsMap for which they are longitude and
   latitude. The second Mapping keeps just those two axes. */
   for ( which = 0; which < 2; which++ ) {
      if ( which == 0 ) {
         for ( i = 0; i < NAX*NAX; i++ ) matrix[ i ] = 0.0;
         for ( i = 0; i < NAX; i++ ) matrix[ i*NAX + i ] = 1.0;
         matrix[ 1*NAX + 2 ] = 0.5;
         matrix[ 2*NAX + 1 ] = -0.5;
         map = (AstMapping *) astMatrixMap( NAX, NAX, 0, matrix, " " );
      } else {
         map = (AstMapping *) astWcsMap( NAX, AST__TAN, 2, 3, " " );
      }
      for ( i = 0; i < NAX; i++ ) inperm[ i ] = 0;
      inperm[ 1 ] = 1;
      inperm[ 2 ] = 2;
      outperm[ 0 ] = 2;
      outperm[ 1 ] = 3;
      cmp = astCmpMap( map, astPermMap( NAX, inperm, 2, outperm, NULL, " " ),
                       1, " " );

/* Neither coupled axis can be split off alone. */
      in[ 0 ] = 2;
      astMapSplit( cmp, 1, in, out, &part );
      if ( part || !astOK ) {
         fprintf( stderr, "testcorruptinput: astMapSplit split one of two "
                  "coupled axes from a %s\n", astGetC( map, "Class" ) );
         nfail++;
         if ( part ) part = astAnnul( part );
      }

/* The two can be split off together, and feed both outputs. */
      in[ 0 ] = 2;
      in[ 1 ] = 3;
      astMapSplit( cmp, 2, in, out, &part );
      if ( !part || !astOK || astGetI( part, "Nout" ) != 2 ||
           out[ 0 ] != 1 || out[ 1 ] != 2 ) {
         fprintf( stderr, "testcorruptinput: astMapSplit did not split two "
                  "coupled axes from a %s\n", astGetC( map, "Class" ) );
         nfail++;
      }
      if ( part ) part = astAnnul( part );
      if ( !astOK ) astClearStatus;
   }
   return nfail;
#undef NAX
}

/* Read an object from a FITS header given as an array of "ncard" cards.
   Returns NULL if nothing could be read. */
static AstObject *read_header( const char *cards[], int ncard ) {
   AstFitsChan *fc = astFitsChan( NULL, NULL, " " );
   AstObject *obj;
   int i;
   for ( i = 0; i < ncard; i++ ) astPutFits( fc, cards[ i ], 0 );
   astClear( fc, "Card" );
   obj = astRead( fc );
   fc = astAnnul( fc );
   return obj;
}

/* Read a FITS header that need only be handled without a memory fault,
   and in a reasonable time.  Transform a position with the result if
   there is one, and write it out again as a dump and as a FITS header.
   Any error reported is cleared.  Always returns 0, since a fault or a
   time limit stops the test. */
static int check_survives( const char *cards[], int ncard ) {
   AstFitsChan *fc;
   AstObject *obj = read_header( cards, ncard );
   char *dump;
   double xin = 10.0;
   double xout;
   double yin = 20.0;
   double yout;
   if ( obj && astOK && astIsAMapping( obj ) &&
        astGetI( obj, "Nin" ) == 2 && astGetI( obj, "Nout" ) == 2 ) {
      astTran2( obj, 1, &xin, &yin, 1, &xout, &yout );
      if ( !astOK ) astClearStatus;
      astTran2( obj, 1, &xin, &yin, 0, &xout, &yout );
   }
   if ( !astOK ) astClearStatus;
   if ( obj ) {
      dump = astToString( obj );
      dump = astFree( dump );
      if ( !astOK ) astClearStatus;
      fc = astFitsChan( NULL, NULL, " " );
      astWrite( fc, obj );
      fc = astAnnul( fc );
      obj = astAnnul( obj );
   }
   if ( !astOK ) astClearStatus;
   return 0;
}

#define NCARD( cards ) ( (int) ( sizeof( cards )/sizeof( cards[ 0 ] ) ) )

/* FITS headers. */
static int test_fits( void ) {
   int nfail = 0;

/* A celestial header to which one card is added. */
   const char *sky[] = {
      "CTYPE1  = 'RA---TAN'",
      "CTYPE2  = 'DEC--TAN'",
      "CRVAL1  =                180.0",
      "CRVAL2  =                 45.0",
      "CRPIX1  =                256.5",
      "CRPIX2  =                256.5",
      "CDELT1  =              -0.0003",
      "CDELT2  =               0.0003",
      NULL };
   const int nsky = NCARD( sky ) - 1;

/* WCSAXES values that are not a number of axes that FITS can describe.
   These were used as the number of axes of several Frames. */
   const char *badaxes[] = {
      "WCSAXES =           2147483647",
      "WCSAXES =                 1000",
      "WCSAXES =              1.0E300",
      "WCSAXES =                    0",
      "WCSAXES =                   -3" };

/* A SIP header that states one WCS axis, where the distortion needs
   two. */
   const char *sip[] = {
      "WCSAXES =                    1",
      "CTYPE1  = 'RA---TAN-SIP'",
      "CTYPE2  = 'DEC--TAN-SIP'",
      "CRVAL1  =                180.0",
      "CRVAL2  =                 45.0",
      "CRPIX1  =                256.5",
      "CRPIX2  =                256.5",
      "CD1_1   =              -0.0003",
      "CD2_2   =               0.0003",
      "A_ORDER =                    2",
      "A_2_0   =               1.0E-6",
      "B_ORDER =                    2",
      "B_0_2   =               1.0E-6" };

/* CTYPE values that are too short to hold a celestial axis type. */
   const char *ctype0[] = {
      "CTYPE1  = ''",
      "CRVAL1  =                  1.0",
      "CRPIX1  =                  1.0",
      "CDELT1  =                  1.0" };
   const char *ctype1[] = {
      "CTYPE1  = 'X'",
      "CRVAL1  =                  1.0",
      "CRPIX1  =                  1.0",
      "CDELT1  =                  1.0" };

/* An IRAF TNX header in which a WAT card has no value, or a value that is
   not a string. */
   const char *tnx[] = {
      "CTYPE1  = 'RA---TNX'",
      "CTYPE2  = 'DEC--TNX'",
      "CRVAL1  =                110.5",
      "CRVAL2  =                 55.8",
      "CRPIX1  =               4189.8",
      "CRPIX2  =                -30.1",
      "CD1_1   =              -7.6E-7",
      "CD2_1   =              -1.2E-4",
      "CD1_2   =               1.2E-4",
      "CD2_2   =              -1.3E-6",
      "WAT0_001= 'system=image'",
      "WAT1_001= 'wtype=tnx axtype=ra lngcor = \"3. 2. 2. 2. 0.005 0.49 '",
      NULL,
      "WAT1_003= '0.25 0.49 -0.2 0.005 -0.001 \"'",
      "WAT2_001= 'wtype=tnx axtype=dec latcor = \"3. 2. 2. 2. 0.005 0.49 '",
      "WAT2_002= '0.25 0.49 0.15 -0.003 6.5E-4 \"'" };

/* A header with two longitude axes, which is rejected. */
   const char *lon2[] = {
      "CTYPE1  = 'RA---NCP'",
      "CTYPE2  = 'DEC--NCP'",
      "CTYPE3  = 'RA---TAN'",
      "CRVAL1  =                180.0",
      "CRVAL2  =                 45.0",
      "CRPIX1  =                256.5",
      "CRPIX2  =                256.5",
      "CDELT1  =              -0.0003",
      "CDELT2  =               0.0003" };

/* A header with a reference latitude that has 302 digits of degrees, and
   one with an equinox that puts the sidereal time far out of range. */
   const char *bigval[] = {
      "CTYPE1  = 'RA---CAR'",
      "CTYPE2  = 'DEC--CAR'",
      "CRVAL1  =                  0.0",
      "CRVAL2  =              1.0E300",
      "CRPIX1  =                  1.0",
      "CRPIX2  =                  1.0",
      "CDELT1  =                 -1.0",
      "CDELT2  =                  1.0" };
/* A header with a pixel size so small that its reciprocal is infinite.
   Writing it never returned, because the simplification of a Mapping
   stops when a step leaves the Mapping equal to what it was, and two
   Mappings that hold the same infinite value did not compare equal. */
   const char *tiny[] = {
      "CTYPE1  = 'WAVE'",
      "CRVAL1  =               5.5E-7",
      "CRPIX1  =                513.0",
      "CDELT1  =               1e-320",
      "CUNIT1  = 'm'" };
/* A logarithmic axis with a negative reference value that has a three
   digit exponent, which needs one more character than was allowed for
   in the text of the expression used to describe the axis. */
   const char *logax[] = {
      "CTYPE1  = 'FREQ-LOG'",
      "CRVAL1  = -1.2345678901234567E-300",
      "CRPIX1  =                512.5",
      "CDELT1  =                  1.0",
      "CUNIT1  = 'Hz'" };

/* A spectral axis with a rest frequency so small that the first step
   used to find a rate of change, a fraction of it, is zero. Reading it
   never returned. */
   const char *tinyrf[] = {
      "CTYPE1  = 'BETA-FRQ'",
      "CRVAL1  =               1.0E-6",
      "CRPIX1  =                512.5",
      "CDELT1  =               1.0E-8",
      "RESTFRQ =               1e-320",
      "SPECSYS = 'BARYCENT'" };

/* A cube with two spectral axes, for which no reference point can be
   found when it is written out as a FITS header. */
   const char *spec2[] = {
      "CTYPE1  = 'FREQ'",
      "CTYPE2  = 'PIXEL'",
      "CTYPE3  = 'FREQ-F2W'",
      "CRVAL1  =             1.4204E9",
      "CRVAL2  =                  0.0",
      "CRVAL3  =                  0.0",
      "CRPIX1  =                513.0",
      "CRPIX2  =                 50.0",
      "CRPIX3  =                 50.0",
      "CDELT1  =           2147483648",
      "CDELT2  =                  1.0",
      "CDELT3  =                  1.0",
      "CUNIT1  = 'Hz'",
      "SPECSYS = 'BARYCENT'",
      "RESTFRQ =             1.4204E9" };
/* More WCS axes than pixel axes, which FITS allows. The two pixel axes
   are rotated, or are celestial axes. Writing either out again used to
   take a time that grew faster than the factorial of the number of WCS
   axes, and was minutes for nine of them. */
   const char *rotated[] = {
      "WCSAXES =                   20",
      "CTYPE1  = 'PIXEL'",
      "CTYPE2  = 'PIXEL'",
      "CRVAL1  =                  0.0",
      "CRVAL2  =                  0.0",
      "CRPIX1  =                 50.0",
      "CRPIX2  =                 50.0",
      "CD1_1   =                  1.0",
      "CD1_2   =                  0.5",
      "CD2_1   =                 -0.5",
      "CD2_2   =                  1.0" };
   const char *skyaxes[] = {
      "WCSAXES =                   20",
      "CTYPE1  = 'RA---TAN'",
      "CTYPE2  = 'DEC--TAN'",
      "CRVAL1  =                180.0",
      "CRVAL2  =                 45.0",
      "CRPIX1  =                256.5",
      "CRPIX2  =                256.5",
      "CDELT1  =              -0.0003",
      "CDELT2  =               0.0003" };
   const char *bigeqx[] = {
      "CTYPE1  = 'RA---TAN'",
      "CTYPE2  = 'DEC--TAN'",
      "CRVAL1  =                180.0",
      "CRVAL2  =                 45.0",
      "CRPIX1  =                256.5",
      "CRPIX2  =                256.5",
      "CDELT1  =              -0.0003",
      "CDELT2  =               0.0003",
      "RADESYS = 'FK5'",
      "EQUINOX =           2147483647" };
   AstFitsChan *fc;
   int i;
   int ival;

   sky[ nsky ] = "WCSAXES =                    2";
   nfail += check_made( "FITS header with WCSAXES = 2",
                        read_header( sky, nsky + 1 ) );
   for ( i = 0; i < NCARD( badaxes ); i++ ) {
      sky[ nsky ] = badaxes[ i ];
      nfail += check_failed( badaxes[ i ], read_header( sky, nsky + 1 ),
                             AST__BDFTS );
   }

   nfail += check_survives( sip, NCARD( sip ) );
   nfail += check_survives( ctype0, NCARD( ctype0 ) );
   nfail += check_survives( ctype1, NCARD( ctype1 ) );
   tnx[ 12 ] = "WAT1_002=";
   nfail += check_survives( tnx, NCARD( tnx ) );
   tnx[ 12 ] = "WAT1_002=                    5";
   nfail += check_survives( tnx, NCARD( tnx ) );
   nfail += check_survives( lon2, NCARD( lon2 ) );
   nfail += check_survives( bigval, NCARD( bigval ) );
   nfail += check_survives( bigeqx, NCARD( bigeqx ) );
   nfail += check_survives( tiny, NCARD( tiny ) );
   nfail += check_survives( rotated, NCARD( rotated ) );
   nfail += check_survives( skyaxes, NCARD( skyaxes ) );
   nfail += check_survives( logax, NCARD( logax ) );
   nfail += check_survives( tinyrf, NCARD( tinyrf ) );
   nfail += check_survives( spec2, NCARD( spec2 ) );

/* A floating point keyword value that an int cannot hold is not
   converted to one. */
   fc = astFitsChan( NULL, NULL, " " );
   astPutFits( fc, "NAXIS1  =              -1.0E300", 0 );
   astPutFits( fc, "NAXIS2  =                100.0", 0 );
   astClear( fc, "Card" );
   ival = 0;
   if ( astGetFitsI( fc, "NAXIS2", &ival ) != 1 || ival != 100 || !astOK ) {
      fprintf( stderr, "testcorruptinput: NAXIS2 = 100.0 was not read as "
               "the integer 100\n" );
      nfail++;
   }
   if ( !astOK ) astClearStatus;
   ival = 0;
   if ( astGetFitsI( fc, "NAXIS1", &ival ) && astOK ) {
      fprintf( stderr, "testcorruptinput: NAXIS1 = -1.0E300 was read as "
               "the integer %d\n", ival );
      nfail++;
   }
   if ( !astOK ) astClearStatus;
   fc = astAnnul( fc );
   return nfail;
}

int main( void ) {
   int status = SAI__OK;
   int nfail = 0;
   astWatch( &status );
   astBegin;

   nfail += test_integers();
   nfail += test_axis_limit();
   nfail += test_components();
   nfail += test_regions();
   nfail += test_framesets();
   nfail += test_frames();
   nfail += test_mappings();
   nfail += test_large_values();
   nfail += test_long_text();
   nfail += test_split();
   nfail += test_fits();

   astEnd;
   if ( !astOK ) astClearStatus;

   if ( nfail == 0 ) {
      printf( " All corrupt input tests passed\n" );
      return 0;
   }
   printf( "Corrupt input tests failed (%d)\n", nfail );
   return 1;
}
