/*
*  Name:
*     testregionhardening

*  Purpose:
*     Check that unusual Regions do not make the library overrun its arrays.

*  Description:
*     A Region is made from positions and Frames supplied by the caller,
*     and can be mapped into a Frame with a different number of axes, so
*     the functions that work with it cannot assume that it has any extent,
*     that it is bounded, that its base and current Frames have the same
*     number of axes, or that a value derived from it will fit in an int.
*     Each case here makes a Region through the public interface and uses
*     it in a way that used to read or write beyond the end of an array,
*     call a function without end, or loop for ever.  The test requires
*     the stated result, or the stated error.  Most of the faults do not
*     change a result, so the test is intended to be run with the address
*     and undefined-behavior sanitizers enabled.  The cases cover:
*
*       - Grids: the grid of points inside a Box that has no width on an
*         axis, or is far longer than it is wide, where the number of
*         points on an axis was unlimited.
*       - Attributes: a string value that contains a "%" character, given
*         to a Region that passes the setting on to its FrameSet.
*       - Masks: masking an array with a Region that has no upper limit,
*         and with a PointList when the grid has the wrong number of
*         dimensions.
*       - Overlaps: a Region with more current Frame axes than base Frame
*         axes, a PointList, two Intervals of which one excludes the range
*         that the other includes, and a Prism that has no boundary mesh.
*       - Axes: the bounds of a PointList and of a NullRegion with fewer
*         current Frame axes than base Frame axes, and the mesh of a Prism
*         with more.
*       - Ellipses: an Ellipse used as the uncertainty of another Region,
*         which moves it, three points that do not define an ellipse, and
*         a Circle mapped so that no point of its boundary can be used to
*         fit an Ellipse.
*       - Polygons: vertices that all coincide, and vertices so far apart
*         that the search for a point inside the Polygon had no end.
*       - Fewer axes: a Box, a Circle and a PointList mapped into a Frame
*         that has fewer axes, after which the Mapping supplied must be
*         unchanged.

*  Notes:
*     - This is a native C test with no Fortran original.
*     - Public C API only.  An expected error is caught through astStatus
*       and cleared with astClearStatus so the remaining cases can run;
*       failures of the test itself are counted and reported at the end.
*/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "sae_par.h"
#include "ast.h"
#include "ast_err.h"

/* Check that the status is "expected" after the operation described by
   "what".  Returns 0 if so, and 1 otherwise.  Any raised status is cleared
   before returning. */
static int check_status( const char *what, int expected ) {
   int code = astStatus;
   if ( !astOK ) astClearStatus;
   if ( code != expected ) {
      fprintf( stderr, "testregionhardening: %s gave status %d but status "
               "%d was expected\n", what, code, expected );
      return 1;
   }
   return 0;
}

/* Check that the integer "value" found by the operation described by
   "what" equals "expected", and that no error was reported.  Returns 0 if
   so, and 1 otherwise. */
static int check_int( const char *what, int value, int expected ) {
   if ( check_status( what, SAI__OK ) ) return 1;
   if ( value != expected ) {
      fprintf( stderr, "testregionhardening: %s gave %d but %d was "
               "expected\n", what, value, expected );
      return 1;
   }
   return 0;
}

/* Check that "value" is within "tol" of "expected".  Returns 0 if so, and
   1 otherwise. */
static int check_double( const char *what, double value, double expected,
                         double tol ) {
   if ( !( fabs( value - expected ) <= tol ) ) {
      fprintf( stderr, "testregionhardening: %s gave %.17g but %.17g was "
               "expected\n", what, value, expected );
      return 1;
   }
   return 0;
}

/* Get the grid of points inside the Box with corners "p1" and "p2" in the
   Frame "frm", which has "naxes" axes, and check that the status is
   "expected".  If the grid was made, check that it has at least "minpoint"
   points and, if "inside" is non-zero, that every point is within the Box.
   The status found is returned in "*code" if "code" is not NULL.  Returns
   the number of failures. */
static int check_grid( const char *what, void *frm, int naxes,
                       const double *p1, const double *p2, int expected,
                       int minpoint, int inside, int *code ) {
   AstRegion *box;
   double *pts;
   double hi;
   double lo;
   double tol;
   int i;
   int j;
   int nfail = 0;
   int np = 0;

   box = (AstRegion *) astBox( frm, 1, p1, p2, NULL, "MeshSize=30" );
   astGetRegionMesh( box, 0, 0, naxes, &np, NULL );
   if ( code ) {
      *code = astStatus;
      if ( !astOK ) astClearStatus;
   } else if ( check_status( what, expected ) ) {
      nfail++;
   } else if ( expected == SAI__OK ) {
      if ( np < minpoint ) {
         fprintf( stderr, "testregionhardening: %s has %d points but at "
                  "least %d were expected\n", what, np, minpoint );
         nfail++;
      }
      pts = astMalloc( sizeof( double )*(size_t)( np*naxes ) );
      astGetRegionMesh( box, 0, np, naxes, &np, pts );
      if ( astOK && inside ) {
         for ( j = 0; j < naxes && !nfail; j++ ) {
            lo = ( p1[ j ] < p2[ j ] ) ? p1[ j ] : p2[ j ];
            hi = ( p1[ j ] < p2[ j ] ) ? p2[ j ] : p1[ j ];
            tol = 1.0E-9*( hi - lo );
            for ( i = 0; i < np; i++ ) {
               if ( !( pts[ j*np + i ] >= lo - tol &&
                       pts[ j*np + i ] <= hi + tol ) ) {
                  fprintf( stderr, "testregionhardening: %s has a point "
                           "at %.17g on axis %d, outside the Box\n", what,
                           pts[ j*np + i ], j + 1 );
                  nfail++;
                  break;
               }
            }
         }
      }
      nfail += check_status( what, SAI__OK );
      pts = astFree( pts );
   }
   box = astAnnul( box );
   return nfail;
}

/* The number of points on each axis of a grid is the nominal number
   scaled by the length of the axis divided by the geometric mean of the
   lengths of the axes that share its Unit.  An axis of zero length made
   that mean zero, so the numbers were infinite or not a number, and were
   converted to an int and multiplied together. */
static int test_grids( void ) {
   double a[ 3 ];
   double b[ 3 ];
   int narrow;
   int nfail = 0;
   int nowidth;

   astBegin;

/* A Box with no width on its last axis.  The other axes are still covered
   by more than the two points that each axis must have. */
   a[ 0 ] = 0.0; a[ 1 ] = 0.0; a[ 2 ] = 0.0;
   b[ 0 ] = 1.0; b[ 1 ] = 0.0; b[ 2 ] = 0.0;
   nfail += check_grid( "The grid of a 2-D Box with no width on axis 2",
                        astFrame( 2, " " ), 2, a, b, SAI__OK, 6, 1, NULL );
   nfail += check_grid( "The grid of a 2-D Box with no width on axis 2, in "
                        "a Frame whose axes have different Units",
                        astFrame( 2, "Unit(1)=m,Unit(2)=s" ), 2, a, b,
                        SAI__OK, 6, 1, NULL );
   b[ 1 ] = 1.0;
   nfail += check_grid( "The grid of a 3-D Box with no width on axis 3",
                        astFrame( 3, " " ), 3, a, b, SAI__OK, 8, 1, NULL );
   nfail += check_grid( "The grid of a 3-D Box with no width on axis 3, in "
                        "a CmpFrame",
                        astCmpFrame( astFrame( 2, " " ), astFrame( 1, " " ),
                                     " " ),
                        3, a, b, SAI__OK, 8, 1, NULL );

/* A Box so much longer than it is wide that it would need more points
   than a PointSet can hold. */
   b[ 0 ] = 1.0E12; b[ 1 ] = 1.0E-12;
   nfail += check_grid( "The grid of a Box 1E12 long and 1E-12 wide",
                        astFrame( 2, " " ), 2, a, b, AST__ATTIN, 0, 0, NULL );

/* On the sky, a Box with no height gets a grid, as before.  A Box with no
   width used to need an unlimited number of parallels, and is now treated
   as a very narrow Box is. */
   a[ 0 ] = 0.5; a[ 1 ] = 0.1;
   b[ 0 ] = 0.6; b[ 1 ] = 0.1;
   nfail += check_grid( "The grid of a sky Box with no height",
                        astSkyFrame( " " ), 2, a, b, SAI__OK, 1, 0, NULL );
   a[ 1 ] = -0.1;
   b[ 0 ] = 0.5 + 1.0E-7;
   nfail += check_grid( "The grid of a narrow sky Box", astSkyFrame( " " ),
                        2, a, b, SAI__OK, 0, 0, &narrow );
   b[ 0 ] = 0.5;
   nfail += check_grid( "The grid of a sky Box with no width",
                        astSkyFrame( " " ), 2, a, b, SAI__OK, 0, 0,
                        &nowidth );
   if ( nowidth != narrow ) {
      fprintf( stderr, "testregionhardening: The grid of a sky Box with no "
               "width gave status %d but a narrow one gave status %d\n",
               nowidth, narrow );
      nfail++;
   }

   astEnd;
   return nfail;
}

/* A Region passes the setting for an attribute that it does not recognise
   to its FrameSet, and used to pass it as the format string for astSet,
   so that a "%" in the value was taken to be a format specifier. */
static int test_attributes( void ) {
   AstRegion *box;
   const char *value;
   double p1[ 2 ] = { 0.0, 0.0 };
   double p2[ 2 ] = { 1.0, 1.0 };
   int adaptive;
   int nfail = 0;

   astBegin;
   for ( adaptive = 0; adaptive < 2; adaptive++ ) {
      box = (AstRegion *) astBox( astFrame( 2, " " ), 1, p1, p2, NULL,
                                  "Adaptive=%d", adaptive );
      astSetC( box, "Title", "50%s done, 100%n to do" );
      value = astGetC( box, "Title" );
      if ( astOK && strcmp( value, "50%s done, 100%n to do" ) ) {
         fprintf( stderr, "testregionhardening: The Title of a Box with "
                  "Adaptive=%d is \"%s\"\n", adaptive, value );
         nfail++;
      }
      astSet( box, "%s", "Format(1)=%10.3f" );
      value = astGetC( box, "Format(1)" );
      if ( astOK && strcmp( value, "%10.3f" ) ) {
         fprintf( stderr, "testregionhardening: Format(1) of a Box with "
                  "Adaptive=%d is \"%s\"\n", adaptive, value );
         nfail++;
      }
      nfail += check_status( "Setting the Title and Format(1) of a Box",
                             SAI__OK );
   }
   astEnd;
   return nfail;
}

/* Mask the pixels 1 to 10 with the Region "reg" and the Mapping "map",
   and check that the number of pixels changed is "count" and that the
   pixels changed are those marked "1" in the string "wanted".  Returns
   the number of failures. */
static int check_mask( const char *what, void *reg, void *map, int inside,
                       int count, const char *wanted ) {
   char got[ 11 ];
   double d[ 10 ];
   int i;
   int lbnd[ 1 ] = { 1 };
   int n;
   int ubnd[ 1 ] = { 10 };

   for ( i = 0; i < 10; i++ ) d[ i ] = 0.0;
   n = astMaskD( reg, map, inside, 1, lbnd, ubnd, d, 9.0 );
   for ( i = 0; i < 10; i++ ) got[ i ] = ( d[ i ] == 9.0 ) ? '1' : '0';
   got[ 10 ] = 0;
   if ( check_int( what, n, count ) ) return 1;
   if ( strcmp( got, wanted ) ) {
      fprintf( stderr, "testregionhardening: %s changed pixels %s but %s "
               "was expected\n", what, got, wanted );
      return 1;
   }
   return 0;
}

/* The bounds of a Region in grid coordinates were converted to an int
   without a check.  The upper bound of a Region with no upper limit is
   DBL_MAX, for which the result was a box containing no pixels, so every
   pixel was taken to be outside the Region. */
static int test_masks( void ) {
   AstRegion *above;
   AstRegion *below;
   AstRegion *list;
   double bad[ 1 ] = { AST__BAD };
   double pts[ 4 ] = { 1.0, 3.0, 2.0, 4.0 };
   double seven[ 1 ] = { 7.0 };
   double three[ 1 ] = { 3.0 };
   int d[ 6 ] = { 0, 0, 0, 0, 0, 0 };
   int lbnd[ 1 ] = { 1 };
   int n;
   int nfail = 0;
   int ubnd[ 1 ] = { 6 };

   astBegin;
   above = (AstRegion *) astInterval( astFrame( 1, " " ), three, bad, NULL,
                                      " " );
   below = (AstRegion *) astInterval( astFrame( 1, " " ), bad, seven, NULL,
                                      " " );

/* A pixel on the edge of a closed Region is in neither the Region nor its
   negation, so it is not changed in either case. */
   nfail += check_mask( "Masking outside an Interval with no upper limit",
                        above, NULL, 0, 2, "1100000000" );
   nfail += check_mask( "Masking inside an Interval with no upper limit",
                        above, NULL, 1, 7, "0001111111" );
   nfail += check_mask( "Masking outside an Interval with no lower limit",
                        below, NULL, 0, 3, "0000000111" );
   nfail += check_mask( "Masking inside an Interval with no lower limit",
                        below, NULL, 1, 6, "1111110000" );

/* A PointList uses its own function, which went on to use the Region for
   the grid after reporting that the grid had the wrong number of
   dimensions, when there is no such Region. */
   list = (AstRegion *) astPointList( astFrame( 2, " " ), 2, 2, 2, pts, NULL,
                                      " " );
   n = astMaskI( list, NULL, 1, 1, lbnd, ubnd, d, 1 );
   nfail += check_status( "Masking a 1-D array with a 2-D PointList",
                          AST__NGDIN );
   if ( n != 0 ) {
      fprintf( stderr, "testregionhardening: Masking a 1-D array with a "
               "2-D PointList changed %d pixels\n", n );
      nfail++;
   }

   astEnd;
   return nfail;
}

/* Check that the overlap of "reg1" and "reg2" is "expected", with no
   error.  Returns the number of failures. */
static int check_overlap( const char *what, void *reg1, void *reg2,
                          int expected ) {
   return check_int( what, astOverlap( reg1, reg2 ), expected );
}

static int test_overlaps( void ) {
   AstRegion *and1;
   AstRegion *and2;
   AstRegion *box;
   AstRegion *circle;
   AstRegion *interval;
   AstRegion *list2;
   AstRegion *list3;
   AstRegion *mapped;
   AstRegion *prism;
   AstRegion *unitbox;
   const char *fwd[ 2 ] = { "a = x", "b = 0*x" };
   const char *inv[ 1 ] = { "x = a" };
   double bad[ 2 ] = { AST__BAD, AST__BAD };
   double c1[ 2 ] = { 0.0, 0.0 };
   double c2[ 2 ] = { 5.0, 5.0 };
   double centre[ 3 ] = { 0.5, 0.5, 0.5 };
   double hi[ 2 ] = { 2.0, 1.0 };
   double in1[ 2 ] = { 1.0, -1.0 };
   double in2[ 2 ] = { 2.0, 1.0 };
   double lo[ 2 ] = { 1.0, 0.0 };
   double pts2[ 4 ] = { 1.0, 3.0, 2.0, 4.0 };
   double pts3[ 6 ] = { 1.0, 3.0, 1.0, 2.0, 4.0, 2.0 };
   double radius[ 1 ] = { 0.3 };
   double u1[ 2 ] = { 0.0, 0.0 };
   double u2[ 2 ] = { 1.0, 1.0 };
   double zero[ 1 ] = { 0.0 };
   int nfail = 0;
   int result;

   astBegin;

/* The bounds of the second Region were held in an array with an element
   for each base Frame axis of the first Region, but were then used as
   bounds in its current Frame.  The first Region here has one base Frame
   axis and two current Frame axes, and no boundary mesh, which is when
   the bounds are used.  A MathMap is used so that the mapped Interval is
   not simplified into an Interval with two base Frame axes.  A position
   is inside the mapped Interval if the inverse of the MathMap puts it in
   the original, whatever its second coordinate, so the Box is inside. */
   interval = (AstRegion *) astInterval( astFrame( 1, " " ), zero, bad, NULL,
                                         " " );
   mapped = astMapRegion( interval, astMathMap( 1, 2, 2, fwd, 1, inv, " " ),
                          astFrame( 2, " " ) );
   box = (AstRegion *) astBox( astFrame( 2, " " ), 1, in1, in2, NULL, " " );
   nfail += check_overlap( "The overlap of a half-line mapped into two "
                           "dimensions and a Box", mapped, box, 3 );
   nfail += check_overlap( "The overlap of a Box and a half-line mapped "
                           "into two dimensions", box, mapped, 2 );

/* A PointList returned a mask with a byte for each point where the caller
   expected an int for each point, and went on to check as many of its own
   points as it had been given, however many it has. */
   list2 = (AstRegion *) astPointList( astFrame( 2, " " ), 2, 2, 2, pts2,
                                       NULL, " " );
   list3 = (AstRegion *) astPointList( astFrame( 2, " " ), 3, 2, 3, pts3,
                                       NULL, " " );
   box = (AstRegion *) astBox( astFrame( 2, " " ), 1, c1, c2, NULL, " " );
   nfail += check_overlap( "The overlap of a PointList and a Box that "
                           "contains it", list2, box, 2 );
   nfail += check_overlap( "The overlap of two PointLists that hold the "
                           "same positions, one of them twice", list2,
                           list3, 5 );
   and1 = (AstRegion *) astCmpRegion( list2, box, AST__AND, " " );
   and2 = (AstRegion *) astCmpRegion( astCopy( list2 ), astCopy( box ),
                                      AST__AND, " " );
   nfail += check_overlap( "The overlap of two CmpRegions that each hold a "
                           "PointList", and1, and2, 5 );

/* Two Intervals are compared an axis at a time, and the results for the
   axes combined with a table that has no entry for an axis on which one
   Interval is the negation of the other.  The Intervals as a whole are
   only each other's negation if no other axis is limited. */
   nfail += check_overlap( "The overlap of two 1-D Intervals, one the "
                           "negation of the other",
                           astInterval( astFrame( 1, " " ), hi, lo, NULL,
                                        " " ),
                           astInterval( astFrame( 1, " " ), lo, hi, NULL,
                                        " " ), 6 );
   lo[ 1 ] = AST__BAD;
   hi[ 1 ] = AST__BAD;
   nfail += check_overlap( "The overlap of two 2-D Intervals, one the "
                           "negation of the other on axis 1, with no limit "
                           "on axis 2",
                           astInterval( astFrame( 2, " " ), hi, lo, NULL,
                                        " " ),
                           astInterval( astFrame( 2, " " ), lo, hi, NULL,
                                        " " ), 6 );
   lo[ 1 ] = 0.0;
   hi[ 1 ] = 1.0;
   nfail += check_overlap( "The overlap of two 2-D Intervals, one the "
                           "negation of the other on axis 1, with the same "
                           "range on axis 2",
                           astInterval( astFrame( 2, " " ), hi, lo, NULL,
                                        " " ),
                           astInterval( astFrame( 2, " " ), lo, hi, NULL,
                                        " " ), 1 );

/* The part of the boundary of a Region that is within a box is found from
   the intersection of the Region and the box, which for a Region with no
   boundary mesh asked for the same thing again, without end.  A Prism has
   no mesh if one of its components has no limit on one side. */
   interval = (AstRegion *) astInterval( astFrame( 1, " " ), zero, bad, NULL,
                                         " " );
   unitbox = (AstRegion *) astBox( astFrame( 2, " " ), 1, u1, u2, NULL,
                                   " " );
   prism = (AstRegion *) astPrism( interval, unitbox, " " );
   circle = (AstRegion *) astCircle( astFrame( 3, " " ), 1, centre, radius,
                                     NULL, " " );
   result = astOverlap( circle, prism );
   nfail += check_status( "The overlap of a Circle and a Prism that has no "
                          "mesh", AST__INTER );
   result += astOverlap( prism, circle );
   nfail += check_status( "The overlap of a Prism that has no mesh and a "
                          "Circle", AST__INTER );
   if ( result != 0 ) {
      fprintf( stderr, "testregionhardening: An overlap that could not be "
               "determined was returned as %d\n", result );
      nfail++;
   }

   astEnd;
   return nfail;
}

/* Functions that fill an array with one element for each axis of the base
   Frame, or of the current Frame, must use the number of axes of that
   Frame.  Mapping a Region into a Frame with a different number of axes
   makes the two numbers differ. */
static int test_axes( void ) {
   AstRegion *mapped;
   AstRegion *prism;
   double centre[ 2 ] = { 0.0, 0.0 };
   double con[ 1 ] = { 0.0 };
   double hi[ 1 ] = { 8.0 };
   double lbnd[ 2 ];
   double lo[ 1 ] = { 2.0 };
   double pts[ 4 ] = { 1.0, 3.0, 2.0, 4.0 };
   double radius[ 1 ] = { 5.0 };
   double ubnd[ 2 ];
   int inperm2[ 2 ] = { 1, -1 };
   int inperm3[ 3 ] = { 1, 2, 3 };
   int nfail = 0;
   int np;
   int outperm1[ 1 ] = { 1 };
   int outperm4[ 4 ] = { 1, 2, 3, -1 };

   astBegin;

/* A PointList and a NullRegion with two base Frame axes and one current
   Frame axis.  The second element of each array is not part of the bounds,
   and must be left alone. */
   mapped = astMapRegion( astPointList( astFrame( 2, " " ), 2, 2, 2, pts,
                                        NULL, " " ),
                          astPermMap( 2, inperm2, 1, outperm1, con, " " ),
                          astFrame( 1, " " ) );
   lbnd[ 1 ] = 123.0;
   ubnd[ 1 ] = 456.0;
   astGetRegionBounds( mapped, lbnd, ubnd );
   nfail += check_status( "The bounds of a PointList mapped into one "
                          "dimension", SAI__OK );
   nfail += check_double( "The lower bound of a PointList mapped into one "
                          "dimension", lbnd[ 0 ], 1.0, 1.0E-6 );
   nfail += check_double( "The upper bound of a PointList mapped into one "
                          "dimension", ubnd[ 0 ], 3.0, 1.0E-6 );
   nfail += check_double( "The value after the lower bound of a PointList",
                          lbnd[ 1 ], 123.0, 0.0 );
   nfail += check_double( "The value after the upper bound of a PointList",
                          ubnd[ 1 ], 456.0, 0.0 );

   mapped = astMapRegion( astNullRegion( astFrame( 2, " " ), NULL, " " ),
                          astPermMap( 2, inperm2, 1, outperm1, con, " " ),
                          astFrame( 1, " " ) );
   astGetRegionBounds( mapped, lbnd, ubnd );
   nfail += check_status( "The bounds of a NullRegion mapped into one "
                          "dimension", SAI__OK );
   if ( !( lbnd[ 0 ] > ubnd[ 0 ] ) ) {
      fprintf( stderr, "testregionhardening: A NullRegion mapped into one "
               "dimension has bounds %g to %g\n", lbnd[ 0 ], ubnd[ 0 ] );
      nfail++;
   }
   nfail += check_double( "The value after the lower bound of a NullRegion",
                          lbnd[ 1 ], 123.0, 0.0 );
   nfail += check_double( "The value after the upper bound of a NullRegion",
                          ubnd[ 1 ], 456.0, 0.0 );

/* A Prism with three base Frame axes and four current Frame axes. */
   prism = (AstRegion *) astPrism( astCircle( astFrame( 2, " " ), 1, centre,
                                              radius, NULL, " " ),
                                   astInterval( astFrame( 1, " " ), lo, hi,
                                                NULL, " " ), " " );
   mapped = astMapRegion( prism,
                          astPermMap( 3, inperm3, 4, outperm4, con, " " ),
                          astFrame( 4, " " ) );
   np = 0;
   astGetRegionMesh( mapped, 1, 0, 4, &np, NULL );
   nfail += check_status( "The mesh of a Prism mapped into four dimensions",
                          SAI__OK );
   if ( np < 1 ) {
      fprintf( stderr, "testregionhardening: The mesh of a Prism mapped "
               "into four dimensions has %d points\n", np );
      nfail++;
   }

   astEnd;
   return nfail;
}

static int test_ellipses( void ) {
   AstFrame *frm;
   AstRegion *box;
   AstRegion *ellipse;
   AstRegion *mapped;
   AstRegion *unc;
   const char *fwd[ 2 ] = { "a = qif( x > 0, 5, <bad> )",
                            "b = qif( x > 0, 7, <bad> )" };
   const char *inv[ 2 ] = { "x = a", "y = b" };
   double angle[ 1 ] = { 0.3 };
   double axes[ 2 ] = { 0.2, 0.1 };
   double c1[ 2 ] = { 0.0, 0.0 };
   double c2[ 2 ] = { 10.0, 10.0 };
   double centre[ 2 ];
   double circen[ 2 ] = { 1.0, 0.0 };
   double elbnd[ 2 ];
   double eubnd[ 2 ];
   double lbnd[ 2 ];
   double p1[ 2 ] = { 2.0, 0.0 };
   double p2[ 2 ] = { 1.0, 0.0 };
   double radius[ 1 ] = { 2.0 };
   double ubnd[ 2 ];
   double x[ 1 ] = { 0.5 };
   double xout[ 1 ];
   double y[ 1 ] = { 0.1 };
   double yout[ 1 ];
   int i;
   int k;
   int nfail = 0;

   astBegin;
   frm = astFrame( 2, " " );

/* A Region centres its uncertainty Region on its own centre.  An Ellipse
   moved its three defining points using the point index where the axis
   index was needed, and so used a pointer from beyond the end of an array
   of two.  Check that the Ellipse is the same size after the move, when
   it starts at the origin and when it is already where it is needed. */
   for ( k = 0; k < 2; k++ ) {
      centre[ 0 ] = centre[ 1 ] = k ? 5.0 : 0.0;
      ellipse = (AstRegion *) astEllipse( frm, 1, centre, axes, angle, NULL,
                                          " " );
      astGetRegionBounds( ellipse, elbnd, eubnd );
      box = (AstRegion *) astBox( frm, 1, c1, c2, ellipse, " " );
      unc = astGetUnc( box, 0 );
      astGetRegionBounds( unc, lbnd, ubnd );
      nfail += check_status( "The uncertainty of a Box given as an Ellipse",
                             SAI__OK );
      for ( i = 0; i < 2; i++ ) {
         nfail += check_double( "The centre of an Ellipse used as the "
                                "uncertainty of a Box",
                                0.5*( lbnd[ i ] + ubnd[ i ] ), 5.0, 1.0E-6 );
         nfail += check_double( "The width of an Ellipse used as the "
                                "uncertainty of a Box",
                                ubnd[ i ] - lbnd[ i ],
                                eubnd[ i ] - elbnd[ i ], 1.0E-6 );
      }
   }

/* Three points in a line do not define an ellipse, which is reported when
   the Ellipse is first used.  The values that could not be found were
   then taken to be up to date, and the second use read them through a
   NULL pointer. */
   ellipse = (AstRegion *) astEllipse( frm, 0, c1, p1, p2, NULL, " " );
   nfail += check_status( "Making an Ellipse from three points in a line",
                          SAI__OK );
   for ( k = 0; k < 2; k++ ) {
      astTran2( ellipse, 1, x, y, 1, xout, yout );
      nfail += check_status( "Using an Ellipse made from three points in a "
                             "line", AST__BADIN );
   }

/* A mapped Circle is simplified by fitting a Circle, and then an Ellipse,
   to its boundary mesh in the new Frame.  This Mapping puts the part of
   the mesh with a positive first coordinate at the new centre and gives
   the rest no position, so no mesh point is away from the centre, and the
   search for one ran off the end of the mesh. */
   mapped = astMapRegion( astCircle( frm, 1, circen, radius, NULL, " " ),
                          astMathMap( 2, 2, 2, fwd, 2, inv, " " ),
                          astFrame( 2, " " ) );
   nfail += check_status( "Mapping a Circle so that no Ellipse fits its "
                          "boundary", SAI__OK );
   if ( !mapped ) {
      fprintf( stderr, "testregionhardening: Mapping a Circle so that no "
               "Ellipse fits its boundary gave no Region\n" );
      nfail++;
   }

   astEnd;
   return nfail;
}

static int test_polygons( void ) {
   AstRegion *polygon;
   double far[ 6 ] = { 0.0, 1.0E308, 0.0, 0.0, 0.0, 1.0E308 };
   double inf[ 6 ] = { 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 };
   double point[ 6 ] = { 1.0, 1.0, 1.0, 2.0, 2.0, 2.0 };
   double x[ 2 ] = { 1.0, 50.0 };
   double xout[ 2 ];
   double xout2[ 2 ];
   double y[ 2 ] = { 2.0, 50.0 };
   double yout[ 2 ];
   double yout2[ 2 ];
   int nfail = 0;

   astBegin;

/* A Polygon keeps a position that is inside it, found from the edges that
   have some length.  If the vertices all coincide there is no such edge,
   and the position was never set.  The first vertex is now used, so two
   uses of such a Polygon must agree. */
   polygon = (AstRegion *) astPolygon( astFrame( 2, " " ), 3, 3, point, NULL,
                                       " " );
   astTran2( polygon, 2, x, y, 1, xout, yout );
   polygon = (AstRegion *) astPolygon( astFrame( 2, " " ), 3, 3, point, NULL,
                                       " " );
   astTran2( polygon, 2, x, y, 1, xout2, yout2 );
   nfail += check_status( "Using a Polygon with three vertices at one "
                          "position", SAI__OK );
   if ( xout[ 0 ] != xout2[ 0 ] || yout[ 0 ] != yout2[ 0 ] ||
        xout[ 1 ] != xout2[ 1 ] || yout[ 1 ] != yout2[ 1 ] ) {
      fprintf( stderr, "testregionhardening: Two Polygons with the same "
               "three coincident vertices gave different results\n" );
      nfail++;
   }

/* The position is found by drawing a line across the Polygon from the
   middle of an edge, starting with a line five times as long as the edge
   and halving it until its end is the right distance from its start.  If
   the first length is infinite, no number of halvings makes it finite, so
   the search had no end.  These calls only have to return. */
   polygon = (AstRegion *) astPolygon( astFrame( 2, " " ), 3, 3, far, NULL,
                                       " " );
   if ( !astOK ) astClearStatus;
   inf[ 1 ] = HUGE_VAL;
   polygon = (AstRegion *) astPolygon( astFrame( 2, " " ), 3, 3, inf, NULL,
                                       " " );
   if ( !astOK ) astClearStatus;

   astEnd;
   return nfail;
}

/* When a Region is mapped into a Frame with fewer axes, the base Frame
   axes that are not needed are removed if the Mapping allows it.  This
   released the same Mapping twice, which is an error, and left the Mapping
   supplied by the caller inverted.  It also took the axes picked from the
   Region to be a Region, which the axes picked from a Circle are not. */
static int test_fewer_axes( void ) {
   AstMapping *copy;
   AstMapping *map;
   AstRegion *mapped;
   AstRegion *regs[ 3 ];
   const char *names[ 3 ] = { "Box", "Circle", "PointList" };
   const char *outside[ 3 ] = { "1001110000", "1111110000", "1000110000" };
   double c1[ 2 ] = { 1.5, 2.5 };
   double c2[ 2 ] = { 3.5, 4.5 };
   double centre[ 2 ] = { 3.0, 3.5 };
   double pts[ 6 ] = { 2.0, 3.0, 4.0, 1.0, 5.0, 9.0 };
   double radius[ 1 ] = { 1.6 };
   int counts[ 3 ] = { 4, 6, 3 };
   int d[ 6 ];
   int i;
   int inperm[ 2 ] = { 1, 0 };
   int j;
   int lbnd[ 1 ] = { 1 };
   int n;
   int nfail = 0;
   int outperm[ 1 ] = { 1 };
   int ubnd[ 1 ] = { 6 };
   char got[ 11 ];

   astBegin;
   regs[ 0 ] = (AstRegion *) astBox( astFrame( 2, " " ), 1, c1, c2, NULL,
                                     " " );
   regs[ 1 ] = (AstRegion *) astCircle( astFrame( 2, " " ), 1, centre,
                                        radius, NULL, " " );
   regs[ 2 ] = (AstRegion *) astPointList( astFrame( 2, " " ), 3, 2, 3, pts,
                                           NULL, " " );

/* The Mapping keeps the first axis, and its inverse gives the second axis
   no value. */
   map = (AstMapping *) astPermMap( 2, inperm, 1, outperm, NULL, " " );
   copy = astCopy( map );
   for ( i = 0; i < 3; i++ ) {
      mapped = astMapRegion( regs[ i ], map, astFrame( 1, " " ) );
      if ( check_status( "Mapping a Region into one dimension", SAI__OK ) ||
           !mapped ) {
         fprintf( stderr, "testregionhardening: A %s could not be mapped "
                  "into one dimension\n", names[ i ] );
         nfail++;
      }
      if ( !astEqual( map, copy ) ) {
         fprintf( stderr, "testregionhardening: The Mapping used to map a "
                  "%s into one dimension was changed\n", names[ i ] );
         nfail++;
      }

/* Masking with the same Mapping maps the Region in the same way. */
      for ( j = 0; j < 6; j++ ) d[ j ] = 0;
      n = astMaskI( regs[ i ], map, 0, 1, lbnd, ubnd, d, 1 );
      for ( j = 0; j < 6; j++ ) got[ j ] = d[ j ] ? '1' : '0';
      strcpy( got + 6, "0000" );
      if ( check_int( "Masking outside a Region mapped into one dimension",
                      n, counts[ i ] ) || strcmp( got, outside[ i ] ) ) {
         fprintf( stderr, "testregionhardening: Masking outside a %s "
                  "mapped into one dimension changed %d pixels, %s\n",
                  names[ i ], n, got );
         nfail++;
      }
      if ( !astEqual( map, copy ) ) {
         fprintf( stderr, "testregionhardening: The Mapping used to mask "
                  "with a %s was changed\n", names[ i ] );
         nfail++;
      }
   }

   astEnd;
   return nfail;
}

int main( void ) {
   int status = SAI__OK;
   int nfail = 0;

   astWatch( &status );

   nfail += test_grids();
   nfail += test_attributes();
   nfail += test_masks();
   nfail += test_overlaps();
   nfail += test_axes();
   nfail += test_ellipses();
   nfail += test_polygons();
   nfail += test_fewer_axes();

   if ( !astOK ) astClearStatus;

   if ( nfail == 0 ) {
      printf( " All Region hardening tests passed\n" );
      return 0;
   }
   printf( "Region hardening tests failed (%d)\n", nfail );
   return 1;
}
