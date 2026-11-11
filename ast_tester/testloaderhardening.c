/*
*  Name:
*     testloaderhardening

*  Purpose:
*     Check that object loaders validate structural data read from a dump.

*  Description:
*     Several class loaders read counts, sizes and indices straight from a
*     dump and used them without the range checks that the corresponding
*     constructor applies.  A hand-crafted dump could therefore build an
*     object whose internal arrays are mutually inconsistent, which led to
*     out-of-bounds access when the object was later used (or, for Moc, to an
*     undefined shift at load time).

*     This test confirms that a crafted dump which violates each invariant is
*     rejected cleanly at load time, and that a normally-constructed object
*     of the same class still loads.  The crafted cases cover:
*
*       - PolyMap: an axis power larger than the stated maximum power.
*       - LutMap:  fewer than two lookup table elements.
*       - WcsMap:  a latitude axis index beyond the number of axes.
*       - SplineMap: too few coefficients for the spline order.
*       - Plot:    a base Frame with more axes than the fixed per-axis arrays.
*       - Moc:     ranges present but no MaxOrder.

*  Notes:
*     - This is a native C test with no Fortran original.
*     - Public C API only.  Errors expected from a deliberately invalid dump
*       are caught through astOK and cleared with astClearStatus so the
*       remaining cases can run; failures of the test itself are counted and
*       reported at the end.
*/

#include <stdio.h>
#include <string.h>
#include "sae_par.h"
#include "ast.h"
#include "ast_err.h"

/* Attempt to load a dump that should be rejected.  Returns 0 if it was
   rejected cleanly (NULL object, error status raised), 1 otherwise.  Any
   raised status is cleared before returning so the caller can continue. */
static int check_rejected( const char *what, const char *dump ) {
   AstObject *obj = astFromString( dump );
   int bad_status = !astOK;
   if ( obj ) obj = astAnnul( obj );
   if ( !astOK ) astClearStatus;
   if ( obj != NULL || !bad_status ) {
      fprintf( stderr, "testloaderhardening: %s was accepted but should "
               "have been rejected\n", what );
      return 1;
   }
   return 0;
}

/* Build a valid object, write it to a dump, and confirm that the dump loads
   again as the expected class.  The supplied object is not annulled.
   Returns 0 on success, 1 on failure. */
static int check_valid( const char *what, AstObject *obj, const char *class ) {
   char *dump;
   AstObject *obj2;
   int nfail = 0;

   if ( !astOK || !obj ) {
      fprintf( stderr, "testloaderhardening: could not build %s\n", what );
      if ( !astOK ) astClearStatus;
      return 1;
   }

   dump = astToString( obj );
   obj2 = astFromString( dump );
   if ( !astOK || !obj2 ) {
      fprintf( stderr, "testloaderhardening: valid %s failed to load\n", what );
      if ( !astOK ) astClearStatus;
      nfail = 1;
   } else if ( strcmp( astGetC( obj2, "Class" ), class ) != 0 ) {
      fprintf( stderr, "testloaderhardening: valid %s loaded as %s\n", what,
               astGetC( obj2, "Class" ) );
      nfail = 1;
   }
   if ( obj2 ) obj2 = astAnnul( obj2 );
   dump = astFree( dump );
   return nfail;
}

int main( void ) {
   int status = SAI__OK;
   int nfail = 0;
   astWatch( &status );

/* PolyMap: a coefficient power must not exceed the maximum power read for
   that axis, or the transformation indexes its work array out of bounds. */
   nfail += check_rejected( "PolyMap power exceeding max power",
      " Begin PolyMap\n Nin = 1\n IsA Mapping\n"
      " MPF1 = 2\n NCF1 = 1\n CF1 = 1\n PF1 = 1000\n IterInv = 0\n"
      " End PolyMap\n" );
   {
      double cf[] = { 1.0, 1.0, 2.0 };
      AstPolyMap *pm = astPolyMap( 1, 1, 1, cf, 0, NULL, "IterInverse=0" );
      nfail += check_valid( "PolyMap", (AstObject *) pm, "PolyMap" );
      if ( pm ) pm = astAnnul( pm );
   }

/* LutMap: at least two table elements are needed for interpolation. */
   nfail += check_rejected( "LutMap with one element",
      " Begin LutMap\n Nin = 1\n IsA Mapping\n"
      " Nlut = 1\n LutInt = 0\n L1 = 1\n End LutMap\n" );
   nfail += check_rejected( "LutMap with zero increment",
      " Begin LutMap\n Nin = 1\n IsA Mapping\n"
      " Nlut = 2\n Incr = 0\n L1 = 1\n L2 = 2\n End LutMap\n" );
   {
      double lut[ 4 ] = { 1.0, 2.0, 4.0, 8.0 };
      AstLutMap *m = astLutMap( 4, lut, 0.0, 1.0, " " );
      nfail += check_valid( "LutMap", (AstObject *) m, "LutMap" );
      if ( m ) m = astAnnul( m );
   }

/* WcsMap: the longitude and latitude axis indices must lie within the
   WcsMap, or the projection parameter code indexes its per-axis arrays out
   of bounds. */
   nfail += check_rejected( "WcsMap with out-of-range latitude axis",
      " Begin WcsMap\n Nin = 2\n IsA Mapping\n"
      " Type = \"SIN\"\n WcsAx2 = 99\n PV1_0 = 0.5\n End WcsMap\n" );
   {
      AstWcsMap *m = astWcsMap( 2, AST__SIN, 1, 2, " " );
      nfail += check_valid( "WcsMap", (AstObject *) m, "WcsMap" );
      if ( m ) m = astAnnul( m );
   }

/* SplineMap: the number of coefficients must be at least the spline order. */
   nfail += check_rejected( "SplineMap with too few X coefficients",
      " Begin SplineMap\n Nin = 2\n IsA Mapping\n"
      " KX = 3\n KY = 1\n NX = 2\n NY = 2\n"
      " TX1 = 0\n TX2 = 0.5\n TX3 = 1\n TX4 = 1\n TX5 = 1\n"
      " TY1 = 0\n TY2 = 0.5\n TY3 = 1\n"
      " CU1 = 0\n CU2 = 1\n CU3 = 2\n CU4 = 3\n"
      " CV1 = 0\n CV2 = 1\n CV3 = 2\n CV4 = 3\n End SplineMap\n" );
   nfail += check_rejected( "SplineMap with non-monotonic X knots",
      " Begin SplineMap\n Nin = 2\n IsA Mapping\n"
      " KX = 1\n KY = 1\n NX = 3\n NY = 3\n"
      " TX1 = 0\n TX2 = 0.9\n TX3 = 0.1\n TX4 = 1\n"
      " TY1 = 0\n TY2 = 0.4\n TY3 = 0.6\n TY4 = 1\n"
      " CU1=0\n CU2=1\n CU3=2\n CU4=3\n CU5=4\n CU6=5\n CU7=6\n CU8=7\n CU9=8\n"
      " CV1=0\n CV2=1\n CV3=2\n CV4=3\n CV5=4\n CV6=5\n CV7=6\n CV8=7\n CV9=8\n"
      " End SplineMap\n" );
   {
      double tx[ 3 ] = { 0.0, 0.5, 1.0 }, ty[ 3 ] = { 0.0, 0.5, 1.0 };
      double cu[ 4 ] = { 0.0, 1.0, 2.0, 3.0 }, cv[ 4 ] = { 0.0, 1.0, 2.0, 3.0 };
      AstSplineMap *m = astSplineMap( 1, 1, 2, 2, tx, ty, cu, cv, " " );
      nfail += check_valid( "SplineMap", (AstObject *) m, "SplineMap" );
      if ( m ) m = astAnnul( m );
   }

/* Plot: the base Frame must not have more axes than the fixed-length
   per-axis arrays in the Plot structure (sized for a Plot3D). */
   nfail += check_rejected( "Plot with a four-axis base Frame",
      " Begin Plot\n Nframe = 2\n Currnt = 2\n Lnk2 = 1\n"
      " Frm1 =\n Begin Frame\n Naxes = 4\n Domain = \"GRAPHICS\"\n End Frame\n"
      " Frm2 =\n Begin Frame\n Naxes = 4\n Domain = \"GRID\"\n End Frame\n"
      " Map2 =\n Begin WinMap\n Nin = 4\n IsA Mapping\n End WinMap\n"
      " IsA FrameSet\n End Plot\n" );

/* Moc: a Moc that contains any ranges must carry a valid MaxOrder. */
   nfail += check_rejected( "Moc with ranges but no MaxOrder",
      " Begin Moc\n Adapt = 1\n"
      " Frm =\n Begin SkyFrame\n Naxes = 2\n System = \"ICRS\"\n IsA Frame\n"
      " End SkyFrame\n RegAxes = 2\n IsA Region\n"
      " NRange = 1\n I0 = 80\n K0 = 95\n End Moc\n" );
   {
      AstMoc *m = astMoc( "MaxOrder=4" );
      astAddCell( m, AST__OR, 2, 5 );
      nfail += check_valid( "Moc", (AstObject *) m, "Moc" );
      if ( m ) m = astAnnul( m );
   }

   if ( !astOK ) astClearStatus;

   if ( nfail == 0 ) {
      printf( " All loader hardening tests passed\n" );
      return 0;
   }
   printf( "Loader hardening tests failed (%d)\n", nfail );
   return 1;
}
