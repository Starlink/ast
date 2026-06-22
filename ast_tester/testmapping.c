/*
 *  Test Mapping class (PolyMap / LinearApprox).
 *  Converted from the Fortran test testmapping.f.
 *  Direct conversion; no material differences from the Fortran original.
 */
#include "ast.h"
#include <stdio.h>
#include <string.h>

static void stopit( int *status, const char *text ) {
   if( *status != 0 ) return;
   *status = 1;
   printf( "%s\n", text );
}

int main( void ) {
   int status_value = 0;
   int *status = &status_value;
   AstPolyMap *pm;
   double coeff[20], fit[6], lbnd[2], ubnd[2];

   /* Coefficients for a 2-in 2-out polynomial:
      out1 = 1.0 + 2.0*x
      out2 = 1.0 + 3.0*y + 3.0*y^2  */
   double c[] = { 1.0, 1, 0, 0,
                  2.0, 1, 1, 0,
                  1.0, 2, 0, 0,
                  3.0, 2, 0, 1,
                  3.0, 1, 0, 2 };

   memcpy( coeff, c, sizeof( coeff ) );

   astWatch( status );
   astBegin;

   pm = astPolyMap( 2, 2, 4, coeff, 0, coeff, " " );

   lbnd[0] = -1.0;
   lbnd[1] = -1.0;
   ubnd[0] = 1.0;
   ubnd[1] = 1.0;
   if( astLinearApprox( pm, lbnd, ubnd, 0.001, fit ) ) {
      if( fit[0] != 1.0 || fit[1] != 1.0 ||
          fit[2] != 2.0 || fit[3] != 0.0 ||
          fit[4] != 0.0 || fit[5] != 3.0 ) {
         stopit( status, "Error 0" );
      }
   } else {
      stopit( status, "Error 1" );
   }

   coeff[12] = AST__BAD;
   pm = astPolyMap( 2, 2, 4, coeff, 0, coeff, " " );

   if( astLinearApprox( pm, lbnd, ubnd, 0.001, fit ) ) {
      if( fit[0] != 1.0 || fit[1] != AST__BAD ||
          fit[2] != 2.0 || fit[3] != 0.0 ||
          fit[4] != AST__BAD || fit[5] != AST__BAD ) {
         stopit( status, "Error 2" );
      }
   } else {
      stopit( status, "Error 3" );
   }

   pm = astPolyMap( 2, 2, 5, coeff, 0, coeff, " " );
   if( astLinearApprox( pm, lbnd, ubnd, 0.001, fit ) ) {
      stopit( status, "Error 4" );
   }

   /* PermMap equality must take the referenced constant *values* into
      account, not just the constant *indices*. Two PermMaps that feed an
      axis from the same constant slot but store different values there
      must compare unequal. */
   {
      int perm1[1] = { -1 };       /* feed axis from constant index 0 */
      int perm2[1] = { -2 };       /* feed axis from constant index 1 */
      double con5[1] = { 5.0 };
      double con7[1] = { 7.0 };
      double con95[2] = { 9.0, 5.0 };  /* index 1 holds 5.0 */
      double con57[2] = { 5.0, 7.0 };
      double con59[2] = { 5.0, 9.0 };
      AstPermMap *a, *b;

      /* Same index (-1), differing constant values -> unequal. */
      a = astPermMap( 1, perm1, 1, perm1, con5, " " );
      b = astPermMap( 1, perm1, 1, perm1, con7, " " );
      if( astEqual( a, b ) ) stopit( status, "Error 5" );
      a = astAnnul( a );
      b = astAnnul( b );

      /* Same index (-1), identical constant values -> equal. */
      a = astPermMap( 1, perm1, 1, perm1, con5, " " );
      b = astPermMap( 1, perm1, 1, perm1, con5, " " );
      if( !astEqual( a, b ) ) stopit( status, "Error 6" );
      a = astAnnul( a );
      b = astAnnul( b );

      /* Different indices (-1 vs -2) referencing equal values (both 5.0)
         -> equal. */
      a = astPermMap( 1, perm1, 1, perm1, con5,  " " );
      b = astPermMap( 1, perm2, 1, perm2, con95, " " );
      if( !astEqual( a, b ) ) stopit( status, "Error 7" );
      a = astAnnul( a );
      b = astAnnul( b );

      /* Differing unreferenced surplus constants -> still equal (only the
         referenced slot, index 0, matters here). */
      a = astPermMap( 1, perm1, 1, perm1, con57, " " );
      b = astPermMap( 1, perm1, 1, perm1, con59, " " );
      if( !astEqual( a, b ) ) stopit( status, "Error 8" );
      a = astAnnul( a );
      b = astAnnul( b );
   }

   astEnd;

   if( *status == 0 ) {
      printf( " All Mapping tests passed\n" );
   } else {
      printf( "Mapping tests failed\n" );
   }
   return *status;
}
