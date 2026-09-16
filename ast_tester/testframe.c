#include <math.h>
#include <string.h>
#include "sae_par.h"
#include "ast.h"
#include "ast_err.h"

int main( void ) {
   int status = SAI__OK;
   astWatch( &status );

   AstFrame* frame = astFrame( 1, " " );
   astSetC( frame, "Unit(1)", "s*(m/s)" );
   const char* result = astGetC( frame, "NormUnit(1)" );

   if( strcmp( result, "m" ) ) {
      astError( AST__INTER, "NormUnit did not give expected result" );
   }

/* Not every Unit value is a units expression that can be parsed. A Frame
   axis has a blank Unit until one is set, and there is then nothing to
   normalise, so NormUnit is blank too rather than an error. The same
   applies to the Frame classes that report a Unit of their own but leave
   the Axis Unit unset. */
   {
      AstFrame *plain = astFrame( 2, " " );
      AstSpecFrame *spec = astSpecFrame( " " );
      AstTimeFrame *time = astTimeFrame( " " );
      const char *blank;
      int iaxis;

      for( iaxis = 1; iaxis <= 2; iaxis++ ) {
         blank = astGetC( plain, iaxis == 1 ? "NormUnit(1)" : "NormUnit(2)" );
         if( astOK && ( !blank || strlen( blank ) ) ) {
            astError( AST__INTER, "NormUnit(%d) of a default Frame is '%s', "
                      "expected a blank string", iaxis, blank ? blank : "<NULL>" );
         }
      }

      blank = astGetC( spec, "NormUnit(1)" );
      if( astOK && ( !blank || strlen( blank ) ) ) {
         astError( AST__INTER, "NormUnit(1) of a default SpecFrame is '%s', "
                   "expected a blank string", blank ? blank : "<NULL>" );
      }

      blank = astGetC( time, "NormUnit(1)" );
      if( astOK && ( !blank || strlen( blank ) ) ) {
         astError( AST__INTER, "NormUnit(1) of a default TimeFrame is '%s', "
                   "expected a blank string", blank ? blank : "<NULL>" );
      }

      plain = astAnnul( plain );
      spec = astAnnul( spec );
      time = astAnnul( time );
   }

/* A SkyFrame axis describes its values with a sexagesimal format such as
   "ddd:mm:ss", which is not a units expression either. Reading NormUnit
   must still succeed; the value follows the Axis Unit. */
   {
      AstSkyFrame *sky = astSkyFrame( " " );
      const char *norm = astGetC( sky, "NormUnit(1)" );

      if( astOK && ( !norm || !strlen( norm ) ) ) {
         astError( AST__INTER, "NormUnit(1) of a default SkyFrame is '%s', "
                   "expected a non-blank string", norm ? norm : "<NULL>" );
      }
      sky = astAnnul( sky );
   }

/* astAxAngle: when the offset position has a zero component on the
   measured axis but a non-zero component on the other axis, the angle is
   still well defined. Previously the nudge applied to break the
   degeneracy was scaled by a stale loop index, collapsing the offset
   point onto the reference point and returning AST__BAD. */
   {
      AstFrame *f2 = astFrame( 2, " " );
      double a1[2] = { 0.0, 0.0 }, b1[2] = { 1.0, 0.0 };
      double a2[2] = { 0.0, 5.0 }, b2[2] = { 3.0, 2.0 };
      double ang1 = astAxAngle( f2, a1, b1, 1 );
      double ang2 = astAxAngle( f2, a2, b2, 1 );

      if( astOK && ang1 == AST__BAD ) {
         astError( AST__INTER, "astAxAngle returned AST__BAD for "
                   "(0,0)->(1,0) on axis 1" );
      }
      if( astOK && fabs( ang2 - M_PI/4.0 ) > 1.0E-10 ) {
         astError( AST__INTER, "astAxAngle gave %g, expected PI/4 for "
                   "(0,5)->(3,2) on axis 1", ang2 );
      }
      f2 = astAnnul( f2 );
   }

   if( astOK ) {
      printf(" All Frame tests passed\n");
   } else {
      printf("Frame tests failed\n");
   }
   return astOK ? 0 : 1;
}
