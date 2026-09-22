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

/* NormUnit is documented to equal Unit when no simplification can be
   performed, and several Unit values cannot be simplified at all: a Frame
   axis has a blank Unit until one is set, and a SkyFrame axis describes its
   values with a sexagesimal format such as "ddd:mm:ss". Check the
   relationship rather than the individual strings, so that this pins the
   documented behaviour and not whatever the library currently returns.

   The classes below are the ones that report a Unit of their own by
   over-riding astGetUnit while leaving the Axis Unit unset, which is where
   NormUnit and Unit are easiest to get out of step. */
   {
      AstFrame *frames[ 4 ];
      int naxes[ 4 ] = { 2, 2, 1, 1 };
      const char *names[ 4 ] = { "Frame", "SkyFrame", "SpecFrame",
                                 "TimeFrame" };
      int iframe;
      int iaxis;

      frames[ 0 ] = astFrame( 2, " " );
      frames[ 1 ] = (AstFrame *) astSkyFrame( " " );
      frames[ 2 ] = (AstFrame *) astSpecFrame( " " );
      frames[ 3 ] = (AstFrame *) astTimeFrame( " " );

      for( iframe = 0; iframe < 4 && astOK; iframe++ ) {
         for( iaxis = 1; iaxis <= naxes[ iframe ] && astOK; iaxis++ ) {
            char attr[ 20 ];
            const char *unit;
            const char *norm;

            snprintf( attr, sizeof( attr ), "Unit(%d)", iaxis );
            unit = astGetC( frames[ iframe ], attr );
            if( !astOK ) break;

/* Copy it, since the next astGetC may re-use the buffer it points into. */
            {
               char saved[ 128 ];
               snprintf( saved, sizeof( saved ), "%s", unit ? unit : "" );

               snprintf( attr, sizeof( attr ), "NormUnit(%d)", iaxis );
               norm = astGetC( frames[ iframe ], attr );
               if( !astOK ) break;

               if( !norm || strcmp( norm, saved ) ) {
                  astError( AST__INTER, "NormUnit(%d) of a default %s is "
                            "'%s', expected '%s' since it cannot be "
                            "simplified", iaxis, names[ iframe ],
                            norm ? norm : "<NULL>", saved );
               }
            }
         }
      }

      for( iframe = 0; iframe < 4; iframe++ ) {
         frames[ iframe ] = astAnnul( frames[ iframe ] );
      }
   }

/* A Unit that can be simplified is, and the value follows the Unit the
   Frame reports even when the axes have been permuted. */
   {
      AstFrame *perm = astFrame( 2, "Unit(1)=km,Unit(2)=s*(m/s)" );
      int outperm[ 2 ] = { 2, 1 };
      const char *norm;

      astPermAxes( perm, outperm );

      norm = astGetC( perm, "NormUnit(1)" );
      if( astOK && ( !norm || strcmp( norm, "m" ) ) ) {
         astError( AST__INTER, "NormUnit(1) of a permuted Frame is '%s', "
                   "expected 'm'", norm ? norm : "<NULL>" );
      }
      norm = astGetC( perm, "NormUnit(2)" );
      if( astOK && ( !norm || strcmp( norm, "km" ) ) ) {
         astError( AST__INTER, "NormUnit(2) of a permuted Frame is '%s', "
                   "expected 'km'", norm ? norm : "<NULL>" );
      }
      perm = astAnnul( perm );
   }

/* A SkyFrame axis reports a description of its values rather than a unit
   for them whenever its Format calls for more than one sexagesimal field,
   or for a single field of a time. Those descriptions contain spaces,
   which the units parser reads as multiplication, so they have to be left
   alone. A Format that names a single field of an angle, and a Unit set
   explicitly, both give units expressions which are simplified as usual. */
   {
      const char *attrs[ 6 ] = { "Format(1)=bhms", "Format(1)=bdms.2",
                                 "Format(1)=btm", "Format(1)=bts",
                                 "Format(1)=bd", "Unit(1)=s*(m/s)" };
      const char *units[ 6 ] = { "hh mm ss", "ddd mm ss.ss",
                                 "minutes of time", "seconds of time",
                                 "degrees", "s*(m/s)" };
      const char *norms[ 6 ] = { "hh mm ss", "ddd mm ss.ss",
                                 "minutes of time", "seconds of time",
                                 "deg", "m" };
      int i;

      for( i = 0; i < 6 && astOK; i++ ) {
         AstSkyFrame *sf = astSkyFrame( attrs[ i ] );
         const char *unit = astGetC( sf, "Unit(1)" );
         char saved[ 128 ];

         if( astOK ) {
            const char *norm;

            snprintf( saved, sizeof( saved ), "%s", unit ? unit : "" );
            if( strcmp( saved, units[ i ] ) ) {
               astError( AST__INTER, "Unit(1) of a SkyFrame with %s is "
                         "'%s', expected '%s'", attrs[ i ], saved,
                         units[ i ] );
            }

            norm = astGetC( sf, "NormUnit(1)" );
            if( astOK && ( !norm || strcmp( norm, norms[ i ] ) ) ) {
               astError( AST__INTER, "NormUnit(1) of a SkyFrame with %s is "
                         "'%s', expected '%s'", attrs[ i ],
                         norm ? norm : "<NULL>", norms[ i ] );
            }
         }

         sf = astAnnul( sf );
      }
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
