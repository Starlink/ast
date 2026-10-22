/*
 * Test astConvert/astFindFrame public behaviour.
 */
#include "ast.h"
#include "ast_err.h"
#include <math.h>
#include <stdio.h>

int main( void ){
   int status_value = 0;
   int *status = &status_value;

   AstFrameSet *fs;
   AstSkyFrame *sf;
   AstSpecFrame *df;
   AstCmpFrame *cf;
   AstFrame *bf;
   AstFrame *target, *template;

   astWatch( status );
   astBegin;

   sf = astSkyFrame( " " );
   df = astSpecFrame( " " );
   cf = astCmpFrame( df, sf, " " );
   bf = astFrame( 2, "Domain=SKY" );

   fs = astConvert( bf, sf, " " );
   if( fs ) {
      if( !astEqual( astGetFrame( fs, AST__BASE ), bf ) && astOK ) {
         astError( AST__INTER, "Error 1\n" );
      } else if( !astEqual( astGetFrame( fs, AST__CURRENT ), sf ) && astOK ) {
         astError( AST__INTER, "Error 2\n" );
      } else if( !astIsAUnitMap( astGetMapping( fs, AST__BASE, AST__CURRENT ) ) ) {
         astError( AST__INTER, "Error 3\n" );
      }
   } else {
      astError( AST__INTER, "Error 4\n" );
   }

   fs = astConvert( sf, bf, " " );
   if( fs ) {
      if( !astEqual( astGetFrame( fs, AST__BASE ), sf ) && astOK ) {
         astError( AST__INTER, "Error 5\n" );
      } else if( !astEqual( astGetFrame( fs, AST__CURRENT ), bf ) && astOK ) {
         astError( AST__INTER, "Error 6\n" );
      } else if( !astIsAUnitMap( astGetMapping( fs, AST__BASE, AST__CURRENT ) ) ) {
         astError( AST__INTER, "Error 7\n" );
      }
   } else {
      astError( AST__INTER, "Error 8\n" );
   }


   astSetC( bf, "Domain", "NOTSKY" );
   fs = astConvert( bf, sf, " " );
   if( fs ) {
      astShow( fs );
      astError( AST__INTER, "Error 9\n" );
   }

   fs = astConvert( sf, bf, " " );
   if( fs ) {
      astShow( fs );
      astError( AST__INTER, "Error 10\n" );
   }

   astClear( bf, "Domain" );

   fs = astConvert( bf, sf, " " );
   if( fs ) {
      if( !astEqual( astGetFrame( fs, AST__BASE ), bf ) && astOK ) {
         astError( AST__INTER, "Error 11\n" );
      } else if( !astEqual( astGetFrame( fs, AST__CURRENT ), sf ) && astOK ) {
         astError( AST__INTER, "Error 12\n" );
      } else if( !astIsAUnitMap( astGetMapping( fs, AST__BASE, AST__CURRENT ) ) ) {
         astError( AST__INTER, "Error 13\n" );
      }
   } else {
      astError( AST__INTER, "Error 14\n" );
   }

   fs = astConvert( sf, bf, " " );
   if( fs ) {
      if( !astEqual( astGetFrame( fs, AST__BASE ), sf ) && astOK ) {
         astError( AST__INTER, "Error 15\n" );
      } else if( !astEqual( astGetFrame( fs, AST__CURRENT ), bf ) && astOK ) {
         astError( AST__INTER, "Error 16\n" );
      } else if( !astIsAUnitMap( astGetMapping( fs, AST__BASE, AST__CURRENT ) ) ) {
         astError( AST__INTER, "Error 17\n" );
      }
   } else {
      astError( AST__INTER, "Error 18\n" );
   }


   fs = astConvert( bf, cf, " " );
   if( fs ) {
      if( !astEqual( astGetFrame( fs, AST__BASE ), bf ) && astOK ) {
         astError( AST__INTER, "Error 19\n" );
      } else if( !astEqual( astGetFrame( fs, AST__CURRENT ), cf ) && astOK ) {
         astError( AST__INTER, "Error 20\n" );
      } else if( !astIsAPermMap( astGetMapping( fs, AST__BASE, AST__CURRENT ) ) ) {
         astError( AST__INTER, "Error 21\n" );
      }
   } else {
      astError( AST__INTER, "Error 22\n" );
   }

   fs = astConvert( cf, bf, " " );
   if( fs ) {
      if( !astEqual( astGetFrame( fs, AST__BASE ), cf ) && astOK ) {
         astError( AST__INTER, "Error 23\n" );
      } else if( !astEqual( astGetFrame( fs, AST__CURRENT ), bf ) && astOK ) {
         astError( AST__INTER, "Error 24\n" );
      } else if( !astIsAPermMap( astGetMapping( fs, AST__BASE, AST__CURRENT ) ) ) {
         astError( AST__INTER, "Error 25\n" );
      }
   } else {
      astError( AST__INTER, "Error 26\n" );
   }


   astSetC( bf, "Domain", "NOTSKY" );
   fs = astConvert( bf, cf, " " );
   if( fs ) {
      astShow( fs );
      astError( AST__INTER, "Error 27\n" );
   }

   fs = astConvert( cf, bf, " " );
   if( fs ) {
      astShow( fs );
      astError( AST__INTER, "Error 28\n" );
   }


   astSetC( bf, "Domain", "SKY" );
   fs = astConvert( bf, cf, " " );
   if( fs ) {
      if( !astEqual( astGetFrame( fs, AST__BASE ), bf ) && astOK ) {
         astError( AST__INTER, "Error 29\n" );
      } else if( !astEqual( astGetFrame( fs, AST__CURRENT ), cf ) && astOK ) {
         astError( AST__INTER, "Error 30\n" );
      } else if( !astIsAPermMap( astGetMapping( fs, AST__BASE, AST__CURRENT ) ) ) {
         astError( AST__INTER, "Error 31\n" );
      }
   } else {
      astError( AST__INTER, "Error 32\n" );
   }

   fs = astConvert( cf, bf, " " );
   if( fs ) {
      if( !astEqual( astGetFrame( fs, AST__BASE ), cf ) && astOK ) {
         astError( AST__INTER, "Error 33\n" );
      } else if( !astEqual( astGetFrame( fs, AST__CURRENT ), bf ) && astOK ) {
         astError( AST__INTER, "Error 34\n" );
      } else if( !astIsAPermMap( astGetMapping( fs, AST__BASE, AST__CURRENT ) ) ) {
         astError( AST__INTER, "Error 35\n" );
      }
   } else {
      astError( AST__INTER, "Error 36\n" );
   }


   fs = astConvert( sf, cf, " " );
   if( fs ) {
      if( !astEqual( astGetFrame( fs, AST__BASE ), sf ) && astOK ) {
         astError( AST__INTER, "Error 37\n" );
      } else if( !astEqual( astGetFrame( fs, AST__CURRENT ), cf ) && astOK ) {
         astError( AST__INTER, "Error 38\n" );
      } else if( !astIsAPermMap( astGetMapping( fs, AST__BASE, AST__CURRENT ) ) ) {
         astError( AST__INTER, "Error 39\n" );
      }
   } else {
      astError( AST__INTER, "Error 40\n" );
   }

   fs = astConvert( cf, sf, " " );
   if( fs ) {
      if( !astEqual( astGetFrame( fs, AST__BASE ), cf ) && astOK ) {
         astError( AST__INTER, "Error 41\n" );
      } else if( !astEqual( astGetFrame( fs, AST__CURRENT ), sf ) && astOK ) {
         astError( AST__INTER, "Error 42\n" );
      } else if( !astIsAPermMap( astGetMapping( fs, AST__BASE, AST__CURRENT ) ) ) {
         astError( AST__INTER, "Error 43\n" );
      }
   } else {
      astError( AST__INTER, "Error 44\n" );
   }


   fs = astFindFrame( sf, cf, " " );
   if( !fs && astOK ) {
      astError( AST__INTER, "Error 45\n" );
   }

   fs = astFindFrame( cf, sf, " " );
   if( fs && astOK ) {
      astError( AST__INTER, "Error 46\n" );
   }

   astSetI( sf, "MaxAxes", 3 );
   astSetI( sf, "MinAxes", 1 );

   fs = astFindFrame( cf, sf, " " );
   if( !fs && astOK ) {
      astError( AST__INTER, "Error 47\n" );
   } else {
      if( !astEqual( astGetFrame( fs, AST__BASE ), cf ) && astOK ) {
         astError( AST__INTER, "Error 48\n" );
      } else if( !astEqual( astGetFrame( fs, AST__CURRENT ), sf ) && astOK ) {
         astError( AST__INTER, "Error 49\n" );
      } else if( !astIsAPermMap( astGetMapping( fs, AST__BASE, AST__CURRENT ) ) ) {
         astError( AST__INTER, "Error 50\n" );
      }
   }

   target = astFrame( 2, "Domain=ARDAPP" );
   template = (AstFrame *) astSkyFrame( "System=GAPPT" );
   fs = astFindFrame( target, template, " " );
   if( fs && astOK ) {
      astError( AST__INTER, "Error 51\n" );
   }

/* A constant in a units string is a coefficient wherever it appears:
   "2/m" means "each unit is 2 per metre", the same as "2*m**-1", so a
   value of 1 per metre is 0.5 in either. */
   {
      const char *units[ 4 ] = { "2/m", "2*m**-1", "0.5/m", "2/s" };
      const char *from[ 4 ] = { "1/m", "1/m", "1/m", "Hz" };
      double expect[ 4 ] = { 0.5, 0.5, 2.0, 0.5 };
      double in, out;
      int i;

      for( i = 0; i < 4 && astOK; i++ ) {
         AstFrame *f1 = astFrame( 1, "Unit(1)=%s", from[ i ] );
         AstFrame *f2 = astFrame( 1, "Unit(1)=%s", units[ i ] );
         astSetActiveUnit( f1, 1 );
         astSetActiveUnit( f2, 1 );
         fs = astConvert( f1, f2, " " );
         if( !fs && astOK ) {
            astError( AST__INTER, "Error 52: no conversion from %s to %s\n",
                      from[ i ], units[ i ] );
         } else if( astOK ) {
            in = 1.0;
            astTran1( fs, 1, &in, 1, &out );
            if( fabs( out - expect[ i ] ) > 1.0E-12 ) {
               astError( AST__INTER, "Error 53: 1 %s is %g %s, expected "
                         "%g\n", from[ i ], out, units[ i ], expect[ i ] );
            }
         }
      }
   }




   astEnd;

   if( astOK ) {
      printf(" All astConvert tests passed\n");
   } else {
      printf("astConvert tests failed\n");
   }

   return *status;
}
