/*
 * testgrismmap: GrismMap tests (C only; no Fortran original).
 *
 * A GrismMap implements the dispersion function of FITS-WCS Paper III
 * (Greisen et al. 2006, A&A 446, 747), section 5.1: its forward
 * transformation takes the grism parameter g to wavelength, and its inverse
 * takes wavelength to g.  The tests check:
 *
 *  - attribute defaults, and setting, testing and clearing attributes;
 *  - that each wavelength produced satisfies the grating equation in the
 *    implicit form of Paper III eqs. (68) and (69), which is independent of
 *    the rearranged form the GrismMap evaluates;
 *  - that g is -tan(theta) at the reference wavelength (eq. (76)), and that
 *    wavelength -> g -> wavelength recovers the wavelength;
 *  - that a bad input, an undispersed wavelength, an undispersed reference
 *    wavelength and a GrismMap whose parameters give no dispersion all
 *    produce AST__BAD;
 *  - that astRate in each direction agrees with a numerical derivative and
 *    with the reciprocal of the other direction, and at the reference
 *    wavelength with Paper III eq. (81);
 *  - that astSimplify cancels a GrismMap with its own inverse, does not
 *    cancel GrismMaps that differ (including only in their order), and
 *    merges a ZoomMap that rescales the wavelength; in every case the
 *    simplified Mapping must transform as the original does;
 *  - astEqual, astCopy and native serialization with astToString and
 *    astFromString.
 *
 * Three parameter sets are used: the KPNO MARS grism of Paper III Fig. 5
 * with a non-zero epsilon and theta, the KPNO Hydra echelle of Fig. 4 in
 * order 11, and a grating used in order -1.  Lengths are in metres and
 * angles in radians.
 */
#include "ast.h"
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>

#define DEG2RAD ( 3.14159265358979323846/180.0 )

static int status = 0;
static int failed = 0;

static void fail( const char *fmt, ... )
   __attribute__((format(printf, 1, 2)));

/* Report a failed check. */
static void fail( const char *fmt, ... ) {
   va_list ap;
   printf( "Error: " );
   va_start( ap, fmt );
   vprintf( fmt, ap );
   va_end( ap );
   printf( "\n" );
   failed = 1;
}

/* The parameters of one disperser. */
typedef struct {
   const char *name;
   double g;          /* Ruling density (m^-1) */
   int m;             /* Order */
   double alpha;      /* Angle of incidence (rad) */
   double nr;         /* Refractive index at the reference wavelength */
   double nrp;        /* dn/dwave at the reference wavelength (m^-1) */
   double eps;        /* Grating tilt (rad) */
   double theta;      /* Detector tilt (rad) */
   double waver;      /* Reference wavelength (m) */
   double dgam;       /* Half-width of the sampled range of diffracted
                         angle about the reference ray (rad) */
} Disperser;

static const Disperser dispersers[] = {
   { "MARS grism", 4.5E5, 1, 27.0*DEG2RAD, 1.765, -1.077E6, 3.0*DEG2RAD,
     5.0*DEG2RAD, 7.2452E-7, 0.3 },
   { "Hydra echelle", 3.16E5, 11, 64.8*DEG2RAD, 1.0, 0.0, 0.0, 0.0,
     5.1368E-7, 0.3 },
   { "order -1 grating", 3.0E5, -1, 0.1, 1.0, 0.0, 0.0, 0.02, 5.0E-7, 0.1 },
};
#define NDISP ( (int) ( sizeof( dispersers )/sizeof( dispersers[ 0 ] ) ) )

/* Create a GrismMap for a disperser. */
static AstGrismMap *make_grismmap( const Disperser *d ) {
   AstGrismMap *map = astGrismMap( " " );
   astSetD( map, "GrismG", d->g );
   astSetI( map, "GrismM", d->m );
   astSetD( map, "GrismAlpha", d->alpha );
   astSetD( map, "GrismNR", d->nr );
   astSetD( map, "GrismNRP", d->nrp );
   astSetD( map, "GrismEps", d->eps );
   astSetD( map, "GrismTheta", d->theta );
   astSetD( map, "GrismWaveR", d->waver );
   return map;
}

/* Grism parameter values sampled for a disperser: those for which the
   diffracted angle lies within d->dgam of the reference ray.  The range is
   chosen for each disperser so that every wavelength is positive. */
#define NG 13
static void sample_g( const Disperser *d, double g[ NG ] ) {
   int i;
   for( i = 0; i < NG; i++ ) {
      g[ i ] = tan( d->dgam*( 2.0*i/( NG - 1.0 ) - 1.0 ) - d->theta );
   }
}

/* Paper III eq. (71): the diffracted angle of the reference ray. */
static double gamma_r( const Disperser *d ) {
   return asin( d->g*d->m*d->waver/cos( d->eps ) - d->nr*sin( d->alpha ) );
}

/* Paper III eq. (81): dg/dwave at the reference wavelength. */
static double eq81( const Disperser *d ) {
   return ( d->g*d->m/cos( d->eps ) - d->nrp*sin( d->alpha ) )/
          ( cos( gamma_r( d ) )*cos( d->theta )*cos( d->theta ) );
}

/* Richardson-extrapolated central difference of a 1-D Mapping. */
static double numerical_rate( AstMapping *map, double x, double h ) {
   double xs[ 4 ], y[ 4 ];
   xs[ 0 ] = x - h;
   xs[ 1 ] = x + h;
   xs[ 2 ] = x - 0.5*h;
   xs[ 3 ] = x + 0.5*h;
   astTran1( map, 4, xs, 1, y );
   return ( 4.0*( y[ 3 ] - y[ 2 ] )/h - ( y[ 1 ] - y[ 0 ] )/( 2.0*h ) )/3.0;
}

static int close_rel( double a, double b, double tol ) {
   return a != AST__BAD && b != AST__BAD && fabs( a - b ) <= tol*fabs( b );
}

static void check_attributes( void ) {
   AstGrismMap *map;
   static const char *names[] = { "GrismNR", "GrismNRP", "GrismWaveR",
      "GrismAlpha", "GrismG", "GrismM", "GrismEps", "GrismTheta" };
   static const double defaults[] = { 1.0, 0.0, 5000.0, 0.0, 0.0, 0.0, 0.0,
      0.0 };
   int i;

   astBegin;
   map = astGrismMap( " " );
   for( i = 0; i < 8; i++ ) {
      if( astTest( map, names[ i ] ) ) {
         fail( "%s is set on a new GrismMap", names[ i ] );
      }
      if( astGetD( map, names[ i ] ) != defaults[ i ] ) {
         fail( "%s defaults to %g, expected %g", names[ i ],
               astGetD( map, names[ i ] ), defaults[ i ] );
      }
   }

   astSet( map, "GrismG=3.0E5,GrismM=-2,GrismTheta=0.05" );
   if( !astTest( map, "GrismG" ) || astGetD( map, "GrismG" ) != 3.0E5 ) {
      fail( "GrismG not set by astSet" );
   }
   if( astGetI( map, "GrismM" ) != -2 ) {
      fail( "GrismM is %d, expected -2", astGetI( map, "GrismM" ) );
   }
   astClear( map, "GrismTheta" );
   if( astTest( map, "GrismTheta" ) || astGetD( map, "GrismTheta" ) != 0.0 ) {
      fail( "GrismTheta not cleared" );
   }
   if( astGetI( map, "Nin" ) != 1 || astGetI( map, "Nout" ) != 1 ) {
      fail( "GrismMap is not 1-D" );
   }
   astEnd;
}

static void check_transform( const Disperser *d ) {
   AstGrismMap *map;
   double g[ NG ], wave[ NG ], gback[ NG ], waveback[ NG ];
   double gr, gref, n, gam, lhs, rhs;
   int i;

   astBegin;
   map = make_grismmap( d );
   sample_g( d, g );
   astTran1( map, NG, g, 1, wave );
   astTran1( map, NG, wave, 0, gback );
   astTran1( map, NG, gback, 1, waveback );

/* Each output must satisfy the grating equation (eq. (68)) with the
   refractive index of eq. (69) and the diffracted angle of eq. (83). */
   gr = gamma_r( d );
   for( i = 0; i < NG; i++ ) {
      if( wave[ i ] == AST__BAD || wave[ i ] <= 0.0 ) {
         fail( "%s: g = %g gives wavelength %g", d->name, g[ i ], wave[ i ] );
         continue;
      }
      gam = atan( g[ i ] ) + gr + d->theta;
      n = d->nr + d->nrp*( wave[ i ] - d->waver );
      lhs = d->g*d->m*wave[ i ]/cos( d->eps );
      rhs = n*sin( d->alpha ) + sin( gam );
      if( fabs( lhs - rhs ) > 1.0E-14 ) {
         fail( "%s: g = %g, wavelength %.17g violates the grating equation"
               " by %g", d->name, g[ i ], wave[ i ], lhs - rhs );
      }
      if( !close_rel( gback[ i ], g[ i ], 1.0E-12 ) &&
          fabs( gback[ i ] - g[ i ] ) > 1.0E-14 ) {
         fail( "%s: g = %.17g round-trips to %.17g", d->name, g[ i ],
               gback[ i ] );
      }
      if( !close_rel( waveback[ i ], wave[ i ], 1.0E-14 ) ) {
         fail( "%s: wavelength %.17g round-trips to %.17g", d->name,
               wave[ i ], waveback[ i ] );
      }
   }

/* Eq. (76): the grism parameter at the reference wavelength. */
   astTran1( map, 1, &d->waver, 0, &gref );
   if( fabs( gref + tan( d->theta ) ) > 1.0E-14 ) {
      fail( "%s: g at the reference wavelength is %.17g, expected %.17g",
            d->name, gref, -tan( d->theta ) );
   }
   astEnd;
}

static void check_bad( void ) {
   const Disperser *d = &dispersers[ 0 ];
   AstGrismMap *map;
   double in, out;

   astBegin;
   map = make_grismmap( d );

   in = AST__BAD;
   astTran1( map, 1, &in, 1, &out );
   if( out != AST__BAD ) fail( "bad grism parameter gives %g", out );
   astTran1( map, 1, &in, 0, &out );
   if( out != AST__BAD ) fail( "bad wavelength gives %g", out );
   if( astRate( map, &in, 1, 1 ) != AST__BAD ) {
      fail( "rate at a bad grism parameter is not bad" );
   }

/* 10 microns is not diffracted by 450 lines/mm in first order at this
   angle of incidence: the sine of the diffracted angle exceeds 1. */
   in = 1.0E-5;
   astTran1( map, 1, &in, 0, &out );
   if( out != AST__BAD ) fail( "undispersed wavelength gives g = %g", out );
   astInvert( map );
   if( astRate( map, &in, 1, 1 ) != AST__BAD ) {
      fail( "rate at an undispersed wavelength is not bad" );
   }

/* Nor is a reference wavelength of 10 microns, and without a reference ray
   no wavelength has a position, so both directions give AST__BAD. */
   map = make_grismmap( d );
   astSetD( map, "GrismWaveR", 1.0E-5 );
   in = 0.1;
   astTran1( map, 1, &in, 1, &out );
   if( out != AST__BAD ) {
      fail( "undispersed reference wavelength gives wavelength %g", out );
   }
   astTran1( map, 1, &d->waver, 0, &out );
   if( out != AST__BAD ) {
      fail( "undispersed reference wavelength gives g = %g", out );
   }

/* With the default GrismG and GrismM there is no dispersion, so no
   wavelength can be found. */
   map = astGrismMap( " " );
   in = 0.1;
   astTran1( map, 1, &in, 1, &out );
   if( out != AST__BAD ) fail( "default GrismMap gives wavelength %g", out );
   if( astRate( map, &in, 1, 1 ) != AST__BAD ) {
      fail( "rate of the default GrismMap is not bad" );
   }

   if( !astOK ) {
      fail( "AST error while checking bad values" );
      astClearStatus;
   }
   astEnd;
}

static void check_rate( const Disperser *d ) {
   AstGrismMap *map;
   double g[ NG ], wave[ NG ];
   double dwdg, dgdw, num, h;
   int i;

   astBegin;
   map = make_grismmap( d );
   sample_g( d, g );
   astTran1( map, NG, g, 1, wave );

   for( i = 0; i < NG; i++ ) {
      dwdg = astRate( map, &g[ i ], 1, 1 );
      num = numerical_rate( (AstMapping *) map, g[ i ], 1.0E-4 );
      if( !close_rel( dwdg, num, 1.0E-9 ) ) {
         fail( "%s: dwave/dg at g = %g is %.17g, numerically %.17g",
               d->name, g[ i ], dwdg, num );
      }

      astInvert( map );
      dgdw = astRate( map, &wave[ i ], 1, 1 );
      h = 1.0E-4/dgdw;
      num = numerical_rate( (AstMapping *) map, wave[ i ], h );
      if( !close_rel( dgdw, num, 1.0E-9 ) ) {
         fail( "%s: dg/dwave at %g m is %.17g, numerically %.17g",
               d->name, wave[ i ], dgdw, num );
      }
      astInvert( map );

      if( !close_rel( dwdg*dgdw, 1.0, 1.0E-14 ) ) {
         fail( "%s: rates at g = %g are not reciprocal (product %.17g)",
               d->name, g[ i ], dwdg*dgdw );
      }
   }

/* Eq. (81) at the reference wavelength. */
   double waver = d->waver;
   astInvert( map );
   dgdw = astRate( map, &waver, 1, 1 );
   if( !close_rel( dgdw, eq81( d ), 1.0E-14 ) ) {
      fail( "%s: dg/dwave at the reference wavelength is %.17g, eq. (81)"
            " gives %.17g", d->name, dgdw, eq81( d ) );
   }
   astEnd;
}

/* Check that "simp" transforms as "map" does at the supplied inputs. */
static void check_same_transform( AstMapping *map, AstMapping *simp,
                                  int n, const double *in,
                                  const char *what ) {
   double out1[ NG ], out2[ NG ];
   int i;

   astTran1( map, n, in, 1, out1 );
   astTran1( simp, n, in, 1, out2 );
   for( i = 0; i < n; i++ ) {
      if( out1[ i ] == AST__BAD && out2[ i ] == AST__BAD ) continue;
      if( !close_rel( out2[ i ], out1[ i ], 1.0E-12 ) &&
          fabs( out2[ i ] - out1[ i ] ) > 1.0E-14 ) {
         fail( "%s: simplified Mapping gives %.17g at %g, original %.17g",
               what, out2[ i ], in[ i ], out1[ i ] );
         return;
      }
   }
}

static void check_simplify( void ) {
   const Disperser *d = &dispersers[ 0 ];
   AstGrismMap *g1, *g2;
   AstMapping *cm, *simp;
   AstZoomMap *zm;
   double g[ NG ], wave[ NG ];
   char text[ 100 ];

   astBegin;
   sample_g( d, g );
   g1 = make_grismmap( d );
   astTran1( g1, NG, g, 1, wave );

/* A GrismMap followed by its own inverse is a UnitMap. */
   g2 = make_grismmap( d );
   astInvert( g2 );
   cm = (AstMapping *) astCmpMap( g1, g2, 1, " " );
   simp = astSimplify( cm );
   if( !astIsAUnitMap( simp ) ) {
      fail( "GrismMap followed by its inverse simplifies to a %s, not a"
            " UnitMap", astGetC( simp, "Class" ) );
   }
   check_same_transform( cm, simp, NG, g, "GrismMap then its inverse" );

/* The same with the inverse first. */
   cm = (AstMapping *) astCmpMap( g2, g1, 1, " " );
   simp = astSimplify( cm );
   if( !astIsAUnitMap( simp ) ) {
      fail( "inverse GrismMap followed by the GrismMap simplifies to a %s,"
            " not a UnitMap", astGetC( simp, "Class" ) );
   }
   check_same_transform( cm, simp, NG, wave, "inverse GrismMap then GrismMap" );

/* GrismMaps that differ only in their order, or only in another
   attribute, are not inverses of each other, so must not cancel. */
   static const char *changes[] = { "GrismM=2", "GrismM=-1", "GrismG=4.6E5",
      "GrismTheta=0.1" };
   int i;
   for( i = 0; i < 4; i++ ) {
      g2 = make_grismmap( d );
      astSet( g2, "%s", changes[ i ] );
      astInvert( g2 );
      cm = (AstMapping *) astCmpMap( g1, g2, 1, " " );
      simp = astSimplify( cm );
      snprintf( text, sizeof( text ), "GrismMap then inverse with %s",
                changes[ i ] );
      if( astIsAUnitMap( simp ) ) {
         fail( "%s simplifies to a UnitMap", text );
      }
      check_same_transform( cm, simp, NG, g, text );
   }

/* A ZoomMap rescaling the wavelength after a forward GrismMap, or before
   an inverse one, is absorbed into the GrismMap's length attributes. */
   zm = astZoomMap( 1, 1.0E10, " " );
   cm = (AstMapping *) astCmpMap( g1, zm, 1, " " );
   simp = astSimplify( cm );
   if( !astIsAGrismMap( simp ) ) {
      fail( "GrismMap then ZoomMap simplifies to a %s, not a GrismMap",
            astGetC( simp, "Class" ) );
   }
   check_same_transform( cm, simp, NG, g, "GrismMap then ZoomMap" );

   g2 = astCopy( g1 );
   astInvert( g2 );
   astInvert( zm );
   cm = (AstMapping *) astCmpMap( zm, g2, 1, " " );
   simp = astSimplify( cm );
   if( !astIsAGrismMap( simp ) ) {
      fail( "ZoomMap then inverse GrismMap simplifies to a %s, not a"
            " GrismMap", astGetC( simp, "Class" ) );
   }
   {
      double ang[ NG ];
      for( i = 0; i < NG; i++ ) ang[ i ] = wave[ i ]*1.0E10;
      check_same_transform( cm, simp, NG, ang, "ZoomMap then inverse GrismMap" );
   }

   if( !astOK ) {
      fail( "AST error while checking simplification" );
      astClearStatus;
   }
   astEnd;
}

static void check_copy_equal_dump( void ) {
   const Disperser *d = &dispersers[ 1 ];
   AstGrismMap *map, *copy;
   AstObject *back;
   char *pickle;
   double g[ NG ];

   astBegin;
   map = make_grismmap( d );
   sample_g( d, g );

   copy = astCopy( map );
   if( !astEqual( map, copy ) ) fail( "copy of a GrismMap is not equal" );
   astSetD( copy, "GrismNR", 1.5 );
   if( astEqual( map, copy ) ) {
      fail( "GrismMaps with different GrismNR are equal" );
   }
   copy = astCopy( map );
   astInvert( copy );
   if( astEqual( map, copy ) ) {
      fail( "a GrismMap equals its inverse" );
   }

   pickle = astToString( map );
   back = pickle ? astFromString( pickle ) : NULL;
   if( !back ) {
      fail( "GrismMap does not survive astToString/astFromString" );
   } else {
      if( !astIsAGrismMap( back ) ) {
         fail( "astFromString returned a %s", astGetC( back, "Class" ) );
      } else if( !astEqual( map, back ) ) {
         fail( "GrismMap read back from its dump is not equal" );
      } else if( astGetI( back, "GrismM" ) != d->m ||
                 astGetD( back, "GrismWaveR" ) != d->waver ) {
         fail( "GrismMap attributes not preserved by its dump" );
      }
      check_same_transform( (AstMapping *) map, (AstMapping *) back, NG, g,
                            "GrismMap read back from its dump" );
   }
   if( pickle ) pickle = astFree( pickle );

   if( !astOK ) {
      fail( "AST error while checking copy, equality and dump" );
      astClearStatus;
   }
   astEnd;
}

int main( void ) {
   int i;

   astWatch( &status );
   astBegin;

   check_attributes();
   for( i = 0; i < NDISP; i++ ) {
      check_transform( &dispersers[ i ] );
      check_rate( &dispersers[ i ] );
   }
   check_bad();
   check_simplify();
   check_copy_equal_dump();

   astEnd;

   if( !failed && astOK ) {
      printf( " All GrismMap tests passed\n" );
      return 0;
   }
   printf( "GrismMap tests failed\n" );
   return 1;
}
