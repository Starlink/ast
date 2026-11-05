/*
 *  wcs_conversion_fixtures.c -- generator for the WCS-conversion merge
 *  fixtures (SpecMap/SlaMap/TimeMap) listed in the "Deferred: WCS-conversion
 *  merge" section of simplify_coverage_gaps.md.
 *
 *  Each WCS-conversion class cancels adjacent conversions by type. Its
 *  MapMerge has two relevant code paths:
 *    - invert-normalisation (SWAP_CODES): swap every conversion code to its
 *      inverse when the map is applied inverted;
 *    - forward cancellation (PAIR_CVT): collapse an adjacent fwd/inv pair.
 *  For each class this program emits two "giant" fixtures that between them
 *  exercise both paths across every enumerated conversion-inverse pair:
 *    <class>_invert_normalize.map  one inverted map, a conversion of each type
 *    <class>_pair_cancel.map       one forward map, adjacent fwd/inv pairs
 *
 *  Build and run (paths relative to the repo root, against a CMake build/):
 *
 *    cc -I build ast_tester/coverage/wcs_conversion_fixtures.c \
 *       -o /tmp/wcsconv build/libast.dylib \
 *       -Wl,-force_load,build/libast_err.a \
 *       -Wl,-force_load,build/libast_grf_5.6.a \
 *       -Wl,-force_load,build/libast_grf3d.a \
 *       -Wl,-rpath,"$PWD/build" -lm
 *    mkdir -p /tmp/sweep && /tmp/wcsconv /tmp/sweep
 *
 *  then install the .map files under ast_tester/fixtures/simplify/, generate
 *  their .simp with build/ast_tester/simplify, and wire them into
 *  fixtures/simplify/cases.txt. Extend the per-class Pair[] tables to cover more
 *  conversion types.
 *
 *  NOTE: this closes the primary merge directions; the astOK/EQUAL
 *  sub-branches of the SWAP_CODES/PAIR_CVT macro expansions are a residual
 *  plateau that no reasonable fixture set removes.
 */

#include "ast.h"
#include <stdio.h>

static void dump( AstMapping *m, const char *path ) {
   AstChannel *ch = astChannel( NULL, NULL, "SinkFile=%s", path );
   astWrite( ch, m ); astAnnul( ch );
}
static double A[6] = { 2000.0, 2000.0, 2000.0, 2000.0, 2000.0, 2000.0 };

/* A conversion-inverse pair, its argument count, and the arguments for the
   forward and inverse conversions.  A NULL argument array means A. */
typedef struct {
   const char *fwd, *inv;
   int narg;
   const double *fargs, *iargs;
} Pair;

#define FARGS(P) ( (P).fargs ? (P).fargs : A )
#define IARGS(P) ( (P).iargs ? (P).iargs : A )

/* Build two fixtures for a class: an inverted map holding one conversion of
   each type (exercises the invert-normalisation SWAP_CODES), and a forward
   map holding adjacent fwd/inv pairs (exercises the PAIR_CVT cancellation). */
#define BUILD(CTOR, ADD, PAIRS, INVNAME, CANCNAME)                       \
   do {                                                                  \
      int n = (int)(sizeof(PAIRS)/sizeof(Pair));                         \
      char path[ 1024 ];                                                 \
      *status = 0;                                                       \
      { void *m = CTOR;                                                  \
        for( int i=0;i<n;i++ ) ADD(m, PAIRS[i].fwd, PAIRS[i].narg,       \
                                   FARGS(PAIRS[i]));                     \
        astInvert(m);                                                    \
        snprintf( path, sizeof path, "%s/%s", dir, INVNAME );            \
        if(astOK) dump((AstMapping*)m, path);                            \
        else fprintf(stderr, "%s invert build failed\n", path); }        \
      *status = 0;                                                       \
      { void *m = CTOR;                                                  \
        for( int i=0;i<n;i++ ){                                          \
           ADD(m, PAIRS[i].fwd, PAIRS[i].narg, FARGS(PAIRS[i]));         \
           ADD(m, PAIRS[i].inv, PAIRS[i].narg, IARGS(PAIRS[i])); }       \
        snprintf( path, sizeof path, "%s/%s", dir, CANCNAME );           \
        if(astOK) dump((AstMapping*)m, path);                            \
        else fprintf(stderr, "%s cancel build failed\n", path); }        \
   } while(0)

int main( int argc, char **argv ) {
   int sv = 0; int *status = &sv; astWatch( status );
   const char *dir = ( argc > 1 ) ? argv[ 1 ] : "/tmp/sweep";

   static Pair spec[] = {
      {"FRTOVL","VLTOFR",1,NULL,NULL},{"ENTOFR","FRTOEN",0,NULL,NULL},
      {"WNTOFR","FRTOWN",0,NULL,NULL},{"WVTOFR","FRTOWV",0,NULL,NULL},
      {"AWTOFR","FRTOAW",0,NULL,NULL},{"VRTOVL","VLTOVR",0,NULL,NULL},
      {"VOTOVL","VLTOVO",0,NULL,NULL},{"ZOTOVL","VLTOZO",0,NULL,NULL},
      {"BTTOVL","VLTOBT",0,NULL,NULL},{"TPF2HL","HLF2TP",6,NULL,NULL},
      {"USF2HL","HLF2US",3,NULL,NULL},{"GEF2HL","HLF2GE",3,NULL,NULL},
      {"BYF2HL","HLF2BY",3,NULL,NULL},{"LKF2HL","HLF2LK",2,NULL,NULL},
      {"LDF2HL","HLF2LD",2,NULL,NULL} };
   BUILD( astSpecMap(1,0,""), astSpecAdd, spec,
          "specmap_invert_normalize.map", "specmap_pair_cancel.map" );

/* SlaMap arguments that are physically meaningful, so that the cancelled
   pairs are inverses and the fixtures can be checked numerically: Besselian
   epoch 1950 for the FK4 conversions; the date MJD 58849 (2020 January 1)
   and the J2000 equinox for AMP and MAP, which take them in opposite order;
   an observer at latitude 0.35 rad with the matching diurnal aberration
   1.4e-6 for H2E and E2H; and a local apparent sidereal time of 1.2 rad. */
   static const double b1950[] = { 1950.0 };
   static const double j2000[] = { 2000.0 };
   static const double date[] = { 58849.0 };
   static const double ampargs[] = { 58849.0, 2000.0 };
   static const double mapargs[] = { 2000.0, 58849.0 };
   static const double obs[] = { 0.35, 1.4E-6 };
   static const double last[] = { 1.2 };
   static Pair sla[] = {
      {"FK54Z","FK45Z",1,b1950,b1950},{"FK5HZ","HFK5Z",1,j2000,j2000},
      {"AMP","MAP",2,ampargs,mapargs},{"ADDET","SUBET",1,b1950,b1950},
      {"ECLEQ","EQECL",1,date,date},{"GALEQ","EQGAL",0,NULL,NULL},
      {"GALSUP","SUPGAL",0,NULL,NULL},{"H2E","E2H",2,obs,obs},
      {"R2H","H2R",1,last,last},{"HEEQ","EQHE",1,date,date} };
   BUILD( astSlaMap(0,""), astSlaAdd, sla,
          "slamap_invert_normalize.map", "slamap_pair_cancel.map" );

   static Pair tim[] = {
      {"TTTOTDB","TDBTOTT",5,NULL,NULL},{"TDBTOTCB","TCBTOTDB",1,NULL,NULL},
      {"TTTOTCG","TCGTOTT",1,NULL,NULL},{"TAITOTT","TTTOTAI",1,NULL,NULL},
      {"UTTOGMST","GMSTTOUT",1,NULL,NULL},
      {"GMSTTOLMST","LMSTTOGMST",3,NULL,NULL},
      {"LASTTOLMST","LMSTTOLAST",3,NULL,NULL},
      {"UTTOUTC","UTCTOUT",1,NULL,NULL} };
   BUILD( astTimeMap(0,""), astTimeAdd, tim,
          "timemap_invert_normalize.map", "timemap_pair_cancel.map" );

   fprintf( stderr, "done\n" );
   return 0;
}
