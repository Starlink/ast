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
 *    mkdir -p /tmp/sweep && /tmp/wcsconv
 *
 *  then install the .map files under ast_tester/simplify_fixtures/, generate
 *  their .simp with build/ast_tester/simplify, and wire them into
 *  simplify_tests.txt. Extend the per-class Pair[] tables to cover more
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

/* A conversion-inverse pair and its argument count. */
typedef struct { const char *fwd, *inv; int narg; } Pair;

/* Build two fixtures for a class: an inverted map holding one conversion of
   each type (exercises the invert-normalisation SWAP_CODES), and a forward
   map holding adjacent fwd/inv pairs (exercises the PAIR_CVT cancellation). */
#define BUILD(CTOR, ADD, PAIRS, INVPATH, CANCPATH)                       \
   do {                                                                  \
      int n = (int)(sizeof(PAIRS)/sizeof(Pair));                         \
      *status = 0;                                                       \
      { void *m = CTOR;                                                  \
        for( int i=0;i<n;i++ ) ADD(m, PAIRS[i].fwd, PAIRS[i].narg, A);   \
        astInvert(m);                                                    \
        if(astOK) dump((AstMapping*)m, INVPATH);                         \
        else fprintf(stderr, "%s invert build failed\n", INVPATH); }     \
      *status = 0;                                                       \
      { void *m = CTOR;                                                  \
        for( int i=0;i<n;i++ ){ ADD(m, PAIRS[i].fwd, PAIRS[i].narg, A);  \
                                ADD(m, PAIRS[i].inv, PAIRS[i].narg, A); }\
        if(astOK) dump((AstMapping*)m, CANCPATH);                        \
        else fprintf(stderr, "%s cancel build failed\n", CANCPATH); }    \
   } while(0)

int main( void ) {
   int sv = 0; int *status = &sv; astWatch( status );

   static Pair spec[] = {
      {"FRTOVL","VLTOFR",1},{"ENTOFR","FRTOEN",0},{"WNTOFR","FRTOWN",0},
      {"WVTOFR","FRTOWV",0},{"AWTOFR","FRTOAW",0},{"VRTOVL","VLTOVR",0},
      {"VOTOVL","VLTOVO",0},{"ZOTOVL","VLTOZO",0},{"BTTOVL","VLTOBT",0},
      {"TPF2HL","HLF2TP",6},{"USF2HL","HLF2US",3},{"GEF2HL","HLF2GE",3},
      {"BYF2HL","HLF2BY",3},{"LKF2HL","HLF2LK",2},{"LDF2HL","HLF2LD",2} };
   BUILD( astSpecMap(1,0,""), astSpecAdd, spec,
          "/tmp/sweep/specmap_invert_normalize.map",
          "/tmp/sweep/specmap_pair_cancel.map" );

   static Pair sla[] = {
      {"FK54Z","FK45Z",1},{"FK5HZ","HFK5Z",1},{"AMP","MAP",2},
      {"ADDET","SUBET",1},{"ECLEQ","EQECL",1},{"GALEQ","EQGAL",0},
      {"GALSUP","SUPGAL",0},{"H2E","E2H",2},{"R2H","H2R",1},
      {"HEEQ","EQHE",1} };
   BUILD( astSlaMap(0,""), astSlaAdd, sla,
          "/tmp/sweep/slamap_invert_normalize.map",
          "/tmp/sweep/slamap_pair_cancel.map" );

   static Pair tim[] = {
      {"TTTOTDB","TDBTOTT",5},{"TDBTOTCB","TCBTOTDB",1},{"TTTOTCG","TCGTOTT",1},
      {"TAITOTT","TTTOTAI",1},{"UTTOGMST","GMSTTOUT",1},
      {"GMSTTOLMST","LMSTTOGMST",3},{"LASTTOLMST","LMSTTOLAST",3},
      {"UTTOUTC","UTCTOUT",1} };
   BUILD( astTimeMap(0,""), astTimeAdd, tim,
          "/tmp/sweep/timemap_invert_normalize.map",
          "/tmp/sweep/timemap_pair_cancel.map" );

   fprintf( stderr, "done\n" );
   return 0;
}
