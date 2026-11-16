/*
*  Name:
*     testfitswrite

*  Purpose:
*     Check that an unusual FrameSet is written to a FITS header without
*     harm.

*  Description:
*     astWrite turns a FrameSet into the keywords of a foreign FITS
*     encoding by taking its Mapping apart.  A FrameSet that cannot be
*     described must simply not be written: astWrite returns zero, reports
*     no error and leaves the FitsChan as it found it.  Each case here is a
*     small FrameSet for which the library used to read or write memory it
*     should not have, call itself without end, report an internal error,
*     or write a header that did not describe the FrameSet.  The test is
*     intended to be run with the address and undefined-behavior
*     sanitizers enabled, which is how the memory faults are detected.
*     The cases cover:
*
*       - A compound Mapping that has an Ident, and so is not taken apart.
*       - Axes that could only be described by a look-up table ("-TAB"),
*         where one world axis depends on two pixel axes, two world axes
*         depend on one pixel axis, or the Mapping has no inverse.
*       - Polynomial distortion ("-SIP"): a polynomial with a linear term
*         after its last term of higher order, and celestial axes that
*         share a FrameSet with axes fed by fewer pixel axes than remain.
*       - Version letters: more Frames than there are letters, and an
*         offset sky coordinate system beside a Frame that claims a letter
*         through its Ident.
*       - FITS-PC: secondary axis descriptions that have a reference
*         frame, projection parameters or a number of axes that the
*         primary description does not have.
*       - FITS-CLASS: a spectral cube with a fourth axis and a single
*         spatial pixel.
*       - Text: a sky axis Symbol too long for a keyword comment, and a
*         Domain that differs from an axis Symbol only in its last
*         character.
*       - A write that fails after it has replaced cards that were already
*         in the FitsChan, which must be put back.

*  Notes:
*     - This is a native C test with no Fortran original.
*     - Public C API only.
*/

#include <stdio.h>
#include <string.h>
#include "sae_par.h"
#include "ast.h"
#include "ast_err.h"

/* The largest number of cards compared by same_cards. */
#define MXCARD 100

/* Read a FrameSet from a FITS-WCS header given as an array of cards that
   ends with a NULL pointer. */
static AstFrameSet *read_cards( const char *cards[] ) {
   AstFitsChan *fc = astFitsChan( NULL, NULL, "Encoding=FITS-WCS" );
   AstFrameSet *fs;
   int i;
   for ( i = 0; cards[ i ]; i++ ) astPutFits( fc, cards[ i ], 0 );
   astClear( fc, "Card" );
   fs = astRead( fc );
   fc = astAnnul( fc );
   return fs;
}

/* A FrameSet that describes a 100 by 200 pixel image of the sky with a
   gnomonic projection. */
static AstFrameSet *tan_frameset( void ) {
   const char *cards[] = { "NAXIS   = 2", "NAXIS1  = 100", "NAXIS2  = 200",
      "CRPIX1  = 50", "CRPIX2  = 100", "CRVAL1  = 10", "CRVAL2  = 20",
      "CTYPE1  = 'RA---TAN'", "CTYPE2  = 'DEC--TAN'", "CDELT1  = -0.001",
      "CDELT2  = 0.001", NULL };
   return read_cards( cards );
}

/* Look for a keyword in a FitsChan.  Returns non-zero if it is there, and
   then "card" holds its card. */
static int find_card( AstFitsChan *fc, const char *name, char card[ 81 ] ) {
   card[ 0 ] = 0;
   astClear( fc, "Card" );
   return astFindFits( fc, name, card, 0 );
}

/* Write a FrameSet to a FitsChan.  No error may be reported, and if
   "expect" is zero or one then astWrite must return that value.  Returns
   0 if so, and 1 otherwise.  Any error reported is cleared. */
static int check_write( const char *what, AstFitsChan *fc, AstFrameSet *fs,
                        int expect ) {
   int nfail = 0;
   int nobj = astWrite( fc, fs );
   if ( !astOK ) {
      fprintf( stderr, "testfitswrite: %s: astWrite reported error %d\n",
               what, astStatus );
      astClearStatus;
      nfail++;
   } else if ( expect >= 0 && nobj != expect ) {
      fprintf( stderr, "testfitswrite: %s: astWrite returned %d, not %d\n",
               what, nobj, expect );
      nfail++;
   }
   return nfail;
}

/* Require a keyword to be present in, or absent from, a FitsChan.  If
   "value" is not NULL and the keyword should be present, its value must
   also start with that text.  Returns 0 if so, and 1 otherwise. */
static int check_card( const char *what, AstFitsChan *fc, const char *name,
                       int present, const char *value ) {
   char card[ 81 ];
   int found = find_card( fc, name, card );
   if ( !astOK ) astClearStatus;
   if ( found && !present ) {
      fprintf( stderr, "testfitswrite: %s: unexpected card \"%s\"\n", what,
               card );
      return 1;
   } else if ( !found && present ) {
      fprintf( stderr, "testfitswrite: %s: no %s card\n", what, name );
      return 1;
   } else if ( found && value &&
               strncmp( card + 10, value, strlen( value ) ) ) {
      fprintf( stderr, "testfitswrite: %s: card \"%s\" does not have a "
               "value starting %s\n", what, card, value );
      return 1;
   }
   return 0;
}

/* Copy the cards of a FitsChan into an array.  Returns their number. */
static int list_cards( AstFitsChan *fc, char cards[ MXCARD ][ 81 ] ) {
   int n = 0;
   astClear( fc, "Card" );
   while ( n < MXCARD && astFindFits( fc, "%f", cards[ n ], 1 ) ) n++;
   return n;
}

/* A compound Mapping with an Ident is left whole by the functions that
   take a Mapping apart, so the search for the projection inside it must
   not go on for ever. */
static int test_ident( void ) {
   AstFitsChan *fc;
   AstFrameSet *fs;
   AstFrameSet *fs2;
   AstMapping *map;
   int nfail = 0;

   astBegin;
   fs = tan_frameset();
   map = (AstMapping *) astCmpMap( astGetMapping( fs, AST__BASE, AST__CURRENT ),
                                   astUnitMap( 2, " " ), 1, "Ident=X" );
   fs2 = astFrameSet( astGetFrame( fs, AST__BASE ), " " );
   astAddFrame( fs2, AST__BASE, map, astGetFrame( fs, AST__CURRENT ) );
   fc = astFitsChan( NULL, NULL, "Encoding=FITS-WCS" );
   nfail += check_write( "a CmpMap with an Ident", fc, fs2, -1 );
   astEnd;
   return nfail;
}

/* Axes that are neither linear nor a known projection can only be written
   as look-up tables, and only when each world axis depends on a single
   pixel axis. */
static int test_tab( void ) {
   AstFitsChan *fc;
   AstFrameSet *fs;
   AstMapping *lm;
   AstMapping *pm;
   double con = 1.0;
   double lut1[ 5 ] = { 0.0, 0.001, 0.004, 0.009, 0.016 };
   double lut2[ 5 ] = { 0.0, 1.0, 4.0, 9.0, 16.0 };
   double lut3[ 3 ] = { 0.0, 1.0, 0.0 };
   double lut4[ 5 ] = { -1.0, 0.0, 1.0, 8.0, 27.0 };
   int inperm1[ 4 ] = { 1, 1, 2, 2 };
   int inperm2[ 3 ] = { 1, 1, 2 };
   int inperm3[ 3 ] = { -1, -1, 1 };
   int outperm1[ 2 ] = { 1, 3 };
   int outperm3[ 2 ] = { 3, 3 };
   int nfail = 0;

/* Four pixel axes and a sky position.  Each sky axis is a table of one
   pixel axis, and the inverse sets two pixel axes from each sky axis. */
   astBegin;
   pm = (AstMapping *) astPermMap( 4, inperm1, 2, outperm1, NULL, " " );
   lm = (AstMapping *) astCmpMap( astLutMap( 5, lut1, 1.0, 1.0, " " ),
                                  astLutMap( 5, lut1, 1.0, 1.0, " " ), 0, " " );
   fs = astFrameSet( astFrame( 4, "Domain=GRID" ), " " );
   astAddFrame( fs, AST__BASE, astCmpMap( pm, lm, 1, " " ), astSkyFrame( " " ) );
   fc = astFitsChan( NULL, NULL,
                     "Encoding=FITS-WCS,TabOK=1,FitsAxisOrder=<copy>" );
   astPutFits( fc, "NAXIS1  = 5", 0 );
   astPutFits( fc, "NAXIS2  = 5", 0 );
   astPutFits( fc, "NAXIS3  = 5", 0 );
   astPutFits( fc, "NAXIS4  = 5", 0 );
   nfail += check_write( "sky axes tabulated against two pixel axes each",
                         fc, fs, -1 );
   astEnd;

/* Three pixel axes and two other axes, the first a table of one pixel
   axis whose inverse sets two pixel axes. */
   astBegin;
   pm = (AstMapping *) astPermMap( 3, inperm2, 2, outperm1, NULL, " " );
   lm = (AstMapping *) astCmpMap( astLutMap( 5, lut2, 1.0, 1.0, " " ),
                                  astUnitMap( 1, " " ), 0, " " );
   fs = astFrameSet( astFrame( 3, "Domain=GRID" ), " " );
   astAddFrame( fs, AST__BASE, astCmpMap( pm, lm, 1, " " ),
                astFrame( 2, "Domain=WORLD" ) );
   fc = astFitsChan( NULL, NULL,
                     "Encoding=FITS-WCS,TabOK=1,FitsAxisOrder=<copy>" );
   astPutFits( fc, "NAXIS1  = 5", 0 );
   astPutFits( fc, "NAXIS2  = 5", 0 );
   astPutFits( fc, "NAXIS3  = 5", 0 );
   nfail += check_write( "an axis tabulated against two pixel axes", fc, fs,
                         -1 );
   astEnd;

/* Three pixel axes and two world axes, both of which depend on the third
   pixel axis, through a table that has no inverse. */
   astBegin;
   pm = (AstMapping *) astPermMap( 3, inperm3, 2, outperm3, &con, " " );
   lm = (AstMapping *) astCmpMap( astLutMap( 3, lut3, 1.0, 1.0, " " ),
                                  astUnitMap( 1, " " ), 0, " " );
   fs = astFrameSet( astFrame( 3, "Domain=GRID" ), " " );
   astAddFrame( fs, AST__BASE, astCmpMap( pm, lm, 1, " " ),
                astFrame( 2, "Domain=WORLD" ) );
   fc = astFitsChan( NULL, NULL, "Encoding=FITS-WCS,TabOK=1" );
   nfail += check_write( "two axes that depend on one pixel axis", fc, fs,
                         -1 );
   astEnd;

/* A single axis that is not linear, with look-up tables not allowed, so
   that no pixel scale can be found for it. */
   astBegin;
   fs = astFrameSet( astFrame( 1, "Domain=GRID" ), " " );
   astAddFrame( fs, AST__BASE, astLutMap( 5, lut4, 1.0, 1.0, " " ),
                astFrame( 1, "Domain=WORLD" ) );
   fc = astFitsChan( NULL, NULL, "Encoding=FITS-WCS" );
   astPutFits( fc, "NAXIS1  = 5", 0 );
   nfail += check_write( "an axis that is not linear", fc, fs, 0 );
   astEnd;
   return nfail;
}

/* A sky Mapping that starts with the polynomial "coeffs" (with "ncoeff"
   terms) about the reference pixel of the gnomonic projection in "fs". */
static AstMapping *poly_sky( AstFrameSet *fs, int ncoeff, double *coeffs ) {
   double shift[ 2 ] = { -50.0, -100.0 };
   double unshift[ 2 ] = { 50.0, 100.0 };
   AstMapping *poly = (AstMapping *) astPolyMap( 2, 2, ncoeff, coeffs, 0, NULL,
                                                 "IterInverse=1" );
   return (AstMapping *) astCmpMap(
             astCmpMap( astShiftMap( 2, shift, " " ), poly, 1, " " ),
             astCmpMap( astShiftMap( 2, unshift, " " ),
                        astGetMapping( fs, AST__BASE, AST__CURRENT ), 1, " " ),
             1, " " );
}

/* Polynomial distortion of the celestial axes. */
static int test_sip( void ) {
   AstFitsChan *fc;
   AstFrameSet *fs;
   AstFrameSet *fs2;
   AstMapping *map;
   AstMapping *pm;
   double coeffs1[] = { 2.0, 1, 1, 0,   1.0E-6, 1, 2, 0,   2.0, 2, 0, 1 };
   double coeffs2[] = { 1.0, 1, 1, 0,   1.0E-3, 1, 2, 0,   1.0, 2, 0, 1 };
   double con = 1.0;
   int inperm[ 4 ] = { 1, 2, 3, -1 };
   int outperm[ 4 ] = { 1, 2, 3, 3 };
   int nfail = 0;

/* A polynomial of the pixel coordinates in which the last term is linear
   and follows the only term of higher order. */
   astBegin;
   fs = tan_frameset();
   map = (AstMapping *) astCmpMap( astPolyMap( 2, 2, 3, coeffs1, 0, NULL,
                                               "IterInverse=1" ),
                                   astGetMapping( fs, AST__BASE, AST__CURRENT ),
                                   1, " " );
   fs2 = astFrameSet( astGetFrame( fs, AST__BASE ), " " );
   astAddFrame( fs2, AST__BASE, map, astGetFrame( fs, AST__CURRENT ) );
   fc = astFitsChan( NULL, NULL, "Encoding=FITS-WCS" );
   astPutFits( fc, "NAXIS1  = 100", 0 );
   astPutFits( fc, "NAXIS2  = 200", 0 );
   nfail += check_write( "a polynomial that ends with a linear term", fc, fs2,
                         1 );
   astEnd;

/* Distorted celestial axes with a third axis that has its own pixel axis.
   This can be described, and is. */
   astBegin;
   fs = tan_frameset();
   map = (AstMapping *) astCmpMap( poly_sky( fs, 3, coeffs2 ),
                                   astUnitMap( 1, " " ), 0, " " );
   fs2 = astFrameSet( astFrame( 3, "Domain=GRID" ), " " );
   astAddFrame( fs2, AST__BASE, map,
                astCmpFrame( astGetFrame( fs, AST__CURRENT ),
                             astFrame( 1, "Domain=OTHER" ), " " ) );
   fc = astFitsChan( NULL, NULL, "Encoding=FITS-WCS" );
   astPutFits( fc, "NAXIS1  = 100", 0 );
   astPutFits( fc, "NAXIS2  = 200", 0 );
   astPutFits( fc, "NAXIS3  = 10", 0 );
   nfail += check_write( "distorted sky axes and a third axis", fc, fs2, 1 );
   nfail += check_card( "distorted sky axes and a third axis", fc, "A_ORDER",
                        1, NULL );
   nfail += check_card( "distorted sky axes and a third axis", fc, "CTYPE1",
                        1, "'RA---TAN-SIP'" );
   astEnd;

/* Distorted celestial axes with two more world axes that both copy the
   third of four pixel axes, so that the other axes are fed by one pixel
   axis where two remain. */
   astBegin;
   fs = tan_frameset();
   pm = (AstMapping *) astPermMap( 4, inperm, 4, outperm, &con, " " );
   map = (AstMapping *) astCmpMap( pm,
                                   astCmpMap( poly_sky( fs, 3, coeffs2 ),
                                              astUnitMap( 2, " " ), 0, " " ),
                                   1, " " );
   fs2 = astFrameSet( astFrame( 4, "Domain=GRID" ), " " );
   astAddFrame( fs2, AST__BASE, map,
                astCmpFrame( astGetFrame( fs, AST__CURRENT ),
                             astFrame( 2, "Domain=OTHER" ), " " ) );
   fc = astFitsChan( NULL, NULL, "Encoding=FITS-WCS" );
   astPutFits( fc, "NAXIS1  = 100", 0 );
   astPutFits( fc, "NAXIS2  = 200", 0 );
   astPutFits( fc, "NAXIS3  = 10", 0 );
   astPutFits( fc, "NAXIS4  = 1", 0 );
   nfail += check_write( "distorted sky axes and two axes that share a "
                         "pixel axis", fc, fs2, -1 );
   astEnd;
   return nfail;
}

/* The letters that distinguish secondary axis descriptions. */
static int test_versions( void ) {
   AstFitsChan *fc;
   AstFrameSet *fs;
   int icur;
   int i;
   int nfail = 0;

/* Twenty seven Frames that need a letter, where there are twenty six
   letters.  The first twenty six are described and the last is not. */
   astBegin;
   fs = tan_frameset();
   icur = astGetI( fs, "Current" );
   for ( i = 0; i < 27; i++ ) {
      astAddFrame( fs, AST__BASE, astZoomMap( 2, 2.0 + i, " " ),
                   astFrame( 2, "Domain=EXTRA" ) );
   }
   astSetI( fs, "Current", icur );
   fc = astFitsChan( NULL, NULL, "Encoding=FITS-WCS" );
   nfail += check_write( "more Frames than version letters", fc, fs, 1 );
   nfail += check_card( "more Frames than version letters", fc, "CDELT1A", 1,
                        "                 2.0" );
   nfail += check_card( "more Frames than version letters", fc, "CDELT1Z", 1,
                        "                27.0" );
   astEnd;

/* An offset sky coordinate system is also described as an absolute one,
   which needs a letter.  It must not be given the letter that another
   Frame has claimed with its Ident. */
   astBegin;
   fs = tan_frameset();
   astSetC( fs, "SkyRef", "0:40:00, 20:00:00" );
   astSetC( fs, "SkyRefIs", "Origin" );
   icur = astGetI( fs, "Current" );
   astAddFrame( fs, AST__BASE, astZoomMap( 2, 3.0, " " ),
                astFrame( 2, "Domain=FOO,Ident=A" ) );
   astSetI( fs, "Current", icur );
   fc = astFitsChan( NULL, NULL, "Encoding=FITS-WCS" );
   nfail += check_write( "an offset sky system and an Ident", fc, fs, 1 );
   nfail += check_card( "an offset sky system and an Ident", fc, "CTYPE1A", 1,
                        "'FOO1" );
   nfail += check_card( "an offset sky system and an Ident", fc, "RADESYSA",
                        0, NULL );
   nfail += check_card( "an offset sky system and an Ident", fc, "CTYPE1B", 1,
                        "'RA---TAN'" );
   astEnd;
   return nfail;
}

/* FITS-PC shares the reference frame, the projection parameters and the
   PC matrix of the primary axis description with every secondary one, so a
   secondary description that needs different ones cannot be written. */
static int test_pc( void ) {
   const char *sin_cards[] = { "NAXIS   = 2", "NAXIS1  = 100",
      "NAXIS2  = 200", "CRPIX1  = 50", "CRPIX2  = 100", "CRVAL1  = 10",
      "CRVAL2  = 20", "CTYPE1  = 'RA---TAN'", "CTYPE2  = 'DEC--TAN'",
      "CDELT1  = -0.001", "CDELT2  = 0.001", "CRPIX1A = 50", "CRPIX2A = 100",
      "CRVAL1A = 10", "CRVAL2A = 20", "CTYPE1A = 'RA---SIN'",
      "CTYPE2A = 'DEC--SIN'", "CDELT1A = -0.001", "CDELT2A = 0.001",
      "PV2_1A  = 0.1", "PV2_2A  = 0.2", NULL };
   AstFitsChan *fc;
   AstFrameSet *fs;
   AstFrameSet *sky;
   int nfail = 0;

/* A primary description with no celestial axes, and so no reference
   frame, and a secondary one that is a sky coordinate system. */
   astBegin;
   sky = tan_frameset();
   fs = astFrameSet( astGetFrame( sky, AST__BASE ), " " );
   astAddFrame( fs, AST__BASE, astGetMapping( sky, AST__BASE, AST__CURRENT ),
                astGetFrame( sky, AST__CURRENT ) );
   astAddFrame( fs, AST__BASE, astUnitMap( 2, " " ),
                astFrame( 2, "Domain=PRIM" ) );
   fc = astFitsChan( NULL, NULL, "Encoding=FITS-PC" );
   astPutFits( fc, "NAXIS1  = 100", 0 );
   astPutFits( fc, "NAXIS2  = 200", 0 );
   nfail += check_write( "a secondary sky system with no primary one", fc,
                         fs, 1 );
   nfail += check_card( "a secondary sky system with no primary one", fc,
                        "C1YPE1", 0, NULL );
   astEnd;

/* A primary description that has no projection parameters, and a
   secondary one that has. */
   astBegin;
   fs = read_cards( sin_cards );
   fc = astFitsChan( NULL, NULL, "Encoding=FITS-PC" );
   astPutFits( fc, "NAXIS1  = 100", 0 );
   astPutFits( fc, "NAXIS2  = 200", 0 );
   nfail += check_write( "secondary projection parameters with no primary "
                         "ones", fc, fs, 1 );
   nfail += check_card( "secondary projection parameters with no primary "
                        "ones", fc, "C1YPE1", 0, NULL );
   astEnd;

/* A primary description with two axes and a secondary one with three. */
   astBegin;
   fs = astFrameSet( astFrame( 2, "Domain=GRID" ), " " );
   astAddFrame( fs, AST__BASE, astPermMap( 2, (int[]){ 1, 2 }, 3,
                                           (int[]){ 1, 2, -1 },
                                           (double[]){ 4.0 }, " " ),
                astFrame( 3, "Domain=ALT" ) );
   astAddFrame( fs, AST__BASE, astUnitMap( 2, " " ),
                astFrame( 2, "Domain=PRIM" ) );
   fc = astFitsChan( NULL, NULL, "Encoding=FITS-PC" );
   astPutFits( fc, "NAXIS1  = 10", 0 );
   astPutFits( fc, "NAXIS2  = 10", 0 );
   nfail += check_write( "a secondary description with an extra axis", fc,
                         fs, 1 );
   nfail += check_card( "a secondary description with an extra axis", fc,
                        "C1YPE1", 0, NULL );
   astEnd;
   return nfail;
}

/* FITS-CLASS describes three axes.  A FrameSet with a fourth axis, and a
   single spatial pixel, cannot be written but must not report that the
   wrong number of coordinates was transformed. */
static int test_class( void ) {
   AstFitsChan *fc;
   AstFrame *freq;
   AstFrame *wcs;
   AstFrameSet *fs;
   AstFrameSet *sky;
   AstMapping *map;
   const char *sky_cards[] = { "NAXIS   = 2", "CRPIX1  = 5", "CRPIX2  = 5",
      "CRVAL1  = 10", "CRVAL2  = 20", "CDELT1  = -0.001", "CDELT2  = 0.001",
      "CTYPE1  = 'RA---TAN'", "CTYPE2  = 'DEC--TAN'", "RADESYS = 'FK5'",
      "EQUINOX = 2000", NULL };
   double ina = 0.0;
   double inb = 1.0;
   double outa = 1.0E11;
   double outb = 1.0E11 + 1.0E6;
   int nfail = 0;

   astBegin;
   sky = read_cards( sky_cards );
   freq = astFrame( 1, "Symbol(1)=FREQ,Unit(1)=Hz,Label(1)=Frequency" );
   map = (AstMapping *) astCmpMap( astWinMap( 1, &ina, &inb, &outa, &outb, " " ),
                                   astGetMapping( sky, AST__BASE, AST__CURRENT ),
                                   0, " " );
   map = (AstMapping *) astCmpMap( map, astUnitMap( 1, " " ), 0, " " );
   wcs = (AstFrame *) astCmpFrame( astCmpFrame( freq,
                                                astGetFrame( sky, AST__CURRENT ),
                                                " " ),
                                   astFrame( 1, "Symbol(1)=STOKES" ), " " );
   fs = astFrameSet( astFrame( 4, "Domain=GRID" ), " " );
   astAddFrame( fs, AST__BASE, map, wcs );
   fc = astFitsChan( NULL, NULL, "Encoding=FITS-CLASS" );
   astPutFits( fc, "NAXIS   = 4", 0 );
   astPutFits( fc, "NAXIS1  = 10", 0 );
   astPutFits( fc, "NAXIS2  = 1", 0 );
   astPutFits( fc, "NAXIS3  = 1", 0 );
   astPutFits( fc, "NAXIS4  = 1", 0 );
   nfail += check_write( "a four axis cube with one spatial pixel", fc, fs,
                         -1 );
   astEnd;
   return nfail;
}

/* Text taken from the FrameSet. */
static int test_text( void ) {
   AstFitsChan *fc;
   AstFrameSet *fs;
   int nfail = 0;

/* The comment of an offset sky axis includes its Symbol, which can be
   longer than the comment. */
   astBegin;
   fs = tan_frameset();
   astSetC( fs, "SkyRef", "0:40:00, 20:00:00" );
   astSetC( fs, "SkyRefIs", "Origin" );
   astSetC( fs, "Symbol(1)", "AVERYLONGSYMBOLNAMEAVERYLONGSYMBOLNAME"
                             "AVERYLONGSYMBOLNAMEAVERYLONGSYMBOLNAME"
                             "AVERYLONGSYMBOLNAME" );
   fc = astFitsChan( NULL, NULL, "Encoding=FITS-WCS" );
   nfail += check_write( "a long Symbol for an offset sky axis", fc, fs, 1 );
   nfail += check_card( "a long Symbol for an offset sky axis", fc, "CTYPE1",
                        1, "'OFLN-TAN'" );
   astEnd;

/* WCSNAME is left out when it says the same as CTYPE1.  A Domain that
   differs from the axis Symbol in its last character does not. */
   astBegin;
   fs = astFrameSet( astFrame( 1, "Domain=GRID" ), " " );
   astAddFrame( fs, AST__BASE, astUnitMap( 1, " " ),
                astFrame( 1, "Symbol(1)=ABCD,Domain=ABCX" ) );
   fc = astFitsChan( NULL, NULL, "Encoding=FITS-WCS" );
   nfail += check_write( "a Domain that is nearly the axis Symbol", fc, fs,
                         1 );
   nfail += check_card( "a Domain that is nearly the axis Symbol", fc,
                        "WCSNAME", 1, "'ABCX" );
   astEnd;

/* A Domain that is the axis Symbol, apart from its case and its spaces,
   does say the same. */
   astBegin;
   fs = astFrameSet( astFrame( 1, "Domain=GRID" ), " " );
   astAddFrame( fs, AST__BASE, astUnitMap( 1, " " ),
                astFrame( 1, "Symbol(1)=ABCD,Domain=abcd" ) );
   fc = astFitsChan( NULL, NULL, "Encoding=FITS-WCS" );
   nfail += check_write( "a Domain that is the axis Symbol", fc, fs, 1 );
   nfail += check_card( "a Domain that is the axis Symbol", fc, "WCSNAME", 0,
                        NULL );
   astEnd;

/* A blank Symbol written over a card that holds a blank value. */
   astBegin;
   fs = astFrameSet( astFrame( 1, "Domain=GRID" ), " " );
   astAddFrame( fs, AST__BASE, astUnitMap( 1, " " ),
                astFrame( 1, "Symbol(1)= " ) );
   fc = astFitsChan( NULL, NULL, "Encoding=FITS-WCS" );
   astPutFits( fc, "CTYPE1  = ' '                  / old", 0 );
   nfail += check_write( "a blank Symbol over a blank CTYPE", fc, fs, -1 );
   astEnd;
   return nfail;
}

/* A FrameSet that cannot be written must leave the FitsChan as it was,
   including cards for the keywords that were stored before the write was
   abandoned. */
static int test_failed_write( void ) {
   const char *tpn_cards[] = { "NAXIS   = 2", "NAXIS1  = 100",
      "NAXIS2  = 200", "CRPIX1  = 50", "CRPIX2  = 100", "CRVAL1  = 10",
      "CRVAL2  = 20", "CTYPE1  = 'RA---TPN'", "CTYPE2  = 'DEC--TPN'",
      "CDELT1  = -0.001", "CDELT2  = 0.001", "PV1_1   = 1.0",
      "PV2_1   = 1.0", NULL };
   AstFitsChan *fc;
   AstFrameSet *fs;
   AstFrameSet *sky;
   char after[ MXCARD ][ 81 ];
   char before[ MXCARD ][ 81 ];
   int i;
   int nafter;
   int nbefore;
   int nfail = 0;

/* FITS-PC cannot describe the TPN projection, and finds that out after
   it has stored the reference pixel and reference value. */
   astBegin;
   fs = read_cards( tpn_cards );
   fc = astFitsChan( NULL, NULL, "Encoding=FITS-PC" );
   astPutFits( fc, "NAXIS1  = 100", 0 );
   astPutFits( fc, "CRVAL1  = 5 / mine", 0 );
   astPutFits( fc, "NAXIS2  = 200", 0 );
   astPutFits( fc, "CRPIX2  = 7.5 / also mine", 0 );
   nbefore = list_cards( fc, before );
   astClear( fc, "Card" );
   nfail += check_write( "a projection that FITS-PC cannot describe", fc, fs,
                         0 );
   nafter = list_cards( fc, after );
   if ( nafter != nbefore ) {
      fprintf( stderr, "testfitswrite: a failed write left %d cards where "
               "there were %d\n", nafter, nbefore );
      nfail++;
   } else {
      for ( i = 0; i < nbefore; i++ ) {
         if ( strcmp( before[ i ], after[ i ] ) ) {
            fprintf( stderr, "testfitswrite: a failed write changed the "
                     "card \"%s\" to \"%s\"\n", before[ i ], after[ i ] );
            nfail++;
         }
      }
   }
   astEnd;

/* The same for a secondary axis description that is abandoned while the
   primary one is written.  The primary description has no celestial axes
   and the secondary one is a sky coordinate system, which FITS-PC cannot
   describe without a reference frame for the primary one. */
   astBegin;
   sky = tan_frameset();
   fs = astFrameSet( astGetFrame( sky, AST__BASE ), " " );
   astAddFrame( fs, AST__BASE, astGetMapping( sky, AST__BASE, AST__CURRENT ),
                astGetFrame( sky, AST__CURRENT ) );
   astAddFrame( fs, AST__BASE, astUnitMap( 2, " " ),
                astFrame( 2, "Domain=PRIM" ) );
   fc = astFitsChan( NULL, NULL, "Encoding=FITS-PC" );
   astPutFits( fc, "NAXIS1  = 100", 0 );
   astPutFits( fc, "NAXIS2  = 200", 0 );
   astPutFits( fc, "C1VAL1  = 7 / mine", 0 );
   astClear( fc, "Card" );
   nfail += check_write( "a secondary description that FITS-PC cannot "
                         "describe", fc, fs, 1 );
   nfail += check_card( "a secondary description that FITS-PC cannot "
                        "describe", fc, "C1VAL1", 1, "                   7 " );
   nfail += check_card( "a secondary description that FITS-PC cannot "
                        "describe", fc, "CRVAL1", 1, NULL );
   astEnd;
   return nfail;
}

int main( void ) {
   int status_value = 0;
   int *status = &status_value;
   int nfail = 0;

   astWatch( status );
   astBegin;

   nfail += test_ident();
   nfail += test_tab();
   nfail += test_sip();
   nfail += test_versions();
   nfail += test_pc();
   nfail += test_class();
   nfail += test_text();
   nfail += test_failed_write();

   astEnd;

   if ( !astOK ) {
      fprintf( stderr, "testfitswrite: unexpected error status %d\n",
               astStatus );
      astClearStatus;
      nfail++;
   }

   if ( nfail == 0 ) {
      printf( " All FITS writing tests passed\n" );
   } else {
      printf( "testfitswrite: %d check(s) failed\n", nfail );
   }
   return nfail ? 1 : 0;
}
