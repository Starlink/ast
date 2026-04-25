/*
 * testfitschanerrors.c — Table-driven tests for FitsChan warning and error
 * paths. Each test uses a minimal malformed FITS header embedded in C, reads
 * it through a FitsChan, and verifies that the expected warning or error was
 * produced.
 *
 * Two mechanisms:
 *   Warnings: Enabled via Warnings attribute, retrieved via astWarnings()
 *   Errors:   Set astStatus, checked via astOK / astStatus / astClearStatus
 */

#include "ast.h"
#include "ast_err.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

typedef struct {
   const char *name;
   const char *cards;       /* Newline-separated FITS cards (no END needed) */
   const char *warn_cat;    /* Warnings attribute, or NULL for error test */
   int expect_error;        /* Expected astStatus, or 0 for warning test */
   const char *expect_text; /* Substring in warning text, or NULL */
   int skip_read;           /* 1 = check warnings after PutFits, skip astRead */
} BadHeaderTest;

/* Minimal valid 2D TAN header used as a base for many tests. */
#define TAN_BASE \
   "NAXIS1  =                  100\n" \
   "NAXIS2  =                  100\n" \
   "CTYPE1  = 'RA---TAN'\n" \
   "CTYPE2  = 'DEC--TAN'\n" \
   "CRVAL1  =              180.000\n" \
   "CRVAL2  =               45.000\n" \
   "CRPIX1  =               50.000\n" \
   "CRPIX2  =               50.000\n" \
   "CDELT1  =              -0.0100\n" \
   "CDELT2  =               0.0100\n" \
   "RADESYS = 'FK5'\n" \
   "EQUINOX =               2000.0"

static const BadHeaderTest bad_headers[] = {

   /* --- distortion-sip-noorder (line 8815) ---
      -SIP distortion on axis 3 (only valid on axes 1-2) → warning. */
   { "distortion-sip-on-axis3",
     "NAXIS1  =                  256\n"
     "NAXIS2  =                  256\n"
     "NAXIS3  =                  100\n"
     "CTYPE1  = 'RA---TAN-SIP'\n"
     "CTYPE2  = 'DEC--TAN-SIP'\n"
     "CTYPE3  = 'FREQ-SIP'\n"
     "CRVAL1  =              180.000\n"
     "CRVAL2  =               45.000\n"
     "CRVAL3  =          1.4204E+09\n"
     "CRPIX1  =              128.000\n"
     "CRPIX2  =              128.000\n"
     "CRPIX3  =               50.000\n"
     "CDELT1  =              -0.0003\n"
     "CDELT2  =               0.0003\n"
     "CDELT3  =           1.0E+06\n"
     "RADESYS = 'FK5'\n"
     "EQUINOX =               2000.0\n"
     "A_ORDER =                    2\n"
     "A_1_1   =             1.0E-05\n"
     "B_ORDER =                    2\n"
     "B_1_1   =             1.0E-05",
     "distortion", 0, "distortion will be ignored", 0 },

   /* --- distortion-unknown (line 8901) ---
      Unknown distortion suffix -XXX on TAN. */
   { "distortion-unknown",
     "NAXIS1  =                  100\n"
     "NAXIS2  =                  100\n"
     "CTYPE1  = 'RA---TAN-XXX'\n"
     "CTYPE2  = 'DEC--TAN-XXX'\n"
     "CRVAL1  =              180.000\n"
     "CRVAL2  =               45.000\n"
     "CRPIX1  =               50.000\n"
     "CRPIX2  =               50.000\n"
     "CDELT1  =              -0.0100\n"
     "CDELT2  =               0.0100\n"
     "RADESYS = 'FK5'\n"
     "EQUINOX =               2000.0",
     "distortion", 0, "ignores this distortion", 0 },

   /* --- badpv-tan-allzero (line 13980) ---
      TAN with PV on latitude axis but all zero → simple TAN warning. */
   { "badpv-tan-allzero",
     TAN_BASE "\n"
     "PV2_0   =                  0.0\n"
     "PV2_1   =                  0.0",
     "badpv", 0, "distortion coefficients", 0 },

   /* --- noctype-missing (line 35625) ---
      2-axis header with CRPIX/CRVAL for 2 axes but CTYPE2 missing entirely.
      Use non-celestial CTYPE1 so the missing CTYPE2 produces a warning
      rather than a "longitude without latitude" error. */
   { "noctype-missing",
     "NAXIS1  =                  100\n"
     "NAXIS2  =                  100\n"
     "CTYPE1  = 'FREQ'\n"
     "CRVAL1  =          1.4204E+09\n"
     "CRVAL2  =                  0.0\n"
     "CRPIX1  =               50.000\n"
     "CRPIX2  =               50.000\n"
     "CDELT1  =           1.0E+06\n"
     "CDELT2  =               1.0",
     "noctype", 0, "not found for one or more", 0 },

   /* TODO: badmat-singular (line 37577) and zpx-unsupported (line 41570)
      need further investigation into what header structure triggers the
      warning path vs. an error path. */

   /* --- badpv-noncelestial (line 35914) ---
      PV1_5 on longitude axis exceeds mxpar_lon for TAN projection. */
   { "badpv-lonaxis",
     TAN_BASE "\n"
     "PV1_5   =                  1.0",
     "badpv", 0, "not used by", 0 },

   /* --- badkeyvalue-unparseable (line 32548) ---
      A keyword with a value that cannot be parsed as any FITS type. */
   { "badkeyvalue-unparseable",
     TAN_BASE "\n"
     "BADKEY  = not_a_valid_value",
     "badkeyvalue", 0, "keyword value is illegal", 1 },

   /* --- sao-allzero (line 26711) ---
      SAO polynomial TAN with all distortion coefficients zero. */
   { "badpv-sao-allzero",
     "NAXIS1  =                  100\n"
     "NAXIS2  =                  100\n"
     "CTYPE1  = 'RA---TAN'\n"
     "CTYPE2  = 'DEC--TAN'\n"
     "CRVAL1  =              180.000\n"
     "CRVAL2  =               45.000\n"
     "CRPIX1  =               50.000\n"
     "CRPIX2  =               50.000\n"
     "CDELT1  =              -0.0100\n"
     "CDELT2  =               0.0100\n"
     "RADESYS = 'FK5'\n"
     "EQUINOX =               2000.0\n"
     "CO1_1   =                  0.0\n"
     "CO1_2   =                  0.0\n"
     "CO1_3   =                  0.0\n"
     "CO1_4   =                  0.0\n"
     "CO1_5   =                  0.0\n"
     "CO1_6   =                  0.0\n"
     "CO1_7   =                  0.0\n"
     "CO1_8   =                  0.0\n"
     "CO1_9   =                  0.0\n"
     "CO1_10  =                  0.0\n"
     "CO1_11  =                  0.0\n"
     "CO1_12  =                  0.0\n"
     "CO1_13  =                  0.0\n"
     "CO2_1   =                  0.0\n"
     "CO2_2   =                  0.0\n"
     "CO2_3   =                  0.0\n"
     "CO2_4   =                  0.0\n"
     "CO2_5   =                  0.0\n"
     "CO2_6   =                  0.0\n"
     "CO2_7   =                  0.0\n"
     "CO2_8   =                  0.0\n"
     "CO2_9   =                  0.0\n"
     "CO2_10  =                  0.0\n"
     "CO2_11  =                  0.0\n"
     "CO2_12  =                  0.0\n"
     "CO2_13  =                  0.0",
     "badpv", 0, "SAO encoded", 0 },

   /* --- tnx-unsupported (line 31913) ---
      TNX header with cross-term type 1 (full cross-terms, unsupported).
      The 4th number in the lngcor/latcor string is xterms: must != 2. */
   { "tnx-unsupported",
     "NAXIS1  =                  100\n"
     "NAXIS2  =                  100\n"
     "CTYPE1  = 'RA---TNX'\n"
     "CTYPE2  = 'DEC--TNX'\n"
     "CRVAL1  =              180.000\n"
     "CRVAL2  =               45.000\n"
     "CRPIX1  =               50.000\n"
     "CRPIX2  =               50.000\n"
     "CD1_1   =              -0.0100\n"
     "CD2_2   =               0.0100\n"
     "WAT0_001= 'system=image'\n"
     "WAT1_001= 'wtype=tnx axtype=ra lngcor = \"3. 3. 3. 1. 0. 1. 0. 1. 0. 0. 0. 0. 0.\"'\n"
     "WAT2_001= 'wtype=tnx axtype=dec latcor = \"3. 3. 3. 1. 0. 1. 0. 1. 0. 0. 0. 0. 0.\"'",
     "tnx", 0, "unsupported IRAF", 0 },

};

#define NTESTS (sizeof(bad_headers)/sizeof(bad_headers[0]))

static int test_bad_header( const BadHeaderTest *t, int *status ) {
   int ok = 1;
   if( *status != 0 ) return 0;

   astBegin;
   AstFitsChan *fc = astFitsChan( NULL, NULL, " " );

   /* Enable warnings before inserting cards — some warnings (badkeyvalue,
      badkeyname) are issued during astPutFits, not during astRead. */
   if( t->warn_cat ) {
      astSetC( fc, "Warnings", t->warn_cat );
   }

   /* Split cards on \n and insert each via astPutFits. */
   {
      char *buf = astStore( NULL, t->cards, strlen( t->cards ) + 1 );
      char *line = buf;
      while( line && *line ) {
         char *nl = strchr( line, '\n' );
         if( nl ) *nl = '\0';
         while( *line == ' ' ) line++;
         if( *line ) astPutFits( fc, line, 0 );
         line = nl ? nl + 1 : NULL;
      }
      buf = astFree( buf );
   }
   astPutFits( fc, "END", 0 );

   astClear( fc, "Card" );
   AstObject *obj = NULL;
   if( !t->skip_read ) {
      obj = (AstObject *)astRead( fc );
      /* For warning tests, astRead may also set an error status. */
      if( !t->expect_error && !astOK ) astClearStatus;
   }

   if( t->expect_error ) {
      if( astOK ) {
         printf( "  %s: expected error %d but astRead succeeded\n",
                 t->name, t->expect_error );
         ok = 0;
      } else if( astStatus != t->expect_error ) {
         printf( "  %s: expected error %d got %d\n",
                 t->name, t->expect_error, astStatus );
         astClearStatus;
         ok = 0;
      } else {
         astClearStatus;
      }
   } else if( t->warn_cat ) {
      AstKeyMap *km = (AstKeyMap *)astWarnings( fc );
      if( !km || astMapSize( km ) == 0 ) {
         printf( "  %s: no warnings issued (expected '%s')\n",
                 t->name, t->warn_cat );
         ok = 0;
      }
      if( ok && t->expect_text && km ) {
         int found = 0;
         int nw = astMapSize( km );
         for( int i = 1; i <= nw; i++ ) {
            char key[20];
            const char *val;
            sprintf( key, "Warning_%d", i );
            if( astMapGet0C( km, key, &val ) &&
                strstr( val, t->expect_text ) ) {
               found = 1;
            }
         }
         if( !found ) {
            printf( "  %s: warning text does not contain '%s'\n",
                    t->name, t->expect_text );
            ok = 0;
         }
      }
      if( km ) km = astAnnul( km );
   }

   if( obj ) obj = astAnnul( obj );
   fc = astAnnul( fc );
   astEnd;
   return ok;
}

int main() {
   int status = 0;
   int fails = 0;

   astWatch( &status );

   for( size_t i = 0; i < NTESTS; i++ ) {
      if( !test_bad_header( &bad_headers[i], &status ) ) {
         fails++;
      }
   }

   if( fails ) {
      printf( "%d of %zu bad-header tests failed\n", fails, NTESTS );
      return 1;
   }
   printf( "All %zu bad-header tests passed\n", NTESTS );
   return 0;
}
