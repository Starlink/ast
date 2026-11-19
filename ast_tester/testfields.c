/*
*  Name:
*     testfields

*  Purpose:
*     Check the functions that split formatted axis values into fields.

*  Description:
*     The Plot class formats a value for every numerical label on an axis,
*     and then uses astFields, astAbbrev and astCentre to split each label
*     into fields, to remove the leading fields that it shares with its
*     neighbour, and to choose a position for the major ticks.  The text
*     given to these functions need not match the Format of the axis, and
*     each case here used to read beyond the end of the text, or use a
*     value that had not been set.  Most of the faults do not change a
*     result, so the test is intended to be run with the address and
*     undefined-behavior sanitizers enabled.  The text for each case is held
*     in memory of exactly its own size, so that a read beyond it is seen.
*     The cases cover:
*
*       - A sky axis value with fewer fields than its Format, for which the
*         pointer to the missing separator was used.
*       - A sky axis value that does not match its Format, given to
*         astAbbrev with no value to compare it with, for which the field
*         before the first was returned.
*       - A value with leading spaces on a simple axis, for which the end of
*         the value was looked for beyond the end of the text.
*       - A sky axis value with a leading plus sign and a single field, and
*         one that is blank.
*       - A date with an invalid month or year, which was accepted and gave
*         a value that had not been set.  A day that is not valid for its
*         month is still accepted, and is counted on from the start of the
*         month, since TimeFrame astCentre relies on that.
*       - TimeFrame astCentre with a gap so small that the dates either
*         side of the value are formatted the same, where the search for a
*         field that differs ran beyond the end of both, and with a gap
*         that crosses the end of a month, where the date it made had a
*         month of zero.

*  Notes:
*     - This is a native C test with no Fortran original.
*     - astFields, astAbbrev and astCentre are protected functions, so the
*       test defines astCLASS and uses the protected interface, in which
*       each function takes a pointer to the inherited status.
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sae_par.h"
#include "ast_err.h"

#define astCLASS testfields
#include "object.h"
#include "frame.h"
#include "skyframe.h"
#include "timeframe.h"

void astBegin_( void );
void astEnd_( int * );

/* Report a failure of the test described by "what". */
static int failed( const char *what ) {
   fprintf( stderr, "testfields: %s\n", what );
   return 1;
}

/* Return a copy of "text" in memory that has no room for anything else. */
static char *exact( const char *text ) {
   char *result = malloc( strlen( text ) + 1 );
   if ( result ) strcpy( result, text );
   return result;
}

/* Split "text" into fields using axis "axis" of "frm" and the Format
   "fmt", and check that there are "nfld" fields, that the first starts
   "off" characters into the text and has "len" characters, and that the
   value is "value" (which is not checked if it is AST__BAD).  Returns the
   number of failures. */
static int check_fields( AstFrame *frm, int axis, const char *fmt,
                         const char *text, int nfld, int off, int len,
                         double value, int *status ) {
   char *copy = exact( text );
   char *fields[ 3 ];
   char what[ 200 ];
   double val = 0.0;
   int n;
   int nc[ 3 ];
   int nfail = 0;

   if ( !copy ) return failed( "Out of memory" );
   n = astFields_( frm, axis, fmt, copy, 3, fields, nc, &val, status );
   snprintf( what, sizeof( what ), "astFields with Format \"%s\" and text "
             "\"%s\" gave %d fields%s", fmt, text, n,
             astOK ? "" : " and an error" );
   if ( !astOK ) {
      astClearStatus;
      nfail += failed( what );
   } else if ( n != nfld ) {
      nfail += failed( what );
   } else if ( n > 0 ) {
      if ( fields[ 0 ] - copy != off || nc[ 0 ] != len ) {
         snprintf( what, sizeof( what ), "astFields with Format \"%s\" and "
                   "text \"%s\" put the first field at offset %d with "
                   "length %d", fmt, text, (int)( fields[ 0 ] - copy ),
                   nc[ 0 ] );
         nfail += failed( what );
      }
      if ( value != AST__BAD && !( fabs( val - value ) <= 1.0E-12 ) ) {
         snprintf( what, sizeof( what ), "astFields with Format \"%s\" and "
                   "text \"%s\" gave the value %.17g", fmt, text, val );
         nfail += failed( what );
      }
   }
   free( copy );
   return nfail;
}

/* A formatted sky axis value is split at the separators that its Format
   calls for.  If one was missing, the pointer that should have been to
   the separator was null, and was still moved on and used. */
static int test_sky( int *status ) {
   AstFrame *frm = (AstFrame *) astSkyFrame_( " ", status );
   const char *result;
   char *copy;
   double deg = 3.14159265358979323846/180.0;
   int nfail = 0;

   nfail += check_fields( frm, 0, "hms", "12:34", 0, 0, 0, AST__BAD, status );
   nfail += check_fields( frm, 0, "hms", "12", 0, 0, 0, AST__BAD, status );
   nfail += check_fields( frm, 0, "dms", "1:2", 0, 0, 0, AST__BAD, status );
   nfail += check_fields( frm, 0, "lhms", "12h34", 0, 0, 0, AST__BAD,
                          status );
   nfail += check_fields( frm, 0, "hms", "12:34:56", 3, 0, 2, AST__BAD,
                          status );

/* With no first value, astAbbrev returns the last field of the second
   value.  If that value has no fields, the element before the first of
   an array was returned.  The whole value is now returned. */
   copy = exact( "a:b:c" );
   result = astAbbrev_( frm, 0, "hms", NULL, copy, status );
   if ( !astOK || result != copy ) {
      nfail += failed( "astAbbrev did not return the whole of a sky axis "
                       "value that has no fields" );
      if ( !astOK ) astClearStatus;
   }
   free( copy );

/* If the Format has one field, the field runs to the end of the text.
   Its end was found by counting from the character after a leading sign,
   which was allowed for if the sign was a minus and not if it was a
   plus, so that the text after a plus sign was read from beyond its
   end. */
   nfail += check_fields( frm, 0, "d", "+12.5", 1, 0, 5, 12.5*deg, status );
   nfail += check_fields( frm, 0, "d", "-12.5", 1, 0, 5, -12.5*deg, status );
   nfail += check_fields( frm, 0, "d", "12.5 ", 1, 0, 4, 12.5*deg, status );
   nfail += check_fields( frm, 0, "h", "+3", 1, 0, 2, 45.0*deg, status );

/* No number can be read from a blank value, and the value of the field
   was then one that had not been set.  It is now zero. */
   nfail += check_fields( frm, 0, "d", "", 1, 0, 0, 0.0, status );
   nfail += check_fields( frm, 0, "d", "   ", 1, 3, 0, 0.0, status );

   return nfail;
}

/* A value on a simple axis is one field, without the spaces round it.
   The end of the field was looked for from the end of the text plus the
   number of leading spaces, which is beyond the text. */
static int test_simple( int *status ) {
   AstFrame *frm = astFrame_( 2, " ", status );
   int nfail = 0;

   nfail += check_fields( frm, 0, "%g", "  -7 ", 1, 2, 2, -7.0, status );
   nfail += check_fields( frm, 0, "%10.3f", "     1.500", 1, 5, 5, 1.5,
                          status );
   nfail += check_fields( frm, 0, "%g", "42", 1, 0, 2, 42.0, status );

   return nfail;
}

/* Check that astCentre returns "wanted" for the value and gap given, or
   if "wanted" is AST__BAD, a value that is within the gap of the value
   given.  Returns the number of failures. */
static int check_centre( AstFrame *frm, double value, double gap,
                         double wanted, int *status ) {
   char what[ 200 ];
   double centre = astCentre_( frm, 0, value, gap, status );

   snprintf( what, sizeof( what ), "astCentre gave %.17g for the value "
             "%.17g and the gap %g%s", centre, value, gap,
             astOK ? "" : ", and an error" );
   if ( !astOK ) {
      astClearStatus;
      return failed( what );
   }
   if ( wanted != AST__BAD ) {
      if ( centre != wanted ) return failed( what );
   } else if ( !( fabs( centre - value ) <= gap ) ) {
      return failed( what );
   }
   return 0;
}

/* The year, month and day of a date are converted to a Modified Julian
   Date by a function that reports a year or month that is not valid,
   without returning a date.  The check on what it reported looked for
   the values used by an earlier library, so every date was accepted, and
   for such a year or month the result had not been set. */
static int test_dates( int *status ) {
   AstFrame *frm = (AstFrame *) astTimeFrame_( "Format=iso.0", status );
   AstFrame *plain = astFrame_( 1, " ", status );
   double value;
   int nc;
   int nfail = 0;

/* An ordinary date. */
   value = 0.0;
   nc = astUnformat_( frm, 0, "2012-03-05", &value, status );
   if ( !astOK || nc != 10 || value != 55991.0 ) {
      nfail += failed( "astUnformat did not read \"2012-03-05\" as MJD "
                       "55991" );
      if ( !astOK ) astClearStatus;
   }

/* A day of zero is read as the day before the first of the month. */
   value = 0.0;
   nc = astUnformat_( frm, 0, "2012-03-00", &value, status );
   if ( !astOK || nc != 10 || value != 55986.0 ) {
      nfail += failed( "astUnformat did not read \"2012-03-00\" as MJD "
                       "55986" );
      if ( !astOK ) astClearStatus;
   }

/* A month of zero or thirteen is not a date, so the whole of the text
   cannot be read. */
   nc = astUnformat_( frm, 0, "2012-00-01", &value, status );
   if ( !astOK || nc >= 10 ) {
      nfail += failed( "astUnformat read the whole of \"2012-00-01\"" );
      if ( !astOK ) astClearStatus;
   }
   nc = astUnformat_( frm, 0, "2012-13-01", &value, status );
   if ( !astOK || nc >= 10 ) {
      nfail += failed( "astUnformat read the whole of \"2012-13-01\"" );
      if ( !astOK ) astClearStatus;
   }

/* The same text is not accepted as an Epoch. */
   astSetC_( (AstObject *) plain, "Epoch", "2012-13-01", status );
   if ( astOK ) {
      nfail += failed( "An Epoch of \"2012-13-01\" was accepted" );
   } else {
      astClearStatus;
   }

/* astCentre formats the dates one gap either side of the value, and makes
   a date from the fields they share, with zero for the first field that
   differs and for those after it.  With a gap of a day the day is zero,
   which gives the day before the first of the month. */
   nfail += check_centre( frm, 56000.37, 1.0, 55986.0, status );
   nfail += check_centre( frm, 56000.37, 0.01, 56000.0, status );

/* With a gap of thirty days the month is zero as well, which is not a
   date.  The result was the value that had not been set.  It is now the
   result for an axis that is not formatted as a date. */
   nfail += check_centre( frm, 56000.37, 30.0, AST__BAD, status );
   nfail += check_centre( frm, 56000.37, 400.0, AST__BAD, status );

/* With a gap too small to change the formatted date no field differs, and
   the search for one ran beyond the end of both dates.  The value itself
   is returned.  A bad value is formatted with no fields at all. */
   nfail += check_centre( frm, 56000.37, 1.0E-9, 56000.37, status );
   nfail += check_centre( frm, AST__BAD, 1.0, AST__BAD, status );

   return nfail;
}

int main( void ) {
   int status_value = SAI__OK;
   int *status = &status_value;
   int nfail = 0;

   astWatch( status );
   astBegin_();

   nfail += test_sky( status );
   nfail += test_simple( status );
   nfail += test_dates( status );

   if ( !astOK ) astClearStatus;
   astEnd_( status );

   if ( nfail == 0 ) {
      printf( " All field tests passed\n" );
      return 0;
   }
   printf( "Field tests failed (%d)\n", nfail );
   return 1;
}
