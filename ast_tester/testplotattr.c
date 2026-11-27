/*
*  Name:
*     testplotattr.c
*
*  Purpose:
*     Check the initial state of Plot attributes on a newly built Plot.
*
*  Description:
*     A Plot attribute that has not been assigned must report itself as
*     unset, so that astTest is false and the attribute is absent from a
*     Native dump.  ForceExterior was left out of the initialisation in
*     astInitPlot and so appeared set (to zero) on every new Plot, which
*     made it impossible for code that wraps a Plot to tell whether the
*     user had chosen a value.  The logging GRF module stands in for a
*     graphics library.
*/

#include <stdio.h>
#include <string.h>
#include "ast.h"
#include "grf_log.h"

int main( void ) {
   int status = 0;
   int nfail = 0;
   AstFrame *frame;
   AstPlot *plot;
   float gbox[4] = { 0.0f, 0.0f, 1.0f, 1.0f };
   double bbox[4] = { 0.5, 0.5, 100.5, 100.5 };
   char *dump;

   astWatch( &status );
   astGrfLogInit( NULL );
   astBegin;

   frame = astFrame( 2, "Domain=PIXEL" );
   plot = astPlot( frame, gbox, bbox, " " );

   if( astTest( plot, "ForceExterior" ) ) {
      printf( "FAIL: ForceExterior is reported as set on a new Plot\n" );
      nfail++;
   }
   if( astGetI( plot, "ForceExterior" ) != 0 ) {
      printf( "FAIL: ForceExterior does not default to 0\n" );
      nfail++;
   }
   dump = astToString( plot );
   if( dump && strstr( dump, "FrcExt" ) ) {
      printf( "FAIL: the dump of a new Plot records ForceExterior\n" );
      nfail++;
   }
   dump = astFree( dump );

   /* A value that is set survives, and clearing it restores the default. */
   astSetI( plot, "ForceExterior", 1 );
   if( !astTest( plot, "ForceExterior" ) || astGetI( plot, "ForceExterior" ) != 1 ) {
      printf( "FAIL: ForceExterior=1 was not stored\n" );
      nfail++;
   }
   astClear( plot, "ForceExterior" );
   if( astTest( plot, "ForceExterior" ) ) {
      printf( "FAIL: ForceExterior is still set after astClear\n" );
      nfail++;
   }

   /* The neighbouring attributes behave the same way. */
   if( astTest( plot, "TickAll" ) || astTest( plot, "Border" ) ||
       astTest( plot, "Invisible" ) || astTest( plot, "Clip" ) ) {
      printf( "FAIL: another Plot attribute is reported as set on a new Plot\n" );
      nfail++;
   }

   astEnd;
   if( !astOK ) {
      printf( "FAIL: AST error\n" );
      nfail++;
   }
   if( nfail == 0 ) {
      printf( "All Plot attribute tests passed\n" );
      return 0;
   }
   printf( "%d Plot attribute test(s) failed\n", nfail );
   return 1;
}
