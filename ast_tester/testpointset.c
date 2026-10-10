/*
 *  Test the public PointSet constructor, which takes a printf-style
 *  options string like every other public constructor.
 */
#include "ast.h"
#include <stdio.h>
#include <string.h>

static void stopit( int *status, const char *text ) {
   if( *status != 0 ) return;
   *status = 1;
   printf( "%s\n", text );
}

int main( void ) {
   int status_value = 0;
   int *status = &status_value;
   AstPointSet *ps;
   double **ptr;

   astWatch( status );
   astBegin;

   ps = astPointSet( 3, 2, "ID=%s", "pset" );
   if( !astOK || !ps ) {
      astClearStatus;
      stopit( status, "Error 1" );
   } else {
      if( astGetI( ps, "Npoint" ) != 3 ) stopit( status, "Error 2" );
      if( astGetI( ps, "Ncoord" ) != 2 ) stopit( status, "Error 3" );
      if( strcmp( astGetC( ps, "ID" ), "pset" ) ) stopit( status, "Error 4" );
      ptr = astGetPoints( ps );
      if( !ptr ) {
         stopit( status, "Error 5" );
      } else {
         ptr[ 1 ][ 2 ] = 42.0;
         if( astGetPoints( ps )[ 1 ][ 2 ] != 42.0 ) stopit( status, "Error 6" );
      }
   }

   astEnd;

   if( *status == 0 ) {
      printf( " All PointSet tests passed\n" );
   } else {
      printf( "PointSet tests failed\n" );
   }
   return *status;
}
