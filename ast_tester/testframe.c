#include <string.h>
#include "sae_par.h"
#include "ast.h"
#include "ast_err.h"

int main() {
   int status = SAI__OK;
   astWatch( &status );

   AstFrame* frame = astFrame( 1, " " );
   astSetC( frame, "Unit(1)", "s*(m/s)" );
   const char* result = astGetC( frame, "NormUnit(1)" );

   if( strcmp( result, "m" ) ) {
      astError( AST__INTER, "NormUnit did not give expected result" );
   }

   if( astOK ) {
      printf(" All Frame tests passed\n");
   } else {
      printf("Frame tests failed\n");
   }
}
