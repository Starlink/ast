#include <string.h>
#include "sae_par.h"
#include "ast_err.h"

#define astCLASS testunitnorm
#include "unit.h"

int main() {
   int _status = SAI__OK;
   int* status = &_status;
   astWatch( status );

   const char* result = astUnitNormaliser_("s*(m/s)", status);

   if( strcmp( result, "m" ) ) {
      astError_( AST__INTER, "UnitNormaliser did not give expected result", status );
   }

   if( astOK ) {
      printf(" All UnitNormaliser tests passed\n");
   } else {
      printf("UnitNormaliser tests failed\n");
   }
}
