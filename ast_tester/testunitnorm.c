#include <stdio.h>
#include <string.h>
#include "sae_par.h"
#include "ast_err.h"

#define astCLASS testunitnorm
#include "memory.h"
#include "unit.h"

void astBegin_( void );
void astEnd_( int * );

int main( void ) {
   int _status = SAI__OK;
   int* status = &_status;
   astWatch( status );
   astBegin_();

/* Each units string and its normalised form. A normalised string is only
   blank when it is wholly a constant; one that merely starts with a
   number, such as "2*m", is a units expression like any other. */
   {
      const char *cases[][ 2 ] = {
         { "s*(m/s)", "m" },
         { "2*m", "2*m" },
         { "m/2", "0.5*m" },
         { "1.5 furlong", "1.5*furlong" },
         { "1000 m", "km" },
         { "2", "" },
         { "0.5", "" },

/* A reciprocal is folded into the unit's multiplier prefix when one is
   closer; otherwise it is left alone. Without a closer prefix these used
   to recurse without end. */
         { "2/m", "2/m" },
         { "0.5/m", "0.5/m" },
         { "3/s", "3/s" },
         { "2/km", "2/km" },
         { "10/m", "1/dm" },
         { "1000/m", "1/mm" },
      };
      size_t ncase = sizeof( cases )/sizeof( cases[ 0 ] );
      size_t i;

      for( i = 0; i < ncase && astOK; i++ ) {
         const char *result = astUnitNormaliser_( cases[ i ][ 0 ], status );
         if( astOK && ( !result || strcmp( result, cases[ i ][ 1 ] ) ) ) {
            astError_( AST__INTER, "UnitNormaliser gave '%s' for '%s', "
                       "expected '%s'", status, result ? result : "<NULL>",
                       cases[ i ][ 0 ], cases[ i ][ 1 ] );
         }
         result = astFree_( (void *) result, status );
      }
   }

   astEnd_( status );

   if( astOK ) {
      printf(" All UnitNormaliser tests passed\n");
   } else {
      printf("UnitNormaliser tests failed\n");
   }
   return astOK ? 0 : 1;
}
