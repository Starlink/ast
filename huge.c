#include "ast.h"

#define LEN 2200000000

main(){
   int64_t i, j;
   char *p;
   double shift[2];
   int result;
   int lbnd[2];
   int ubnd[2];

   char *array_in = astMalloc( 2*LEN );
   char *array_out = astMalloc( 2*LEN );
   if( astOK ) {

      printf("Initialising array...\n");
      p = array_in;
      for( j = 1; j <= 2; j++ ) {
         for( i = 1; i <= LEN; i++ ) *(p++) = i*j;
      }

      shift[0] = 1;
      shift[1] = 0;
      AstShiftMap *map = astShiftMap( 2, shift, " " );

      printf("Resampling array...\n");
      lbnd[0] = 1;
      lbnd[1] = 1;
      ubnd[0] = LEN;
      ubnd[1] = 2;
      result = astResampleUB( map, 2, lbnd, ubnd, array_in, NULL,
                              AST__NEAREST, NULL, NULL, 0, 0.0, 0, 0, 2,
                              lbnd, ubnd, lbnd, ubnd, array_out, NULL );

      printf("%d %d\n", array_out[LEN], array_out[LEN+1] );
      printf("%d %d\n", array_out[0], array_out[1] );

   } else {
      printf("Failed to allocate memory\n");
   }

   array_in = astFree( array_in );
   array_out = astFree( array_out );

}
