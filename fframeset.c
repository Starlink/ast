/*
*+
*  Name:
*     fframeset.c

*  Purpose:
*     Define a FORTRAN 77 interface to the AST FrameSet class.

*  Type of Module:
*     C source file.

*  Description:
*     This file defines FORTRAN 77-callable C functions which provide
*     a public FORTRAN 77 interface to the FrameSet class.

*  Routines Defined:
*     AST_ADDFRAME
*     AST_FRAMESET
*     AST_GETFRAME
*     AST_GETMAPPING
*     AST_ISAFRAMESET
*     AST_REMAPFRAME
*     AST_REMOVEFRAME

*  Copyright:
*     <COPYRIGHT_STATEMENT>

*  Authors:
*     RFWS: R.F. Warren-Smith (Starlink)

*  History:
*     1-AUG-1996 (RFWS):
*        Original version.
*/

/* Define the astFORTRAN77 macro which prevents error messages from
   AST C functions from reporting the file and line number where the
   error occurred (since these would refer to this file, they would
   not be useful). */
#define astFORTRAN77

/* Header files. */
/* ============= */
#include "f77.h"                 /* FORTRAN <-> C interface macros (SUN/209) */
#include "c2f77.h"               /* F77 <-> C support functions/macros */
#include "error.h"               /* Error reporting facilities */
#include "memory.h"              /* Memory handling facilities */
#include "mapping.h"             /* C interface to the Mapping class */
#include "frame.h"               /* C interface to the Frame class */
#include "frameset.h"            /* C interface to the FrameSet class */

F77_SUBROUTINE(ast_addframe)( INTEGER(THIS),
                              INTEGER(IFRAME),
                              INTEGER(MAP),
                              INTEGER(FRAME),
                              INTEGER(STATUS) ) {
   GENPTR_INTEGER(THIS)
   GENPTR_INTEGER(IFRAME)
   GENPTR_INTEGER(MAP)
   GENPTR_INTEGER(FRAME)

   astAt( "AST_ADDFRAME", NULL, 0 );
   astWatchSTATUS(
      astAddFrame( astI2P( *THIS ), *IFRAME, astI2P( *MAP ),
                   astI2P( *FRAME ) );
   )
}

F77_INTEGER_FUNCTION(ast_frameset)( INTEGER(FRAME),
                                    CHARACTER(OPTIONS),
                                    INTEGER(STATUS)
                                    TRAIL(OPTIONS) ) {
   GENPTR_INTEGER(FRAME)
   GENPTR_CHARACTER(OPTIONS)
   F77_INTEGER_TYPE(RESULT);
   char *options;
   int i;

   astAt( "AST_FRAMESET", NULL, 0 );
   astWatchSTATUS(
      options = astString( OPTIONS, OPTIONS_length );

/* Change ',' to '\n' (see AST_SET in fobject.c for why). */
      if ( astOK ) {
         for ( i = 0; options[ i ]; i++ ) {
            if ( options[ i ] == ',' ) options[ i ] = '\n';
         }
      }
      RESULT = astP2I( astFrameSet( astI2P( *FRAME ), "%s", options ) );
      astFree( options );
   )
   return RESULT;
}

F77_INTEGER_FUNCTION(ast_getframe)( INTEGER(THIS),
                                    INTEGER(IFRAME),
                                    INTEGER(STATUS) ) {
   GENPTR_INTEGER(THIS)
   GENPTR_INTEGER(IFRAME)
   F77_INTEGER_TYPE(RESULT);

   astAt( "AST_GETFRAME", NULL, 0 );
   astWatchSTATUS(
      RESULT = astP2I( astGetFrame( astI2P( *THIS ), *IFRAME ) );
   )
   return RESULT;
}

F77_INTEGER_FUNCTION(ast_getmapping)( INTEGER(THIS),
                                      INTEGER(IFRAME1),
                                      INTEGER(IFRAME2),
                                      INTEGER(STATUS) ) {
   GENPTR_INTEGER(THIS)
   GENPTR_INTEGER(IFRAME1)
   GENPTR_INTEGER(IFRAME2)
   F77_INTEGER_TYPE(RESULT);

   astAt( "AST_GETMAPPING", NULL, 0 );
   astWatchSTATUS(
      RESULT = astP2I( astGetMapping( astI2P( *THIS ), *IFRAME1, *IFRAME2 ) );
   )
   return RESULT;
}

F77_LOGICAL_FUNCTION(ast_isaframeset)( INTEGER(THIS),
                                       INTEGER(STATUS) ) {
   GENPTR_INTEGER(THIS)
   F77_LOGICAL_TYPE(RESULT);

   astAt( "AST_ISAFRAMESET", NULL, 0 );
   astWatchSTATUS(
      RESULT = astIsAFrameSet( astI2P( *THIS ) ) ? F77_TRUE : F77_FALSE;
   )
   return RESULT;
}

F77_SUBROUTINE(ast_remapframe)( INTEGER(THIS),
                                INTEGER(IFRAME),
                                INTEGER(MAP),
                                INTEGER(STATUS) ) {
   GENPTR_INTEGER(THIS)
   GENPTR_INTEGER(IFRAME)
   GENPTR_INTEGER(MAP)

   astAt( "AST_REMAPFRAME", NULL, 0 );
   astWatchSTATUS(
      astRemapFrame( astI2P( *THIS ), *IFRAME, astI2P( *MAP ) );
   )
}

F77_SUBROUTINE(ast_removeframe)( INTEGER(THIS),
                                 INTEGER(IFRAME),
                                 INTEGER(STATUS) ) {
   GENPTR_INTEGER(THIS)
   GENPTR_INTEGER(IFRAME)

   astAt( "AST_REMOVEFRAME", NULL, 0 );
   astWatchSTATUS(
      astRemoveFrame( astI2P( *THIS ), *IFRAME );
   )
}
