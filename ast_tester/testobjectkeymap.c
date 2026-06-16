/* Test program for the protected astGetKeyMap method of the Object class.
 *
 * astGetKeyMap is a protected method, so (like testaxis.c and
 * testunitnorm.c) this test is compiled with astCLASS defined and uses the
 * protected interface directly.
 *
 * The behaviour covered here includes:
 *   - Lazy creation of an empty KeyMap on first access
 *   - A stable, borrowed reference returned by repeated requests
 *   - A deep-copy of the KeyMap taken by astCopy
 *   - Serialisation that preserves a non-empty KeyMap but omits an empty
 *     one (ensuring that existing Object dumps are unmodified).
 */

#include <stdio.h>
#include <string.h>
#include "sae_par.h"
#include "ast_err.h"

/* Just ensure astCLASS is defined; what it's defined as is arbitrary for this
   purpose */
#define astCLASS testobjectkeymap
#include "memory.h"
#include "object.h"
#include "frame.h"
#include "keymap.h"

void astBegin_( void );
void astEnd_( int * );

static void TestCreate( int *status ) {
   AstFrame *frm = astFrame( 2, " ", status );
   AstKeyMap *km1;
   AstKeyMap *km2;

/* A new Object lazily creates an empty KeyMap on the first request. */
   km1 = astGetKeyMap( frm );
   if( !km1 && astOK ) {
      astError( AST__INTER, "TestCreate: astGetKeyMap returned NULL.\n",
                status );
   } else if( astMapSize( km1 ) != 0 && astOK ) {
      astError( AST__INTER, "TestCreate: New KeyMap not empty (size %d).\n",
                status, astMapSize( km1 ) );
   }

/* Repeated requests return the same borrowed pointer. */
   km2 = astGetKeyMap( frm );
   if( km1 != km2 && astOK ) {
      astError( AST__INTER,
                "TestCreate: astGetKeyMap returned different pointers.\n",
                status );
   }

   frm = astAnnul( frm );
}

static void TestCopy( int *status ) {
   AstFrame *frm = astFrame( 2, " ", status );
   AstFrame *copy;
   AstKeyMap *km1;
   AstKeyMap *km2;
   int v;

/* Store some data in the KeyMap. */
   km1 = astGetKeyMap( frm );
   astMapPut0I( km1, "answer", 42, NULL );

/* astCopy takes a deep-copy of the KeyMap. */
   copy = astCopy( frm );
   km2 = astGetKeyMap( copy );
   if( km2 == km1 && astOK ) {
      astError( AST__INTER,
                "TestCopy: Copy shares its KeyMap with the original.\n",
                status );
   }
   v = 0;
   if( ( !astMapGet0I( km2, "answer", &v ) || v != 42 ) && astOK ) {
      astError( AST__INTER,
                "TestCopy: Copied KeyMap is missing data (got %d).\n",
                status, v );
   }

/* Modifying the copy's KeyMap does not affect the original. */
   astMapPut0I( km2, "answer", 7, NULL );
   v = 0;
   astMapGet0I( km1, "answer", &v );
   if( v != 42 && astOK ) {
      astError( AST__INTER,
                "TestCopy: Modifying the copy affected the original "
                "(got %d).\n", status, v );
   }

   frm = astAnnul( frm );
   copy = astAnnul( copy );
}

static void TestToString( int *status ) {
   AstFrame *frm1 = astFrame( 2, " ", status );
   AstFrame *frm2;
   AstFrame *empty;
   AstKeyMap *km1;
   AstKeyMap *km2;
   char *str;
   int v;

/* Store some data in the KeyMap. */
   km1 = astGetKeyMap( frm1 );
   astMapPut0I( km1, "answer", 42, NULL );

/* A non-empty KeyMap survives a round trip. */
   str = astToString( frm1 );
   if( str && !strstr( str, "KeyMap" ) && astOK ) {
      astError( AST__INTER,
                "TestToString: Non-empty KeyMap was not serialised.\n",
                status );
   }
   frm2 = astFromString( str );
   km2 = astGetKeyMap( frm2 );
   v = 0;
   if( ( !astMapGet0I( km2, "answer", &v ) || v != 42 ) && astOK ) {
      astError( AST__INTER,
                "TestToString: Reloaded KeyMap is missing data (got %d).\n",
                status, v );
   }
   str = astFree( str );

/* An empty KeyMap is not serialised at all, so the dump of an Object whose
   KeyMap exists but is empty must be free of any "KeyMap" item. */
   empty = astFrame( 2, " ", status );
   (void) astGetKeyMap( empty );
   str = astToString( empty );
   if( str && strstr( str, "KeyMap" ) && astOK ) {
      astError( AST__INTER, "TestToString: Empty KeyMap was serialised.\n",
                status );
   }
   str = astFree( str );

   frm1 = astAnnul( frm1 );
   frm2 = astAnnul( frm2 );
   empty = astAnnul( empty );
}

int main( void ) {
   int _status = SAI__OK;
   int *status = &_status;
   astWatch( status );
   astBegin_();

   TestCreate( status );
   TestCopy( status );
   TestToString( status );

   astEnd_( status );

   if( astOK ) {
      printf(" All Object KeyMap tests passed\n");
      return 0;
   } else {
      printf("Object KeyMap tests failed\n");
      return 1;
   }
}
