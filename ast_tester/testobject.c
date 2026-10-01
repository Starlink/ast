/* Test program for the AstObject base class. */

#include "ast.h"
#include <libgen.h>
#include <stdio.h>
#include <string.h>

/* Helper to extract basename from a file path, since __FILE__ may
   include a directory prefix in out-of-source builds. */
static const char *baseName( const char *path ) {
   static char buffer[ 1024 ];
   if( !path ) return path;
   snprintf( buffer, sizeof( buffer ), "%s", path );
   return basename( buffer );
}

/* Verify that an Object survives a checkdump (astToString followed by
   astFromString) round trip, and that the reloaded Object compares equal
   to the original. */
static void TestCheckToString( void ) {
   AstSkyFrame *sf = astSkyFrame( " " );
   AstFrame *bf = astFrame( 2, "Domain=SKY" );
   AstFrameSet *fs = astConvert( bf, sf, " " );
   AstFrameSet *fs2 = NULL;
   char *pickle1 = NULL;
   char *pickle2 = NULL;

   if( !fs ) {
      if( astOK )
         astError( AST__INTER, "TestCheckToString: astConvert failed.\n" );  /* LCOV_EXCL_LINE */
   } else {
      pickle1 = astToString( fs );
      fs2 = astFromString( pickle1 );
      pickle2 = astToString( fs2 );

      if( pickle1 && pickle2 ) {
         if( strcmp( pickle1, pickle2 ) && astOK )
            astError( AST__INTER, "TestCheckToString: round-tripped strings differ.\n" );  /* LCOV_EXCL_LINE */
      } else if( astOK ) {
         astError( AST__INTER, "TestCheckToString: astToString returned NULL.\n" );  /* LCOV_EXCL_LINE */
      }

      if( fs2 && !astEqual( fs, fs2 ) && astOK )
         astError( AST__INTER, "TestCheckToString: reloaded Object is not equal to the original.\n" );  /* LCOV_EXCL_LINE */
   }

   pickle1 = astFree( pickle1 );
   pickle2 = astFree( pickle2 );
   sf = astAnnul( sf );
   bf = astAnnul( bf );

   if( fs )
      fs = astAnnul( fs );

   if( fs2 )
      fs2 = astAnnul( fs2 );
}

/* Verify that astCreatedAt reports the routine, file and line at which an
   Object was created. */
static void TestCreatedAt( void ) {
   const char *routine;
   const char *file;
   int line;
   int bf_line = __LINE__ + 1;
   AstFrame *bf = astFrame( 2, "Domain=SKY" );

   astCreatedAt( bf, &routine, &file, &line );
   if( ( !routine || strcmp( routine, "TestCreatedAt" ) ) && astOK )
      astError( AST__INTER, "TestCreatedAt: routine is '%s'.\n", routine ? routine : "<NULL>" );  /* LCOV_EXCL_LINE */
   if( ( !file || strcmp( baseName( file ), "testobject.c" ) ) && astOK )
      astError( AST__INTER, "TestCreatedAt: file is '%s'.\n", file ? file : "<NULL>" );  /* LCOV_EXCL_LINE */
   if( line != bf_line && astOK )
      astError( AST__INTER, "TestCreatedAt: line is %d, expected %d.\n", line, bf_line );  /* LCOV_EXCL_LINE */

   bf = astAnnul( bf );
}

/* Verify that astActiveObjects reports the Objects active within the
   current context, keyed by class, including the creation information for
   each one. The objects are created within their own AST context so the
   "current" filter yields a deterministic set regardless of any objects
   left alive by other tests. */
static void TestActiveObjects( void ) {
   AstSkyFrame *sf;
   AstFrame *bf;
   AstFrameSet *fs;
   AstFrameSet *fs2;
   AstKeyMap *km;
   const char *routine;
   const char *file;
   const char *key;
   int bf_line;
   int fs2_line;
   int idx;
   int jdx;
   int found;
   int line;
   int nkey;
   void *p;

   astBegin;

   sf = astSkyFrame( " " );
   bf_line = __LINE__ + 1;
   bf = astFrame( 2, "Domain=SKY" );
   fs = astConvert( bf, sf, " " );
   fs2_line = __LINE__ + 1;
   fs2 = astCopy( fs );
   km = astActiveObjects( NULL, 0, 1 );

   if( !km && astOK ) {
      astError( AST__INTER, "TestActiveObjects: astActiveObjects returned NULL.\n" );  /* LCOV_EXCL_LINE */
   } else {
      nkey = astMapSize( km );
      if( nkey != 3 && astOK )
         astError( AST__INTER, "TestActiveObjects: nkey is %d, expected 3.\n", nkey );  /* LCOV_EXCL_LINE */

      astSetC( km, "SortBy", "KeyUp" );
      for( idx = 0; idx < nkey; idx++ ){
         key = astMapKey( km, idx );
         if( idx == 0 ) {
            if( strcmp( key, "Frame" ) && astOK ) {
               astError( AST__INTER, "TestActiveObjects: key 0 is '%s'.\n", key );  /* LCOV_EXCL_LINE */
            } else if( astMapLength( km, key ) != 1 && astOK ) {
               astError( AST__INTER, "TestActiveObjects: Frame count is %d.\n", astMapLength( km, key ) );  /* LCOV_EXCL_LINE */
            } else if( ( !astMapGetElemP( km, key, 0, &p ) || ( p != bf ) ) && astOK ) {
               astError( AST__INTER, "TestActiveObjects: Frame pointer mismatch.\n" );  /* LCOV_EXCL_LINE */
            } else {
               astCreatedAt( p, &routine, &file, &line );
               if( ( !routine || strcmp( routine, "TestActiveObjects" ) ) && astOK )
                  astError( AST__INTER, "TestActiveObjects: Frame routine is '%s'.\n", routine ? routine : "<NULL>" );  /* LCOV_EXCL_LINE */
               if( ( !file || strcmp( baseName( file ), "testobject.c" ) ) && astOK )
                  astError( AST__INTER, "TestActiveObjects: Frame file is '%s'.\n", file ? file : "<NULL>" );  /* LCOV_EXCL_LINE */
               if( line != bf_line && astOK )
                  astError( AST__INTER, "TestActiveObjects: Frame line is %d, expected %d.\n", line, bf_line );  /* LCOV_EXCL_LINE */
            }
         } else if( idx == 1 ) {
            if( strcmp( key, "FrameSet" ) && astOK ) {
               astError( AST__INTER, "TestActiveObjects: key 1 is '%s'.\n", key );  /* LCOV_EXCL_LINE */
            } else if( astMapLength( km, key ) != 2 && astOK ) {
               astError( AST__INTER, "TestActiveObjects: FrameSet count is %d.\n", astMapLength( km, key ) );  /* LCOV_EXCL_LINE */
            } else {

/* The two FrameSets are both active, but their order within the list is
   determined by the order in which the underlying object handles happen
   to be allocated, which is not the order of creation. Locate the entry
   that corresponds to "fs2" (the FrameSet created in this routine by
   astCopy) and verify that its recorded creation information is correct. */
               found = 0;
               for( jdx = 0; jdx < 2; jdx++ ) {
                  if( astMapGetElemP( km, key, jdx, &p ) && p == fs2 ) {
                     found = 1;
                     astCreatedAt( p, &routine, &file, &line );
                     if( ( !routine || strcmp( routine, "TestActiveObjects" ) ) && astOK )
                        astError( AST__INTER, "TestActiveObjects: FrameSet routine is '%s'.\n", routine ? routine : "<NULL>" );  /* LCOV_EXCL_LINE */
                     if( ( !file || strcmp( baseName( file ), "testobject.c" ) ) && astOK )
                        astError( AST__INTER, "TestActiveObjects: FrameSet file is '%s'.\n", file ? file : "<NULL>" );  /* LCOV_EXCL_LINE */
                     if( line != fs2_line && astOK )
                        astError( AST__INTER, "TestActiveObjects: FrameSet line is %d, expected %d.\n", line, fs2_line );  /* LCOV_EXCL_LINE */
                  }
               }
               if( !found && astOK )
                  astError( AST__INTER, "TestActiveObjects: fs2 not found among the active FrameSets.\n" );  /* LCOV_EXCL_LINE */
            }
         } else {
            if( strcmp( key, "SkyFrame" ) && astOK ) {
               astError( AST__INTER, "TestActiveObjects: key 2 is '%s'.\n", key );  /* LCOV_EXCL_LINE */
            } else if( astMapLength( km, key ) != 1 && astOK ) {
               astError( AST__INTER, "TestActiveObjects: SkyFrame count is %d.\n", astMapLength( km, key ) );  /* LCOV_EXCL_LINE */
            }
         }
      }
      km = astAnnul( km );
   }

   astEnd;
}

/* The number of error messages reported while RecordError is installed
   as the error handler using astSetPutErr, and the last of them. */
#define MAX_MESSAGE_LEN 400
static int nmessage = 0;
static char last_message[ MAX_MESSAGE_LEN ];

static void RecordError( int status_value, const char *message ) {
   (void) status_value;
   nmessage++;
   snprintf( last_message, sizeof( last_message ), "%s", message );
}

/* astClone does not require the Object to be locked by the calling
   thread, but it must still reject a pointer that does not identify an
   Object, such as one that has been annulled, returning a null pointer
   and reporting an error. */
static void TestCloneInvalid( void ) {
   static const char *expected = "This pointer has been annulled, or the "
                                 "associated Object deleted.";
   AstFrame *bf;
   AstFrame *stale;
   AstObject *clone;
   int status_value;

   bf = astFrame( 2, " " );
   stale = bf;
   bf = astAnnul( bf );

   nmessage = 0;
   astSetPutErr( RecordError );
   clone = astClone( stale );
   status_value = astStatus;
   astClearStatus;
   astSetPutErr( NULL );

   if( clone && astOK ) {
      astError( AST__INTER, "TestCloneInvalid: astClone returned a pointer.\n" );  /* LCOV_EXCL_LINE */
   }
   if( status_value != AST__OBJIN && astOK ) {
      astError( AST__INTER, "TestCloneInvalid: status is %d, expected AST__OBJIN (%d).\n", status_value, AST__OBJIN );  /* LCOV_EXCL_LINE */
   }
   if( ( !nmessage || strcmp( last_message, expected ) ) && astOK ) {
      astError( AST__INTER, "TestCloneInvalid: last error message is '%s'.\n", nmessage ? last_message : "<none>" );  /* LCOV_EXCL_LINE */
   }
}

/* Like other AST functions, astClone must do nothing if the error status
   is already set: it must return a null pointer, leave the Object's
   reference count and the error status unchanged, and report nothing. */
static void TestCloneStatusSet( void ) {
   AstFrame *bf;
   AstObject *clone;
   int refcount;
   int status_value;

   bf = astFrame( 2, " " );
   refcount = astGetI( bf, "RefCount" );

   nmessage = 0;
   astSetPutErr( RecordError );
   astSetStatus( AST__ATTIN );  /* Arbitrary status not set by astClone itself */
   clone = astClone( bf );
   status_value = astStatus;
   astClearStatus;
   astSetPutErr( NULL );

   if( clone && astOK ) {
      astError( AST__INTER, "TestCloneStatusSet: astClone returned a pointer.\n" );  /* LCOV_EXCL_LINE */
   }
   if( status_value != AST__ATTIN && astOK ) {
      astError( AST__INTER, "TestCloneStatusSet: status is %d, expected it unchanged (%d).\n", status_value, AST__ATTIN );  /* LCOV_EXCL_LINE */
   }
   if( nmessage && astOK ) {  /* Should not be any message */
      astError( AST__INTER, "TestCloneStatusSet: astClone reported '%s'.\n", last_message );  /* LCOV_EXCL_LINE */
   }
   if( astGetI( bf, "RefCount" ) != refcount && astOK ) {
      astError( AST__INTER, "TestCloneStatusSet: RefCount changed from %d to %d.\n", refcount, astGetI( bf, "RefCount" ) );  /* LCOV_EXCL_LINE */
   }

   bf = astAnnul( bf );
}

int main( void ){
   astBegin;

   TestCheckToString();
   TestCreatedAt();
   TestActiveObjects();
   TestCloneInvalid();
   TestCloneStatusSet();

   astEnd;

   if( astOK ) {
      printf(" All Object tests passed\n");
      return 0;
   } else {
      printf("Object tests failed\n");
      return 1;
   }
   return astOK ? 0 : 1;
}
