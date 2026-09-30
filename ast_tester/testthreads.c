/*
 * Tests passing Objects between threads: a thread must lock an Object
 * before it can use it.
 */
#include "sae_par.h"
#include "ast.h"
#include "ast_err.h"
#include <pthread.h>
#include <stdio.h>
#include <string.h>

#define MAX_MESSAGES 4
#define MAX_MESSAGE_LEN 400

typedef struct Worker {
   AstObject *obj;
   int lock;
   int status;
   int nmessage;
   char messages[ MAX_MESSAGES ][ MAX_MESSAGE_LEN ];
} Worker;

/* Each worker records the error messages AST reports in its thread.  The
   error handler set with astSetPutErr is given no argument identifying
   the caller, so it finds the calling thread's Worker through this
   thread-specific data key. */
static pthread_key_t worker_key;

/* Error handler that records each message in the calling thread's
   Worker. */
static void record_error( int status_value, const char *message ) {
   Worker *data = pthread_getspecific( worker_key );

   (void) status_value;
   if( data && data->nmessage < MAX_MESSAGES ) {
      snprintf( data->messages[ data->nmessage ], MAX_MESSAGE_LEN, "%s",
                message );
      data->nmessage++;
   }
}

/* Print the error messages a Worker recorded. */
static void print_messages( const Worker *data ) {
   int i;

   for( i = 0; i < data->nmessage; i++ ) {
      printf( "      %s\n", data->messages[ i ] );
   }
}

/* Transform a point using the Worker's Object, locking it first if "lock"
   is set, and record the AST status the thread ended with, and any error
   messages reported, in the Worker. */
static void *worker( void *ptr ) {
   double xin, xout;
   Worker *data = (Worker *) ptr;
   int status = SAI__OK;

/* AST maintains a separate status value for each thread.  Watch a
   thread-local variable here so the outcome of the calls below can be
   reported back to the main thread through the shared structure. */
   astWatch( &status );

   data->nmessage = 0;
   pthread_setspecific( worker_key, data );

/* Record error messages through the Worker struct so we can check them
   later in the tests's assertions. */
   astSetPutErr( record_error );

   if( data->lock ) {
      astLock( data->obj, 0 );
   }

   xin = 0;
   astTran1( data->obj, 1, &xin, 1, &xout );

   if( data->lock ) {
      astUnlock( data->obj, 1 );
   }

   data->status = status;
   if( !astOK ) {
      astClearStatus;
   }
   return NULL;
}

/* Run "worker" in two threads at once, one for each Worker, and wait for
   both to finish.  Returns zero if the threads could not be run. */
static int run_workers( Worker *data1, Worker *data2 ) {
   pthread_t thread1;
   pthread_t thread2;

   if( pthread_create( &thread1, NULL, worker, data1 ) ) {
      printf( "FAIL: could not create thread 1\n" );
      return 0;
   }

   if( pthread_create( &thread2, NULL, worker, data2 ) ) {
      printf( "FAIL: could not create thread 2\n" );
      pthread_join( thread1, NULL );
      return 0;
   }

   if( pthread_join( thread1, NULL ) || pthread_join( thread2, NULL ) ) {
      printf( "FAIL: could not join threads\n" );
      return 0;
   }

   return 1;
}

/* Check that a Worker failed with AST__LCKERR, reporting "expected" as the
   last of its error messages (the ones before it give the context of the
   error). */
static int check_lock_error( const char *name, const Worker *data,
                             const char *expected ) {
   if( data->status == AST__LCKERR && data->nmessage > 0 &&
       !strcmp( data->messages[ data->nmessage - 1 ], expected ) ) {
      return 1;
   }

   printf( "FAIL: %s: thread status %d, expected AST__LCKERR (%d), and "
           "error messages:\n", name, data->status, AST__LCKERR );
   print_messages( data );
   printf( "   expected the last to be:\n      %s\n", expected );
   return 0;
}

/* Check that a Worker succeeded without reporting any errors. */
static int check_success( const char *name, const Worker *data ) {
   if( data->status == SAI__OK && data->nmessage == 0 ) {
      return 1;
   }

   printf( "FAIL: %s: thread status %d, expected success, and error "
           "messages:\n", name, data->status );
   print_messages( data );
   return 0;
}

/* A thread cannot use an Object it has not locked, even when no other
   thread has it locked.  The main thread unlocks a UnitMap, and two
   workers use it without locking it.  Both must fail with a locking
   error. */
static int test_unlocked_use_fails( void ) {
   static Worker data1, data2;  /* static to avoid ASan stack-use-after-return
                                   false positive when threads access these */
   static const char *expected = "astCheckLock(UnitMap): The supplied "
      "UnitMap cannot be used since it is not locked for use by the "
      "current thread (programming error).";
   AstUnitMap *map;
   int ok = 1;

   map = astUnitMap( 1, " " );
   astUnlock( map, 1 );

   data1.obj = (AstObject *) map;
   data1.lock = 0;
   data2.obj = (AstObject *) map;
   data2.lock = 0;

   if( !run_workers( &data1, &data2 ) ) {
      ok = 0;
   } else {
      ok = check_lock_error( "unlocked-use-fails", &data1, expected ) && ok;
      ok = check_lock_error( "unlocked-use-fails", &data2, expected ) && ok;
   }

   astLock( map, 0 );
   map = astAnnul( map );
   return ok && astOK;
}

/* A thread can use an Object it has locked.  The main thread gives each
   worker its own copy of a UnitMap, unlocked, and the workers lock their
   copies before using them.  Both must succeed. */
static int test_locked_copies( void ) {
   static Worker data1, data2;  /* static to avoid ASan stack-use-after-return
                                   false positive when threads access these */
   AstUnitMap *map1;
   AstUnitMap *map2;
   int ok = 1;

   map1 = astUnitMap( 1, " " );
   map2 = astCopy( map1 );
   astUnlock( map1, 1 );
   astUnlock( map2, 1 );

   data1.obj = (AstObject *) map1;
   data1.lock = 1;
   data2.obj = (AstObject *) map2;
   data2.lock = 1;

   if( !run_workers( &data1, &data2 ) ) {
      ok = 0;
   } else {
      ok = check_success( "locked-copies", &data1 ) && ok;
      ok = check_success( "locked-copies", &data2 ) && ok;
   }

   astLock( map1, 0 );
   map1 = astAnnul( map1 );
   astLock( map2, 0 );
   map2 = astAnnul( map2 );
   return ok && astOK;
}

int main( void ) {
   int status = SAI__OK;
   int fails = 0;

   astWatch( &status );
   if( pthread_key_create( &worker_key, NULL ) ) {
      printf( "FAIL: could not create thread-specific data key\n" );
      return 1;
   }

   if( !test_unlocked_use_fails() ) {
      fails++;
   }
   if( !test_locked_copies() ) {
      fails++;
   }

   if( !astOK ) {
      printf( "FAIL: main thread status %d\n", astStatus );
      fails++;
   }

   if( fails ) {
      printf( "%d thread test(s) failed\n", fails );
   } else {
      printf( " All thread tests passed\n" );
   }
   return fails ? 1 : 0;
}
