/*
 * Tests that astLock does not deadlock when several threads contend for
 * the same Object.
 *
 * Each Object has two mutexes: a primary one, whose ownership is the lock
 * on the Object, and a secondary one guarding the record of which thread
 * owns it.  A thread that waits for another thread's lock blocks on the
 * primary mutex, then takes the secondary one to record itself as owner.
 * If any thread could block on the primary mutex while holding the
 * secondary one, a waiter taking over the lock and a third thread
 * arriving at that moment could each end up waiting for the mutex the
 * other holds.
 *
 * Each thread uses its own pointer (handle, really) to the shared Object,
 * cloned from the original, since a pointer belongs to the thread that last
 * locked it and cannot be used by another thread until it is unlocked.
 *
 * Also tests that astUnlock, given a pointer that belongs to another
 * thread, reports an error without disturbing that thread's pointer.
 *
 * A deadlock hangs rather than fails, so a watchdog alarm aborts the test
 * if it runs for too long.
 */
#include "sae_par.h"
#include "ast.h"
#include "ast_err.h"
#include <pthread.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

/* Three threads are the minimum that can deadlock; more make it likely
   sooner. */
#define NTHREAD 4
#define NITER 20000

/* Seconds to allow each case before assuming it has deadlocked.  The
   cases take a fraction of a second even in a sanitizer build. */
#define TIMEOUT 10

/* What a worker does with the shared Object. */
enum { WAIT, NOWAIT, UNLOCK_ONLY };

typedef struct Worker {
   AstObject *shared;
   int id;    /* Used in debugging so workers can ID themselves */
   int role;
   int status;
   int nbusy; /* Used in debugging; count of how many times a Worker failed to
                 obtain the object lock without waiting */
} Worker;

static const char *current_case = "";

static void watchdog( int sig ) {
   (void) sig;
   fprintf( stderr, "FAIL: %s timed out, probable deadlock in astLock\n",
            current_case );
   abort();
}

/* An error handler that discards the messages for the lock failures that
   non-waiting workers expect. */
static void discard_error( int status_value, const char *message ) {
   (void) status_value;
   (void) message;
}

/* Repeatedly lock and unlock the shared Object.  A worker that does not
   wait counts the attempts that failed because another thread had the
   Object locked, and stops at any other error. */
static void *lock_worker( void *arg ) {
   Worker *worker = arg;
   int i;
   int status = SAI__OK;
   int wait = worker->role == WAIT;

   astWatch( &status );
   if( !wait ) {
      astSetPutErr( discard_error );
   }

   for( i = 0; i < NITER; i++ ) {
      astLock( worker->shared, wait );
      if( !astOK ) {
         if( wait || astStatus != AST__LCKERR ) {
            break;
         }
         worker->nbusy++;
         astClearStatus;
         continue;
      }
      astUnlock( worker->shared, 1 );
      if( !astOK ) {
         break;
      }
   }

   worker->status = status;
   if( !astOK ) {
      astClearStatus;
   }
   return NULL;
}

/* Repeatedly try to unlock the shared Object, without ever locking it.
   Whenever another thread has the Object locked, astUnlock must fail with
   AST__LCKERR; otherwise it does nothing.  Stops at any other error. */
static void *unlock_worker( void *arg ) {
   Worker *worker = arg;
   int i;
   int status = SAI__OK;

   astWatch( &status );
   astSetPutErr( discard_error );

   for( i = 0; i < NITER; i++ ) {
      astUnlock( worker->shared, 1 );
      if( !astOK ) {
         if( astStatus != AST__LCKERR ) {
            break;
         }
         worker->nbusy++;
         astClearStatus;
      }
   }

   worker->status = status;
   if( !astOK ) {
      astClearStatus;
   }
   return NULL;
}

/* Run NTHREAD workers on "object", each with its own pointer to it. The
   first "nwait" of them lock it, waiting for it when another thread has
   it locked, and the next "nunlock" only try to unlock it. The rest lock
   it without waiting. Annuls the supplied pointer. */
static int run_workers( const char *name, AstObject *object, int nwait,
                        int nunlock ) {
   static const char *role_names[] = { "waiting", "not waiting",
                                       "unlocking" };
   static Worker workers[ NTHREAD ];
   pthread_t threads[ NTHREAD ];
   int i;
   int ok = 1;

   current_case = name;
   alarm( TIMEOUT );

   for( i = 0; i < NTHREAD; i++ ) {
      workers[ i ].id = i;
      workers[ i ].shared = astClone( object );
   }
   object = astAnnul( object );
   for( i = 0; i < NTHREAD; i++ ) {
      astUnlock( workers[ i ].shared, 1 );
   }

   for( i = 0; i < NTHREAD; i++ ) {
      if( i < nwait ) {
         workers[ i ].role = WAIT;
      } else if( i < nwait + nunlock ) {
         workers[ i ].role = UNLOCK_ONLY;
      } else {
         workers[ i ].role = NOWAIT;
      }
      workers[ i ].status = SAI__OK;
      workers[ i ].nbusy = 0;
      if( pthread_create( &threads[ i ], NULL,
                          workers[ i ].role == UNLOCK_ONLY ? unlock_worker
                                                           : lock_worker,
                          &workers[ i ] ) ) {
         printf( "FAIL: %s: could not create thread %d\n", name, i );
         exit( 1 );
      }
   }
   for( i = 0; i < NTHREAD; i++ ) {
      if( pthread_join( threads[ i ], NULL ) ) {
         printf( "FAIL: %s: could not join thread %d\n", name, i );
         exit( 1 );
      }
   }

   alarm( 0 );

   /* Check that all workers exited successfully */
   for( i = 0; i < NTHREAD; i++ ) {
      if( workers[ i ].status != SAI__OK ) {
         printf( "FAIL: %s: thread %d (%s) failed with status %d\n", name,
                 i, role_names[ workers[ i ].role ],
                 workers[ i ].status );
         ok = 0;
      }
   }

   /* The object should not be locked to any worker, so obtaining the lock
    * in this thread in order to annul the object should work */
   for( i = 0; i < NTHREAD; i++ ) {
      astLock( workers[ i ].shared, 0 );
      workers[ i ].shared = astAnnul( workers[ i ].shared );
   }
   return ok && astOK;
}

/* Threads that all wait for the Object. */
static int test_all_wait( void ) {
   return run_workers( "all-wait", (AstObject *) astFrame( 2, " " ),
                       NTHREAD, 0 );
}

/* Threads that wait for the Object, together with threads that do not
   and so must fail immediately, with AST__LCKERR, whenever another thread
   has it locked or is taking the lock over.

   This uses a ZoomMap, which contains no other Objects.  Locking an
   Object that does, such as a Frame and its Axes, locks each in turn, and
   a non-waiting astLock that fails part way through leaves the Objects
   already locked still locked by the calling thread. */
static int test_mixed_wait( void ) {
   return run_workers( "mixed-wait", (AstObject *) astZoomMap( 2, 2.0, " " ),
                       NTHREAD / 2, 0 );
}

/* Threads that wait for the Object, together with threads that try to
   unlock it while other threads have it locked.  An astUnlock that does
   not own the Object must report an error rather than touch the data of
   the thread that does, which changes as that thread locks and unlocks
   the Object. */
static int test_unlock_other_thread( void ) {
   return run_workers( "unlock-other-thread",
                       (AstObject *) astFrame( 2, " " ), NTHREAD / 2,
                       NTHREAD / 2 );
}

int main( void ) {
   int status = SAI__OK;
   int fails = 0;

   astWatch( &status );
   signal( SIGALRM, watchdog );

   if( !test_all_wait() ) {
      fails++;
   }
   if( !test_mixed_wait() ) {
      fails++;
   }
   if( !test_unlock_other_thread() ) {
      fails++;
   }

   if( fails ) {
      printf( "%d lock contention test(s) failed\n", fails );
   } else {
      printf( " All lock contention tests passed\n" );
   }
   return fails ? 1 : 0;
}
