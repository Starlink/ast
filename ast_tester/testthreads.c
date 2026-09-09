#include "sae_par.h"
#include "ast.h"
#include "ast_err.h"
#include <pthread.h>
#include <math.h>

void *worker( void *ptr );

typedef struct MyData {
   AstObject *obj;
   int lock;
   int status;
   double expected;
} MyData;

int main( void ){
   pthread_t thread1, thread2;
   static MyData data1, data2;  /* static to avoid ASan stack-use-after-return
                                   false positive when threads access these */
   int status = SAI__OK;

   astWatch( &status );


/* Create a UnitMap, unlock it, and then create two threads and tell them
   to use the UnitMap to transform a point. This should fail in the threads
   do not lock the UnitMap pointer */
   data1.obj = (AstObject *) astUnitMap( 1, " " );
   data2.obj = data1.obj;
   astUnlock( data1.obj, 1 );
   data1.lock = 0;
   data2.lock = 0;

   if( pthread_create( &thread1, NULL, worker, &data1 ) ) {
      astError( AST__INTER, "Error creating thread1");
   } else if( pthread_create( &thread2, NULL, worker, &data2 ) ) {
      astError( AST__INTER, "Error creating thread2");
   }

/* Wait for both threads to finish. */
   if( astOK ) {
      if( pthread_join( thread1, NULL) ) {
         astError( AST__INTER, "Error joining thread1\n");
      } else if( pthread_join( thread2, NULL) ) {
         astError( AST__INTER, "Error joining thread2\n");
      }
   }

/* Both worker threads should have failed with a locking error, because
   they used the UnitMap without locking it for their own use.  Each thread
   reports its own status back through the shared structure (this build has
   no Starlink EMS error stack to carry status between threads). */
   if( astOK && ( data1.status != AST__LCKERR ||
                  data2.status != AST__LCKERR ) ) {
      astError( AST__INTER, "Error 1 (thread status %d and %d)",
                data1.status, data2.status );
   }


/* Do the same again but this time, send unlocked independent copies of the
   UnitMap to the workers and tell the workers to lock the pointer before
   using it. This should work. */
   astLock( data1.obj, 0 );
   data2.obj = astCopy( data1.obj );
   astUnlock( data1.obj, 1 );
   astUnlock( data2.obj, 1 );
   data1.lock = 1;
   data2.lock = 1;

   if( pthread_create( &thread1, NULL, worker, &data1 ) ) {
      astError( AST__INTER, "Error creating thread1");
   } else if( pthread_create( &thread2, NULL, worker, &data2 ) ) {
      astError( AST__INTER, "Error creating thread2");
   }

   if( astOK ) {
      if( pthread_join( thread1, NULL) ) {
         astError( AST__INTER, "Error joining thread1\n");
      } else if( pthread_join( thread2, NULL) ) {
         astError( AST__INTER, "Error joining thread2\n");
      }
   }

/* Both worker threads should have succeeded this time, because they locked
   the UnitMap before using it. */
   if( astOK && ( data1.status != SAI__OK || data2.status != SAI__OK ) ) {
      astError( AST__INTER, "Error 2 (thread status %d and %d)",
                data1.status, data2.status );
   }









/* Warm both ChebyMap inverse caches, then transfer the inverted Mapping
   to another thread. All cached child Mappings must transfer their locks. */
   if( astOK ) {
      double coeffs[] = { 1, 1, 1, .125, 1, 2 };
      double lo = 0, hi = 10, target = 0, value;
      AstChebyMap *cm = astChebyMap( 1, 1, 2, coeffs, 0, NULL,
                                     &lo, &hi, NULL, NULL,
                                     "NiterInverse=20,TolInverse=1e-12" );
      astTran1( cm, 1, &target, 0, &value );
      astInvert( cm );
      data1.obj = (AstObject *) cm;
      data1.expected = 5*(1 + .25/(1+sqrt(1.125)));
      astUnlock( cm, 1 );
      if( pthread_create( &thread1, NULL, worker, &data1 ) ) {
         astError( AST__INTER, "Error creating ChebyMap worker" );
      } else if( pthread_join( thread1, NULL ) ) {
         astError( AST__INTER, "Error joining ChebyMap worker" );
      }
      if( astOK && data1.status != SAI__OK ) {
         astError( AST__INTER, "ChebyMap cache transfer failed: %d", data1.status );
      }
      astLock( cm, 0 );
      cm = astAnnul( cm );
   }

   if( astOK ) {
      printf(" All thread tests passed\n");
   } else {
      printf("Thread tests failed\n");
   }
   return astOK ? 0 : 1;

}

void *worker( void *ptr ) {
   double xin, xout;
   MyData *data = (MyData *) ptr;
   int status = SAI__OK;

/* AST maintains a separate status value for each thread.  Watch a
   thread-local variable here so the outcome of the calls below can be
   reported back to the main thread through the shared structure. */
   astWatch( &status );

   if( data->lock ) astLock( data->obj, 0 );

   xin = 0;
   astTran1( data->obj, 1, &xin, 1, &xout );
   if( astOK && (xout == AST__BAD || !isfinite(xout) ||
                 fabs(xout-data->expected) > 1e-10) ) {
      astError( AST__INTER, "Incorrect worker transformation: %.17g", xout );
   }

   if( data->lock ) astUnlock( data->obj, 1 );

   data->status = status;
   if( !astOK ) astClearStatus;
   return NULL;
}
