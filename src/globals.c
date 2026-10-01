#if defined( THREAD_SAFE )

#define astCLASS

#include "globals.h"
#include "error.h"
#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>

/* Configuration results. */
/* ---------------------- */
#if HAVE_CONFIG_H
#include <config.h>
#endif

/* Select the appropriate memory management functions. These will be the
   system's malloc, free and realloc unless AST was configured with the
   "--with-starmem" option, in which case they will be the starmem
   malloc, free and realloc. */
#ifdef HAVE_STAR_MEM_H
#  include <star/mem.h>
#  define MALLOC starMalloc
#  define FREE starFree
#  define REALLOC starRealloc
#else
#  define MALLOC malloc
#  define FREE free
#  define REALLOC realloc
#endif

/* Module variables */
/* ================ */

/* A count of the number of thread-specific data structures created so
   far. Create a mutex to serialise access to this static variable. */
static int nthread = 0;
static pthread_mutex_t nthread_mutex = PTHREAD_MUTEX_INITIALIZER;

/* External variables visible throughout AST */
/* ========================================= */

/* Set a flag indicating that the thread-specific data key has not yet
   been created. */
pthread_once_t starlink_ast_globals_initialised = PTHREAD_ONCE_INIT;

/* Declare the pthreads key that will be associated with the thread-specific
   data for each thread. */
pthread_key_t starlink_ast_globals_key;

/* Declare the pthreads key that will be associated with the thread-specific
   status value for each thread. */
pthread_key_t starlink_ast_status_key;


/* Function prototypes: */
/* ==================== */

static void FreeGlobals( AstGlobals * );
static void ThreadExit( void * );


/* Function definitions: */
/* ===================== */


void astGlobalsCreateKey_( void ) {
/*
*+
*  Name:
*     astGlobalsCreateKey_

*  Purpose:
*     Create the thread specific data key used for accessing global data.

*  Type:
*     Protected function.

*  Synopsis:
*     #include "globals.h"
*     astGlobalsCreateKey_()

*  Description:
*     This function creates the thread-specific data key. It is called
*     once only by the pthread_once function, which is invoked via the
*     astGET_GLOBALS(this) macro by each AST function that requires access to
*     global data.

*  Returned Value:
*     Zero for success.

*-
*/

/* Create the key used to access thread-specific global data values,
   arranging for the data to be released when each thread exits. Report
   an error if it fails. */
   if( pthread_key_create( &starlink_ast_globals_key, ThreadExit ) ) {
      fprintf( stderr, "ast: Failed to create Thread-Specific Data key" );

/* If succesful, create the key used to access the thread-specific status
   value. Report an error if it fails. */
   } else if( pthread_key_create( &starlink_ast_status_key, NULL ) ) {
      fprintf( stderr, "ast: Failed to create Thread-Specific Status key" );

   }

}

AstGlobals *astGlobalsInit_( void ) {
/*
*+
*  Name:
*     astGlobalsInit

*  Purpose:
*     Create and initialise a structure holding thread-specific global
*     data values.

*  Type:
*     Protected function.

*  Synopsis:
*     #include "globals.h"
*     AstGlobals *astGlobalsInit;

*  Description:
*     This function allocates memory to hold thread-specific global data
*     for use throughout AST, and initialises it.

*  Returned Value:
*     Pointer to the structure holding global data values for the
*     currently executing thread.

*-
*/

/* Local Variables: */
   AstGlobals *globals;
   AstStatusBlock *status;

/* Allocate memory to hold the global data values for the currently
   executing thread. Use malloc rather than astMalloc (the AST memory
   module uses global data managed by this module and so using astMalloc
   could put us into an infinite loop). */
   globals = MALLOC( sizeof( AstGlobals ) );

   if ( !globals ){
      fprintf( stderr, "ast: Failed to allocate memory to hold AST "
               "global data values" );

/* Initialise the global data values. */
   } else {

/* Each thread has a unique integer identifier. */
      pthread_mutex_lock( &nthread_mutex );
      globals->thread_identifier = nthread++;
      pthread_mutex_unlock( &nthread_mutex );

/* The owning thread holds a reference until it exits. */
      globals->nref = 1;
      pthread_mutex_init( &globals->ref_mutex, NULL );
      globals->status_block = NULL;

#define INIT(class) astInit##class##Globals_( &(globals->class) );
      INIT( Error );
      INIT( Memory );
      INIT( Object );
      INIT( YamlChan );
      INIT( Axis );
      INIT( Mapping );
      INIT( Frame );
      INIT( Channel );
      INIT( CmpMap );
      INIT( KeyMap );
      INIT( FitsChan );
      INIT( FitsTable );
      INIT( CmpFrame );
      INIT( DSBSpecFrame );
      INIT( FrameSet );
      INIT( LutMap );
      INIT( MathMap );
      INIT( PcdMap );
      INIT( PointSet );
      INIT( SkyAxis );
      INIT( SkyFrame );
      INIT( SlaMap );
      INIT( SpecFrame );
      INIT( SphMap );
      INIT( TimeFrame );
      INIT( WcsMap );
      INIT( ZoomMap );
      INIT( FluxFrame );
      INIT( SpecFluxFrame );
      INIT( GrismMap );
      INIT( IntraMap );
      INIT( Plot );
      INIT( Plot3D );
      INIT( Region );
      INIT( Xml );
      INIT( XmlChan );
      INIT( Box );
      INIT( Circle );
      INIT( CmpRegion );
      INIT( DssMap );
      INIT( Ellipse );
      INIT( Interval );
      INIT( MatrixMap );
      INIT( Moc );
      INIT( MocChan );
      INIT( NormMap );
      INIT( NullRegion );
      INIT( PermMap );
      INIT( PointList );
      INIT( PolyMap );
      INIT( ChebyMap );
      INIT( Polygon );
      INIT( Prism );
      INIT( RateMap );
      INIT( SelectorMap );
      INIT( ShiftMap );
      INIT( SpecMap );
      INIT( Stc );
      INIT( SplineMap );
      INIT( StcCatalogEntryLocation );
      INIT( StcObsDataLocation );
      INIT( SwitchMap );
      INIT( Table );
      INIT( TimeMap );
      INIT( TranMap );
      INIT( UnitMap );
      INIT( UnitNormMap );
      INIT( WinMap );
      INIT( StcResourceProfile );
      INIT( StcSearchLocation );
      INIT( StcsChan );
      INIT( XphMap );
#undef INIT

/* Save the pointer as the value of the starlink_ast_globals_key
   thread-specific data key. */
      if( pthread_setspecific( starlink_ast_globals_key, globals ) ) {
         fprintf( stderr, "ast: Failed to store Thread-Specific Data pointer." );

/* We also take this opportunity to allocate and initialise the
   thread-specific status value. */
      } else {
         status = MALLOC( sizeof( AstStatusBlock ) );
         globals->status_block = status;
         if( status ) {
            status->internal_status = 0;
            status->status_ptr = &( status->internal_status );

/* If succesful, store the pointer to this memory as the value of the
   status key for the currently executing thread. Report an error if
   this fails. */
            if( pthread_setspecific( starlink_ast_status_key, status ) ) {
               fprintf( stderr, "ast: Failed to store Thread-Specific Status pointer." );
            }

         } else {
            fprintf( stderr, "ast: Failed to allocate memory for Thread-Specific Status pointer." );
         }
      }
   }

/* Return a pointer to the data structure holding the global data values. */
   return globals;
}

void astGlobalsRef_( AstGlobals *globals ) {
/*
*+
*  Name:
*     astGlobalsRef

*  Purpose:
*     Add a reference to a structure holding thread-specific global data.

*  Type:
*     Protected function.

*  Synopsis:
*     #include "globals.h"
*     void astGlobalsRef( AstGlobals *globals )

*  Description:
*     This function increments the reference count of a structure
*     holding thread-specific global data, preventing it from being freed
*     until a matching call to astGlobalsUnref. Each Object holds a
*     reference to the structure containing its vtab. May be called by
*     any thread.

*  Parameters:
*     globals
*        Pointer to the structure.
*-
*/
   if( !globals ) {
      return;
   }

   pthread_mutex_lock( &globals->ref_mutex );
   globals->nref++;
   pthread_mutex_unlock( &globals->ref_mutex );
}

void astGlobalsUnref_( AstGlobals *globals ) {
/*
*+
*  Name:
*     astGlobalsUnref

*  Purpose:
*     Remove a reference to a structure holding thread-specific global
*     data.

*  Type:
*     Protected function.

*  Synopsis:
*     #include "globals.h"
*     void astGlobalsUnref( AstGlobals *globals )

*  Description:
*     This function decrements the reference count of a structure holding
*     thread-specific global data, freeing the structure if no references
*     remain. This happens once the owning thread has exited and every
*     Object using one of the structure's vtabs has been deleted, so the
*     structure may be freed by any thread. May be called by any thread.

*  Parameters:
*     globals
*        Pointer to the structure.
*-
*/

/* Local Variables: */
   int nref;

   if( !globals ) {
      return;
   }

   pthread_mutex_lock( &globals->ref_mutex );
   nref = --globals->nref;
   pthread_mutex_unlock( &globals->ref_mutex );

   if( nref == 0 ) {
      FreeGlobals( globals );
   }
}

static void FreeGlobals( AstGlobals *globals ) {
/*
*  Name:
*     FreeGlobals

*  Purpose:
*     Free a structure holding thread-specific global data.

*  Type:
*     Private function.

*  Synopsis:
*     void FreeGlobals( AstGlobals *globals )

*  Description:
*     This function frees a structure holding thread-specific global
*     data, together with the memory used by the vtabs it contains. The
*     per-thread resources referred to by the structure must already have
*     been released by ThreadExit.

*  Parameters:
*     globals
*        Pointer to the structure.
*/
   astFreeObjectVtabs_( &globals->Object, astGetStatusPtr );
   pthread_mutex_destroy( &globals->ref_mutex );
   FREE( globals );
}

static void ThreadExit( void *data ) {
/*
*  Name:
*     ThreadExit

*  Purpose:
*     Release a thread's thread-specific global data when it exits.

*  Type:
*     Private function.

*  Synopsis:
*     void ThreadExit( void *data )

*  Description:
*     This function is the destructor for the starlink_ast_globals_key
*     thread-specific data key, invoked by pthreads when a thread that has
*     used AST exits. It releases the resources used only by the thread,
*     then drops the thread's reference to its global data structure. The
*     structure itself survives until no Object refers to any of the
*     vtabs it contains, since Objects created by the thread may have been
*     handed to other threads.

*  Parameters:
*     data
*        Pointer to the AstGlobals structure for the exiting thread.
*/

/* Local Variables: */
   AstGlobals *globals;
   AstStatusBlock *status_block;
   AstStatusBlock temp_status_block;
   int *status;

   globals = (AstGlobals *) data;
   status_block = globals->status_block;

/* astGlobalsInit_ leaves the thread without a status block if it could
   not allocate one. Use a temporary one in that case, since the AST
   functions used below need a status variable. */
   if( !status_block ) {
      status_block = &temp_status_block;
   }

/* pthreads clears the thread-specific data before invoking destructors.
   Re-instate it so that the AST functions used below find this thread's
   data rather than creating new data for the thread. */
   pthread_setspecific( starlink_ast_globals_key, globals );
   pthread_setspecific( starlink_ast_status_key, status_block );

/* Any status variable registered with astWatch may have gone out of scope
   by now, so revert to the internal status variable. */
   status_block->internal_status = 0;
   status_block->status_ptr = &( status_block->internal_status );
   status = status_block->status_ptr;

/* Empty the memory cache, and disable it so that memory freed from here
   on is returned to the system. */
   astMemCaching_( 0, status );

/* Release the per-thread resources held by each class that has any.
   Classes whose resources are Objects come first, since deleting an
   Object may use the resources of the classes below. */
#define FREE_GLOBALS(class) astFree##class##Globals_( &(globals->class), status );
   FREE_GLOBALS( FitsChan );
   FREE_GLOBALS( Plot3D );
   FREE_GLOBALS( SpecFrame );
   FREE_GLOBALS( Frame );
   FREE_GLOBALS( SkyFrame );
   FREE_GLOBALS( Plot );
   FREE_GLOBALS( KeyMap );
   FREE_GLOBALS( Object );
   FREE_GLOBALS( Error );
#undef FREE_GLOBALS

/* Drop the thread's own reference. */
   astGlobalsUnref_( globals );

   pthread_setspecific( starlink_ast_globals_key, NULL );
   pthread_setspecific( starlink_ast_status_key, NULL );
   if( status_block != &temp_status_block ) {
      FREE( status_block );
   }
}

#endif

