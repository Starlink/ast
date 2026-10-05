/*
 * Tests that AST releases a thread's thread-specific data when the thread
 * exits, and that objects created by a thread remain usable after it
 * exits.
 *
 * AST keeps a block of global data for every thread that calls it,
 * including that thread's copy of every class vtab.  An object refers to
 * the vtab of the thread that created it, and keeps doing so after being
 * unlocked and handed to another thread, so the block must outlive the
 * thread for as long as any such object exists.
 *
 * Leaks are detected by LeakSanitizer, which the CMake build enables for
 * this test when sanitizers are on.  Each case does its AST work in a
 * worker thread named after the case, so a leak report's allocation stack
 * identifies the case that leaked.  Without sanitizers only checks that the
 * code does not crash.
 */
#include "sae_par.h"
#include "ast.h"
#include "ast_err.h"
#include <pthread.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

typedef void *(*Worker)( void * );

/* Run "worker" in a new thread with argument "arg" and wait for it to
   exit.  Returns the AST status the worker ended with. */
static int run_thread( Worker worker, void *arg ) {
   pthread_t thread;
   void *result;

   if( pthread_create( &thread, NULL, worker, arg ) ) {
      printf( "FAIL: could not create thread\n" );
      return AST__INTER;
   }
   if( pthread_join( thread, &result ) ) {
      printf( "FAIL: could not join thread\n" );
      return AST__INTER;
   }
   return (int) (size_t) result;
}

/* Return from a worker thread, passing its AST status to run_thread. */
static void *worker_exit( int *status ) {
   int result = *status;

   if( result != SAI__OK ) {
      astClearStatus;
   }
   return (void *) (size_t) result;
}


/* A thread that uses AST and exits must not leave its thread-specific data
   behind. */
static void *exit_worker( void *arg ) {
   AstFrame *frame;
   int status = SAI__OK;

   (void) arg;
   astWatch( &status );
   frame = astFrame( 2, " " );
   frame = astAnnul( frame );
   return worker_exit( &status );
}

static int test_exit( void ) {
   if( run_thread( exit_worker, NULL ) != SAI__OK ) {
      printf( "FAIL: exit: worker failed\n" );
      return 0;
   }
   return 1;
}


/* Create a Galactic SkyFrame, unlock it and exit, leaving the SkyFrame
   referring to this thread's vtabs. */
static void *create_skyframe_worker( void *arg ) {
   AstSkyFrame **sky = arg;
   int status = SAI__OK;

   astWatch( &status );
   *sky = astSkyFrame( "System=Galactic" );
   astExempt( *sky );
   /* astUnlock leaves the SkyFrame free to associate with a different thread
    * but it does not move its vtab reference until it is locked again by a
    * different thread with its own globals. */
   astUnlock( *sky, 1 );
   return worker_exit( &status );
}


/* Check that a SkyFrame created by create_skyframe_worker is usable. */
static int check_skyframe( AstSkyFrame *sky, const char *test_name ) {
   const char *system = astGetC( sky, "System" );

   if( astOK && strcmp( system, "GALACTIC" ) ) {
      printf( "FAIL: %s: System is %s, expected GALACTIC\n", test_name, system );
      return 0;
   }
   return astOK;
}

/* An object handed to a thread that has no vtab of its own for the
   object's class keeps using the vtab of the thread that created it,
   after that thread has exited. */
static void *handoff_foreign_vtab_worker( void *arg ) {
   AstSkyFrame *sky;
   int status = SAI__OK;
   int *ok = arg;

   astWatch( &status );
   if( run_thread( create_skyframe_worker, &sky ) == SAI__OK ) {
      astLock( sky, 0 );
      *ok = check_skyframe( sky, "handoff-foreign-vtab" );
      sky = astAnnul( sky );
   }
   return worker_exit( &status );
}

/* Test the handoff of an object to a another thread while it still
 * holds a reference to a foreign thread's vtab */
static int test_handoff_foreign_vtab( void ) {
   static int ok = 0;

/* Run in a new thread that is guaranteed to have no SkyFrame vtab. */
   if( run_thread( handoff_foreign_vtab_worker, &ok ) != SAI__OK ) {
      printf( "FAIL: handoff-foreign-vtab: worker failed\n" );
      return 0;
   }
   return ok;
}

/* An object handed to a thread that has its own vtab for the object's
   class switches to that vtab when locked. */
static void *handoff_own_vtab_worker( void *arg ) {
   AstSkyFrame *own;
   AstSkyFrame *sky;
   int status = SAI__OK;
   int *ok = arg;

   astWatch( &status );
   own = astSkyFrame( " " );
   if( run_thread( create_skyframe_worker, &sky ) == SAI__OK ) {
      astLock( sky, 0 );
      *ok = check_skyframe( sky, "handoff-own-vtab" );
      sky = astAnnul( sky );
   }
   own = astAnnul( own );
   return worker_exit( &status );
}


/* Test that an object created on a foreign thread can use the
 * current thread's already initialized global vtab (in this case
 * from an existing SkyFrame already created in that thread) */
static int test_handoff_own_vtab( void ) {
   static int ok = 0;

   if( run_thread( handoff_own_vtab_worker, &ok ) != SAI__OK ) {
      printf( "FAIL: handoff-own-vtab: worker failed\n" );
      return 0;
   }
   return ok;
}

/* A copy of an object that refers to an exited thread's vtab refers to
   the same vtab, and stays usable after the original is annulled. */
static void *copy_worker( void *arg ) {
   AstSkyFrame *copy;
   AstSkyFrame *sky;
   int status = SAI__OK;
   int *ok = arg;

   astWatch( &status );
   if( run_thread( create_skyframe_worker, &sky ) == SAI__OK ) {
      astLock( sky, 0 );
      copy = astCopy( sky );
      sky = astAnnul( sky );
      *ok = check_skyframe( copy, "copy" );
      copy = astAnnul( copy );
   }
   return worker_exit( &status );
}

static int test_copy( void ) {
   static int ok = 0;

/* Run in a new thread that is guaranteed to have no SkyFrame vtab. */
   if( run_thread( copy_worker, &ok ) != SAI__OK ) {
      printf( "FAIL: copy: worker failed\n" );
      return 0;
   }
   return ok;
}


/* With memory caching switched on, astFree keeps small blocks in a cache
   held for the thread rather than freeing them.  Both the blocks cached
   while the thread runs and those freed as its resources are released on
   exit must be returned to the system. */
static void *memory_cache_worker( void *arg ) {
   AstFrame *frame;
   int status = SAI__OK;

   (void) arg;
   astWatch( &status );
   (void) astTune( "MemoryCaching", 1 );

/* Deleting the Frame caches the small blocks it used. */
   frame = astFrame( 2, " " );
   (void) astGetC( frame, "Title" );
   frame = astAnnul( frame );

   return worker_exit( &status );
}

static int test_memory_cache( void ) {
   if( run_thread( memory_cache_worker, NULL ) != SAI__OK ) {
      printf( "FAIL: memory-cache: worker failed\n" );
      return 0;
   }
   return 1;
}


/* Read a FITS-WCS header with an MJD-OBS value.  The FitsChan converts the
   MJD-OBS value to the TDB timescale, using a pair of TimeFrames that the
   FitsChan class creates for the calling thread.  Returns the FrameSet
   read from the header, or NULL.

   This is used as a test case in test_class_globals_freed and in
   test_fitschan_other_thread specifically targeting the case of the FitsChan
   class.
   */
static AstObject *read_fits_header( void ) {
   static const char *cards[] = {
      "NAXIS   = 2",
      "NAXIS1  = 100",
      "NAXIS2  = 100",
      "CTYPE1  = 'RA---TAN'",
      "CTYPE2  = 'DEC--TAN'",
      "CRPIX1  = 50.0",
      "CRPIX2  = 50.0",
      "CRVAL1  = 10.0",
      "CRVAL2  = 20.0",
      "CDELT1  = -0.001",
      "CDELT2  = 0.001",
      "MJD-OBS = 58849.5",  /* This is needed to invoke the TDB conversion */
      NULL
   };
   AstFitsChan *fc;
   AstObject *obj;
   int i;

   fc = astFitsChan( NULL, NULL, " " );
   for( i = 0; cards[ i ]; i++ ) {
      astPutFits( fc, cards[ i ], 0 );
   }
   astClear( fc, "Card" );
   obj = astRead( fc );
   fc = astAnnul( fc );
   return obj;
}

/* Some classes keep resources in their thread-specific data for the life
   of the thread, which their astFree<Class>Globals function releases when
   the thread exits.  This worker makes each such class allocate them.
   (Those of the Error module and the Plot class are not exercised.)

   This test exercises all such classes, including:

   - FitsChan
   - Frame
   - FrameSet
   - KeyMap
   - Plot3D
   - SkyFrame
   - SpecFrame */
static void *class_globals_freed_worker( void *arg ) {
   static const float gbox[ 6 ] = { 0.0f, 0.0f, 0.0f, 1.0f, 1.0f, 1.0f };
   static const double bbox[ 6 ] = { 0.0, 0.0, 0.0, 1.0, 1.0, 1.0 };
   AstFrame *frame3d;
   AstFrameSet *cvt;
   AstKeyMap *km;
   AstPlot3D *plot3d;
   AstSkyFrame *azel;
   AstSkyFrame *icrs;
   AstSpecFrame *spec;
   const char *cval;
   double xin[ 1 ] = { 0.1 };
   double yin[ 1 ] = { 0.2 };
   double xout[ 1 ];
   double yout[ 1 ];
   int status = SAI__OK;

   (void) arg;  /* Not used */
   astWatch( &status );
   astBegin;

/* astFreeSkyFrameGlobals: converting from AzEl needs the local sidereal
   time, which the SkyFrame class computes using a pair of TimeFrames it
   creates for the thread. */
   azel = astSkyFrame( "System=AzEl,Epoch=2020.0" );
   icrs = astSkyFrame( "System=ICRS,Epoch=2020.0" );
   cvt = astConvert( azel, icrs, " " );
   if( cvt ) {
      astTran2( cvt, 1, xin, yin, 1, xout, yout );
   }

/* astFreeObjectGlobals: astGetC returns strings from a buffer held for
   the thread. */
   (void) astGetC( azel, "System" );

/* astFreeFrameGlobals: likewise for astFormat.  The Frame class also
   formats and unformats ObsLon and ObsLat values using a SkyFrame it
   creates for the thread. */
   (void) astFormat( azel, 1, 0.5 );
   astSetC( azel, "ObsLon", "E10:00:00" );
   (void) astGetC( azel, "ObsLon" );

/* astFreeSpecFrameGlobals: the SpecFrame class formats and unformats
   RefRA and RefDec values using a SkyFrame it creates for the thread. */
   spec = astSpecFrame( "RefRA=10:00:00,RefDec=20:00:00" );
   (void) astGetC( spec, "RefRA" );

/* astFreeFitsChanGlobals: the TimeFrames used by read_fits_header. */
   (void) read_fits_header();

/* astFreePlot3DGlobals: changing an attribute of a Plot3D updates its
   component Plots, using a FrameSet that the Plot3D class creates for the
   thread. */
   frame3d = astFrame( 3, " " );
   plot3d = astPlot3D( frame3d, gbox, bbox, " " );
   astSet( plot3d, "Title=Test" );

/* astFreeKeyMapGlobals: likewise for astMapGet0C and astMapKey, which use
   separate buffers. */
   km = astKeyMap( " " );
   astMapPut0I( km, "Key", 1, NULL );
   (void) astMapGet0C( km, "Key", &cval );
   (void) astMapKey( km, 0 );

   astEnd;
   return worker_exit( &status );
}

static int test_class_globals_freed( void ) {
   if( run_thread( class_globals_freed_worker, NULL ) != SAI__OK ) {
      printf( "FAIL: class-globals-freed: worker failed\n" );
      return 0;
   }
   return 1;
}

/* Seconds to allow test_fitschan_other_thread before assuming it is
   waiting for ever. */
#define TIMEOUT 30

static void watchdog( int sig ) {
   static const char msg[] = "FAIL: fitschan-other-thread: timed out, "
                             "probably waiting for an Object locked by "
                             "another thread\n";
   ssize_t n;
   (void) sig;  /* Not used */

   n = write( STDERR_FILENO, msg, sizeof( msg ) - 1 );
   (void) n;
   abort();
}

/* Read a FITS header with an MJD-OBS value, and check that it gives a
   FrameSet. */
static void *fitschan_worker( void *arg ) {
   AstObject *obj;
   int status = SAI__OK;

   (void) arg;
   astWatch( &status );

   obj = read_fits_header();
   if( astOK && !( obj && astIsAFrameSet( obj ) ) ) {
      printf( "FAIL: fitschan-other-thread: no FrameSet read from the "
              "header\n" );
      astSetStatus( AST__INTER );
   }

   if( obj ) {
      obj = astAnnul( obj );
   }

   return worker_exit( &status );
}

/* The TimeFrames that the FitsChan class uses for converting MJD values
   belong to the thread using them, so a thread can read a FITS header
   with an MJD-OBS value even after another thread has used the FitsChan
   class first, and still holds any Objects it created while doing so.
   This only tests anything if the main thread is the first in the process
   to use the FitsChan class, and does not itself convert an MJD value. */
static int test_fitschan_other_thread( void ) {
   AstFitsChan *fc;
   int ok = 1;

   fc = astFitsChan( NULL, NULL, " " );

   alarm( TIMEOUT );
   if( run_thread( fitschan_worker, NULL ) != SAI__OK ) {
      printf( "FAIL: fitschan-other-thread: worker failed\n" );
      ok = 0;
   }
   alarm( 0 );

   fc = astAnnul( fc );
   return ok;
}

int main( void ) {
   int status = SAI__OK;
   int fails = 0;

   astWatch( &status );
   signal( SIGALRM, watchdog );

/* This must be the first use of the FitsChan class in the process, so
   that the main thread is the first to use it. */
   if( !test_fitschan_other_thread() ) {
      fails++;
   }
   if( !test_exit() ) {
      fails++;
   }
   if( !test_handoff_foreign_vtab() ) {
      fails++;
   }
   if( !test_handoff_own_vtab() ) {
      fails++;
   }
   if( !test_copy() ) {
      fails++;
   }
   if( !test_memory_cache() ) {
      fails++;
   }
   if( !test_class_globals_freed() ) {
      fails++;
   }

   if( !astOK ) {
      printf( "FAIL: main thread status %d\n", astStatus );
      fails++;
   }

   if( fails ) {
      printf( "%d thread exit test(s) failed\n", fails );
   } else {
      printf( " All thread exit tests passed\n" );
   }
   return fails ? 1 : 0;
}
