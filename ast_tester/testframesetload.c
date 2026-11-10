/*
*  Name:
*     testframesetload

*  Purpose:
*     Check that astLoadFrameSet validates the node graph it reads.

*  Description:
*     A FrameSet dump describes a graph: each Frame sits at a node, each
*     node but the first is derived from another node, Frames may mirror
*     another Frame's variants, and the Base and Current attributes are
*     Frame indices.  astLoadFrameSet reads all of these as plain integers
*     with no range check, and they are later used as indices into the
*     Frame and node arrays (for example by astGetFrame and astGetMapping).
*     An out-of-range value therefore has to be rejected at load time, or
*     it becomes an out-of-bounds read or write when the FrameSet is used.

*     This test builds a valid two-Frame FrameSet, confirms it loads and
*     remains usable, and then corrupts individual graph values to confirm
*     each is rejected.  A cycle in the variants chain is in range, so it
*     loads, but must then be caught when the chain is walked rather than
*     looping for ever.

*  Notes:
*     - This is a native C test with no Fortran original.
*     - Public C API only.  Errors expected from a deliberately invalid
*       dump are caught through astOK and cleared with astClearStatus so the
*       remaining cases can run; failures of the test itself are counted and
*       reported at the end.
*/

#include <stdio.h>
#include <string.h>
#include "sae_par.h"
#include "ast.h"

/* Return a newly allocated copy of "src" with the first occurrence of
   "find" replaced by "repl", or NULL if "find" does not occur. */
static char *dup_replace( const char *src, const char *find,
                          const char *repl ) {
   const char *hit = strstr( src, find );
   char *out = NULL;
   if ( hit ) {
      size_t pre = (size_t) ( hit - src );
      size_t flen = strlen( find );
      size_t rlen = strlen( repl );
      out = astMalloc( strlen( src ) - flen + rlen + 1 );
      if ( out ) {
         memcpy( out, src, pre );
         memcpy( out + pre, repl, rlen );
         strcpy( out + pre + rlen, hit + flen );
      }
   }
   return out;
}

/* Return a newly allocated copy of "src" with "text" inserted immediately
   after the first occurrence of "marker", or NULL if "marker" is absent. */
static char *dup_insert( const char *src, const char *marker,
                         const char *text ) {
   const char *hit = strstr( src, marker );
   char *out = NULL;
   if ( hit ) {
      size_t pos = (size_t) ( hit - src ) + strlen( marker );
      size_t tlen = strlen( text );
      out = astMalloc( strlen( src ) + tlen + 1 );
      if ( out ) {
         memcpy( out, src, pos );
         memcpy( out + pos, text, tlen );
         strcpy( out + pos + tlen, src + pos );
      }
   }
   return out;
}

/* Attempt to load a dump that should be rejected.  Returns 0 if it was
   rejected cleanly, 1 otherwise.  Frees "dump" and clears any status. */
static int check_rejected( const char *what, char *dump ) {
   int result = 0;
   if ( !dump ) {
      fprintf( stderr, "testframesetload: could not build the %s dump\n",
               what );
      return 1;
   }
   {
      AstObject *obj = astFromString( dump );
      int bad_status = !astOK;
      if ( obj ) obj = astAnnul( obj );
      if ( !astOK ) astClearStatus;
      if ( obj != NULL || !bad_status ) {
         fprintf( stderr, "testframesetload: %s was accepted but should "
                  "have been rejected\n", what );
         result = 1;
      }
   }
   dump = astFree( dump );
   return result;
}

int main( void ) {
   int status = SAI__OK;
   int nfail = 0;
   char *dump;
   astWatch( &status );

/* Build a valid two-Frame FrameSet (two nodes, one link) and dump it. */
   {
      AstFrame *f1 = astFrame( 2, "Domain=BASE" );
      AstFrameSet *fs = astFrameSet( f1, " " );
      AstFrame *f2 = astFrame( 2, "Domain=CURR" );
      AstMapping *um = (AstMapping *) astUnitMap( 2, " " );
      astAddFrame( fs, AST__BASE, um, f2 );
      dump = astToString( fs );
      fs = astAnnul( fs );
   }

/* The valid FrameSet must load, round-trip unchanged and be usable through
   the paths that index the Frame and node arrays. */
   {
      AstObject *obj = astFromString( dump );
      if ( astOK && obj && !strcmp( astGetC( obj, "Class" ), "FrameSet" ) ) {
         char *dump2 = astToString( obj );
         AstFrame *fb, *fc;
         AstMapping *map;
         if ( astOK && strcmp( dump, dump2 ) != 0 ) {
            fprintf( stderr, "testframesetload: valid FrameSet did not "
                     "round-trip unchanged\n" );
            nfail++;
         }
         fb = astGetFrame( (AstFrameSet *) obj, AST__BASE );
         fc = astGetFrame( (AstFrameSet *) obj, AST__CURRENT );
         map = astGetMapping( (AstFrameSet *) obj, AST__BASE, AST__CURRENT );
         if ( !astOK || !fb || !fc || !map ) {
            fprintf( stderr, "testframesetload: valid FrameSet was not "
                     "usable after loading\n" );
            nfail++;
         }
         if ( fb ) fb = astAnnul( fb );
         if ( fc ) fc = astAnnul( fc );
         if ( map ) map = astAnnul( map );
         dump2 = astFree( dump2 );
      } else {
         fprintf( stderr, "testframesetload: valid FrameSet failed to "
                  "load\n" );
         nfail++;
      }
      if ( obj ) obj = astAnnul( obj );
      if ( !astOK ) astClearStatus;
   }

/* Each out-of-range graph value must be rejected at load time. The keyword
   case does not matter to the loader; the Dump writes capitalised forms. */
   nfail += check_rejected( "out-of-range Current",
                            dup_replace( dump, "Currnt = 2", "Currnt = 99" ) );
   nfail += check_rejected( "out-of-range Base",
                            dup_insert( dump, "Nframe = 2\n", " Base = 99\n" ) );
   nfail += check_rejected( "out-of-range node",
                            dup_insert( dump, "Nframe = 2\n", " Nod1 = 99\n" ) );
   nfail += check_rejected( "out-of-range link",
                            dup_replace( dump, "Lnk2 = 1", "Lnk2 = 99" ) );
   nfail += check_rejected( "out-of-range variants",
                            dup_insert( dump, "Nframe = 2\n", " Vfr1 = 99\n" ) );

/* A cycle in the variants chain is in range, so it loads, but walking the
   chain (here through the Variant attribute) must report corruption rather
   than loop for ever. */
   {
      char *cyc = dup_insert( dump, "Nframe = 2\n", " Vfr1 = 2\n Vfr2 = 1\n" );
      if ( cyc ) {
         AstObject *obj = astFromString( cyc );
         if ( !astOK || !obj ) {
            fprintf( stderr, "testframesetload: cyclic-variants FrameSet "
                     "failed to load\n" );
            nfail++;
            if ( !astOK ) astClearStatus;
         } else {
            (void) astGetC( obj, "Variant" );
            if ( astOK ) {
               fprintf( stderr, "testframesetload: cyclic variants chain "
                        "was not detected\n" );
               nfail++;
            } else {
               astClearStatus;
            }
         }
         if ( obj ) obj = astAnnul( obj );
         cyc = astFree( cyc );
      } else {
         fprintf( stderr, "testframesetload: could not build the cyclic "
                  "variants dump\n" );
         nfail++;
      }
   }

   dump = astFree( dump );

   if ( !astOK ) astClearStatus;

   if ( nfail == 0 ) {
      printf( " All FrameSet load tests passed\n" );
      return 0;
   }
   printf( "FrameSet load tests failed (%d)\n", nfail );
   return 1;
}
