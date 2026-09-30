/*
 *  Test the UseSIMD attribute shared by the PolyMap, SphMap and MatrixMap
 *  classes (and others to follow).  This attribute selects between the
 *  SIMD-vectorised and scalar transform code paths at run time.
 *
 *  For each class we exercise the attribute accessors (set, get, test, clear)
 *  and run a transform with UseSIMD=0 to drive the scalar fall-back path,
 *  checking it produces correct results.  We deliberately do not assert any
 *  numerical difference between the SIMD and scalar paths: by design they
 *  agree to within a few ULP, so such an assertion would be both meaningless
 *  and fragile.
 */
#include "ast.h"
#include <stdio.h>
#include <math.h>

static void stopit( int i, int *status ) {
   if( *status != 0 ) return;
   printf( "Error %d\n", i );
   *status = 1;
}

/* Exercise the UseSIMD attribute accessors on the supplied Mapping.  Error
   codes are offset by "base" so failures identify the class under test. */
static void check_usesimd( AstMapping *map, int base, int *status ) {
   int have_simd;

   if( !astOK )
      return;

/* The value reported when the attribute is unset reflects whether SIMD
   support was compiled in: 1 if so, 0 otherwise.  Read it first so the
   remaining checks work in both kinds of build. */
   have_simd = astGetI( map, "UseSIMD" );
   if( astTest( map, "UseSIMD" ) ) {
      stopit( base + 1, status ); /* LCOV_EXCL_LINE */
   }

/* Disabling is always permitted and persists. */
   astSetI( map, "UseSIMD", 0 );
   if( astGetI( map, "UseSIMD" ) != 0 ) {
      stopit( base + 2, status ); /* LCOV_EXCL_LINE */
   }
   if( !astTest( map, "UseSIMD" ) ) {
      stopit( base + 3, status ); /* LCOV_EXCL_LINE */
   }

/* In a SIMD build, any non-zero value normalises to 1.  In a non-SIMD build,
   enabling SIMD is not possible and should report an error, which we check
   for and then clear so the remaining tests can continue. */
   if( have_simd ) {
      astSetI( map, "UseSIMD", 5 );
      if( astGetI( map, "UseSIMD" ) != 1 ) {
         stopit( base + 4, status ); /* LCOV_EXCL_LINE */
      }
   } else {
      astSetI( map, "UseSIMD", 1 );
      if( astOK ) {
         stopit( base + 4, status ); /* LCOV_EXCL_LINE */
      } else {
         astClearStatus;
      }
   }

/* Clearing restores the compiled-in default. */
   astClear( map, "UseSIMD" );
   if( astTest( map, "UseSIMD" ) ) {
      stopit( base + 5, status ); /* LCOV_EXCL_LINE */
   }
   if( astGetI( map, "UseSIMD" ) != have_simd ) {
      stopit( base + 6, status ); /* LCOV_EXCL_LINE */
   }
}

/* Compare a transform result against an expected value, treating AST__BAD as
   a distinct value that must match exactly. */
static int approx( double got, double want ) {
   if( got == AST__BAD || want == AST__BAD ) {
      return got == want;
   }
   return fabs( got - want ) <= 1.0e-10;
}

/* Build a 2-in/2-out MatrixMap from the given form and matrix, transform the
   "np" supplied points through it both with SIMD enabled (pass 0) and disabled
   (pass 1), and check the results against the expected outputs (AST__BAD
   permitted).  This exercises the matching branches in both TransformLoopSIMD
   and TransformLoopScalar.  Error codes are offset by "base".
   If compiled without SIMD support both passes are identical/redundant. */
static void check_tran2( int form, const double *matrix, int np,
                         const double *xin, const double *yin,
                         const double *xexp, const double *yexp,
                         int base, int *status ) {
   double xout[ 8 ];
   double yout[ 8 ];
   AstMatrixMap *mm;
   int pass;
   int i;

   if( !astOK )
      return;

   for( pass = 0; pass < 2; pass++ ) {
      mm = astMatrixMap( 2, 2, form, matrix, " " );
      if( pass == 1 ) {
         astSetI( mm, "UseSIMD", 0 );
      }
      astTran2( mm, np, xin, yin, 1, xout, yout );
      for( i = 0; i < np; i++ ) {
         if( !approx( xout[ i ], xexp[ i ] ) ) {
            stopit( base + 10 * pass + i, status ); /* LCOV_EXCL_LINE */
         }
         if( !approx( yout[ i ], yexp[ i ] ) ) {
            stopit( base + 10 * pass + 100 + i, status ); /* LCOV_EXCL_LINE */
         }
      }
      mm = astAnnul( mm );
   }
}

/* A simple linear 2-in/2-out forward transformation: out0 = 2*in0,
   out1 = 3*in1 (no inverse). */
static void TestPolyMap( int *status ) {
   double coeff[] = { 2.0, 1, 1, 0,
                      3.0, 2, 0, 1 };
   double xin[ 16 ];
   double yin[ 16 ];
   double xout[ 16 ];
   double yout[ 16 ];
   double tol = 1.0e-10;
   AstPolyMap *pm;
   int idx;

   if( !astOK )
      return;

   pm = astPolyMap( 2, 2, 2, coeff, 0, NULL, " " );
   check_usesimd( (AstMapping *) pm, 1000, status );

/* Force the scalar path and check the transform is still correct. */
   astSetI( pm, "UseSIMD", 0 );

   for( idx = 0; idx < 16; idx++ ) {
      xin[ idx ] = 0.1 * idx - 0.5;
      yin[ idx ] = 0.2 * idx + 1.0;
   }

   astTran2( pm, 16, xin, yin, 1, xout, yout );

   for( idx = 0; idx < 16; idx++ ) {
      if( fabs( xout[ idx ] - 2.0 * xin[ idx ] ) > tol ) {
         stopit( 1100 + idx, status ); /* LCOV_EXCL_LINE */
      }
      if( fabs( yout[ idx ] - 3.0 * yin[ idx ] ) > tol ) {
         stopit( 1200 + idx, status ); /* LCOV_EXCL_LINE */
      }
   }
}

/* A 2x2 diagonal scaling: out0 = 2*in0, out1 = 0.5*in1. */
static void TestMatrixMap( int *status ) {
   double diag[] = { 2.0, 0.5 };
   double xin[ 16 ];
   double yin[ 16 ];
   double xout[ 16 ];
   double yout[ 16 ];
   double tol = 1.0e-10;
   AstMatrixMap *mm;
   int idx;

   if( !astOK )
      return;

   mm = astMatrixMap( 2, 2, 1, diag, " " );
   check_usesimd( (AstMapping *) mm, 2000, status );

   astSetI( mm, "UseSIMD", 0 );

   for( idx = 0; idx < 16; idx++ ) {
      xin[ idx ] = 0.1 * idx - 0.5;
      yin[ idx ] = 0.2 * idx + 1.0;
   }

   astTran2( mm, 16, xin, yin, 1, xout, yout );

   for( idx = 0; idx < 16; idx++ ) {
      if( fabs( xout[ idx ] - 2.0 * xin[ idx ] ) > tol ) {
         stopit( 2100 + idx, status ); /* LCOV_EXCL_LINE */
      }
      if( fabs( yout[ idx ] - 0.5 * yin[ idx ] ) > tol ) {
         stopit( 2200 + idx, status ); /* LCOV_EXCL_LINE */
      }
   }
}

/* Exercise the FULL, DIAGONAL and UNIT MatrixMap forms together with the
   AST__BAD-handling and non-square branches of the transform loops, on both
   the SIMD and scalar paths. */
static void TestMatrixForms( int *status ) {

/* FULL matrix [[2,1],[0,3]]
   This covers the zero-element skip and the bad-input propagation: a bad x
   feeds out0 (coefficient 2) but not out1 (coefficient 0), while a bad y feeds
   both. */
   double full[] = { 2.0, 1.0, 0.0, 3.0 };
   double full_x[] = { 1.0, AST__BAD, 1.0, 3.0 };
   double full_y[] = { 2.0, 2.0, AST__BAD, -1.0 };
   double full_x_exp[] = { 4.0, AST__BAD, AST__BAD, 5.0 };
   double full_y_exp[] = { 6.0, 6.0, AST__BAD, -3.0 };

/* FULL matrix whose second output row is entirely zero, covering the
   all-zero-row branch (output initialised to 0 with no contributing term). */
   double zrow[] = { 2.0, 1.0, 0.0, 0.0 };
   double zrow_x[] = { 1.0, 3.0 };
   double zrow_y[] = { 2.0, -1.0 };
   double zrow_x_exp[] = { 4.0, 5.0 };
   double zrow_y_exp[] = { 0.0, 0.0 };

/* FULL matrix containing a bad element, forcing the SIMD path to fall back to
   the scalar loop, which propagates AST__BAD through the affected output. */
   double badm[] = { 2.0, 1.0, AST__BAD, 3.0 };
   double badm_x[] = { 1.0 };
   double badm_y[] = { 2.0 };
   double badm_x_exp[] = { 4.0 };
   double badm_y_exp[] = { AST__BAD };

/* DIAGONAL matrix with a bad diagonal term: y output is all AST__BAD. */
   double dbad[] = { 2.0, AST__BAD };
   double dbad_x[] = { 1.0, 3.0 };
   double dbad_y[] = { 2.0, -1.0 };
   double dbad_x_exp[] = { 2.0, 6.0 };
   double dbad_y_exp[] = { AST__BAD, AST__BAD };

/* UNIT matrix: output is a copy of the input. */
   double unit_x[] = { 1.0, 3.0 };
   double unit_y[] = { 2.0, -1.0 };

/* Non-square diagonal (2 inputs, 3 outputs): the extra output row is filled
   with zeros (the nax < ncoord_out branch). */
   double diag3[] = { 2.0, 3.0 };
   double in3[] = { 1.0, 3.0, 2.0, -1.0 };  /* 2 coords x 2 points */
   double out3[ 6 ];  /* 3 coords x 2 points */
   AstMatrixMap *mm;
   int pass;

   if( !astOK )
      return;

   check_tran2( 0, full, 4, full_x, full_y, full_x_exp, full_y_exp, 2300,
                status );
   check_tran2( 0, zrow, 2, zrow_x, zrow_y, zrow_x_exp, zrow_y_exp, 2400,
                status );
   check_tran2( 0, badm, 1, badm_x, badm_y, badm_x_exp, badm_y_exp, 2500,
                status );
   check_tran2( 1, dbad, 2, dbad_x, dbad_y, dbad_x_exp, dbad_y_exp, 2600,
                status );
   check_tran2( 2, NULL, 2, unit_x, unit_y, unit_x, unit_y, 2700, status );

/* Non-square case via astTranN, on both the SIMD and scalar paths. */
   for( pass = 0; pass < 2; pass++ ) {
      mm = astMatrixMap( 2, 3, 1, diag3, " " );
      if( pass == 1 ) {
         astSetI( mm, "UseSIMD", 0 );
      }
      astTranN( mm, 2, 2, 2, in3, 1, 3, 2, out3 );
      if( !approx( out3[ 0 ], 2.0 ) || !approx( out3[ 1 ], 6.0 ) ) {
         stopit( 2800 + pass, status ); /* LCOV_EXCL_LINE */
      }
      if( !approx( out3[ 2 ], 6.0 ) || !approx( out3[ 3 ], -3.0 ) ) {
         stopit( 2810 + pass, status ); /* LCOV_EXCL_LINE */
      }
      if( !approx( out3[ 4 ], 0.0 ) || !approx( out3[ 5 ], 0.0 ) ) {
         stopit( 2820 + pass, status ); /* LCOV_EXCL_LINE */
      }
      mm = astAnnul( mm );
   }
}

/* Forward maps Cartesian (x,y,z) to spherical (longitude, latitude).  Build
   8 points (>= the SIMD batch threshold) from known angles, force the scalar
   path, and check the recovered angles. */
static void TestSphMap( int *status ) {
   double lon[ 8 ] = { 0.1, 0.7, -0.4, 2.0, -1.5, 0.3, 1.2, -2.6 };
   double lat[ 8 ] = { 0.2, -0.5, 0.9, -0.1, 0.6, -1.0, 0.4, -0.3 };
   double in[ 24 ];
   double out[ 16 ];
   double tol = 1.0e-9;
   AstSphMap *sm;
   int idx;

   if( !astOK )
      return;

   sm = astSphMap( " " );
   check_usesimd( (AstMapping *) sm, 3000, status );

   for( idx = 0; idx < 8; idx++ ) {
      in[ idx ] = cos( lat[ idx ] ) * cos( lon[ idx ] );
      in[ 8 + idx ] = cos( lat[ idx ] ) * sin( lon[ idx ] );
      in[ 16 + idx ] = sin( lat[ idx ] );
   }

   astSetI( sm, "UseSIMD", 0 );
   astTranN( sm, 8, 3, 8, in, 1, 2, 8, out );

   for( idx = 0; idx < 8; idx++ ) {
      if( fabs( out[ idx ] - lon[ idx ] ) > tol ) {
         stopit( 3100 + idx, status ); /* LCOV_EXCL_LINE */
      }
      if( fabs( out[ 8 + idx ] - lat[ idx ] ) > tol ) {
         stopit( 3200 + idx, status ); /* LCOV_EXCL_LINE */
      }
   }
}

int main( void ) {
   int status_value = 0;
   int *status = &status_value;

   astWatch( status );
   astBegin;

   TestPolyMap( status );
   TestMatrixMap( status );
   TestMatrixForms( status );
   TestSphMap( status );

   astEnd;
   astFlushMemory( 1 );

   if( *status == 0 ) {
      printf( " All UseSIMD tests passed\n" );
   } else {
      printf( "UseSIMD tests failed\n" ); /* LCOV_EXCL_LINE */
   }
   return *status;
}
