#include "ast.h"
#include "mers.h"
#include "sae_par.h"
#include <math.h>
#include <string.h>


static void stopit( int i, int *status );
static void check_imaging_wcs( AstObject *obj, const char *text, int *status );
static void check_tansip_wcs( AstObject *obj, const char *text, int *status );
static void check_equal_transforms( AstObject *obj, AstObject *obj2, const char *text, int *status );
static void check_sphmap_mappings( AstMapping *map, AstMapping *map2, const char *text, int *status );

static void test_yamlencoding_attribute( int *status );
static void test_imaging_wcs_roundtrip( int *status );
static void test_tansip_wcs_roundtrip( int *status );
static void test_lsst_wcs_roundtrip( int *status );
static void test_native_encoding_roundtrip( int *status );
static void test_sphmap_roundtrip( int *status );

static int chrMatch( const char *a, const char *b ){
   int result = 0;
   if( a && b ) result = !strcmp( a, b );
   return result;
}

int main(){

   int status_value;
   int *status = &status_value;

   status_value = SAI__OK;

   astBegin;

   test_yamlencoding_attribute( status );
   test_imaging_wcs_roundtrip( status );
   test_tansip_wcs_roundtrip( status );
   test_lsst_wcs_roundtrip( status );
   test_native_encoding_roundtrip( status );
   test_sphmap_roundtrip( status );

   astEnd;

   if( *status == SAI__OK ) {
      printf( " All YamlChan tests passed\n" );
   } else {
      printf( "YamlChan tests failed\n" );
   }
   return *status == SAI__OK ? 0 : 1;

}



void stopit( int i, int *status ){
   if( *status == SAI__OK ) {
      printf("Error %d\n", i );
      *status = SAI__ERROR;
   }
}



/* Test getting, setting, and clearing the YamlEncoding attribute. */
void test_yamlencoding_attribute( int *status ){
   AstYamlChan *ch;

   if( *status != SAI__OK ) return;

   ch = astYamlChan( NULL, NULL, " " );

   if( !chrMatch( astGetC( ch, "YAMLENCODING" ), "ASDF" ) ) stopit( -1, status );

   if( astTest( ch, "YAMLENCODING" ) ) stopit( -2, status );

   astSetC( ch, "YamlEncoding", "ASDF" );

   if( !astTest( ch, "YAMLENCODING" ) ) stopit( -3, status );

   if( !chrMatch( astGetC( ch, "YAMLENCODING" ), "ASDF" ) ) stopit( -4, status );

   astClear( ch, "YamlEncoding" );

   if( astTest( ch, "YAMLENCODING" ) ) stopit( -5, status );

   if( !chrMatch(astGetC( ch, "YAMLENCODING" ),"ASDF") ) stopit( -6, status );

   astAnnul( ch );
}



/* Test round-trip of imaging_wcs.asdf: write to yamltest.asdf, read back,
   verify numerical transform results are unchanged. */
void test_imaging_wcs_roundtrip( int *status ){
   AstYamlChan *ch;
   AstObject *obj;
   AstObject *obj2;

   if( *status != SAI__OK ) return;

   ch = astYamlChan( NULL, NULL, " " );
   astSet( ch, "SourceFile=imaging_wcs.asdf,SinkFile=yamltest.asdf" );

   obj = astRead( ch );
   check_imaging_wcs( obj, "Read tests failed for imaging_wcs.asdf", status );

   if( astWrite( ch, obj ) != 1 ) stopit( 13, status );

   astClear( ch, "SourceFile,SinkFile" );
   astSet( ch, "SourceFile=yamltest.asdf" );
   astClear( ch, "YamlEncoding" );

   obj2 = astRead( ch );
   if( !chrMatch( astGetC( ch, "YAMLENCODING" ), "ASDF" ) ) stopit( 131, status );
   check_imaging_wcs( obj2, "Read tests failed for yamltest.asdf", status );

   astAnnul( obj );
   astAnnul( obj2 );
   astAnnul( ch );
}



/* Test round-trip of tanSipWcs: read from tanSipWcs.txt (native AST channel),
   write to tanSipWcs.asdf, read back from ASDF, verify transform values. */
void test_tansip_wcs_roundtrip( int *status ){
   AstYamlChan *ch;
   AstChannel *ch2;
   AstObject *obj;
   AstObject *obj2;

   if( *status != SAI__OK ) return;

   ch  = astYamlChan( NULL, NULL, " " );
   ch2 = astChannel( NULL, NULL, " " );

   astSet( ch2, "SourceFile=tanSipWcs.txt" );
   obj = astRead( ch2 );
   check_tansip_wcs( obj, "Read tests failed for tanSipWcs.txt", status );

   astSet( ch, "SinkFile=tanSipWcs.asdf" );
   if( astWrite( ch, obj ) != 1 ) stopit( 14, status );

   astClear( ch, "SourceFile,SinkFile" );
   astSet( ch, "SourceFile=tanSipWcs.asdf" );
   obj2 = astRead( ch );
   check_tansip_wcs( obj2, "Read tests failed for tanSipWcs.asdf", status );

   astAnnul( obj );
   astAnnul( obj2 );
   astAnnul( ch );
   astAnnul( ch2 );
}



/* Test round-trip of lsst_wcs: read from lsst_wcs.txt, write to
   lsst_wcs.asdf, read back, verify transforms agree numerically. */
void test_lsst_wcs_roundtrip( int *status ){
   AstYamlChan *ch;
   AstChannel *ch2;
   AstObject *obj;
   AstObject *obj2;

   if( *status != SAI__OK ) return;

   ch  = astYamlChan( NULL, NULL, " " );
   ch2 = astChannel( NULL, NULL, " " );

   astSet( ch2, "SourceFile=lsst_wcs.txt" );
   obj = astRead( ch2 );

   astSet( ch, "SinkFile=lsst_wcs.asdf" );
   if( astWrite( ch, obj ) != 1 ) stopit( 15, status );

   astClear( ch, "SinkFile" );
   astSet( ch, "SourceFile=lsst_wcs.asdf" );
   obj2 = astRead( ch );
   if( !obj2 ) stopit( 16, status );

   check_equal_transforms( obj, obj2, "Tests failed for lsst_wcs.txt", status );

   astAnnul( obj );
   astAnnul( obj2 );
   astAnnul( ch );
   astAnnul( ch2 );
}



/* Test NATIVE encoding round-trip: write and read back using YamlEncoding=NATIVE,
   verify the recovered object is equal to the original. */
void test_native_encoding_roundtrip( int *status ){
   AstYamlChan *ch;
   AstChannel *ch2;
   AstObject *obj;
   AstObject *obj2;

   if( *status != SAI__OK ) return;

   ch2 = astChannel( NULL, NULL, " " );
   astSet( ch2, "SourceFile=lsst_wcs.txt" );
   obj = astRead( ch2 );
   astAnnul( ch2 );

   ch = astYamlChan( NULL, NULL, "YamlEncoding=NATIVE " );

   if( !chrMatch( astGetC( ch, "YAMLENCODING" ), "NATIVE" ) ) stopit( 17, status );

   if( !astTest( ch, "YAMLENCODING" ) ) stopit( 18, status );

   astSet( ch, "SinkFile=nativetest.yaml" );
   if( astWrite( ch, obj ) != 1 ) stopit( 19, status );

   astClear( ch, "YamlEncoding" );
   astClear( ch, "SinkFile" );
   astSet( ch, "SourceFile=nativetest.yaml" );
   obj2 = astRead( ch );
   if( !obj2 ) stopit( 20, status );

   if( !astEqual( obj2, obj ) ) stopit( 21, status );

   if( !chrMatch( astGetC( ch, "YAMLENCODING" ), "NATIVE") ) stopit( 22, status );

   astAnnul( obj );
   astAnnul( obj2 );
   astAnnul( ch );
}



/* Test SphMap round-trip: write a FrameSet containing an AST SphMap to
   ASDF, read it back, and verify the recovered mapping gives the same
   numerical results as the original (i.e. keeps the correct units in
   radians). */
void test_sphmap_roundtrip( int *status ){
   AstYamlChan *ch;
   AstFrame *frm3d;
   AstFrame *frm2d;
   AstFrameSet *sphfs;
   AstObject *sphfs2;
   AstMapping *sphmap;
   AstMapping *sphmap2;

   if( *status != SAI__OK ) return;

   frm3d = astFrame( 3, "Domain=CART3D" );
   frm2d = astFrame( 2, "Domain=SPH" );
   sphmap = (AstMapping *) astSphMap( " " );
   sphfs = astFrameSet( frm3d, " " );
   astAddFrame( sphfs, AST__BASE, sphmap, frm2d );
   astAnnul( sphmap );
   astAnnul( frm3d );
   astAnnul( frm2d );

   ch = astYamlChan( NULL, NULL, " " );
   astSet( ch, "SinkFile=sphmap_roundtrip.asdf" );
   if( astWrite( ch, sphfs ) != 1 ) stopit( 30, status );

   astClear( ch, "SinkFile" );
   astSet( ch, "SourceFile=sphmap_roundtrip.asdf" );
   sphfs2 = astRead( ch );
   if( !sphfs2 ) stopit( 31, status );

   sphmap  = astGetMapping( sphfs,  AST__BASE, AST__CURRENT );
   sphmap2 = astGetMapping( (AstFrameSet *) sphfs2, AST__BASE, AST__CURRENT );

   check_sphmap_mappings( sphmap, sphmap2, "SphMap round-trip test failed", status );

   astAnnul( sphmap );
   astAnnul( sphmap2 );
   astAnnul( sphfs );
   astAnnul( sphfs2 );
   astAnnul( ch );
}



void check_imaging_wcs( AstObject *obj, const char *text, int *status ){
   double xout[ 6 ];
   double yout[ 6 ];
   double xin[ 6 ] = {-0.25, -2.0, -1.0,  0.1, 1.5, 1.0 };
   double yin[ 6 ] = { 0.0, 0.0, -2.5, -0.2, 2.5, 2.5};

   if( *status != SAI__OK ) return;

   astTran2( obj, 6, xin, yin, 1, xout, yout );

   if( fabs( xout[ 0 ] - 0.0964300052 ) > 1.0E-10 ) {
      stopit( 1, status );
   } else if( fabs( xout[ 1 ] - 0.0964301088 ) > 1.0E-10 ) {
      stopit( 2, status );
   } else if( fabs( xout[ 2 ] - 0.0964301532 ) > 1.0E-10 ) {
      stopit( 3, status );
   } else if( fabs( xout[ 3 ] - 0.0964299927 ) > 1.0E-10 ) {
      stopit( 4, status );
   } else if( fabs( xout[ 4 ] - 0.0964297983 ) > 1.0E-10 ) {
      stopit( 5, status );
   } else if( fabs( xout[ 5 ] - 0.0964298276 ) > 1.0E-10 ) {
      stopit( 6, status );
   } else if( fabs( yout[ 0 ] + 1.2575438634 ) > 1.0E-10 ) {
      stopit( 7, status );
   } else if( fabs( yout[ 1 ] + 1.25754399404 ) > 1.0E-10 ) {
      stopit( 8, status );
   } else if( fabs( yout[ 2 ] + 1.25754387354 ) > 1.0E-10 ) {
      stopit( 9, status );
   } else if( fabs( yout[ 3 ] + 1.25754383357 ) > 1.0E-10 ) {
      stopit( 10, status );
   } else if( fabs( yout[ 4 ] + 1.25754377796 ) > 1.0E-10 ) {
      stopit( 11, status );
   } else if( fabs( yout[ 5 ] + 1.25754381562 ) > 1.0E-10 ) {
      stopit( 12, status );
   }

   if( *status != SAI__OK ) printf( "%s\n", text );

}



void check_tansip_wcs( AstObject *obj, const char *text, int *status ){
   int i;
   double xout[ 6 ];
   double xin[ 6 ] = {-25.0, -200.0, -100.0,  10.0, 150.0, 100.0 };
   double xrec[ 6 ];
   double xoutT[ 6 ] = { 0.7408749950008373, 0.7397105172925298, 0.7404083214098485, 0.7411104249728614, 0.7420090246546116, 0.741675898352057 };
   double yout[ 6 ];
   double yin[ 6 ] = {0.0,  0.0, -250.0, -20.0, 250.0, 250.0 };
   double yrec[ 6 ];
   double youtT[ 6 ] = { 0.7658201704375505, 0.7658040776691578, 0.7646142651379874, 0.7657273765608887, 0.7670347808637439, 0.7670304352978341 };

   if( *status != SAI__OK ) return;

   astTran2( obj, 6, xin, yin, 1, xout, yout );
   for( i = 0; i < 6; i++ ){
      if( fabs( xout[ i ] - xoutT[ i ] ) > 1.0E-12 ) {
         if( *status == SAI__OK ) printf("%g %g %g\n", xout[ i ], xoutT[ i ], fabs( xout[ i ] - xoutT[ i ] ));
         stopit( i, status );
      } else if( fabs( yout[ i ] - youtT[ i ] ) > 1.0E-12 ) {
         if( *status == SAI__OK ) printf("%g %g %g\n", yout[ i ], youtT[ i ], fabs( yout[ i ] - youtT[ i ] ));
         stopit( i + 6, status );
      }
   }

   astTran2( obj, 6, xout, yout, 0, xrec, yrec );
   for( i = 0; i < 6; i++ ){
      if( fabs( xin[ i ] - xrec[ i ] ) > 1.0E-8 ) {
         if( *status == SAI__OK ) printf("%g %g %g\n", xin[ i ], xrec[ i ], fabs( xin[ i ] - xrec[ i ] ));
         stopit( i + 12, status );
      } else if( fabs( yin[ i ] - yrec[ i ] ) > 1.0E-8 ) {
         if( *status == SAI__OK ) printf("%g %g %g\n", yin[ i ], yrec[ i ], fabs( yin[ i ] - yrec[ i ] ));
         stopit( i + 18, status );
      }
   }

   if( *status != SAI__OK ) printf( "%s\n", text );

}



void check_equal_transforms( AstObject *obj, AstObject *obj2, const char *text, int *status ){

   int i;
   double xout[ 6 ];
   double xout2[ 6 ];
   double xin[ 6 ] = {-25.0, -200.0, -100.0,  10.0, 150.0, 100.0};
   double xrec[ 6 ];
   double xrec2[ 6 ];
   double yout[ 6 ];
   double yout2[ 6 ];
   double yin[ 6 ] = {0.0,  0.0, -250.0, -20.0, 250.0, 250.0};
   double yrec[ 6 ];
   double yrec2[ 6 ];

   if( *status != SAI__OK ) return;

   astTran2( obj, 6, xin, yin, 1, xout, yout );
   astTran2( obj2, 6, xin, yin, 1, xout2, yout2 );

   for( i = 0; i < 6; i++ ){
      if( fabs( xout[ i ] - xout2[ i ] ) > 1.0E-12 ) {
         if( *status == SAI__OK ) printf("%g %g %g\n", xout[ i ], xout2[ i ], fabs( xout[ i ] - xout2[ i ] ));
         stopit( i, status );
      } else if( fabs( yout[ i ] - yout2[ i ] ) > 1.0E-12 ) {
         if( *status == SAI__OK ) printf("%g %g %g\n", yout[ i ], yout2[ i ], fabs( yout[ i ] - yout2[ i ] ));
         stopit( i + 6, status );
      }
   }

   astTran2( obj, 6, xout, yout, 0, xrec, yrec );
   astTran2( obj, 6, xout2, yout2, 0, xrec2, yrec2 );

   for( i = 0; i < 6; i++ ){
      if( fabs( xrec[ i ] - xrec2[ i ] ) > 1.0E-8 ) {
         if( *status == SAI__OK ) printf("%g %g %g\n", xrec[ i ], xrec2[ i ], fabs( xrec[ i ] - xrec2[ i ] ) );
         stopit( i + 12, status );
      } else if( fabs( yrec[ i ] - yrec2[ i ] ) > 1.0E-8 ) {
         if( *status == SAI__OK ) printf("%g %g %g\n", yrec[ i ], yrec2[ i ], fabs( yrec[ i ] - yrec2[ i ] ) );
         stopit( i + 18, status );
      }
   }

   if( *status != SAI__OK ) printf( "%s\n", text );

}



void check_sphmap_mappings( AstMapping *map, AstMapping *map2, const char *text, int *status ){
/* Three test points as unit 3-D Cartesian vectors.
   Layout for astTranN: in[coord * npoint + point], npoint=3, indim=3. */
   double s = 1.0 / sqrt( 3.0 );
   double xin[9] = {
      1.0, 0.0, s,    /* x coords of the 3 points */
      0.0, 1.0, s,    /* y coords */
      0.0, 0.0, s     /* z coords */
   };
   double xout1[6], xout2[6];  /* 2 output coords x 3 points */
   int i;

   if( *status != SAI__OK ) return;

   astTranN( map,  3, 3, 3, xin, 1, 2, 3, xout1 );
   astTranN( map2, 3, 3, 3, xin, 1, 2, 3, xout2 );

   for( i = 0; i < 6; i++ ) {
      if( fabs( xout1[i] - xout2[i] ) > 1.0E-10 ) {
         if( *status == SAI__OK )
            printf( "SphMap round-trip mismatch at output[%d]: %g vs %g\n",
                    i, xout1[i], xout2[i] );
         stopit( 32 + i, status );
         break;
      }
   }

   if( *status != SAI__OK ) printf( "%s\n", text );
}
