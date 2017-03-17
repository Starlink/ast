      program testchebymap
      implicit none

      include 'SAE_PAR'
      include 'AST_PAR'
      include 'PRM_PAR'

      integer status, lstat, cm, cm2, cm3, i, j, nco
      double precision lbnd( 2 ), ubnd( 2 ), dval
      double precision tlbnd( 2 ), tubnd( 2 )

      double precision coeffs_1( 4*3 ), xin( 5 ), xout( 5 ), xrec( 5 ),
     :                 yrec( 5 ),
     :                 work( 5 ), coeffs_2( 2*3 ), coeffs_3( 4*5 ),
     :                 yin(5), yout(5), xi, yi, xv, yv, y, a, x,
     :                 cofs( 100 )


C  f(x) = 1.5*T0(x') - 1.0*T2(x') + 2.0*T3(x') - 1.3*T4(x')
      data coeffs_1 /1.5D0, 1.0D0, 0.0D0,
     :              -1.0D0, 1.0D0, 2.0D0,
     :               2.0D0, 1.0D0, 3.0D0,
     :               1.3D0, 1.0D0, 4.0D0 /

C  f(x) = 1.0*T0(x') - 2.0*T1(x')
      data coeffs_2 /1.0D0, 1.0D0, 0.0D0,
     :              -2.0D0, 1.0D0, 1.0D0 /


C  fx(x,y) = 1.0*T0(x')*T0(y') - 2.0*T1(x')*T2(y') + T1(y')
C  fy(x,y) = 1.5*T0(x')*T0(y') - 2.5*T1(x')*T2(y')
      data coeffs_3 /1.0D0, 1.0D0, 0.0D0, 0.0D0,
     :              -2.0D0, 1.0D0, 1.0D0, 2.0D0,
     :               1.0D0, 1.0D0, 0.0D0, 1.0D0,
     :               1.5D0, 2.0D0, 0.0D0, 0.0D0,
     :              -2.5D0, 2.0D0, 1.0D0, 2.0D0/




      status = sai__ok
      call ast_begin( status )

C  One-dimensional ChebyMaps, order 1: a constant equal to 1.5
      lbnd( 1 ) = -1.0D0
      lbnd( 2 ) = -1.0D0
      ubnd( 1 ) = 1.0D0
      ubnd( 2 ) = 1.0D0

      cm = ast_chebymap( 1, 1, 1, coeffs_1, 0, 0.0D0, lbnd, ubnd, 1.0D0,
     :                   1.0D0, ' ', status )

      xin( 1 ) = -1.0D0
      xin( 2 ) = -0.5D0
      xin( 3 ) =  0.0D0
      xin( 4 ) =  0.5D0
      xin( 5 ) =  1.0D0

      call ast_tran1( cm, 5, xin, .true., xout, status )
      do i = 1, 5
         if( xout( i ) .ne. 1.5D0 ) call stopit( 0, status )
      end do

C  One-dimensional ChebyMaps, order 3: 2.5 - 2*x*x
      cm = ast_chebymap( 1, 1, 2, coeffs_1, 0, 0.0D0, lbnd, ubnd, 1.0D0,
     :                   1.0D0, ' ', status )
      call ast_tran1( cm, 5, xin, .true., xout, status )
      do i = 1, 5
         if( xout( i ) .ne. 2.5D0 - 2.0D0*xin( i )**2 )
     :        call stopit( 1, status )
      end do

C  One-dimensional ChebyMaps, order 4: 2.5 - 6*x - 2*x*x + 8*x*x*x
      cm = ast_chebymap( 1, 1, 3, coeffs_1, 0, 0.0D0, lbnd, ubnd, 1.0D0,
     :                   1.0D0, ' ', status )
      call ast_tran1( cm, 5, xin, .true., xout, status )
      do i = 1, 5
         dval = 2.5D0-6.0D0*xin(i)-2.0D0*xin(i)**2+8.0D0*xin(i)**3
         if( xout( i ) .ne. dval ) call stopit( 2, status )
      end do

C  One-dimensional ChebyMaps, order 5
      cm = ast_chebymap( 1, 1, 4, coeffs_1, 0, 0.0D0, lbnd, ubnd, 1.0D0,
     :                   1.0D0, ' ', status )
      call ast_tran1( cm, 5, xin, .true., xout, status )

      do i = 1, 5
         work( 1 ) = 1.0D0
         work( 2 ) = xin( i )
         do j = 3, 5
            work( j ) =  2.0D0 * xin( i ) * work( j - 1 )
     :                                    - work( j - 2 )
         end do

         if( 1.5D0*work(1) - 1.0D0*work(3) +
     :       2.0D0*work(4) + 1.3D0*work(5) .ne.
     :       xout( i ) ) call stopit( 3, status )

      end do

C  Check the IterInverse attribute is zero.
      if( ast_getl( cm, 'IterInverse', status ) ) then
         call stopit( 4, status )
      end if

c  The astPolyTran method on a 1-dimensional ChebyMaps, order 2: 1 - 2*x
      cm = ast_chebymap( 1, 1, 2, coeffs_2, 0, 0.0D0, lbnd, ubnd, 1.0D0,
     :                   1.0D0, ' ', status )
      cm2 = ast_polytran( cm, .false., 0.01D0, 0.01D0, 5, lbnd,
     :                    ubnd, status )
      if( cm2 .eq. AST__NULL ) then
         call stopit( 5, status )
      else
         call ast_tran1( cm2, 5, xout, .false., xrec, status )
         do i = 1, 5
            dval = ( 1.0D0 - xout(i) )/2.0D0
            if( abs( xrec( i ) - dval ) .gt. 1.0D-3*abs( dval ) )
     :             call stopit( 6, status )
         end do
      end if

c  The astPolyTran method on a 1-dimensional ChebyMaps, order 5.
      lbnd( 1 ) = -100.0D0
      ubnd( 1 ) = 100.0D0
      cm = ast_chebymap( 1, 1, 4, coeffs_1, 0, 0.0D0, lbnd, ubnd, 1.0D0,
     :                   1.0D0, ' ', status )
      cm2 = ast_polytran( cm, .false., 0.01D0, 0.01D0, 10, -5.0D0,
     :                    50.0D0, status )

      if( cm2 .eq. AST__NULL ) then
         call stopit( 7, status )
      else
         xin(1) = 0.0D0;
         xin(2) = 10.0D0;
         xin(3) = 20.0D0;
         xin(4) = 30.0D0;
         xin(5) = 40.0D0;
         call ast_tran1( cm2, 5, xin, .true., xout, status )
         call ast_tran1( cm2, 5, xout, .false., xrec, status )

         do i = 1, 5
            if( abs( xrec( i ) - xin( i ) ) .gt. 0.01D0 )
     :             call stopit( 8, status )
         end do
      end if

c ast_equal and ast_copy
      cm3 = ast_copy( cm2, status )
      if( .not. ast_equal( cm2, cm3, status ) ) then
         call stopit( 9, status )
      end if

c astDump and astLoadChebyMap
      call checkdump( cm2, status )

*  2-dimensional ChebyMaps: forward transformation
      lbnd(1) = 0.0D0
      lbnd(2) = 0.0D0
      ubnd(1) = 10.0D0
      ubnd(2) = 10.0D0

      cm = ast_chebymap( 2, 2, 5, coeffs_3, 0, 0.0D0, lbnd, ubnd,
     :                   1.0D0, 1.0D0, ' ', status )

      xin(1) = 0.0D0
      xin(2) = 2.0D0
      xin(3) = 6.0D0
      xin(4) = 10.0D0

      yin(1) = 2.0D0
      yin(2) = 5.0D0
      yin(3) = 8.0D0
      yin(4) = 0.0D0

      call ast_tran2( cm, 4, xin, yin, .true., xout, yout, status )
      do i = 1, 4
         xi = 2.0D0*( xin(i) - lbnd(1) )/( ubnd(1) - lbnd(1) ) - 1.0D0
         yi = 2.0D0*( yin(i) - lbnd(2) )/( ubnd(2) - lbnd(2) ) - 1.0D0

         xv = 1 - 2*xi*(2*yi**2 - 1) + yi
         yv = 1.5 - 2.5*xi*(2*yi**2 - 1)

         if( abs( xout(i) - xv ) .gt. 1.0D-6*abs(xv) .or.
     :       abs( yout(i) - yv ) .gt. 1.0D-6*abs(yv) ) then
            call stopit( 10, status )
         end if
      end do


*  2-dimensional ChebyMaps: fitted inverse transformation
      tlbnd(1) = 4.0D0
      tlbnd(2) = 4.0D0
      tubnd(1) = 6.0D0
      tubnd(2) = 6.0D0

      cm2 = ast_polytran( cm, .false., 0.01D0, 0.01D0, 10, tlbnd,
     :                    tubnd, status )

      if( cm2 .eq. AST__NULL ) then
         call stopit( 11, status )
      else
         xin(1) = 4.0D0
         xin(2) = 4.5D0
         xin(3) = 5.0D0
         xin(4) = 5.5D0

         yin(1) = 6.0D0
         yin(2) = 5.5D0
         yin(3) = 5.0D0
         yin(4) = 4.5D0

         call ast_tran2( cm2, 4, xin, yin, .true., xout, yout, status )
         call ast_tran2( cm2, 4, xout, yout, .false., xrec, yrec,
     :                   status )
         do i = 1, 4
            if( abs( xrec(i) - xin(i) ) .gt. 0.01D0 .or.
     :          abs( yrec(i) - yin(i) ) .gt. 0.01D0 ) then
               call stopit( 12, status )
            end if
         end do
      end if



      call ast_polycoeffs( cm2, .true., 0, 0.0D0, nco, status )
      if( nco .ne. 5 ) then
         call stopit( 13, status )
      endif

      call ast_polycoeffs( cm2, .true., 100, cofs, nco, status )
      if( nco .ne. 5 ) then
         call stopit( 14, status )
      else
         do i = 1, 20
            if( cofs( i ) .ne. coeffs_3( i ) ) then
               call stopit( 15, status )
            end if
         end do
      endif


      call ast_polycoeffs( cm2, .false., 0, 0.0D0, nco, status )
      if( nco .ne. 12 ) then
         call stopit( 16, status )
      endif

      call ast_polycoeffs( cm2, .false., 100, cofs, nco, status )
      if( nco .ne. 12 ) then
         call stopit( 17, status )
      else if( cofs( 2 ) .ne. 1.0D0 ) then
         call stopit( 18, status )
      else if( abs( cofs( 13 ) - (-6.5376876905037529) ) .gt.
     :        1.0E-6 ) then
         call stopit( 19, status )
      else if( cofs( 15 ) .ne. 2.0D0 ) then
         call stopit( 20, status )
      end if







      call ast_end( status )
      call ast_activememory( 'testchebymap' );
      call ast_flushmemory( 1 )

      if( status .eq. sai__ok ) then
         write(*,*) 'All ChebyMap tests passed'
      else
         write(*,*) 'ChebyMap tests failed'
      end if

      end


      subroutine stopit( i, status )
      implicit none
      include 'SAE_PAR'
      integer i, status
      if( status .eq. sai__ok ) then
         write( *,* ) 'Error ',i
         status = sai__error
      end if
      end


      subroutine checkdump( obj, status )
      implicit none
      include 'SAE_PAR'
      include 'AST_PAR'
      integer obj, status, ch, result

      if( status .ne. sai__ok ) return

      ch = ast_channel( AST_NULL, AST_NULL, ' ', status )

      call ast_set( ch, 'SinkFile=fred.tmp', status )
      if( ast_write( ch, obj, status ) .ne. 1 ) then
         call stopit( -1, status )
      end if
      call ast_clear( ch, 'SinkFile', status )

      call ast_set( ch, 'SourceFile=fred.tmp', status )
      result = ast_read( ch, status )
      if( result .eq. ast__null ) then
         call stopit( -2, status )
      end if
      call ast_clear( ch, 'SourceFile', status )

      if( .not. ast_equal( result, obj, status ) ) then
         call ast_show( obj, status )
         call ast_show( result, status )
         call stopit( -3, status )
      end if

      end
