      program testmapping
      implicit none

      include 'AST_PAR'
      include 'SAE_PAR'

      integer status, pm, m1, m2, m3, m4, m5, m6
      logical series, inv1, inv2
      double precision  coeff(20), fit(6), lbnd(2), ubnd(2), shift(1)

      data coeff / 1.0, 1.0, 0.0, 0.0,
     :             2.0, 1.0, 1.0, 0.0,
     :             1.0, 2.0, 0.0, 0.0,
     :             3.0, 2.0, 0.0, 1.0,
     :             3.0, 1.0, 0.0, 2.0 /



      status = sai__ok
      call err_mark( status )
      call ast_begin( status )

      pm = ast_polymap( 2, 2, 4, coeff, 0, coeff, ' ', status )

      lbnd( 1 ) = -1.0D0
      lbnd( 2 ) = -1.0D0
      ubnd( 1 ) = 1.0D0
      ubnd( 2 ) = 1.0D0
      if( ast_linearapprox(pm, lbnd, ubnd, 0.001D0, fit, status) ) then
         if( fit(1) .ne. 1.0D0 .or. fit(2) .ne. 1.0D0 .or.
     :       fit(3) .ne. 2.0D0 .or. fit(4) .ne. 0.0D0 .or.
     :       fit(5) .ne. 0.0D0 .or. fit(6) .ne. 3.0D0 ) then
            call stopit( status, 'Error 0' )
         end if
      else
         call stopit( status, 'Error 1' )
      end if

      coeff( 13 ) = AST__BAD
      pm = ast_polymap( 2, 2, 4, coeff, 0, coeff, ' ', status )

      if( ast_linearapprox(pm, lbnd, ubnd, 0.001D0, fit, status) ) then
         if( fit(1) .ne. 1.0D0 .or. fit(2) .ne. AST__BAD .or.
     :       fit(3) .ne. 2.0D0 .or. fit(4) .ne. 0.0D0 .or.
     :       fit(5) .ne. AST__BAD .or. fit(6) .ne. AST__BAD ) then
            call stopit( status, 'Error 2' )
         end if
      else
         call stopit( status, 'Error 3' )
      end if

      pm = ast_polymap( 2, 2, 5, coeff, 0, coeff, ' ', status )

      if( ast_linearapprox(pm, lbnd, ubnd, 0.001D0, fit, status) ) then
         write(*,*) fit
         call stopit( status, 'Error 4' )
      end if


*  Test Protect attribute
      shift( 1 ) = 1.0D0
      m1 = ast_shiftmap( 1, shift, ' ', status )
      m2 = ast_zoommap( 1, 2.0D0, ' ', status )
      m3 = ast_cmpmap( m1, m2, .true., ' ', status )
      m4 = ast_copy( m3, status )
      call ast_invert( m4, status )
      m5 = ast_cmpmap( m3, m4, .true., ' ', status )
      m6 = ast_simplify( m5, status )
      if( .not. ast_isaunitmap( m6, status ) )then
         call stopit( status, 'Error 5' )
      end if

      call ast_setl( m1, 'Protect', .true., status )
      if( .not. ast_getl( m1, 'Protect', status ) ) then
         call stopit( status, 'Error 5a' )
      else if( ast_getl( m3, 'Protect', status ) ) then
         call stopit( status, 'Error 5b' )
      endif

      m4 = ast_copy( m3, status )
      call ast_invert( m4, status )
      m5 = ast_cmpmap( m3, m4, .true., ' ', status )
      m6 = ast_simplify( m5, status )
      if( ast_isaunitmap( m6, status ) )then
         call stopit( status, 'Error 6' )
      end if
      call ast_decompose( m6, m1, m2, series, inv1, inv2, status )
      if( .not. ast_isashiftmap( m1, status ) ) then
         call stopit( status, 'Error 7' )
      end if
      if( .not. ast_isashiftmap( m2, status ) ) then
         call stopit( status, 'Error 8' )
      end if
      if( .not. series .or. inv1 .or. .not. inv2 ) then
         call stopit( status, 'Error 9' )
      end if

      call ast_setl( m5, 'Protect', .TRUE., status )
      m6 = ast_simplify( m5, status )
      call ast_decompose( m6, m1, m2, series, inv1, inv2, status )
      if( .not. ast_isacmpmap( m1, status ) ) then
         call stopit( status, 'Error 10' )
      end if
      if( .not. ast_isacmpmap( m2, status ) ) then
         call stopit( status, 'Error 11' )
      end if
      if( .not. series .or. inv1 .or. .not. inv2 ) then
         call stopit( status, 'Error 12' )
      end if

      if( .not. ast_getl( m1, 'Protect', status ) ) then
         call stopit( status, 'Error 13' )
      else if( .not. ast_getl( m2, 'Protect', status ) ) then
         call stopit( status, 'Error 14' )
      endif

      call ast_end( status )
      call err_rlse( status )

      if( status .eq. sai__ok ) then
         write(*,*) 'All Mapping tests passed'
      else
         write(*,*) 'Mapping tests failed'
      end if

      end



      subroutine stopit( status, text )
      implicit none
      include 'SAE_PAR'
      integer status
      character text*(*)
      if( status .ne. sai__ok ) return
      status = sai__error
      write(*,*) text
      end



