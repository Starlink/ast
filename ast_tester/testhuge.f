      program testhuge
      implicit none

      include 'AST_PAR'
      include 'SAE_PAR'
      include 'CNF_PAR'

      integer*8 NX
      parameter( NX = 100 )

      integer*8 NY
      parameter( NY = 100 )

      integer status, map, ip1, ip2
      integer*8 ubnd(2), lbnd(2), lbnd2(2), ubnd2(2), nbad
      double precision shifts(2)

c      call ast_watchmemory( 2276565 )

      status = sai__ok
      call err_mark( status )
      call ast_begin( status )

      call psx_calloc8( NX*NY, '_BYTE', ip1, status )
      call psx_calloc8( NX*NY, '_BYTE', ip2, status )


      call fill( NX, NY, %val(cnf_pval(ip1)), status )

      shifts(1) = 0.5D0
      shifts(2) = 0.5D0
      map = ast_shiftmap( 2, shifts, ' ', status )

      lbnd(1) = -10
      ubnd(1) = lbnd(1) - 1 + NX
      lbnd(2) = -10
      ubnd(2) = lbnd(2) - 1 + NY

      lbnd2(1) = lbnd(1) + 100
      ubnd2(1) = ubnd(1) - 100
      lbnd2(2) = lbnd(2) + 100
      ubnd2(2) = ubnd(2) - 100

      write(*,*) 'Resampling...'
      nbad = ast_resampleB( map, 2, lbnd, ubnd, %val(cnf_pval(ip1)),
     :                      shifts, AST__LINEAR, AST_NULL, shifts,
     :                      AST__USEBAD, 0.1, 50, 0, 2, lbnd, ubnd,
     :                      lbnd2, ubnd2, %val(cnf_pval(ip2)), shifts,
     :                      status )

      write(*,*) 'nbad = ',nbad


      call dump( NX, NY, %val(cnf_pval(ip1)), 'input.asc', status )
      call dump( NX, NY, %val(cnf_pval(ip2)), 'output.asc', status )




      call psx_free( ip1, status )
      call psx_free( ip2, status )

      call ast_end( status )
      call err_rlse( status )

c      call ast_activememory( 'testhuge' )
      call ast_flushmemory( 1 )

      if( status .eq. sai__ok ) then
         write(*,*) 'All Huge tests passed'
      else
         write(*,*) 'Huge tests failed'
      end if

      end



      subroutine stopit( text, status )
      implicit none
      include 'SAE_PAR'
      integer status
      character text*(*)
      if( status .ne. sai__ok ) return
      status = sai__error
      write(*,*) text
      end


      subroutine fill( nx, ny, array, status )
      implicit none
      include 'SAE_PAR'
      integer status
      integer*8 nx, ny, ix, iy, cx, cy, dy, jy
      byte array(nx,ny)
      real factor

      if( status .ne. sai__ok ) return

      cx = nx/2
      cy = ny/2

      factor = -250.0/sqrt( real( cx**2 + cy**2 ) )

      dy = ny/10
      jy = dy
      do iy = 1, ny

         if( iy .eq. jy ) then
            write(*,*) 'Filling row ',iy,' of ',ny
            jy = jy + dy
         end if

         do ix = 1, nx
            array( ix, iy ) = sqrt( real((ix - cx)**2 + (iy - cy)**2) )*
     :                        factor + 128.0
         end do

      end do

      end

      subroutine dump( nx, ny, array, fname, status )
      implicit none
      include 'SAE_PAR'
      integer status
      integer*8 nx, ny, ix, iy
      byte array(nx,ny)
      character fname*(*)

      if( status .ne. sai__ok ) return

      write(*,*) 'Writing ',fname
      open( 10, file=fname, status='new' )

      write(10,*) '# ix iy array'
      do iy = 1, ny
         do ix = 1, nx
            write(10,*) ix, iy, array( ix, iy )
         end do
      end do

      close(10)

      end



