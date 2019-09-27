      program testhuge
      implicit none

      include 'AST_PAR'
      include 'SAE_PAR'
      include 'CNF_PAR'

      integer*8 NX
      parameter( NX = 60000 )

      integer*8 NY
      parameter( NY = 60000 )

      integer status
      integer*8 npix, size

c      call ast_watchmemory( 2276565 )

      status = sai__ok
      call err_mark( status )
      call ast_begin( status )

      call psx_calloc8( NX*NY, '_BYTE', ip, status )


      call fill( NX, NY, %val(cnf_pval(ip)), status )


xxx



      call ast_end( status )
      call err_rlse( status )

c      call ast_activememory( 'testmoc' )
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


