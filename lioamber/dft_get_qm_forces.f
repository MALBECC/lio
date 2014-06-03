      subroutine dft_get_qm_forces(dxyzqm)
      use garcha_mod
c      use qmmm_module, only : qmmm_struct
      implicit real*8 (a-h,o-z)

       REAL*8 , intent(inout) :: dxyzqm(3,natom)
       real*8, dimension (:,:), ALLOCATABLE :: ff
       real*8 ftot(3)


       allocate(ff(natom,3))
         ff=0
        call g2g_timer_start('int1G'//CHAR(0))
        call int1G(ff)
        call g2g_timer_stop('int1G'//CHAR(0))
        call g2g_timer_start('intSG'//CHAR(0))
        call intSG(ff)
        call g2g_timer_stop('intSG'//CHAR(0))
c         write(77,*) ff
        call g2g_timer_start('int3G'//CHAR(0))
       call int3G(ff,.true.)
        call g2g_timer_stop('int3G'//CHAR(0))

c        factor=627.509391D0/0.5291772108D0
        factor=1.D0
       do i=1,natom 
        do j=1,3
       dxyzqm(j,i)=ff(i,j)*factor
        enddo
         enddo
      
         deallocate (ff)

      end subroutine dft_get_qm_forces
c---------------------------------------------------------------------
