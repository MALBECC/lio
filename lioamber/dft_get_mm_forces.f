!½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½!
       subroutine dft_get_mm_forces(dxyzcl,dxyzqm)
!½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½!
       use garcha_mod
       implicit real*8 (a-h,o-z)
       REAL*8 , intent(inout) :: dxyzqm(3,natom)
       REAL*8 , intent(inout) :: dxyzcl(3,nsol)
       real*8, dimension (:,:), ALLOCATABLE :: ff,ffcl
c       real*8, dimension (:,:), ALLOCATABLE :: ffs,ffcls
!
!
!--------------------------------------------------------------------!
       allocate(ff(natom,3), ffcl(ntatom,3))
c       allocate(ffs(natom,3), ffcls(ntatom,3))
c       real*8 ftot(3)
       ffcl=0
       ff=0

       call g2g_timer_start('intsolG')
       call intsolG(ff,ffcl)
       call g2g_timer_stop('intsolG')

       factor=1.D0
       do i=1,natom 
       do j=1,3
         dxyzqm(j,i)=ff(i,j)*factor+dxyzqm(j,i)
       enddo
       enddo

       do jj=1,nsol
       do j=1,3
         dxyzcl(j,jj)=ffcl(natom+jj,j)*factor!+dxyzcl(j,jj)
       enddo
       enddo         
c       write(*,*) 'aca toyyyy'
!
!
!--------------------------------------------------------------------!       
       deallocate (ff,ffcl) 
c       deallocate (ffs,ffcls) 
897    format (F17.11)
       return;end subroutine
!½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½!
