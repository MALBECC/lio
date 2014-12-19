!½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½!
       subroutine dft_get_mm_forces(dxyzcl,dxyzqm)
!½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½!
       use garcha_mod
       implicit real*8 (a-h,o-z)
       REAL*8 , intent(inout) :: dxyzqm(3,natom)
       REAL*8 , intent(inout) :: dxyzcl(3,nsol)
       real*8, dimension (:,:), ALLOCATABLE :: ff,ffcl!,g2gff,g2gffcl
       !real*8 diff,rms,rmscl,mx,mxcl,s
c       real*8, dimension (:,:), ALLOCATABLE :: ffs,ffcls
!
!
!--------------------------------------------------------------------!
       allocate(ff(natom,3), ffcl(nsol,3))!ntatom,3))
       !allocate(g2gff(natom,3), g2gffcl(ntatom,3))
c       allocate(ffs(natom,3), ffcls(ntatom,3))
c       real*8 ftot(3)
       ffcl=0
       ff=0

       !call g2g_timer_start('intsolG')
       !call intsolG(ff,ffcl)
       !call g2g_timer_stop('intsolG')

       !g2gff=0
       !g2gffcl=0
       call g2g_timer_start('g2g_qmmm_forces')
       call g2g_qmmm_forces(ff,ffcl)
       call g2g_timer_stop('g2g_qmmm_forces')

       !mx = 0
       !s = 0
       !do i=1,natom
       !  do j=1,3
       !    diff = abs(ff(i,j)-g2gff(i,j))
       !    mx = max(diff,mx)
       !    s = s + diff**2
       !  enddo
       !enddo
       !rms = sqrt(s/(natom*3))

       !clmx = 0
       !s = 0
       !do i=1,nsol
       !  do j=1,3
       !    diff = abs(ffcl(i+natom,j)-g2gffcl(i,j))
       !    mxcl = max(diff,mxcl)
       !    s = s + diff**2
       !  enddo
       !enddo
       !rmscl = sqrt(s/(nsol*3))
       !write(*,*) "MAX DIFF QM: ",mx
       !write(*,*) "MAX DIFF MM: ",mxcl
       !write(*,*) "RMS DIFF QM: ",rms
       !write(*,*) "RMS DIFF MM: ",rmscl

       factor=1.D0
       do i=1,natom 
       do j=1,3
         dxyzqm(j,i)=ff(i,j)*factor+dxyzqm(j,i)
       enddo
       enddo

       do jj=1,nsol
       do j=1,3
         dxyzcl(j,jj)=ffcl(jj,j)*factor!natom+jj,j)*factor!+dxyzcl(j,jj)
       enddo
       enddo         
c       write(*,*) 'aca toyyyy'
!
!
!--------------------------------------------------------------------!       
       deallocate (ff,ffcl)!,g2gff,g2gffcl) 
c       deallocate (ffs,ffcls) 
897    format (F17.11)
       return;end subroutine
!½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½!
