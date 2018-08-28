!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       subroutine dft_get_mm_forces(dxyzcl,dxyzqm)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       use garcha_mod
       use faint_cpu, only: int1G, intsolG
       implicit real*8 (a-h,o-z)
       REAL*8 , intent(inout) :: dxyzqm(3,natom)
       REAL*8 , intent(inout) :: dxyzcl(3,nsol)
       real*8, dimension (:,:), ALLOCATABLE :: ff,ffcl!,g2gff,g2gffcl
       integer :: MM
       !real*8 diff,rms,rmscl,mx,mxcl,s
c       real*8, dimension (:,:), ALLOCATABLE :: ffs,ffcls
!
!
!--------------------------------------------------------------------!
       !allocate(ff(natom,3), ffcl(nsol,3))!ntatom,3))
       !allocate(g2gff(natom,3), g2gffcl(ntatom,3))
c       allocate(ffs(natom,3), ffcls(ntatom,3))
c       real*8 ftot(3)
       MM = M*(M+1)/2
       if (nsol.le.0.or.cubegen_only) return
       factor=1.D0

       call aint_query_gpu_level(igpu)
       call g2g_timer_sum_start('QM/MM gradients')
       if (igpu.lt.2) then
         ! The old version of intsolG expected the MM force array to be
         ! padded in front with # QM atoms spots for some reason
         allocate(ff(natom,3), ffcl(ntatom,3))
         ffcl=0
         ff=0

         call g2g_timer_start('intsolG')
         call intsolG(ff,ffcl,natom, ntatom, RMM(1:MM), d, r, pc, Iz)
         call g2g_timer_stop('intsolG')

         do jj=1,nsol
         do j=1,3
           dxyzcl(j,jj)=ffcl(natom+jj,j)*factor!+dxyzcl(j,jj)
         enddo
         enddo
       else
         ! The GPU version of the QM/MM gradients only uses space for the MM
         ! forces in the MM force array
         allocate(ff(natom,3), ffcl(nsol,3))
         ffcl=0
         ff=0

         if (igpu.gt.3) call int1G(ff, RMM(1:MM), d, r, Iz, natom,
     >    ntatom)
         call g2g_timer_start('aint_qmmm_forces')
         call aint_qmmm_forces(ff,ffcl)
         call g2g_timer_stop('aint_qmmm_forces')

         do jj=1,nsol
         do j=1,3
           dxyzcl(j,jj)=ffcl(jj,j)*factor!+dxyzcl(j,jj)
         enddo
         enddo
       endif

       !mx = 0
       !s = 0
       !do i=1,natom
         !do j=1,3
           !diff = abs(ff(i,j)-g2gff(i,j))
           !mx = max(diff,mx)
           !s = s + diff**2
         !enddo
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

       do i=1,natom
       do j=1,3
         dxyzqm(j,i)=ff(i,j)*factor+dxyzqm(j,i)
       enddo
       enddo

       call g2g_timer_sum_stop('QM/MM gradients')
       call g2g_timer_sum_stop('Forces')

!
!
!--------------------------------------------------------------------!
       deallocate (ff,ffcl)!,g2gff,g2gffcl)
c       deallocate (ffs,ffcls)
897    format (F17.11)
       return;end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
