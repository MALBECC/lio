!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       subroutine dft_get_qm_forces(dxyzqm)
!--------------------------------------------------------------------!
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       use garcha_mod, only: natom, nsol, cubegen_only, number_restr,
     &                       first_step, doing_ehrenfest, r,
     &                       qm_forces_ds, qm_forces_total,
     &                       RMM, d, Iz, natom, ntatom, M, Md

       use ehrendata, only: nullify_forces
       use faint_cpu, only: int1G, intSG
       use faint_cpu77, only: int3G
       implicit none
       real*8,intent(out) :: dxyzqm(3,natom)
       real*8,allocatable :: ff1G(:,:),ffSG(:,:),ff3G(:,:)
       real*8             :: factor
       integer            :: fileunit,kk,ii,igpu, MM, M15, MMd
       logical            :: print_forces
!variables for restrain calculations
       real*8 :: f_r
       integer :: i

       MM  = M  * (M  +1) / 2
       MMd = Md * (Md +1) / 2
       M15 = 1 + M + 4*MM + 2*MMd
!--------------------------------------------------------------------!
       if(cubegen_only) return
       call g2g_timer_sum_start('Forces')
       allocate(ff1G(natom,3),ffSG(natom,3),ff3G(natom,3))

       call g2g_timer_start('int1G')
       ff1G=0.0d0
       call aint_query_gpu_level(igpu)
       if (igpu.lt.4) then
         call g2g_timer_sum_start('Nuclear attraction gradients')
         call int1G(ff1G, RMM(1:MM), d, r, Iz, natom, ntatom)
         call g2g_timer_sum_stop('Nuclear attraction gradients')
       elseif (nsol.le.0) then
         call g2g_timer_sum_start('Nuclear attraction gradients')
         call int1G(ff1G, RMM(1:MM), d, r, Iz, natom, ntatom)
         call aint_qmmm_forces(ff1G,0)
         call g2g_timer_sum_stop('Nuclear attraction gradients')
       endif
       call g2g_timer_stop('int1G')

       call g2g_timer_start('intSG')
       call g2g_timer_sum_start('Overlap gradients')
       ffSG=0.0d0
!       if ( (doing_ehrenfest) .and. (.not.first_step) ) then
       if (doing_ehrenfest) then
          ffSG=-transpose(qm_forces_ds)
       else
          call intSG(ffSG, RMM(M15:M15+MM), r, d, natom, ntatom)
       endif
       call g2g_timer_stop('intSG')
       call g2g_timer_sum_stop('Overlap gradients')

       call g2g_timer_start('int3G')
       call g2g_timer_sum_start('Coulomb+Exchange-correlation')
       ff3G=0.0d0
       call int3G(ff3G,.true.)
       call g2g_timer_stop('int3G')
       call g2g_timer_sum_stop('Coulomb+Exchange-correlation')

       factor=1.D0
c       factor=627.509391D0/0.5291772108D0
       do kk=1,natom
       do ii=1,3
         dxyzqm(ii,kk)=ff1G(kk,ii)+ffSG(kk,ii)+ff3G(kk,ii)
         dxyzqm(ii,kk)=dxyzqm(ii,kk)*factor
       enddo
       enddo

!
! FFR - Ehrenfest needs to keep track of forces
!--------------------------------------------------------------------!
       if ( nullify_forces ) dxyzqm(:,:)=0.0d0
       if ( doing_ehrenfest ) then
         qm_forces_total=qm_forces_ds
         qm_forces_total=qm_forces_total-transpose(ff1G)
         qm_forces_total=qm_forces_total-transpose(ff3G)
       endif


!--------------------------------------------------------------------!
        IF (number_restr.GT.0) THEN
! distance restrain case
          call get_restrain_forces(dxyzqm, f_r)
	  WRITE(*,*) "DISTANCE RESTRAIN ADDED TO FORCES"
        END IF


! FFR: force calculation should be separated from force passing and
!      force writing. All can be in the same module, but different
!      subroutines.
!--------------------------------------------------------------------!
       print_forces=.false.
       if (print_forces) then
         fileunit=3242
         if (first_step) then
            open(unit=fileunit,file='Forces.log')
         else
            open(unit=fileunit,file='Forces.log',access='APPEND')
         endif

         write(fileunit,'(A)')
     >   '------------------------------------------------------------'
         do kk=1,natom
            write(fileunit,200) 'TOTS', kk,
     >         dxyzqm(1,kk), dxyzqm(2,kk), dxyzqm(3,kk)
!           write(fileunit,200) 'TOTS',kk,
!     >       ff1G(kk,1)+ffSG(kk,1)+ff3G(kk,1),
!     >       ff1G(kk,2)+ffSG(kk,2)+ff3G(kk,2),
!     >       ff1G(kk,3)+ffSG(kk,3)+ff3G(kk,3)
         enddo
         write(fileunit,'(A)')
     >   '------------------------------------------------------------'
         write(fileunit,*)

         do kk=1,natom
           write(fileunit,200)
     >       'FF1G',kk,ff1G(kk,1),ff1G(kk,2),ff1G(kk,3)
           write(fileunit,200)
     >       'FFSG',kk,ffSG(kk,1),ffSG(kk,2),ffSG(kk,3)
           write(fileunit,200)
     >       'FF3G',kk,ff3G(kk,1),ff3G(kk,2),ff3G(kk,3)
           write(fileunit,*)
         enddo
         write(fileunit,*)

         close(fileunit)
       endif

       if (nsol.le.0) then
         call g2g_timer_sum_stop('Forces')
       endif

! FFR: No other place for this to go right now.
       if ( first_step ) first_step = .false.
!--------------------------------------------------------------------!
       deallocate(ff1G,ffSG,ff3G)
 200   format(1X,A4,1X,I4,3(2X,E14.7))
       return;end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
