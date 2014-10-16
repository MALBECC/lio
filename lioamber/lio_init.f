!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       subroutine lio_init()
!--------------------------------------------------------------------!
!
! (*) Allocation of globar variables.
! (*) Open files to write.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       use garcha_mod
       implicit none
       call g2g_timer_start('lio_init')
       write(6,'(A)') ''
!
!
! VARIABLE INITIALIZATION
!----------------------------------------------------------!
       write(6,'(A)') ' (*) Initializing lio variables...'
       if (.not.allocated(Smat))    allocate(Smat(M,M))
       if (.not.allocated(RealRho)) allocate(RealRho(M,M))
!
!
!
!----------------------------------------------------------!
       call g2g_timer_stop('lio_init')
       return;end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
