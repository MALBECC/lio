!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% RESTART_TD.F90 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! * read_td_restart_verlet_d (reads TD restart for Verlet propagator -double)  !
! * read_td_restart_verlet_s (reads TD restart for Verlet propagator -single)  !
! * read_td_restart_magnus_d (reads TD restart for Magnus propagator -double)  !
! * read_td_restart_magnus_s (reads TD restart for Magnus propagator -single)  !
! * write_td_restart_verlet_d (writes TD restart for Verlet propagator -double)!
! * write_td_restart_verlet_s (writes TD restart for Verlet propagator -single)!
! * write_td_restart_magnus_d (writes TD restart for Magnus propagator -double)!
! * write_td_restart_magnus_s (writes TD restart for Magnus propagator -single)!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!


!% READ_TD_RESTART_XX_YY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Reads TD-DFT restarts in single and double precision. Verlet restarts need   !
! only electronic density (Rho), while Magnus needs two Fock matrices (F1a and !
! F1b.                                                                         !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine read_td_restart_verlet_d(rho, M, file_name)
   ! rho       : Complex density matrix.
   ! M         : Density matrix size.
   ! file_name : Name of the restart file.
   implicit none
   character(len=20), intent(in)  :: file_name
   integer          , intent(in)  :: M
   complex*16       , intent(out) :: rho(M,M)
   logical                        :: exists
   integer                        :: UID
   UID = 1550

   inquire(file = file_name,  exist = exists)
   if (.not.exists) then
      write(*,*) 'ERROR: TD Verlet restart not found.'
      write(*,*) 'Check file existence or use tdrestart=false.'
      stop
   endif

   open(unit = UID, file = file_name, status = 'old')
   call read_rho_restart(rho, M, UID)
   close(UID)

   return
end subroutine read_td_restart_verlet_d

subroutine read_td_restart_verlet_s(rho, M, file_name)
   ! rho       : Complex density matrix.
   ! M         : Density matrix size.
   ! file_name : Name of the restart file.
   implicit none
   character(len=20), intent(in)  :: file_name
   integer          , intent(in)  :: M
   complex*8        , intent(out) :: rho(M,M)
   logical                        :: exists
   integer                        :: UID
   UID = 1550

   inquire(file = file_name,  exist = exists)
   if (.not.exists) then
      write(*,*) 'ERROR: TD Verlet restart not found.'
      write(*,*) 'Check file existence or use tdrestart=false.'
      stop
   endif

   open(unit = UID, file = file_name, status = 'old')
   call read_rho_restart(rho, M, UID)
   close(UID)

   return
end subroutine read_td_restart_verlet_s

subroutine read_td_restart_magnus_d(rho, fock_a, fock_b, M, file_name, is_fock)
   ! rho       : Complex density matrix.
   ! M         : Density matrix size.
   ! file_name : Name of the restart file.
   ! UID       : Inputfile unit.
   ! fock_a/b  : Fock matrices for Magnus.

   implicit none
   character(len=20), intent(in)  :: file_name
   integer          , intent(in)  :: M
   complex*16       , intent(out) :: rho(M,M)
   logical          , intent(out) :: is_fock
   real*8           , intent(out) :: fock_a(M,M), fock_b(M,M)
   logical                        :: exists
   integer                        :: UID, file_stat
   UID = 1550

   inquire(file = file_name,  exist = exists)
   if (.not.exists) then
      write(*,*) 'ERROR: TD Magnus restart file ', file_name, ' not found.'
      write(*,*) 'Check file existence or use tdrestart=false.'
      stop
   endif

   open(unit = UID, file = file_name, status = 'old')
   rewind(UID)
   call read_sqmatrix(rho, M, UID)
   read(UID,*,iostat=file_stat)
   if (file_stat .eq. 0) then
      call read_sqmatrix(fock_a, M, UID)
      call read_sqmatrix(fock_b, M, UID)
      is_fock = .true.
   else
      if (IS_IOSTAT_END(file_stat)) then
         write(*,*) '  Fock not found in TD Magnus restart. Making leapfrog.'
         is_fock = .false.
      endif
   endif
   close(UID)

   return
end subroutine read_td_restart_magnus_d

subroutine read_td_restart_magnus_s(rho, fock_a, fock_b, M, file_name, is_fock)
   ! rho       : Complex density matrix.
   ! M         : Density matrix size.
   ! file_name : Name of the restart file.
   ! UID       : Inputfile unit.
   ! fock_a/b  : Fock matrices for Magnus.

   implicit none
   character(len=20), intent(in)  :: file_name
   integer          , intent(in)  :: M
   logical          , intent(out) :: is_fock
   complex*8        , intent(out) :: rho(M,M)
   real*8           , intent(out) :: fock_a(M,M), fock_b(M,M)
   logical                        :: exists
   integer                        :: UID, file_stat
   UID = 1550

   inquire(file = file_name,  exist = exists)
   if (.not.exists) then
      write(*,*) 'ERROR: TD Magnus restart file ', file_name, ' not found.'
      write(*,*) 'Check file existence or use tdrestart=false.'
      stop
   endif

   open(unit = UID, file = file_name, status = 'old')
   rewind(UID)
   call read_sqmatrix(rho, M, UID)
   read(UID,*,iostat=file_stat)
   if (file_stat .eq. 0) then
      call read_sqmatrix(fock_a, M, UID)
      call read_sqmatrix(fock_b, M, UID)
      is_fock = .true.
   else
      if (IS_IOSTAT_END(file_stat)) then
         write(*,*) '  Fock not found in TD Magnus restart. Making leapfrog.'
         is_fock = .false.
      endif
   endif
   close(UID)

   return
end subroutine read_td_restart_magnus_s
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!


!% WRITE_TD_RESTART_XX_YY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Writes TD-DFT restarts in single and double precision. Verlet restarts need  !
! only electronic density (Rho), while Magnus needs two Fock matrices (F1a and !
! F1b.                                                                         !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine write_td_restart_verlet_d(rho, M, file_name)
   ! rho     : Complex density matrix.
   ! M       : Density matrix size.
   ! file_name : Name of the file containing the density matrix.
   implicit none
   character(len=20), intent(in) :: file_name
   integer          , intent(in) :: M
   complex*16       , intent(in) :: rho(M,M)
   integer                       :: UID
   UID = 1550

   open(unit = UID, file = file_name)
   call write_rho_restart(rho, M, UID)
   close(UID)

   return
end subroutine write_td_restart_verlet_d

subroutine write_td_restart_verlet_s(rho, M, file_name)
   ! rho       : Complex density matrix.
   ! M         : Density matrix size.
   ! file_name : Name of the file containing the density matrix.
   implicit none
   character(len=20), intent(in) :: file_name
   integer          , intent(in) :: M
   complex*8        , intent(in) :: rho(M,M)
   integer                       :: UID
   UID = 1550

   open(unit = UID, file = file_name)
   call write_rho_restart(rho, M, UID)
   close(UID)

   return
end subroutine write_td_restart_verlet_s

subroutine write_td_restart_magnus_d(rho, fock_a, fock_b, M, file_name)
   ! rho       : Complex density matrix.
   ! M         : Density matrix size.
   ! file_name : Name of the restart file.
   ! UID       : Inputfile unit.
   ! fock_a/b  : Fock matrices for Magnus.
   implicit none
   character(len=20), intent(in) :: file_name
   integer          , intent(in) :: M
   complex*16       , intent(in) :: rho(M,M)
   real*8           , intent(in) :: fock_a(M,M), fock_b(M,M)
   integer                       :: UID
   UID = 1550

   open(unit = UID, file = file_name)
   rewind(UID)
   call write_sqmatrix(rho, M, UID)
   call write_sqmatrix(fock_a, M, UID)
   call write_sqmatrix(fock_b, M, UID)
   close(UID)

   return
end subroutine write_td_restart_magnus_d

subroutine write_td_restart_magnus_s(rho, fock_a, fock_b, M, file_name)
   ! rho       : Complex density matrix.
   ! M         : Density matrix size.
   ! file_name : Name of the restart file.
   ! UID       : Inputfile unit.
   ! fock_a/b  : Fock matrices for Magnus.
   implicit none
   character(len=20), intent(in) :: file_name
   integer          , intent(in) :: M
   complex*8        , intent(in) :: rho(M,M)
   real*8           , intent(in) :: fock_a(M,M), fock_b(M,M)
   integer                       :: UID
   UID = 1550

   open(unit = UID, file = file_name)
   rewind(UID)
   call write_sqmatrix(rho, M, UID)
   call write_sqmatrix(fock_a, M, UID)
   call write_sqmatrix(fock_b, M, UID)
   close(UID)

   return
end subroutine write_td_restart_magnus_s
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
