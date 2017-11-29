!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% RESTART_FOCK.F90 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! * read_fock_restart_cd   (reads Fock restart for closed shell - double)      !
! * read_fock_restart_cs   (reads Fock restart for closed shell - single)      !
! * read_fock_restart_od   (reads Fock restart for open shell - double)        !
! * read_fock_restart_os   (reads Fock restart for open shell - single)        !
! * write_fock_restart_cd  (writes Fock restart for closed shell - double)     !
! * write_fock_restart_cs  (writes Fock restart for closed shell - single)     !
! * write_fock_restart_od  (writes Fock restart for open shell - double)       !
! * write_fock_restart_os  (writes Fock restart for open shell - single)       !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!


!% READ_FOCK_XX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Reads the Fock matrix from restart, rearranging it if needed, for open and   !
! closed shell cases in both single and double precision.                      !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine read_fock_restart_cd(fock, M, UID)
   ! fock : Fock matrix.
   ! M    : Fock matrix size.
   ! UID  : Input file unit ID.
   implicit none
   integer, intent(in)  :: M, UID
   real*8 , intent(out) :: fock(M,M)
 
   rewind(UID)
   call read_sqmatrix(fock, M, UID)
   return
end subroutine read_fock_restart_cd

subroutine read_fock_restart_cs(fock, M, UID)
   ! fock : Fock matrix.
   ! M    : Fock matrix size.
   ! UID  : Input file unit ID.
   implicit none
   integer, intent(in)  :: M, UID
   real*4 , intent(out) :: fock(M,M)

   rewind(UID)
   call read_sqmatrix(fock, M, UID)
   return
end subroutine read_fock_restart_cs

subroutine read_fock_restart_od(fock_a, fock_b, M, UID)
   ! fock_a : Fock matrix Alpha.
   ! fock_b : Fock matrix Beta.
   ! M      : Fock matrix size.
   ! UID    : Input file unit ID.
   implicit none
   integer, intent(in)  :: M, UID
   real*8 , intent(out) :: fock_a(M,M), fock_b(M,M)

   rewind(UID)
   call read_sqmatrix(fock_a, M, UID)
   call read_sqmatrix(fock_b, M, UID)
   return
end subroutine read_fock_restart_od

subroutine read_fock_restart_os(fock_a, fock_b, M, UID)
   ! fock_a : Fock matrix Alpha.
   ! fock_b : Fock matrix Beta.
   ! M      : Fock matrix size.
   ! UID    : Input file unit ID.
   implicit none
   integer, intent(in)  :: M, UID
   real*4 , intent(out) :: fock_a(M,M), fock_b(M,M)

   rewind(UID)
   call read_sqmatrix(fock_a, M, UID)
   call read_sqmatrix(fock_b, M, UID)
   return
end subroutine read_fock_restart_os
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!


!% WRITE_FOCK_XX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Rearranges and prints Fock matrix for open and closed shell cases in both    !
! single and double precision.                                                 !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine write_fock_restart_cd(fock, M, UID)
   ! fock : Fock matrix.
   ! M    : Fock matrix size.
   ! UID  : Output file unit ID.
   implicit none
   integer, intent(in) :: M, UID
   real*8 , intent(in) :: fock(M,M)

   rewind(UID)
   call write_sqmatrix(fock, M, UID)
   return
end subroutine write_fock_restart_cd

subroutine write_fock_restart_cs(fock, M, UID)
   ! fock : Fock matrix.
   ! M    : Fock matrix size.
   ! UID  : Output file unit ID.
   implicit none
   integer, intent(in) :: M, UID
   real*4 , intent(in) :: fock(M,M)

   rewind(UID)
   call write_sqmatrix(fock, M, UID)
   return
end subroutine write_fock_restart_cs

subroutine write_fock_restart_od(fock_a, fock_b, M, UID)
   ! fock_a : Fock matrix Alpha.
   ! fock_b : Fock matrix Beta.
   ! M      : Fock matrix size.
   ! UID    : Output file unit ID.
   implicit none
   integer, intent(in) :: M, UID
   real*8 , intent(in) :: fock_a(M,M), fock_b(M,M)

   rewind(UID)
   call write_sqmatrix(fock_a, M, UID)
   call write_sqmatrix(fock_b, M, UID)
   return
end subroutine write_fock_restart_od

subroutine write_fock_restart_os(fock_a, fock_b, M, UID)
   ! fock_a : Fock matrix Alpha.
   ! fock_b : Fock matrix Beta.
   ! M      : Fock matrix size.
   ! UID    : Output file unit ID.
   implicit none
   integer, intent(in) :: M, UID
   real*4 , intent(in) :: fock_a(M,M), fock_b(M,M)

   rewind(UID)
   call write_sqmatrix(fock_a, M, UID)
   call write_sqmatrix(fock_b, M, UID)
   return
end subroutine write_fock_restart_os
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
