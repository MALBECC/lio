!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% RESTART_RHO.F90 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! * read_rho_restart_ccd  (reads Rho restart for closed shell - cmplx double)  !
! * read_rho_restart_ccs  (reads Rho restart for closed shell - cmplx single)  !
! * read_rho_restart_cd   (reads Rho restart for closed shell - double)        !
! * read_rho_restart_cs   (reads Rho restart for closed shell - single)        !
! * read_rho_restart_ocd  (reads Rho restart for open shell - cmplx double)    !
! * read_rho_restart_ocs  (reads Rho restart for open shell - cmplx single)    !
! * read_rho_restart_od   (reads Rho restart for open shell - double)          !
! * read_rho_restart_os   (reads Rho restart for open shell - single)          !
! * write_rho_restart_ccd (writes Rho restart for closed shell - cmplx double) !
! * write_rho_restart_ccs (writes Rho restart for closed shell - cmplx single) !
! * write_rho_restart_cd  (writes Rho restart for closed shell - double)       !
! * write_rho_restart_cs  (writes Rho restart for closed shell - single)       !
! * write_rho_restart_ocd (writes Rho restart for open shell - cmplx double)   !
! * write_rho_restart_ocs (writes Rho restart for open shell - cmplx single)   !
! * write_rho_restart_od  (writes Rho restart for open shell - double)         !
! * write_rho_restart_os  (writes Rho restart for open shell - single)         !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!


!% READ_RHO_XX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Reads the Rho matrix from restart, rearranging it if needed, for open and    !
! closed shell cases in both single and double precision, real and complex.    !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine read_rho_restart_ccd(rho, M, UID)
   ! rho : Rho matrix.
   ! M   : Rho matrix size.
   ! UID : Input file unit ID.
   implicit none
   integer   , intent(in)  :: M, UID
   complex*16, intent(out) :: rho(M,M)

   rewind(UID)
   call read_sqmatrix(rho, M, UID)
   return
end subroutine read_rho_restart_ccd

subroutine read_rho_restart_ccs(rho, M, UID)
   ! rho : Rho matrix.
   ! M   : Rho matrix size.
   ! UID : Input file unit ID.
   implicit none
   integer  , intent(in)  :: M, UID
   complex*8, intent(out) :: rho(M,M)

   rewind(UID)
   call read_sqmatrix(rho, M, UID)
   return
end subroutine read_rho_restart_ccs

subroutine read_rho_restart_cd(rho, M, UID)
   ! rho : Rho matrix.
   ! M   : Rho matrix size.
   ! UID : Input file unit ID.
   implicit none
   integer, intent(in)  :: M, UID
   real*8 , intent(out) :: rho(M,M)
 
   rewind(UID)
   call read_sqmatrix(rho, M, UID)
   return
end subroutine read_rho_restart_cd

subroutine read_rho_restart_cs(rho, M, UID)
   ! rho : Rho matrix.
   ! M   : Rho matrix size.
   ! UID : Input file unit ID.
   implicit none
   integer, intent(in)  :: M, UID
   real*4 , intent(out) :: rho(M,M)

   rewind(UID)
   call read_sqmatrix(rho, M, UID)
   return
end subroutine read_rho_restart_cs

subroutine read_rho_restart_ocd(rho_a, rho_b, M, UID)
   ! rho_a : Rho matrix Alpha.
   ! rho_b : Rho matrix Beta.
   ! M     : Rho matrix size.
   ! UID   : Input file unit ID.
   implicit none
   integer   , intent(in)  :: M, UID
   complex*16, intent(out) :: rho_a(M,M), rho_b(M,M)

   rewind(UID)
   call read_sqmatrix(rho_a, M, UID)
   call read_sqmatrix(rho_b, M, UID)
   return
end subroutine read_rho_restart_ocd

subroutine read_rho_restart_ocs(rho_a, rho_b, M, UID)
   ! rho_a : Rho matrix Alpha.
   ! rho_b : Rho matrix Beta.
   ! M     : Rho matrix size.
   ! UID   : Input file unit ID.
   implicit none
   integer  , intent(in)  :: M, UID
   complex*8, intent(out) :: rho_a(M,M), rho_b(M,M)

   rewind(UID)
   call read_sqmatrix(rho_a, M, UID)
   call read_sqmatrix(rho_b, M, UID)
   return
end subroutine read_rho_restart_ocs

subroutine read_rho_restart_od(rho_a, rho_b, M, UID)
   ! rho_a : Rho matrix Alpha.
   ! rho_b : Rho matrix Beta.
   ! M     : Rho matrix size.
   ! UID   : Input file unit ID.
   implicit none
   integer, intent(in)  :: M, UID
   real*8 , intent(out) :: rho_a(M,M), rho_b(M,M)

   rewind(UID)
   call read_sqmatrix(rho_a, M, UID)
   call read_sqmatrix(rho_b, M, UID)
   return
end subroutine read_rho_restart_od

subroutine read_rho_restart_os(rho_a, rho_b, M, UID)
   ! rho_a : Rho matrix Alpha.
   ! rho_b : Rho matrix Beta.
   ! M     : Rho matrix size.
   ! UID   : Input file unit ID.
   implicit none
   integer, intent(in)  :: M, UID
   real*4 , intent(out) :: rho_a(M,M), rho_b(M,M)

   rewind(UID)
   call read_sqmatrix(rho_a, M, UID)
   call read_sqmatrix(rho_b, M, UID)
   return
end subroutine read_rho_restart_os
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!


!% WRITE_RHO_XX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Rearranges and prints Rho matrix for open and closed shell cases in both     !
! single and double precision, real and complex                                !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine write_rho_restart_ccd(rho, M, UID)
   ! rho : Rho matrix.
   ! M   : Rho matrix size.
   ! UID : Output file unit ID.
   implicit none
   integer   , intent(in) :: M, UID
   complex*16, intent(in) :: rho(M,M)

   rewind(UID)
   call write_sqmatrix(rho, M, UID)
   return
end subroutine write_rho_restart_ccd

subroutine write_rho_restart_ccs(rho, M, UID)
   ! rho : Rho matrix.
   ! M   : Rho matrix size.
   ! UID : Output file unit ID.
   implicit none
   integer  , intent(in) :: M, UID 
   complex*8, intent(in) :: rho(M,M)

   rewind(UID)
   call write_sqmatrix(rho, M, UID)
   return
end subroutine write_rho_restart_ccs

subroutine write_rho_restart_cd(rho, M, UID)
   ! rho : Rho matrix.
   ! M   : Rho matrix size.
   ! UID : Output file unit ID.
   implicit none
   integer, intent(in) :: M, UID
   real*8 , intent(in) :: rho(M,M)

   rewind(UID)
   call write_sqmatrix(rho, M, UID)
   return
end subroutine write_rho_restart_cd

subroutine write_rho_restart_cs(rho, M, UID)
   ! rho : Rho matrix.
   ! M   : Rho matrix size.
   ! UID : Output file unit ID.
   implicit none
   integer, intent(in) :: M, UID
   real*4 , intent(in) :: rho(M,M)

   rewind(UID)
   call write_sqmatrix(rho, M, UID)
   return
end subroutine write_rho_restart_cs

subroutine write_rho_restart_ocd(rho_a, rho_b, M, UID)
   ! rho_a : Rho matrix Alpha.
   ! rho_b : Rho matrix Beta.
   ! M     : Rho matrix size.
   ! UID   : Output file unit ID.
   implicit none
   integer   , intent(in) :: M, UID
   complex*16, intent(in) :: rho_a(M,M), rho_b(M,M)

   rewind(UID)
   call write_sqmatrix(rho_a, M, UID)
   call write_sqmatrix(rho_b, M, UID)
   return
end subroutine write_rho_restart_ocd

subroutine write_rho_restart_ocs(rho_a, rho_b, M, UID)
   ! rho_a : Rho matrix Alpha.
   ! rho_b : Rho matrix Beta.
   ! M     : Rho matrix size.
   ! UID   : Output file unit ID.
   implicit none
   integer  , intent(in) :: M, UID
   complex*8, intent(in) :: rho_a(M,M), rho_b(M,M)

   rewind(UID)
   call write_sqmatrix(rho_a, M, UID)
   call write_sqmatrix(rho_b, M, UID)
   return
end subroutine write_rho_restart_ocs

subroutine write_rho_restart_od(rho_a, rho_b, M, UID)
   ! rho_a : Rho matrix Alpha.
   ! rho_b : Rho matrix Beta.
   ! M     : Rho matrix size.
   ! UID   : Output file unit ID.
   implicit none
   integer, intent(in) :: M, UID
   real*8 , intent(in) :: rho_a(M,M), rho_b(M,M)

   rewind(UID)
   call write_sqmatrix(rho_a, M, UID)
   call write_sqmatrix(rho_b, M, UID)
   return
end subroutine write_rho_restart_od

subroutine write_rho_restart_os(rho_a, rho_b, M, UID)
   ! rho_a : Rho matrix Alpha.
   ! rho_b : Rho matrix Beta.
   ! M     : Rho matrix size.
   ! UID   : Output file unit ID.
   implicit none
   integer, intent(in) :: M, UID
   real*4 , intent(in) :: rho_a(M,M), rho_b(M,M)

   rewind(UID)
   call write_sqmatrix(rho_a, M, UID)
   call write_sqmatrix(rho_b, M, UID)
   return
end subroutine write_rho_restart_os
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
