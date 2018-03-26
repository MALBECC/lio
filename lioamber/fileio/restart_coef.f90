!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% RESTART_COEF.F90 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! * read_coef_restart_cd  (reads coefficient restart for closed shell - double)!
! * read_coef_restart_cs  (reads coefficient restart for closed shell - single)!
! * read_coef_restart_od  (reads coefficient restart for open shell - double)  !
! * read_coef_restart_os  (reads coefficient restart for open shell - single)  !
! * write_coef_restart_cd (writes coefficient restart for closed shell -double)!
! * write_coef_restart_cs (writes coefficient restart for closed shell -single)!
! * write_coef_restart_od (writes coefficient restart for open shell - double) !
! * write_coef_restart_os (writes coefficient restart for open shell - single) !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!


!% READ_COEF_XX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Reads the Coefficient matrix from restart, rearranging it if needed, for     !
! open and closed shell cases in both single and double precision.             !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine read_coef_restart_cd(coef, dens, M, NCO, UID)
   ! coef : Coefficient matrix.
   ! dens : Density matrix.
   ! M    : Coefficient matrix size.
   ! NCO  : Number of occupied orbitals.
   ! UID  : Input file unit ID.
   implicit none
   integer, intent(in)  :: M, NCO, UID
   real*8 , intent(out) :: coef(M,NCO), dens(M,M)
   integer :: i, j, k

   coef = 0.0D0
   dens = 0.0D0
   rewind(UID)
   call read_matrix(coef, M, NCO, 1, 1, UID)

   do i=1, M
   do j=1, M
   do k=1, NCO
      dens(i,j) = dens(i,j) + 2.0D0*coef(i,k)*coef(j,k)
   enddo
   enddo
   enddo

   return
end subroutine read_coef_restart_cd

subroutine read_coef_restart_cs(coef, dens, M, NCO, UID)
   ! coef : Coefficient matrix.
   ! dens : Density matrix.
   ! M    : Coefficient matrix size.
   ! NCO  : Number of occupied orbitals.
   ! UID  : Input file unit ID.
   implicit none
   integer, intent(in)  :: M, UID, NCO
   real*4 , intent(out) :: coef(M,NCO), dens(M,M)
   integer :: i, j, k

   coef = 0.0D0
   dens = 0.0D0
   rewind(UID)
   call read_matrix(coef, M, NCO, 1, 1, UID)

   do i=1, M
   do j=1, M
   do k=1, NCO
      dens(i,j) = dens(i,j) + 2.0D0*coef(i,k)*coef(j,k)
   enddo
   enddo
   enddo

   return
end subroutine read_coef_restart_cs

subroutine read_coef_restart_od(coef_a, coef_b, dens_t, dens_a, dens_b,        &
                                M, NCOa, NCOb, UID)
   ! coef_a : Coefficient matrix Alpha.
   ! coef_b : Coefficient matrix Beta.
   ! dens_t : Density matrix.
   ! dens_a : Alpha density matrix.
   ! dens_b : Beta density matrix.
   ! M      : Coefficient matrix size.
   ! NCOa   : Number of occupied Alpha orbitals.
   ! NCOb   : Number of occupied Beta orbitals.
   ! UID    : Input file unit ID.
   implicit none
   integer, intent(in)  :: M, NCOa, NCOb, UID
   real*8 , intent(out) :: coef_a(M,NCOa), coef_b(M,NCOb), dens_t(M,M),        &
                           dens_a(M,M), dens_b(M,M)
   integer :: i, j, k

   coef_a = 0.0D0
   coef_b = 0.0D0
   dens_t = 0.0D0
   rewind(UID)
   call read_matrix(coef_a, M, NCOa, 1, 1, UID)
   call read_matrix(coef_b, M, NCOb, 1, 1, UID)
   do i=1, M
   do j=1, M
      do k=1, NCOa
         dens_t(i,j) = dens_t(i,j) + coef_a(i,k)*coef_a(j,k)
         dens_a(i,j) = dens_a(i,j) + coef_a(i,k)*coef_a(j,k)
      enddo
      do k=1, NCOb
         dens_t(i,j) = dens_t(i,j) + coef_b(i,k)*coef_b(j,k)
         dens_b(i,j) = dens_b(i,j) + coef_b(i,k)*coef_b(j,k)
      enddo
   enddo
   enddo

   return
end subroutine read_coef_restart_od

subroutine read_coef_restart_os(coef_a, coef_b, dens_t, dens_a, dens_b,        &
                                M, NCOa, NCOb, UID)
   ! coef_a : Coefficient matrix Alpha.
   ! coef_b : Coefficient matrix Beta.
   ! dens_t : Density matrix.
   ! dens_a : Alpha density matrix.
   ! dens_b : Beta density matrix.
   ! M      : Coefficient matrix size.
   ! NCOa   : Number of occupied Alpha orbitals.
   ! NCOb   : Number of occupied Beta orbitals.
   ! UID    : Input file unit ID.
   implicit none
   integer, intent(in)  :: M, NCOa, NCOb, UID
   real*4 , intent(out) :: coef_a(M,NCOa), coef_b(M,NCOb), dens_t(M,M),        &
                           dens_a(M,M), dens_b(M,M)
   integer :: i, j, k

   coef_a = 0.0D0
   coef_b = 0.0D0
   dens_t = 0.0D0
   rewind(UID)
   call read_matrix(coef_a, M, NCOa, 1, 1, UID)
   call read_matrix(coef_b, M, NCOb, 1, 1, UID)

   do i=1, M
   do j=1, M
      do k=1, NCOa
         dens_t(i,j) = dens_t(i,j) + coef_a(i,k)*coef_a(j,k)
         dens_a(i,j) = dens_a(i,j) + coef_a(i,k)*coef_a(j,k)
      enddo
      do k=1, NCOb
         dens_t(i,j) = dens_t(i,j) + coef_b(i,k)*coef_b(j,k)
         dens_b(i,j) = dens_b(i,j) + coef_b(i,k)*coef_b(j,k)
      enddo
   enddo
   enddo

   return
end subroutine read_coef_restart_os
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!


!% WRITE_COEF_XX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Rearranges and prints Coefficient matrix for open and closed shell cases in  !
! both single and double precision.                                            !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine write_coef_restart_cd(coef, M, NCO, UID)
   ! coef : Coefficient matrix.
   ! M    : Coefficient matrix size.
   ! NCO  : Number of occupied orbitals.
   ! UID  : Output file unit ID.
   implicit none
   integer, intent(in) :: M, NCO, UID
   real*8 , intent(in) :: coef(M,NCO)

   rewind(UID)
   call write_matrix(coef, M, NCO, 1, 1, UID)
   return
end subroutine write_coef_restart_cd

subroutine write_coef_restart_cs(coef, M, NCO, UID)
   ! coef : Coefficient matrix.
   ! M    : Coefficient matrix size.
   ! NCO  : Number of occupied orbitals.
   ! UID  : Output file unit ID.
   implicit none
   integer, intent(in) :: M, NCO, UID
   real*4 , intent(in) :: coef(M,NCO)

   rewind(UID)
   call write_matrix(coef, M, NCO, 1, 1, UID)
   return
end subroutine write_coef_restart_cs

subroutine write_coef_restart_od(coef_a, coef_b, M, NCOa, NCOb, UID)
   ! coef_a : Coefficient matrix Alpha.
   ! coef_b : Coefficient matrix Beta.
   ! M      : Coefficient matrix size.
   ! NCOa   : Number of occupied Alpha orbitals.
   ! NCOb   : Number of occupied Beta orbitals.
   ! UID    : Output file unit ID.
   implicit none
   integer, intent(in) :: M, NCOa, NCOb, UID
   real*8 , intent(in) :: coef_a(M,NCOa), coef_b(M,NCOb)

   rewind(UID)
   call write_matrix(coef_a, M, NCOa, 1, 1, UID)
   call write_matrix(coef_b, M, NCOb, 1, 1, UID)
   return
end subroutine write_coef_restart_od

subroutine write_coef_restart_os(coef_a, coef_b, M, NCOa, NCOb, UID)
   ! coef_a : Coefficient matrix Alpha.
   ! coef_b : Coefficient matrix Beta.
   ! M      : Coefficient matrix size.
   ! NCOa   : Number of occupied Alpha orbitals.
   ! NCOb   : Number of occupied Beta orbitals.
   ! UID    : Output file unit ID.
   implicit none
   integer, intent(in) :: M, NCOa, NCOb, UID
   real*4 , intent(in) :: coef_a(M,NCOa), coef_b(M,NCOb)

   rewind(UID)
   call write_matrix(coef_a, M, NCOa, 1, 1, UID)
   call write_matrix(coef_b, M, NCOb, 1, 1, UID)
   return
end subroutine write_coef_restart_os
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
