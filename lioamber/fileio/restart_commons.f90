!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% RESTART_COMMONS.F90 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! This file contains common read write routines for restart calls.             !
! File read routines:                                                    [RDR] !
! * read_matrix_cd                                                             !
! * read_matrix_cs                                                             !
! * read_matrix_d                                                              !
! * read_matrix_s                                                              !
! * read_sqmatrix_cd                                                           !
! * read_sqmatrix_cs                                                           !
! * read_sqmatrix_d                                                            !
! * read_sqmatrix_s                                                            !
! File write routines:                                                   [WRR] !
! * write_matrix_cd                                                            !
! * write_matrix_cs                                                            !
! * write_matrix_d                                                             !
! * write_matrix_s                                                             !
! * write_sqmatrix_cd                                                          !
! * write_sqmatrix_cs                                                          !
! * write_sqmatrix_d                                                           !
! * write_sqmatrix_s                                                           !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!


!%% READ_MATRIX_XX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Reads a matrix from an input file, where XX indicates the type: single       !
! precision real (s), double precision real (d), single precision complex (cs) !
! and double precision complex (cd).                                           !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine read_matrix_cd(matrix, size_a, size_b, start_a, start_b, UID)
   implicit none
   integer   , intent(in)  :: size_a, size_b, start_a, start_b, UID
   complex*16, intent(out) :: matrix(size_a, size_b)
   real*8                  :: temp_a, temp_b
   integer                 :: icount, jcount

   do icount = start_a, size_a
   do jcount = start_b, size_b
      read(UID,*) temp_a, temp_b
      matrix(icount, jcount) = cmplx(temp_a, temp_b)
   enddo
   enddo

   return
end subroutine read_matrix_cd

subroutine read_matrix_cs(matrix, size_a, size_b, start_a, start_b, UID)
   implicit none
   integer   , intent(in)  :: size_a, size_b, start_a, start_b, UID
   complex*8 , intent(out) :: matrix(size_a, size_b)
   real*4                  :: temp_a, temp_b
   integer                 :: icount, jcount

   do icount = start_a, size_a
   do jcount = start_b, size_b
      read(UID,*) temp_a, temp_b
      matrix(icount, jcount) = cmplx(temp_a, temp_b)
   enddo
   enddo

   return
end subroutine read_matrix_cs

subroutine read_matrix_d(matrix, size_a, size_b, start_a, start_b, UID)
   implicit none
   integer, intent(in)  :: size_a, size_b, start_a, start_b, UID
   real*8 , intent(out) :: matrix(size_a, size_b)
   integer              :: icount, jcount

   do icount = start_a, size_a
      read(UID,*) (matrix(icount, jcount), jcount = start_b, size_b)
   enddo

   return
end subroutine read_matrix_d

subroutine read_matrix_s(matrix, size_a, size_b, start_a, start_b, UID)
   implicit none
   integer, intent(in)  :: size_a, size_b, start_a, start_b, UID
   real*4 , intent(out) :: matrix(size_a, size_b)
   integer              :: icount, jcount

   do icount = start_a, size_a
      read(UID,*) (matrix(icount, jcount), jcount = start_b, size_b)
   enddo

   return
end subroutine read_matrix_s
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!


!%% READ_SQMATRIX_XX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Reads a square matrix from an input file, where XX indicates the type: single!
! precision real (s), double precision real (d), single precision complex (cs) !
! and double precision complex (cd).                                           !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine read_sqmatrix_cd(sqmatrix, M, UID)
   ! sqmatrix : Readable matrix.
   ! M        : Matrix dimension.
   ! UID      : Output file unit.
   implicit none
   integer   , intent(in)  :: M, UID
   complex*16, intent(out) :: sqmatrix(M, M)

   call read_matrix_cd(sqmatrix, M, M, 1, 1, UID)
   return
end subroutine read_sqmatrix_cd

subroutine read_sqmatrix_cs(sqmatrix, M, UID)
   ! sqmatrix : Readable matrix.
   ! M        : Matrix dimension.
   ! UID      : Output file unit.
   implicit none
   integer  , intent(in)  :: M, UID
   complex*8, intent(out) :: sqmatrix(M, M)

   call read_matrix_cs(sqmatrix, M, M, 1, 1, UID)
   return
end subroutine read_sqmatrix_cs

subroutine read_sqmatrix_d(sqmatrix, M, UID)
   ! sqmatrix : Readable matrix.
   ! M        : Matrix dimension.
   ! UID      : Output file unit.
   implicit none
   integer, intent(in)  :: M, UID
   real*8 , intent(out) :: sqmatrix(M, M)

   call read_matrix_d(sqmatrix, M, M, 1, 1, UID)
   return
end subroutine read_sqmatrix_d

subroutine read_sqmatrix_s(sqmatrix, M, UID)
   ! sqmatrix : Readable matrix.
   ! M        : Matrix dimension.
   ! UID      : Output file unit.
   implicit none
   integer, intent(in)  :: M, UID
   real*4 , intent(out) :: sqmatrix(M, M)

   call read_matrix_s(sqmatrix, M, M, 1, 1, UID)
   return
end subroutine read_sqmatrix_s
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!


!%% WRITE_MATRIX_XX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Prints a matrix to an output file, where XX indicates the type: single       !
! precision real (s), double precision real (d), single precision complex (cs) !
! and double precision complex (cd).                                           !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine write_matrix_cd(matrix, size_a, size_b, start_a, start_b, UID)
   implicit none
   integer   , intent(in) :: size_a, size_b, start_a, start_b, UID
   complex*16, intent(in) :: matrix(size_a, size_b)
   integer                :: icount, jcount

   do icount = start_a, size_a
   do jcount = start_b, size_b
      write(UID, 400) dble(matrix(icount, jcount)), dimag(matrix(icount, jcount))
   enddo
   enddo

   return
400 format(4(E14.7E2, 2x))
end subroutine write_matrix_cd

subroutine write_matrix_cs(matrix, size_a, size_b, start_a, start_b, UID)
   implicit none
   integer  , intent(in) :: size_a, size_b, start_a, start_b, UID
   complex*8, intent(in) :: matrix(size_a, size_b)
   integer               :: icount, jcount

   do icount = start_a, size_a
   do jcount = start_b, size_b
      write(UID, 400) real(matrix(icount, jcount)), aimag(matrix(icount, jcount))
   enddo
   enddo

   return
400 format(4(E14.7E2, 2x))
end subroutine write_matrix_cs

subroutine write_matrix_d(matrix, size_a, size_b, start_a, start_b, UID)
   implicit none
   integer, intent(in) :: size_a, size_b, start_a, start_b, UID
   real*8 , intent(in) :: matrix(size_a, size_b)
   integer             :: icount, jcount

   do icount = start_a, size_a
      write(UID,*) (matrix(icount, jcount), jcount = start_b, size_b)
   enddo

   return
400 format(4(E14.7E2, 2x))
end subroutine write_matrix_d

subroutine write_matrix_s(matrix, size_a, size_b, start_a, start_b, UID)
   implicit none
   integer, intent(in) :: size_a, size_b, start_a, start_b, UID
   real*4 , intent(in) :: matrix(size_a, size_b)
   integer             :: icount, jcount

   do icount = start_a, size_a
      write(UID,*) (matrix(icount, jcount), jcount = start_b, size_b)
   enddo

   return
400 format(4(E14.7E2, 2x))
end subroutine write_matrix_s
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!


!%% WRITE_SQMATRIX_XX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Prints a square matrix to an output file, where XX indicates the type: single!
! precision real (s), double precision real (d), single precision complex (cs) !
! and double precision complex (cd).                                           !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine write_sqmatrix_cd(sqmatrix, M, UID)
   ! sqmatrix : Printable matrix.
   ! M        : Matrix dimension.
   ! UID      : Output file unit.
   implicit none
   integer   , intent(in) :: M, UID
   complex*16, intent(in) :: sqmatrix(M, M)
   integer :: icount, jcount

   call write_matrix_cd(sqmatrix, M, M, 1, 1, UID)
   return
end subroutine write_sqmatrix_cd

subroutine write_sqmatrix_cs(sqmatrix, M, UID)
   ! sqmatrix : Printable matrix.
   ! M        : Matrix dimension.
   ! UID      : Output file unit.
   implicit none
   integer  , intent(in) :: M, UID
   complex*8, intent(in) :: sqmatrix(M,M)
   integer :: icount, jcount

   call write_matrix_cs(sqmatrix, M, M, 1, 1, UID)
   return
end subroutine write_sqmatrix_cs

subroutine write_sqmatrix_d(sqmatrix, M, UID)
   ! sqmatrix : Printable matrix.
   ! M        : Matrix dimension.
   ! UID      : Output file unit.
   implicit none
   integer, intent(in) :: M, UID
   real*8 , intent(in) :: sqmatrix(M, M)
   integer :: icount, jcount

   call write_matrix_d(sqmatrix, M, M, 1, 1, UID)
   return
end subroutine write_sqmatrix_d

subroutine write_sqmatrix_s(sqmatrix, M, UID)
   ! sqmatrix : Printable matrix.
   ! M        : Matrix dimension.
   ! UID      : Output file unit.
   implicit none
   integer, intent(in) :: M, UID
   real*4 , intent(in) :: sqmatrix(M,M)
   integer :: icount, jcount

   call write_matrix_s(sqmatrix, M, M, 1, 1, UID)
   return
end subroutine write_sqmatrix_s
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
