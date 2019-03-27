!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Build Density Tool                                                           !
! This tool takes two open-shell density restart files and adds them together  !
! for either a new joint restart file (joint_dens) or restart files for time-  !
! dependent DFT.                                                               !
!                                                                              !
! The options file should contain the following input in a &bdens namelist:    !
!   M (integer)         : The number of basis functions.                       !
!   file_rho_a (char*40): File containing the density matrix of the first      !
!                         system.                                              !
!   file_rho_b (char*40): File containing the density matrix of the second     !
!                         system.                                              !
!   invert (logical)    : For the second system (b), exchange Alpha with Beta  !
!                         densities before adding them to the first system (a).!
!                                                                              !
! Usage:                                                                       !
!   ./bdens options_file                                                       !
!                                                                              !
! March 2019 - F. Pedron                                                       !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

program build_dens

   implicit none
   integer           :: M = 0
   logical           :: invert_b = .false.
   character(len=40) :: file_rho_a, file_rho_b, file_input
   real(kind=8), allocatable :: rho_alpha_a(:,:), rho_beta_a(:,:), &
                                rho_alpha_b(:,:), rho_beta_b(:,:), &
                                rho_alpha_t(:,:), rho_beta_t(:,:)

   call get_command_argument(1, file_input)
   call read_input(M, file_rho_a, file_rho_b, invert_b, trim(file_input))

   allocate(rho_alpha_a(M, M), rho_beta_a(M, M))
   allocate(rho_alpha_b(M, M), rho_beta_b(M, M))

   call read_dens(M,rho_alpha_a, rho_beta_a, rho_alpha_b, &
                  rho_beta_b, file_rho_a, file_rho_b)
   call calc_trace(rho_alpha_a, M, "rho_alpha_a")
   call calc_trace(rho_beta_a , M, "rho_beta_a")
   call calc_trace(rho_alpha_b, M, "rho_alpha_b")
   call calc_trace(rho_beta_b , M, "rho_beta_b")

   allocate(rho_alpha_t(M, M), rho_beta_t(M, M))

   call join_dens(rho_alpha_a, rho_beta_a, rho_alpha_b, rho_beta_b, &
                  rho_alpha_t, rho_beta_t,  M, invert_b)
   call calc_trace(rho_alpha_t, M, "rho_alpha_t")
   call calc_trace(rho_beta_t , M, "rho_beta_t")

   call write_dens(rho_alpha_t, rho_beta_t, M)

   deallocate(rho_alpha_a, rho_beta_a, rho_alpha_b, rho_beta_b)
   deallocate(rho_alpha_t, rho_beta_t)
   return
contains

   ! Reads input options and parameters.
   subroutine read_input(M, file_rho_a, file_rho_b, invert, file_in)
      implicit none
      character(len=*) , intent(in)    :: file_in
      integer          , intent(inout) :: M
      logical          , intent(inout) :: invert
      character(len=40), intent(inout) :: file_rho_a, file_rho_b

      namelist /bdens/ M, file_rho_a, file_rho_b, invert
      open(unit = 100, file = file_in)
      read(unit = 100, nml = bdens)
      close(100)
   end subroutine read_input

   ! Reads input densities.
   subroutine read_dens(M, rho_a_a, rho_b_a, rho_a_b, &
                        rho_b_b, file_rho_a, file_rho_b)
      implicit none
      integer          , intent(in)  :: M
      character(len=40), intent(in)  :: file_rho_a, file_rho_b
      real(kind=8)     , intent(out) :: rho_a_a(M, M), rho_b_a(M, M), &
                                        rho_a_b(M, M), rho_b_b(M, M)
      integer :: icount, jcount

      open(unit = 100, file = file_rho_a)
      do icount = 1, M
         read(100,*) (rho_a_a(icount, jcount), jcount = 1, M)
      enddo
      do icount = 1, M
         read(100,*) (rho_b_a(icount, jcount), jcount = 1, M)
      enddo
      close(100)

      open(unit = 100, file = file_rho_b)
      do icount = 1, M
         read(100,*) (rho_a_b(icount, jcount), jcount = 1, M)
      enddo
      do icount = 1, M
         read(100,*) (rho_b_b(icount, jcount), jcount = 1, M)
      enddo

      close(100)

      ! Small correction for non-diagonal terms in LIO
      do icount = 1       , M
      do jcount = icount+1, M
         rho_a_a(icount, jcount) = rho_a_a(icount, jcount) / 2.0D0
         rho_b_a(icount, jcount) = rho_b_a(icount, jcount) / 2.0D0
         rho_a_b(icount, jcount) = rho_a_b(icount, jcount) / 2.0D0
         rho_b_b(icount, jcount) = rho_b_b(icount, jcount) / 2.0D0

         rho_a_a(jcount, icount) = rho_a_a(icount, jcount)
         rho_b_a(jcount, icount) = rho_b_a(icount, jcount)
         rho_a_b(jcount, icount) = rho_a_b(icount, jcount)
         rho_b_b(jcount, icount) = rho_b_b(icount, jcount)
      enddo
      enddo
   end subroutine read_dens

   subroutine write_dens(rho_a_t, rho_b_t, M)
      implicit none
      integer         , intent(in)    :: M
      double precision, intent(inout) :: rho_a_t(M, M),&
                                      rho_b_t(M, M)
      integer :: icount, jcount

      open(unit = 100, file = "td_a_in.restart")
      do icount = 1, M
      do jcount = 1, M
         write(100,*) rho_a_t(icount, jcount), "0.0D0"
      enddo
      enddo
      close(100)

      open(unit = 100, file = "td_b_in.restart")
      do icount = 1, M
      do jcount = 1, M
         write(100,*) rho_b_t(icount, jcount), "0.0D0"
      enddo
      enddo
      close(100)


      ! Fixes non-diagonal elements for output density matrix.
      do icount = 1       , M
      do jcount = icount+1, M
         rho_a_t(icount, jcount) = rho_a_t(icount, jcount) * 2.0D0
         rho_b_t(icount, jcount) = rho_b_t(icount, jcount) * 2.0D0

         rho_a_t(jcount, icount) = rho_a_t(icount, jcount)
         rho_b_t(jcount, icount) = rho_b_t(icount, jcount)
      enddo
      enddo

      open(unit = 100, file = "joint-dens.out")
      do icount = 1, M
         write(100,*) (rho_a_t(icount, jcount), jcount = 1, M)
      enddo
      do icount = 1, M
         write(100,*) (rho_b_t(icount, jcount), jcount = 1, M)
      enddo
      close(100)
   end subroutine write_dens

   ! Joins densities in a new one.
   subroutine join_dens(rho_aa, rho_ba, rho_ab, rho_bb, rho_at, rho_bt, M, &
                        invert)
      implicit none
      integer         , intent(in)  :: M
      logical         , intent(in)  :: invert
      double precision, intent(in)  :: rho_aa(M,M), rho_ba(M,M), &
                                       rho_ab(M,M), rho_bb(M,M)
      double precision, intent(out) :: rho_at(M,M), rho_bt(M,M)
      integer :: Mt, icount, jcount
      
      rho_at = 0.0D0
      rho_bt = 0.0D0

      if (invert) then
         do icount = 1, M
         do jcount = 1, M
            rho_at(icount, jcount) = rho_aa(icount, jcount) +&
                                     rho_bb(icount, jcount)
            rho_bt(icount, jcount) = rho_ba(icount, jcount) +&
                                     rho_ab(icount, jcount)
         enddo
         enddo
      else
         do icount = 1, M
         do jcount = 1, M
            rho_at(icount, jcount) = rho_aa(icount, jcount) +&
                                     rho_ab(icount, jcount)
            rho_bt(icount, jcount) = rho_ba(icount, jcount) +&
                                     rho_bb(icount, jcount)
         enddo
         enddo
      endif
   end subroutine join_dens

   subroutine calc_trace(matrix, msize, message)

      implicit none
      integer         , intent(in) :: msize
      character(len=*), intent(in) :: message
      double precision, intent(in) :: matrix(msize, msize)

      integer          :: icount
      double precision :: trace

      trace = 0.0D0
      do icount = 1, msize
         trace = trace + matrix(icount, icount)
      enddo

      write(*,*) "Trace of ", trim(message), " equals to ", trace
   end subroutine calc_trace


end program build_dens
