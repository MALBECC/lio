!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module fockbias_subs
!------------------------------------------------------------------------------!
!
!    This module controls the aplication of charge biases that are directly
! applied through the fock matrix. It is implementation specific and consists
! of the following subroutines:
!
! fockbias_setorb: This subroutine sets up the atomic charges.
!
! fockbias_setmat: This subroutine sets up the bias matrix (full amplitude).
!                ( changes if atomic positions change )
!
! fockbias_loads: subroutine that reads the atomic biases as input list.
!               ( and calls setorb )
!
! fockbias_apply: Takes the fock matrix and adds the fockbias term (with the
!                 corresponding time shape).
!
! REFERENCE: Physical Review B 74, 155112, 2006
!
!------------------------------------------------------------------------------!
   implicit none

   interface fockbias_apply
      module procedure fockbias_apply_d
      module procedure fockbias_apply_c
      module procedure fockbias_apply_z
   end interface fockbias_apply

   contains
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine fockbias_setorb( qweight_of_atom, atom_of_orb )

   use fockbias_data, only: fockbias_is_active, fockbias_orbqw

   real*8 , intent(in) :: qweight_of_atom(:)
   integer, intent(in) :: atom_of_orb(:)
   integer             :: Nbasis
   integer             :: nn

   if ( .not. fockbias_is_active ) return

   Nbasis = size(atom_of_orb)
   if ( allocated(fockbias_orbqw) ) deallocate(fockbias_orbqw)
   allocate( fockbias_orbqw(Nbasis) )

   do nn = 1, Nbasis
      fockbias_orbqw(nn) = qweight_of_atom( atom_of_orb(nn) )
   end do

end subroutine fockbias_setorb
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine fockbias_setmat( sqsmat )
   use fockbias_data, only: fockbias_is_active, fockbias_orbqw, fockbias_matrix

   real*8 , intent(in) :: sqsmat(:,:)
   real*8              :: newterm
   integer             :: Nbasis
   integer             :: ii, jj, kk

   if (.not.fockbias_is_active) return

   if (.not.allocated(fockbias_orbqw)) then
      print*, "Error inside fockbias_setmat: setorb never performed."
      print*; stop
   end if

   Nbasis = size( fockbias_orbqw )
   if (( Nbasis /= size(sqsmat,1) ).or.( Nbasis /= size(sqsmat,2) )) then
      print*, "Error inside fockbias_setmat: bad sqsmat input size."
      print*, "   Nbasis    = ", Nbasis
      print*, "   sqsmat(1) = ", size(sqsmat,1)
      print*, "   sqsmat(2) = ", size(sqsmat,2)
      print*; stop
   end if

   if ( allocated(fockbias_matrix) ) deallocate(fockbias_matrix)
   allocate( fockbias_matrix(Nbasis,Nbasis) )
   do jj = 1, Nbasis
   do ii = 1, Nbasis
      fockbias_matrix(ii,jj) = 0.0d0
      do kk = 1, Nbasis
         newterm = fockbias_orbqw(kk) * sqsmat(ii,kk) * sqsmat(kk,jj)
         fockbias_matrix(ii,jj) = fockbias_matrix(ii,jj) + newterm
      end do
   end do
   end do

end subroutine fockbias_setmat
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine fockbias_loads( Natom, atom_of_orb, file_unit_in, file_name_in )
   use fockbias_data, only: fockbias_is_active, fockbias_readfile

   integer         , intent(in)           :: Natom
   integer         , intent(in)           :: atom_of_orb(:)
   integer         , intent(in), optional :: file_unit_in
   character(len=*), intent(in), optional :: file_name_in
   integer                                :: file_unit
   real*8          , allocatable          :: qweight_of_atom(:)
   integer                                :: nn, ios

   if ( .not. fockbias_is_active ) return

   if ( .not. present(file_unit_in) ) then
   if ( .not. present(file_name_in) ) then
   if ( fockbias_readfile == "" )     then
      print*, "Error inside fockbias_loads: without a file_name or a file_unit"
      print*, "I can't do anything..."
      print*; stop
   end if
   end if
   end if

!  If there is input fileunit, we will read from there directly.
!  Else, we will open the file in unit=2427 (BIAS)
   if ( present(file_unit_in) ) then
      file_unit = file_unit_in

   else
      file_unit = 2427

      open( file=fockbias_readfile, unit=file_unit, iostat=ios )
      if ( ios /= 0 ) then
         print*, "Error inside fockbias_loads while opening input."
         print*, "  file_name = ", fockbias_readfile
         print*, "  file_unit = ", file_unit
         print*, "  iostatus  = ", ios
         print*; stop
      end if

   end if

!  Now read all atomic weights and set it up.
   allocate( qweight_of_atom(Natom) )
   do nn = 1, Natom
      read( unit=file_unit, fmt=*, iostat=ios ) qweight_of_atom(nn)
      if ( ios /= 0 ) then
         print*, "Error inside fockbias_loads while reading input."
         print*, "  file_name = ", fockbias_readfile
         print*, "  file_unit = ", file_unit
         print*, "  iostatus  = ", ios
         print*; stop
      end if
   end do
   call fockbias_setorb( qweight_of_atom, atom_of_orb )

!  If there was no input unit, we had to open the file. Now we close it.
   if ( .not. present(file_unit_in) ) then
      close( unit=file_unit, iostat=ios )
      if ( ios /= 0 ) then
         print*, "Error inside fockbias_loads while closing input."
         print*, "  file_name = ", fockbias_readfile
         print*, "  file_unit = ", file_unit
         print*, "  iostatus  = ", ios
         print*; stop
      end if
   end if

end subroutine fockbias_loads
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine fockbias_apply_d( timepos, fockmat)
   use fockbias_data, only: fockbias_is_active, fockbias_is_shaped &
                         &, fockbias_timegrow , fockbias_timefall  &
                         &, fockbias_timeamp0 , fockbias_matrix

   implicit none
   real*8 , intent(in)    :: timepos
   real*8 , intent(inout) :: fockmat(:,:)

   real*8  :: time_shape, exparg
   integer :: Nbasis
   integer :: ii, jj

   if (.not.fockbias_is_active) return

   if (.not.allocated(fockbias_matrix)) then
      print*, "Error inside fockbias_apply: setmat never performed."
      print*; stop
   end if

   Nbasis = size( fockbias_matrix, 1 )
   if (( Nbasis /= size(fockmat,1) ).or.( Nbasis /= size(fockmat,2) )) then
      print*, "Error inside fockbias_apply: bad fockmat input size."
      print*, "   Nbasis     = ", Nbasis
      print*, "   fockmat(1) = ", size(fockmat,1)
      print*, "   fockmat(2) = ", size(fockmat,2)
      print*; stop
   end if

   time_shape = 1.0d0

   if ( fockbias_is_shaped) then

      if ( timepos <= fockbias_timegrow ) then
         exparg = (timepos - fockbias_timegrow) / fockbias_timeamp0
         exparg = (-1.0d0) * exparg * exparg
         time_shape = time_shape * dexp( exparg )
      end if

      if ( timepos >= fockbias_timefall ) then
         exparg = (timepos - fockbias_timefall) / fockbias_timeamp0
         exparg = (-1.0d0) * exparg * exparg
         time_shape = time_shape * dexp( exparg )
      end if

   end if

   if (time_shape < 1.00d-16) time_shape = 0.0d0

   do jj = 1, Nbasis
   do ii = 1, Nbasis
      fockmat(ii,jj) = fockmat(ii,jj) + time_shape * fockbias_matrix(ii,jj)
   enddo
   enddo

end subroutine fockbias_apply_d
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine fockbias_apply_c( timepos, fockmat)
   implicit none
   real*8    , intent(in)    :: timepos
   complex*8 , intent(inout) :: fockmat(:,:)
   real*8    , allocatable   :: fockmat_r(:,:)
   integer                   :: N1, N2

   N1 = size(fockmat,1)
   N2 = size(fockmat,2)
   allocate( fockmat_r( N1, N2 ) )
   call fockbias_apply_d( timepos, fockmat_r)
   fockmat(:,:) = fockmat(:,:) + CMPLX( fockmat_r(:,:), 0.0d0 )
   deallocate( fockmat_r )

end subroutine fockbias_apply_c
!
!
!------------------------------------------------------------------------------!
subroutine fockbias_apply_z( timepos, fockmat)
   implicit none
   real*8    , intent(in)    :: timepos
   complex*16, intent(inout) :: fockmat(:,:)
   real*8    , allocatable   :: fockmat_r(:,:)
   integer                   :: N1, N2

   N1 = size(fockmat,1)
   N2 = size(fockmat,2)
   allocate( fockmat_r( N1, N2 ) )
   call fockbias_apply_d( timepos, fockmat_r)
   fockmat(:,:) = fockmat(:,:) + DCMPLX( fockmat_r(:,:), 0.0d0 )
   deallocate( fockmat_r )

end subroutine fockbias_apply_z
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end module fockbias_subs
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
