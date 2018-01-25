!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module atompot_subs
!------------------------------------------------------------------------------!
!
!    This module controls the aplication of charge biases that are directly
! applied through the fock matrix.
!
! REFERENCE: Physical Review B 74, 155112 ͑2006͒
!
! atompot_check_ready: internal subroutine to check if setup was performed.
! atompot_check_msize: internal subroutine to check if size is consistent.
! atompot_setup_basics: setups the size and the application of the potential.
! atompot_setup_shaper: setups the parameters of the shaper.
! atompot_reads_charges: reads the charges to apply from given file.
! atompot_fockadd: adds the potential to a fiven fock matrix, provided the
!                  shape 
!
!------------------------------------------------------------------------------!
   implicit none
   contains
#  include "atompot_oldinit.f90"
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine atompot_check_ready()
   use atompot_data, only: atompot_apply, atompot_ready
   implicit none
   if (.not.atompot_ready) then
      print*,'FATAL ERROR:'
      print*,'  A module subroutine is being used without proper setup.'
      stop
   endif
end subroutine atompot_check_ready
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine atompot_check_msize( msize_inp )
   use atompot_data, only: atompot_apply, atompot_msize
   implicit none
   integer, intent(in) :: msize_inp
   if (msize_inp.ne.atompot_msize) then
      print*,'FATAL ERROR:'
      print*,'  Wrong setup of atompot, inconsistent sizes.'
      print*,'  * msize_inp (input):   ', msize_inp
      print*,'  * atompot_msize (mod): ', atompot_msize
      stop
   endif
end subroutine atompot_check_msize
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine atompot_setup_basics( activate, msize )

   use atompot_data, &
   &   only: atompot_apply, atompot_ready, atompot_msize, qweight_of_orb

   implicit none
   logical, intent(in) :: activate
   integer, intent(in) :: msize
   integer             :: ii, jj

   atompot_apply = activate
   atompot_msize = msize
   atompot_ready = .true.
   if ( allocated(qweight_of_orb) ) deallocate(qweight_of_orb)
   allocate( qweight_of_orb(msize) )
   qweight_of_orb(:) = 0.0d0

end subroutine atompot_setup_basics
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine atompot_setup_shaper( timegrow, timefall, timepos0, timeamp1 )
   
   use atompot_data, only: atompot_timegrow, atompot_timefall &
                        &, atompot_timepos0, atompot_timeamp1

   implicit none
   logical, intent(in) :: timegrow
   logical, intent(in) :: timefall
   real*8 , intent(in) :: timepos0
   real*8 , intent(in) :: timeamp1

   atompot_timegrow = timegrow
   atompot_timefall = timefall
   atompot_timepos0 = timepos0
   atompot_timeamp1 = timeamp1

end subroutine atompot_setup_shaper
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine atompot_reads_charges( nsize, msize, atom_of_orb, source_fname )

   use liosubs       , only: read_list, atmvec_to_orbvec
   use atompot_data  , only: atompot_apply, qweight_of_orb

   implicit none
   integer          , intent(in) :: nsize
   integer          , intent(in) :: msize
   integer          , intent(in) :: atom_of_orb(msize)
   character(len=80), intent(in) :: source_fname

   real*8, allocatable :: qweight_of_atom(:)

   if (.not.atompot_apply) return
   call atompot_check_ready()
   call atompot_check_msize( msize )

   allocate( qweight_of_atom(nsize) )
   call read_list( source_fname, qweight_of_atom )
   call atmvec_to_orbvec( qweight_of_atom, atom_of_orb, qweight_of_orb )
   deallocate( qweight_of_atom )
   
end subroutine atompot_reads_charges
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine atompot_fockadd( msize, timepos, sqsmat, fockao )

   use liosubs     , only: gaussian_shaper
   use atompot_data, only: atompot_apply,    qweight_of_orb   &
                        &, atompot_timegrow, atompot_timefall &
                        &, atompot_timepos0, atompot_timeamp1

   implicit none
   integer, intent(in)    :: msize
   real*8 , intent(in)    :: timepos
   real*8 , intent(in)    :: sqsmat( msize, msize )
   real*8 , intent(inout) :: fockao( msize, msize )

   real*8  :: newterm
   real*8  :: shape_factor
   integer :: ii, jj, kk

   if (.not.atompot_apply) return
   call atompot_check_ready()
   call atompot_check_msize( msize )

   shape_factor = 1.0d0
   call gaussian_shaper( atompot_timegrow, atompot_timefall, timepos &
                      &, atompot_timepos0, atompot_timeamp1, shape_factor )

   do ii = 1, msize
   do jj = 1, msize
      do kk = 1, msize
         newterm       = shape_factor * qweight_of_orb(kk)
         newterm       = newterm * sqsmat(ii,kk) * sqsmat(kk,jj)
         fockao(ii,jj) = fockao(ii,jj) + newterm
      enddo
   enddo
   enddo

end subroutine atompot_fockadd
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
end module atompot_subs
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
