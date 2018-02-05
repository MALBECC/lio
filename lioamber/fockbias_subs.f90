!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module fockbias_subs
!------------------------------------------------------------------------------!
!
!    This module controls the aplication of charge biases that are directly
! applied through the fock matrix.
!
! REFERENCE: Physical Review B 74, 155112 ͑2006͒
!
! fockbias_check_ready: internal subroutine to check if setup was performed.
! fockbias_check_msize: internal subroutine to check if size is consistent.
! fockbias_setup_basics: setups the size and the application of the potential.
! fockbias_setup_shaper: setups the parameters of the shaper.
! fockbias_reads_charges: reads the charges to apply from given file.
! fockbias_fockadd: adds the potential to a fiven fock matrix, provided the
!                  shape 
!
!------------------------------------------------------------------------------!
   implicit none
   contains
#  include "fockbias_oldinit.f90"
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine fockbias_check_ready()
   use fockbias_data, only: fockbias_apply, fockbias_ready
   implicit none
   if (.not.fockbias_ready) then
      print*,'FATAL ERROR:'
      print*,'  A module subroutine is being used without proper setup.'
      stop
   endif
end subroutine fockbias_check_ready
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine fockbias_check_msize( msize_inp )
   use fockbias_data, only: fockbias_apply, fockbias_msize
   implicit none
   integer, intent(in) :: msize_inp
   if (msize_inp.ne.fockbias_msize) then
      print*,'FATAL ERROR:'
      print*,'  Wrong setup of fockbias, inconsistent sizes.'
      print*,'  * msize_inp (input):   ', msize_inp
      print*,'  * fockbias_msize (mod): ', fockbias_msize
      stop
   endif
end subroutine fockbias_check_msize
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine fockbias_setup_basics( activate, msize )

   use fockbias_data, &
   &   only: fockbias_apply, fockbias_ready, fockbias_msize, qweight_of_orb

   implicit none
   logical, intent(in) :: activate
   integer, intent(in) :: msize
   integer             :: ii, jj

   fockbias_apply = activate
   fockbias_msize = msize
   fockbias_ready = .true.
   if ( allocated(qweight_of_orb) ) deallocate(qweight_of_orb)
   allocate( qweight_of_orb(msize) )
   qweight_of_orb(:) = 0.0d0

end subroutine fockbias_setup_basics
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine fockbias_setup_shaper( timegrow, timefall, timepos0, timeamp1 )
   
   use fockbias_data, only: fockbias_timegrow, fockbias_timefall &
                        &, fockbias_timepos0, fockbias_timeamp1

   implicit none
   logical, intent(in) :: timegrow
   logical, intent(in) :: timefall
   real*8 , intent(in) :: timepos0
   real*8 , intent(in) :: timeamp1

   fockbias_timegrow = timegrow
   fockbias_timefall = timefall
   fockbias_timepos0 = timepos0
   fockbias_timeamp1 = timeamp1

end subroutine fockbias_setup_shaper
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine fockbias_reads_charges( nsize, msize, atom_of_orb, source_fname )

   use liosubs       , only: read_list, atmvec_to_orbvec
   use fockbias_data  , only: fockbias_apply, qweight_of_orb

   implicit none
   integer          , intent(in) :: nsize
   integer          , intent(in) :: msize
   integer          , intent(in) :: atom_of_orb(msize)
   character(len=80), intent(in) :: source_fname

   real*8, allocatable :: qweight_of_atom(:)

   if (.not.fockbias_apply) return
   call fockbias_check_ready()
   call fockbias_check_msize( msize )

   allocate( qweight_of_atom(nsize) )
   call read_list( source_fname, qweight_of_atom )
   call atmvec_to_orbvec( qweight_of_atom, atom_of_orb, qweight_of_orb )
   deallocate( qweight_of_atom )
   
end subroutine fockbias_reads_charges
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine fockbias_fockadd( msize, timepos, sqsmat, fockao )

   use liosubs     , only: gaussian_shaper
   use fockbias_data, only: fockbias_apply,    qweight_of_orb   &
                        &, fockbias_timegrow, fockbias_timefall &
                        &, fockbias_timepos0, fockbias_timeamp1

   implicit none
   integer, intent(in)    :: msize
   real*8 , intent(in)    :: timepos
   real*8 , intent(in)    :: sqsmat( msize, msize )
   real*8 , intent(inout) :: fockao( msize, msize )

   real*8  :: newterm
   real*8  :: shape_factor
   integer :: ii, jj, kk

   if (.not.fockbias_apply) return
   call fockbias_check_ready()
   call fockbias_check_msize( msize )

   shape_factor = 1.0d0
   call gaussian_shaper( fockbias_timegrow, fockbias_timefall, timepos &
                      &, fockbias_timepos0, fockbias_timeamp1, shape_factor )

   do ii = 1, msize
   do jj = 1, msize
      do kk = 1, msize
         newterm       = shape_factor * qweight_of_orb(kk)
         newterm       = newterm * sqsmat(ii,kk) * sqsmat(kk,jj)
         fockao(ii,jj) = fockao(ii,jj) + newterm
      enddo
   enddo
   enddo

end subroutine fockbias_fockadd
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
end module fockbias_subs
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
