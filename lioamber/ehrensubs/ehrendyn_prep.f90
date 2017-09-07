!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine ehrendyn_prep( Nbasis, Natoms, time, recalc_forces, Xmat, Xtrp,     &
                        & Sinv, Rmon, Fmat, Bmat, Dmat, Tmat, nucpos, nucvel,  &
                        & nucfor, dipmom, energy )
!------------------------------------------------------------------------------!
!
! Multi-purpose subroutine that does the following:
!
! (1) Adds to the input Fmat (in AO) the part of fock that depends on the
!     introduced Rmon (in ON) and returns it in the ON base.
!
! (2) If recalc_forces, it will recalculate the part of the forces that depend
!     on the Rmon and will also recalculate the matrices Bmat and Dmat (note 
!     that these two do not depend on Rmon, but only on the velocity and the
!     position of the nuclei/basis).
!
! (3) It will recalculate Tmat  ( Fmat + i * Dmat ) with the new Fmat and
!     either the newly calculated Dmat, or the one introduced as input.
!
! NOTE: this subroutine leaves the Rmon and the corresponding Fock matrix
!       inside of RMM, both in AO.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   implicit none
   integer   , intent(in)    :: Nbasis
   integer   , intent(in)    :: Natoms
   real*8    , intent(in)    :: time
   logical   , intent(in)    :: recalc_forces

   real*8    , intent(in)    :: Xmat(Nbasis, Nbasis)
   real*8    , intent(in)    :: Xtrp(Nbasis, Nbasis)
   real*8    , intent(in)    :: Sinv(Nbasis, Nbasis)
   complex*16, intent(in)    :: Rmon(Nbasis, Nbasis)
   real*8    , intent(inout) :: Fmat(Nbasis, Nbasis)
   real*8    , intent(inout) :: Bmat(Nbasis, Nbasis)
   real*8    , intent(inout) :: Dmat(Nbasis, Nbasis)
   complex*16, intent(inout) :: Tmat(Nbasis, Nbasis)

   real*8    , intent(in)    :: nucpos(3, Natoms)
   real*8    , intent(in)    :: nucvel(3, Natoms)
   real*8    , intent(inout) :: nucfor(3, natoms)
   real*8    , intent(inout) :: dipmom(3)
   real*8    , intent(inout) :: energy

   complex*16, allocatable   :: Rmao(:,:)
   real*8    , allocatable   :: nucfor_add(:,:)
   real*8                    :: elec_field(3)
!
!
!------------------------------------------------------------------------------!
   allocate( Rmao(Nbasis, Nbasis) )

   Rmao = matmul( Rmon, Xtrp )
   Rmao = matmul( Xmat, Rmao )

   call ehrenaux_setfld(  time, elec_field )
   call RMMcalc3_FockMao( Rmao, elec_field, Fmat, dipmom, energy)

   if (recalc_forces) then
      allocate( nucfor_add(3,natoms) )
      nucfor_add(:,:) = 0.0d0
      call calc_forceDS( natoms, nbasis, nucpos, nucvel, Rmao, Fmat, Sinv, &
                       & Bmat, nucfor_add )
      Dmat = calc_Dmat( nbasis, Xtrp, Xmat, Bmat )
      nucfor(:,:) = nucfor(:,:) + nucfor_add(:,:)
      deallocate( nucfor_add )
   endif

   Fmat = matmul( Fmat, Xmat )
   Fmat = matmul( Xtrp, Fmat )
   Tmat = DCMPLX(Fmat) + DCMPLX(0.0d0,1.0d0) * DCMPLX(Dmat)

   deallocate( Rmao )
end subroutine ehrendyn_prep
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
