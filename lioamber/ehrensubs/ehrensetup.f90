!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
! Save the forces and density matrix of the last position for restart purposes
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine ehrensetup( Natoms, Nbasis, RealRho )
   use ehrendata,  only: RhoSaveA, RhoSaveB
   use garcha_mod, only: qm_forces_ds, qm_forces_total
   implicit none
   integer, intent(in) :: Natoms
   integer, intent(in) :: Nbasis
   real*8,  intent(in) :: RealRho( Nbasis, Nbasis )

   RhoSaveA = DCMPLX( 0.0d0, 0.0d0 )
   RhoSaveB = DCMPLX( RealRho )

   if (.not.allocated(qm_forces_total)) then
      allocate( qm_forces_total(3, Natoms) )
      qm_forces_total = 0.0d0
   endif

   if (.not.allocated(qm_forces_ds)) then
      allocate( qm_forces_ds(3, Natoms) )
      qm_forces_ds = 0.0d0
   endif

end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
