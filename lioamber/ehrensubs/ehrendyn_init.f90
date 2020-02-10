!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine ehrendyn_init( Natoms, Nbasis, RealRho )
!------------------------------------------------------------------------------!
!
! DESCRIPTION
! Save the forces and density matrix of the last position for restart purposes
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   use ehrendata,  only: stored_densM1, stored_densM2
   use garcha_mod, only: qm_forces_ds, qm_forces_total

   implicit none
   integer, intent(in) :: Natoms
   integer, intent(in) :: Nbasis
   real*8,  intent(in) :: RealRho( Nbasis, Nbasis )

   if (allocated(stored_densM1)) deallocate(stored_densM1)
   allocate(stored_densM1( Nbasis, Nbasis ))
   stored_densM1 = DCMPLX( 0.0d0, 0.0d0 )

   if (allocated(stored_densM2)) deallocate(stored_densM2)
   allocate(stored_densM2( Nbasis, Nbasis ))
   stored_densM2 = DCMPLX( RealRho )

   if (allocated(qm_forces_total)) deallocate(qm_forces_total)
   allocate( qm_forces_total(3, Natoms) )
   qm_forces_total = 0.0d0

   if (allocated(qm_forces_ds)) deallocate(qm_forces_ds)
   allocate( qm_forces_ds(3, Natoms) )
   qm_forces_ds = 0.0d0

end subroutine ehrendyn_init
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
