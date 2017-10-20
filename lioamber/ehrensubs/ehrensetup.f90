!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
! Save the forces and density matrix of the last position for restart purposes
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine ehrensetup( Nbasis, RealRho )
   use ehrendata, only: RhoSaveA, RhoSaveB
   implicit none
   integer, intent(in) :: Nbasis
   real*8,  intent(in) :: RealRho( Nbasis, Nbasis )

   RhoSaveA = DCMPLX( RealRho )
   RhoSaveB = DCMPLX( RealRho )

end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
