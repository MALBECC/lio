!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% DIP.F90 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Subroutine used for the calculation of the dipole moment in NEUTRAL          !
! (non-ionic) systems. Integrals are evaluated using the Obara Saika method.   !
! Inputs the density basis and outputs the dipole moment components.           !
! Original file: 19-1-1993                                                     !
!                                                                              !
! A loop is performed over all basis functions. Basis are supposed to be       !
! ordered according to type: first all s, then all p, then all d, etc.; and    !
! inside each type, they are ordered in shells: px, py, pz, dx2, dxy, dyy, dzx,!
! dzy, dzz, and so on.                                                         !
!                                                                              !
! ns, np, nd are markers for the end of s, p and d sections respectively.      !
! r(Nuc(i),j) is j component of position of nucleus i, j = 1..3.               !
!                                                                              !
! The output is a Matrix with dipole moment integral in AO basis: uDip         !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine get_DipMatrix(uDip)

    use garcha_mod   , only: NCO, Nunp, r, d
    use basis_data   , only: a, c, Nuc, ncont, M, nshell, norm
    use constants_mod, only: pi32
    use faint_cpu    , only: intdip
    
    implicit none
    LIODBLE, intent(inout) :: uDip(3,M*(M+1)/2)

    call intdip(uDip, r, d)

end subroutine get_DipMatrix
