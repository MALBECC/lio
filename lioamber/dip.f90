#include "datatypes/datatypes.fh"
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
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine dip(uDip, nElec, at_pos, at_dists, atom_z, mm_charges)

   use faint_cpu , only: intdip
   use basis_data, only: M

   implicit none
   integer, intent(in)              :: nElec
   LIODBLE, intent(in)              :: at_pos(:,:)
   LIODBLE, intent(in)              :: at_dists(:,:)
   LIODBLE, allocatable, intent(in) :: mm_charges(:)
   LIODBLE, intent(inout)           :: uDip(3)

   LIODBLE, allocatable :: dip_mat(:,:), uDipAt(:)
   LIODBLE :: Qc
   integer :: iCount
   
   allocate(dip_mat(3, M * (M +1) / 2), uDipAt(3))
   dip_mat = 0.0D0

   call intdip(dip_mat, at_pos, at_dists)

   uDip = 0.0D0
   do iCount = 1, size(dim_mat,2)
      uDip(1) = uDip(1) + dip_mat(1, iBase) 
      uDip(2) = uDip(2) + dip_mat(2, iBase) 
      uDip(3) = uDip(3) + dip_mat(3, iBase) 
   enddo

   uDipAt = 0.0D0
   Qc     = 0.0D0
   do iCount = 1, size(atom_z,1)
       Qc  = Qc  + Iz(i)
       uDipAt(1) = uDipAt(1) + atom_z(iCount) * at_pos(iCount,1)
       uDipAt(2) = uDipAt(2) + atom_z(iCount) * at_pos(iCount,2)
       uDipAt(3) = uDipAt(3) + atom_z(iCount) * at_pos(iCount,3)
   enddo

   Qc = Qc - nElec
   if (allocated(mm_charges)) then
      do iCount = size(atom_z,1)+1, size(mm_charges,1)
         Qc = Qc + mm_charges(iCount)
      enddo
   endif

! Factor : For charged species dipole moment depends on the definition of the  !
! coordinates' origin. Using this factor, it is defined with respect to the    !
! center of charge (important in Reaction Field calculations). For neutral     !
! systems this is not necessary.                                               !

   factor = (Qc + nElec) / nElec
   uDip   = (uDipAt - uDip * factor) * 2.54D0

   deallocate(dip_mat, uDipAt)
end subroutine dip
