#include "datatypes/datatypes.fh"
module DOS_data

!Important input variables:                                                    !
! * dos_calc  : Logical, indicates if a DOS calculation is performed           !
! * pdos_calc : Logical, indicates if a PDOS calculation is performed          !
! * pdos_allb : Logical, indicates if a PDOS calculation is performed for each !
!               the basis.                                                     !
!                                                                              !
!Input variables readed from the file PDOS_dat.in in the next order:           !
! * min_level : Integer, indicates the minimum energy level to represent.      !
! * dos_nsteps: Integer, indicates the number of points for the graphic.       !
! * dos_sigma : Integer, indicates the with of the gaussian function centered  !
!               in each state.                                                 !
! * dos_Eref    : Double precision real, Fermi energy for reference.           !
! * pdos_natoms : Integer, number of atoms where project the PDOS.             !
! * pdos_nbases : Integer, number of basis where project the PDOS.             !
! * pdos_nuc    : Integer array, atoms where project the PDOS.                 !
! * pdos_b      : Integer array, basis where project the PDOS.                 !

   implicit none
   logical               :: dos_calc     = .false.
   logical               :: pdos_calc    = .false.
   logical               :: pdos_allb    = .false.
   integer               :: min_level    = 1
   integer               :: dos_nsteps   = 200000
   LIODBLE          :: dos_sigma    = 0.0004d0
   LIODBLE          :: dos_Eref     = 0.0d0
   integer               :: pdos_nbases  = 0     
   integer               :: pdos_natoms  = 0

   integer     , allocatable :: pdos_nuc(:)
   integer     , allocatable :: pdos_base(:)
   LIODBLE, allocatable :: pdos(:)
   LIODBLE, allocatable :: pdos_b(:,:)

end module DOS_data
