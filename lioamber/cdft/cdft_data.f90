#include "../datatypes/datatypes.fh"
module cdft_data
   implicit none
   logical      :: cdft_chrg  = .false.
   logical      :: cdft_spin  = .false.
   logical      :: doing_cdft = .false.

   LIODBLE, allocatable :: at_chrg(:)       ! List of atomic charges.
   LIODBLE, allocatable :: at_spin(:)       ! List of atomic spin charges.
   LIODBLE, allocatable :: jacob(:,:)       ! Jacobian for advancing Vi
   integer                   :: sp_idx  = 0      ! Starting index for spin.

   type cdft_region_data
      integer      :: n_regions = 0 ! Number of regions.
      integer      :: max_nat   = 0 ! Maximum number of atoms in a region.

      ! Main data for regions.
      integer     , allocatable :: natom(:)   ! Number of atoms.
      integer     , allocatable :: atoms(:,:) ! Atom indexes.
      LIODBLE, allocatable :: chrg(:)    ! Charge targets.
      LIODBLE, allocatable :: spin(:)    ! Spin targets.
      LIODBLE, allocatable :: Vc(:)      ! Charge potentials.
      LIODBLE, allocatable :: Vs(:)      ! Spin potentials.
      LIODBLE, allocatable :: Vmix(:)    ! Array with both potentials.
      LIODBLE, allocatable :: cst(:)     ! Charge/Spin constraints.
      
      ! Arrays for Jacobian calculation and propagation
      LIODBLE, allocatable :: Vc_old(:)  ! Charge potentials.
      LIODBLE, allocatable :: Vs_old(:)  ! Spin potentials.
      LIODBLE, allocatable :: Vm_old(:)  ! Array with both potentials.
      LIODBLE, allocatable :: cst_old(:) ! Charge/Spin constraint value.

   end type cdft_region_data

   type(cdft_region_data) :: cdft_reg
end module cdft_data