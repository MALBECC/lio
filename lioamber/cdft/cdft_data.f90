#include "../datatypes/datatypes.fh"
module cdft_data
   implicit none
   logical      :: doing_cdft = .false.

   ! This data structure contains all CDFT internal variables, except
   ! for CDFT region-specific data.
   type cdft_control_base
      logical :: do_chrg  = .false.     ! If applying charge constraints.
      logical :: do_spin  = .false.     ! If applying spin constraints.

      LIODBLE, allocatable :: at_chrg(:)     ! List of atomic charges.
      LIODBLE, allocatable :: at_spin(:)     ! List of atomic spin charges.
      LIODBLE, allocatable :: jacob(:,:)     ! Jacobian for advancing Vi
      integer              :: sp_idx    = 0  ! Starting index for spin.
      integer              :: n_regions = 0  ! Number of regions.
      integer              :: max_nat   = 0  ! Max number of atoms per region.
   end type cdft_control_base

   ! Even though it would be much clearer to have a single cdft_regions class
   ! that is allocated a an array of regions (instead of having N region-sized
   ! arrays), this is done this way to avoid potential conflicts with g2g.
   type cdft_region_data
      ! Main data for regions.
      integer, allocatable :: natom(:)   ! Number of atoms.
      integer, allocatable :: atoms(:,:) ! Atom indexes.
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

   type(cdft_control_base) :: cdft_c
   type(cdft_region_data)  :: cdft_reg
   type(cdft_region_data)  :: cdft_reg2
end module cdft_data