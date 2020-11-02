#include "../datatypes/datatypes.fh"
module cdft_data
   implicit none
   logical      :: doing_cdft = .false.

   ! This data structure contains all CDFT internal variables, except
   ! for CDFT region-specific data.
   type cdft_control_base
      logical :: do_chrg  = .false.     ! If applying charge constraints.
      logical :: do_spin  = .false.     ! If applying spin constraints.
      logical :: mixed    = .false.     ! Mixed (Hab) CDFT calculation.
      logical :: dual     = .false. 
      ! Dual is a special treatment if there are only two regions that
      ! contain all of the atoms.

      LIODBLE, allocatable :: at_chrg(:)     ! List of atomic charges.
      LIODBLE, allocatable :: at_spin(:)     ! List of atomic spin charges.
      LIODBLE, allocatable :: jacob(:,:)     ! Jacobian for advancing Vi
      integer              :: sp_idx    = 0  ! Starting index for spin.
      integer              :: n_regions = 0  ! Number of regions.
      integer              :: max_nat   = 0  ! Max number of atoms per region.
   end type cdft_control_base

   ! Contains internal variables for mixed CDFT calculations.
   type cdft_mcontrol_base
      ! Atomic coeficients alpha and beta for state 1. In closed shell
      ! beta is not allocated.
      LIODBLE, allocatable :: coefs_a1(:,:)
      LIODBLE, allocatable :: coefs_b1(:,:)

      ! Atomic coeficients alpha and beta for state 2.
      LIODBLE, allocatable :: coefs_a2(:,:)
      LIODBLE, allocatable :: coefs_b2(:,:)
   end type cdft_mcontrol_base

   ! Even though it would be much clearer to have a single cdft_regions class
   ! that is allocated a an array of regions (instead of having N region-sized
   ! arrays), this is done this way to avoid potential conflicts with g2g.
   type cdft_region_data
      ! Main data for regions.
      integer, allocatable :: natom(:)   ! Number of atoms.
      integer, allocatable :: nelecs(:)  ! Number of electrons (w/o charge)
      integer, allocatable :: atoms(:,:) ! Atom indexes.
      LIODBLE, allocatable :: chrg(:)    ! Charge targets.
      LIODBLE, allocatable :: spin(:)    ! Spin targets.
      LIODBLE, allocatable :: Vc(:)      ! Charge potentials.
      LIODBLE, allocatable :: Vs(:)      ! Spin potentials.
      LIODBLE, allocatable :: Vmix(:)    ! Array with both potentials.
      LIODBLE, allocatable :: cst(:)     ! Charge/Spin constraints.
      
      ! Second region constraints for mixed CDFT calculations.
      LIODBLE, allocatable :: chrg2(:)   ! Charge targets.
      LIODBLE, allocatable :: spin2(:)   ! Spin targets.
      LIODBLE, allocatable :: Vc2(:)     ! Charge potentials.
      LIODBLE, allocatable :: Vs2(:)     ! Spin potentials.
      
      ! Arrays for Jacobian calculation and propagation
      LIODBLE, allocatable :: Vc_old(:)  ! Charge potentials.
      LIODBLE, allocatable :: Vs_old(:)  ! Spin potentials.
      LIODBLE, allocatable :: Vm_old(:)  ! Array with both potentials.
      LIODBLE, allocatable :: cst_old(:) ! Charge/Spin constraint value.
   end type cdft_region_data

   type(cdft_control_base)  :: cdft_c
   type(cdft_mcontrol_base) :: cdft_mc
   type(cdft_region_data)   :: cdft_reg
end module cdft_data