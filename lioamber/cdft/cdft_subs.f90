!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% CONSTRAINED DFT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! These files contains routines (cdft_subs) and variables (cdft_data) for      !
! Constrained DFT (CDFT) calculations. This implementation uses the Becke      !
! partitioning scheme for the constraints, as used by Holmberg and Laasoen     !
! (JCTC 2016, DOI: 10.1021/acs.jctc.6b01085).                                  !
!                                                                              !
! The following subroutines are called externally:                             !
!   * cdft_options_check: checks internal consistency in input options.        !
!   * cdft_input_read:    reads CDFT block input.                              !
!   * CDFT:               main CDFT loop outside SCF.                          !
!                                                                              !
! The following subroutines are only called internally:                        !
!   * cdft_get_constraints: gets sum(atomic_spin/charges) - spin/charge_target.!
!   * cdft_add_energy:      adds CDFT terms to energy.                         !
!   * cdft_initialise:      allocates arrays for CDFT.                         !
!   * cdft_finalise:        deallocates arrays for CDFT.                       !
!   * cdft_get_deltaV:      gets the increment for bias potentials.            !
!   * cdft_set_potential:   sets new potential for grid integration.           !
!   * cdft_check_conver:    checks CDFT convergence.                           !
!                                                                              !
!------------------------------------------------------------------------------!
! How to input CDFT options:                                                   !
!   CDFT options must be provided in the LIO input file (lio.in) but outside   !
!   namelists. Said file should contain the following sections:                !
!                                                                              !
!{CDFT}                                                                        !
!  N_REG CONST_CHARGE CONST_SPIN                                               !
!  REGION1_NATOM REGION1_CHARGE REGION1_SPIN                                   !
!  REGION2_NATOM REGION2_CHARGE REGION2_SPIN                                   !
!  REGIONN_NATOM REGIONN_CHARGE REGIONN_SPIN                                   !
!  REGION1_ATOMS                                                               !
!  REGION2_ATOMS                                                               !
!  REGIONN_ATOMS                                                               !
!{END}                                                                         !
!                                                                              !
! The {CDFT} and {END} terms indicate the beginning and end of CDFT input.     !
! N_REG (integer) is the number of regions to constrain, while CONST_CHARGE    !
! and CONST_SPIN (both integers too) indicate whether either/both charge       !
! (CONST_CHARGE=1) or/and spin (CONST_SPIN=1) are constrained.                 !
! After that first line, the following lines contain region information; each  !
! line belongs to a specific region of the molecule. REGION_CHARGE indicates   !
! the target charge for a region (in LIODBLE), REGION_SPIN indicates  !
! the target spin of a region, and REGION_NATOM indicates the number of atoms  !
! in a region. Finally the last lines contain the atom indexes belonging to a  !
! region, in the same order as specified in the above lines. These should be   !
! written as a space-separated integer list. See the CDFT test in /tests for   !
! an example.                                                                  !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
#include "../datatypes/datatypes.fh"
module cdft_subs
   implicit none
   private
   public :: cdft_input_read
   public :: cdft_options_check
   public :: CDFT
contains

! Public subroutines.
#include "cdft_main.f90"
#include "cdft_input.f90"

! Private subroutines.
#include "cdft_init_fin.f90"
#include "cdft_utils.f90"
#include "cdft_mixed_utils.f90"
#include "cdft_mixed_hab.f90"
#include "cdft_Vk.f90"

end module cdft_subs