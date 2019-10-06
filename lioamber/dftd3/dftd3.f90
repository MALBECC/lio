!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%% DFTD3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! This is a module which adds DFT-D3 dispersion corrections. It is based on    !
! Grimme's group DFTD3 code, which is freely available at                      ! 
! https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/dft-d3/ . It    !
! includes a much more detailed implementation.                                !
! Also, see doi.org/10.1063/1.3382344 for references.                          !
!                                                                              !
! The module has four subroutines which are called externally, and are all     !
! present in dtd3_main.f90:                                                    !
!   * dftd3_setup allocates the arrays needed.                                 !
!   * dftd3_finalise deallocates said arrays.                                  !
!   * dftd3_energy first calculates C6 and C8 coefficients for current geometry!
!                  and then adds two- and three-body dispersion corrections to !
!                  energy.                                                     !
!   * dftd3_gradients calculates atomic gradients (not forces). It needs to be !
!                     called after a call to dftd3_energy since C6/C8          !
!                     coefficients are not recomputed.                         !
!                                                                              !
! NOTE: On gradient calculation, the dependence of C6(R) is discarded as a     !
!       small term. Future implementations may include said corrections.       !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

module dftd3_data
   implicit none
   logical      :: dftd3     = .false.
   real(kind=8) :: dftd3_cut = 100.0D0

   ! Parameters dependent on the XC functional. Defaults are for PBE.
   real(kind=8) :: dftd3_s6  = 1.0D0
   real(kind=8) :: dftd3_sr6 = 1.217D0
   real(kind=8) :: dftd3_s8  = 0.722D0

   ! Variables only used internally.
   real(kind=8), allocatable :: c6_ab(:,:), c8_ab(:,:), r0_ab(:,:)
   real(kind=8), allocatable :: c6_cn(:,:,:,:,:), c8_coef(:), r_cov(:)
contains
end module dftd3_data

module dftd3

   contains

# include "dftd3_c6c8.f90"
# include "dftd3_2.f90"
# include "dftd3_3.f90"
# include "dftd3_read_params.f90"
# include "dftd3_main.f90"

end module dftd3