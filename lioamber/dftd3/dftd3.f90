! A module for DFTD3 Grimme corrections.
! See doi.org/10.1063/1.3382344 for references.
! Also, visit https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/dft-d3/ for
! a more detailed implementation.

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

# include "dftd3_2.f90"
# include "dftd3_3.f90"
# include "dftd3_read_params.f90"
# include "dftd3_main.f90"


end module dftd3