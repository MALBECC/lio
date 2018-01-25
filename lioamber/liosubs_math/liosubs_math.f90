!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module liosubs_math
!--------------------------------------------------------------------!
   implicit none

   interface matmul3
      module procedure matmul3_ddd
      module procedure matmul3_dcd
   end interface matmul3

   interface purge_zeros
      module procedure purge_zeros_v
      module procedure purge_zeros_m
   end interface purge_zeros

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   contains
#     include "matdcmp_cholesky.f90"
#     include "matdcmp_svd.f90"
#     include "matmul3_head.f90"
#     include "purge_zeros.f90"

end module liosubs_math
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
