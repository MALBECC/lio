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

   interface transform
      module procedure transform_r
      module procedure transform_c
      module procedure transform_z
   end interface transform

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   contains
#     include "matdcmp_cholesky.f90"
#     include "matdcmp_svd.f90"
#     include "matmul3_head.f90"
#     include "purge_zeros.f90"

#     define transform_gen transform_r
#     define GEN_TYPE real*8
#     include "transform_gen.f90"
#     undef transform_gen
#     undef GEN_TYPE

#     define transform_gen transform_c
#     define GEN_TYPE complex*8
#     include "transform_gen.f90"
#     undef transform_gen
#     undef GEN_TYPE

#     define transform_gen transform_z
#     define GEN_TYPE complex*16
#     include "transform_gen.f90"
#     undef transform_gen
#     undef GEN_TYPE

end module liosubs_math
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
