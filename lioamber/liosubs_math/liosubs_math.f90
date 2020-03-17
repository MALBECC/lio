!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
#include "../datatypes/datatypes.fh"
module math_data
   implicit none
   ! Arrays used for Boys function.
   LIODBLE :: STR(880,0:21) = 0.0D0, FAC(0:16) = 0.0D0
end module math_data

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
#     include "funct.f90"
#     include "matdcmp_cholesky.f90"
#     include "matdcmp_svd.f90"
#     include "matmul3_head.f90"
#     include "purge_zeros.f90"

#     define transform_gen transform_r
#     define GEN_TYPE LIODBLE
#     include "transform_gen.f90"
#     undef transform_gen
#     undef GEN_TYPE

#     define CONVERT_R

#     define transform_gen transform_c
#     define GEN_TYPE complex(kind=4)
#     define CPSIZE 4
#     include "transform_gen.f90"
#     undef transform_gen
#     undef CPSIZE
#     undef GEN_TYPE

#     define transform_gen transform_z
#     define GEN_TYPE complex(kind=8)
#     define CPSIZE 8
#     include "transform_gen.f90"
#     undef transform_gen
#     undef GEN_TYPE
#     undef CPSIZE

#     undef CONVERT_R

end module liosubs_math
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
