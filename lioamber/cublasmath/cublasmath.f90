!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
#include "../datatypes/datatypes.fh"
module cublasmath
!--------------------------------------------------------------------!
! INCLUDE FILES WITH HEADERS:
!--------------------------------------------------------------------!
   implicit none
#   include "cumpx_h.f90"
#   include "cumpxt_h.f90"
#   include "cumxp_h.f90"
#   include "cumxtp_h.f90"
#   include "basechange_cublas_h.f90"
#   include "commutator_cublas_h.f90"

!--------------------------------------------------------------------!
! INCLUDE FILES WITH PROCEDURES:
!--------------------------------------------------------------------!
contains

#   include "cumfx.f90"
#   include "cumpxt.f90"
#   include "cumpxt_r.f90"
#   include "cumxtf.f90"
#   include "cumpx.f90"
#   include "cumpx_r.f90"
#   include "cumxp.f90"
#   include "cumxp_r.f90"
#   include "cumxtp.f90"
#   include "basechange_cublas.f90"
#   include "commutator_cublas.f90"
#   include "cu_fock_commuts.f90"
end module
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
