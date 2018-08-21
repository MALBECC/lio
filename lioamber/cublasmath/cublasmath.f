!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       module cublasmath
!--------------------------------------------------------------------!
! INCLUDE FILES WITH HEADERS:
!--------------------------------------------------------------------!
c       implicit none
       use faint_cpu, only: int3lu
       include 'cuconmut_h.f'
       include 'cumagnusfac_h.f'
       include 'cumpx_h.f'
       include 'cumpxt_h.f'
       include 'cumxp_h.f'
       include 'cumxtp_h.f'
       include 'cupredictor_h.f'
       include 'basechange_cublas_h.f'
       include 'commutator_cublas_h.f'
       include 'magnus_cublas_h.f'

!--------------------------------------------------------------------!
! INCLUDE FILES WITH PROCEDURES:
!--------------------------------------------------------------------!
       CONTAINS
       include 'cuconmut_r.f'
       include 'cumfx.f'
       include 'cumpxt.f'
       include 'cumpxt_r.f'
       include 'cumxtf.f'
       include 'cuconmut.f'
       include 'cumagnusfac.f'
       include 'cumatmul_r.f'
       include 'cumpx.f'
       include 'cumpx_r.f'
       include 'cumxp.f'
       include 'cumxp_r.f'
       include 'cumxtp.f'
       include 'cupredictor.f'
       include 'basechange_cublas.f'
       include 'commutator_cublas.f'
       include 'magnus_cublas.f'
       include 'cu_fock_commuts.f'
       end module
!
!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
