!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       module mathsubs
!--------------------------------------------------------------------!
!
!
! INCLUDE FILES WITH HEADERS:
!--------------------------------------------------------------------!
       implicit none
       include 'gaussbell_h.f'
       include 'commutator_h.f'
       include 'basechange_h.f'
       include 'basechange_gemm_h.f'
       contains
!
!
! INCLUDE FILES WITH PROCEDURES:
!--------------------------------------------------------------------!
       include 'gaussbell.f'
       include 'commutator.f'
       include 'basechange.f'
       include 'basechange_gemm.f'
       end module
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
