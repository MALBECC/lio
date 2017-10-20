!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       module general_module
!--------------------------------------------------------------------!
!
!
! INCLUDE FILES WITH HEADERS:
!--------------------------------------------------------------------!
       implicit none
       include 'read_list_h.f'
       include 'atmorb_h.f'
       contains
!
!
! INCLUDE FILES WITH PROCEDURES:
!--------------------------------------------------------------------!
       include 'read_list.f'
       include 'vector_selection.f'
       include 'atmorb.f'
       include 'sdcmp_cholesky.f'
       include 'sdiag_canonical.f'
       end module
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
