!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  module maskrmm
!
!
!--------------------------------------------------------------------!
! INCLUDE FILES WITH HEADERS:
!--------------------------------------------------------------------!
  implicit none
  include 'rmmput_dens_h.f90'
  include 'rmmget_dens_h.f90'
!  include 'rmmput_fock_h.f90'
  include 'rmmget_fock_h.f90'
  contains
!
!--------------------------------------------------------------------!
! INCLUDE FILES WITH PROCEDURES:
!--------------------------------------------------------------------!
  include 'rmmput_dens_all.f90'
  include 'rmmget_dens_all.f90'
!  include 'rmmput_fock_all.f90'
  include 'rmmget_fock_all.f90'
  end module
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
