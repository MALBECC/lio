module fileio_data
   implicit none
   logical :: style   = .false.
   integer :: verbose = 3
end module fileio_data

module fileio
   implicit none
   include 'restart_commons_h.f90'
   include 'restart_coef_h.f90'
   include 'restart_fock_h.f90'
   include 'restart_rho_h.f90'
   include 'restart_td_h.f90'
contains
   include 'restart_commons.f90'
   include 'restart_coef.f90'
   include 'restart_fock.f90'
   include 'restart_rho.f90'
   include 'restart_td.f90'
   include 'output_init.f90'
   include 'output_scf.f90'
end module fileio
