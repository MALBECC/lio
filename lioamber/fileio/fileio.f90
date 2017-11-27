module fileio
   implicit none
   include 'restart_commons_h.f90'
   include 'restart_rho_h.f90'
   include 'restart_fock_h.f90'
contains
   include 'restart_commons.f90'
   include 'restart_rho.f90'
   include 'restart_fock.f90'
end module fileio
