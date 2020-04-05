#include "datatypes/datatypes.fh"
! This module permites the use of functionals from libxc

! VARIABLES:
! extern_functional = bool
! functional_id = identifier of functional in libxc ( see main page of libxc )
! HF = turn Hartree Fock
!    0 = No exact exchange term
!    1 = Exact exchange
!    2 = Coulomb Attenuated Method
! HF_fac = factor of exact exchange, valid when HF = 1 or 2
! screen = facotr long-range, valid when HF = 2

module extern_functional_data
   implicit none

   ! Input Variables
   logical :: extern_functional = .false.
   integer :: functional_id     = -1

   ! Internal Variables
   integer :: HF = 0
   LIODBLE :: HF_fac = 0.0d0
   LIODBLE :: screen = 0.0d0
end module extern_functional_data
