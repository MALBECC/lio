#include "datatypes/datatypes.fh"
module ceed_data
   use typedef_operator, only: operator
   use typedef_cumat   , only: cumat_r

   implicit none

!Important input variables:
! * ceed_calc    : logical, indicates if a CEED calculation will be performed.
! * ceed_td_step : integer, indicates what step CEED will be start.
! * k_ceed       : Double precision, acceleration factor for CEED term.
   logical                     :: ceed_calc      = .false.
   integer                     :: ceed_td_step = 100
   LIODBLE                     :: k_ceed         = 1.0d0
   LIODBLE       , parameter   :: A_ceed         = 2.592668788247536d-7
                                                  !2/3 \hbar \alpha/c^2
   LIODBLE       , allocatable :: d2dip_ceed(:,:) !Matrix to store the total 
                                                  !second derivatives of the
                                                  !dipole moment in the 3
                                                  !coordinates.
   type(cumat_r)               :: Xmat_ceed !Base change operator
   type(operator), allocatable :: dip_ceed_op(:) !Dipole matrix operator for
                                                 !each coordinate.
   type(operator), allocatable :: d2ip_ceed_op(:,:)!Second derivative dipole
                                                   !matrix operator for each
                                                   !coordinate and each spin.
   type(operator), allocatable :: fock_ceed_op(:)!Auxiliar fock matrix operator.

end module ceed_data
