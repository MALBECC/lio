#include "datatypes/datatypes.fh"
module excited_data
   implicit none

   ! Save Integrals in Memory
   integer :: libint_recalc = 0 ! Recalculated Integrals, 1 = Save in Memory

   ! Linear Response
   logical :: lresp = .false.
   integer :: nstates = 4 ! poner en el input
   logical :: fittExcited = .false. ! poner de input
   LIODBLE :: tolv = 1.0d-7 ! conv. crit. to vectors ! pon inp
   LIODBLE :: tole = 1.0d-7 ! conv. crit. to energy !pon inp
   integer :: root = 0 ! Relaxed Density Matrix of Excited State root

   ! Frozen Core and Valence Approximation
   logical :: FCA = .false.
   integer :: nfo = 3
   integer :: nfv = 3

   ! Basis Change
   LIODBLE, dimension(:,:), allocatable :: Coef_trans, Cocc 
   LIODBLE, dimension(:,:), allocatable :: Cocc_trans, Cvir
   LIODBLE, dimension(:,:), allocatable :: Cvir_trans

   ! Excited States Forces
   logical :: excited_forces = .false.
   LIODBLE, dimension(:,:), allocatable :: for_exc

   ! Save Relaxed Density matrix in vector form in order to obtain
   ! differents properties
   LIODBLE, dimension(:), allocatable   :: pack_dens_exc

   ! Trajectory Surface Hopping
   logical :: TSH = .false.
   integer :: tsh_Jstate, tsh_Kstate
   LIODBLE :: dE_accum, lambda, tsh_time_dt
   LIODBLE, dimension(:,:), allocatable :: gamma_old
   complex(kind=8) :: B_old = cmplx(0.0d0,0.0d0,8)
   complex(kind=8), dimension(:), allocatable :: tsh_coef
   

end module excited_data
