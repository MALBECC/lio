#include "datatypes/datatypes.fh"
module fstsh_data
   implicit none

   ! Trajectory Surface Hopping Variable
   logical :: FSTSH = .false. ! input
   integer :: call_number = 0
   integer :: tsh_file = 456
   logical :: after_hopp = .false.
 
   ! Nuclear Time Step
   LIODBLE :: tsh_time_dt = 0.0d0

   ! Overlap in differents times
   LIODBLE, dimension(:,:), allocatable :: Sovl_old, Sovl_now

   ! Old Variables needed in order to obtain overlap at differents times
   LIODBLE, allocatable :: a_old(:,:)
   LIODBLE, allocatable :: c_old(:,:)
   LIODBLE, allocatable :: r_old(:,:)

   ! Molecular Orbital variables
   LIODBLE, allocatable :: C_scf_old(:,:)
   LIODBLE, allocatable :: WFcis_old(:,:)

   ! Potential Energy Surface
   LIODBLE, allocatable :: Nesup_now(:)
   LIODBLE, allocatable :: Nesup_old(:)

   ! Coupling vectors at differents times
   LIODBLE, allocatable :: sigma_old(:,:)
   LIODBLE, allocatable :: sigma_now(:,:)
   LIODBLE, allocatable :: sigma_0(:,:)
   LIODBLE, allocatable :: sigma_1(:,:)
  
   ! Nuclear Velocities of QM part
   LIODBLE, allocatable :: vel_old(:,:)

   ! Old phases in sigma
   LIODBLE, allocatable :: phases_old(:)
  
   ! Steps of dynamics
   integer :: tsh_nucStep = 0
   logical :: first_interp= .false.
   integer :: tsh_Enstep  = 20 ! input

   ! Kind of coupling to be performed
   integer :: type_coupling = 1 ! input

   ! Current State at Nuclear Step
   integer :: current_state = 0

   ! All States: GS + ES
   integer :: all_states = 0

   ! Electronic sub-time step variables
   TDCOMPLEX, allocatable :: coef_Stat(:,:), dot_Stat(:,:)
   LIODBLE, allocatable   :: elec_Coup(:,:,:), elec_Ene(:,:)
   LIODBLE, allocatable   :: elec_Pha(:,:,:)
   



end module fstsh_data
