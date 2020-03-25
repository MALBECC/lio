#include "datatypes/datatypes.fh"
module fstsh_data
   implicit none

   ! Trajectory Surface Hopping Variable
   logical :: FSTSH = .true.
   integer :: call_number = 0
   integer :: tsh_file = 456

   ! Nuclear Velocities
   LIODBLE, dimension(:,:), allocatable :: vel_old

   ! Overlap in differents times
   LIODBLE, dimension(:,:), allocatable :: Sovl_old, Sovl_now

   ! Old Variables needed in order to obtain overlap at differents times
   LIODBLE, allocatable :: a_old(:,:)
   LIODBLE, allocatable :: c_old(:,:)
   LIODBLE, allocatable :: r_old(:,:)

   ! Molecular Orbital variables
   LIODBLE, allocatable :: C_scf_old(:,:)
   LIODBLE, allocatable :: WFcis_old(:,:)
   LIODBLE, allocatable :: sigma_old(:,:)

   ! Old phases in sigma
   LIODBLE, allocatable :: phases_old(:)
  
   ! Steps of dynamics
   integer :: tsh_nucStep = 0
   integer :: tsh_Enstep  = 20

   ! Kind of coupling to be performed
   integer :: type_coupling = 0 

   ! Current State at Nuclear Step
   integer :: current_state = 0

   ! All States: GS + ES
   integer :: all_states = 0



end module fstsh_data
