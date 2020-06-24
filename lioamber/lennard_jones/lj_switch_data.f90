#include "../datatypes/datatypes.fh"
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%% LJ SWITCH DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! This module contains necesary data for the LJ Switch module, see             !
! lj_switch.f90 for further details.                                           !
!                                                                              !
! This module also defines the datatypes lj_atom (used for those QM atoms with !
! LJ parameters depending on atomic charge) and mm_atoms (used for the QM/MM   !
! LJ interactions). These types both have a %destroy (aliased as %kill) method !
! which deallocates all data contained within said datatypes. In addition, the !
! method %set_eps_sig calculates the LJ parameters (and their derivatives) for !
! a given QM atom, receiving the atomic partial charge as an input.            !
!                                                                              !
! First written by: Federico Pedron, Jun/2020                                  !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module LJ_switch_data
   implicit none

   type lj_atom
      integer :: idx    = 0         ! QM atom index in LIO arrays (r, Iz, etc).
      integer :: Z      = 0         ! Atomic number of the QM atom.
      integer :: mmtype = 0         ! MM LJ type, imported from the MM software.
      LIODBLE :: q1     = 0.0D0     ! Reference charge for atom state 1.
      LIODBLE :: q2     = 0.0D0     ! Reference charge for atom state 2.
      LIODBLE :: s1     = 0.0D0     ! Reference sigma for atom state 1.
      LIODBLE :: s2     = 0.0D0     ! Reference sigma for atom state 2.
      LIODBLE :: e1     = 0.0D0     ! Reference epsilon (*4) for atom state 1.
      LIODBLE :: e2     = 0.0D0     ! Reference epsilon (*4) for atom state 2.
      LIODBLE :: eps    = 0.0D0     ! Epsilon (*4) for current atomic charge.
      LIODBLE :: sig    = 0.0D0     ! Sigma for current atomic charge.
      LIODBLE :: deps   = 0.0D0     ! 4*depsilon/dQ for current atomic charge.
      LIODBLE :: dsig   = 0.0D0     ! dSigma/dQ for current atomic charge.

      ! This array contains the indeces of the basis functions belonging to
      ! this atom. This eases the calculations of partial atomic charges.
      integer, allocatable :: basis_id(:)

      contains
         ! See below for these pŕocedures.
         procedure, pass :: set_eps_sig
         procedure, pass :: kill => destroy_lj
   end type lj_atom

   type mm_atom
      integer :: mmtype = 0            ! MM LJ type, imported from the MM software.
      LIODBLE, allocatable :: dist(:)  ! Distance to the nth QM atom with
                                       ! variable LJ parameters (i.e. lj_atom).

      contains
         ! See below for these pŕocedures.
         procedure, pass :: kill => destroy_mm
   end type mm_atom

   LIODBLE :: k_fermi = 10.0D0                 ! Exponent constant for the 
                                               ! Fermi interpolation function.
   integer :: n_lj_atoms                       ! Number of QM atoms with variable LJ.
   type(lj_atom), allocatable :: lj_atoms(:)   ! List of LJ atoms (size n_lj_atoms)
   type(mm_atom), allocatable :: mm_atoms(:)   ! List of MM atoms in a given step.

   ! AMBER type sigma/eps, these arrays are sized ntypes.
   LIODBLE, allocatable :: mmlj_eps(:)  ! Epsilon for a given MM type. Imported
                                        ! from MM software.
   LIODBLE, allocatable :: mmlj_sig(:)  ! Epsilon for a given MM type. Imported
                                        ! from MM software.

contains
   ! set_eps_sig calculates, epsilon, sigma and their derivatives with respect
   ! to atomic charge, using a sigmoid interpolation between two states.
   subroutine set_eps_sig(this, chrg)
      implicit none
      LIODBLE       , intent(in)    :: chrg
      class(lj_atom), intent(inout) :: this

      LIODBLE :: exp_term, inv_exp

      exp_term = exp( -k_fermi * ( chrg - 0.5D0 * (this%q2 + this%q1) ) )
      inv_exp  = 1.0D0 / (1.0D0 + exp_term)

      ! Base values
      this%sig = this%s1 + inv_exp * (this%s2 - this%s1)
      this%eps = this%e1 + inv_exp * (this%e2 - this%e1)

      ! Derivatives dSig/dQ and dEps/dQ
      inv_exp = k_fermi * exp_term * inv_exp * inv_exp

      this%dsig = inv_exp * (this%s2 - this%s1)
      this%deps = inv_exp * (this%e2 - this%e1)      

   end subroutine set_eps_sig

   ! The folowing destroy_xx procedures deallocate arrays in 
   ! the custom atom types.
   subroutine destroy_lj(this)
      implicit none
      class(lj_atom), intent(inout) :: this
      
      if (allocated(this%basis_id)) deallocate(this%basis_id)
   end subroutine destroy_lj

   subroutine destroy_mm(this)
      implicit none
      class(mm_atom), intent(inout) :: this
      
      if (allocated(this%dist)) deallocate(this%dist)
   end subroutine destroy_mm

end module LJ_switch_data

