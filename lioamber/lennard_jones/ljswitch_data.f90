!! VER EN AMBER: Opnq_LJ_atom_pair en SQM
! En params (nttyp, cn1, cn2)
! nttyp = ntypes*(ntypes+1)/2  (number of LJ types)
! ntypes es el total de combinaciones de LJ, por lo que
! solo se necesitan los elementos diagonales (i=i+1) de
! cn1 (A) y cn2 (B) de modo de sacar epsilon y sigma para
! cada átomo.

! En qm module, qmmm_struct: (iqmatoms, qm_mmpairs)
! iqmatoms tiene los indices de los átomos. 
! qm_mm_pairs tiene todos los pares mm-qm


!! Los MMTYPE van en las listas de epsilon y sigma. Estos vienen de
! ix, que es un array en memory_module.F90 aunque en ese mismo modulo
! tambien está atom type index en ese mismo modulo (importarlo)
! qmType=qmmm_struct%qm_atom_type(iqm)
! mmtype_for_iqm=qmmm_opnq%MM_atomType( qmmm_struct%iqmatoms(iqm) )
! jmm_index=qmmm_struct%qm_mm_pair_list(jmm)
! mmtype=qmmm_opnq%MM_atomType(jmm_index)


!! Number of pairs per QM atom. - length of pair_list.  
  ! integer :: qm_mm_pairs

  !! Non bond pair list for each QM atom
  ! integer, dimension(:), pointer :: qm_mm_pair_list => null()

#include "../datatypes/datatypes.fh"
module LJ_switch_data
   implicit none

   type lj_atom
      integer :: idx    = 0
      integer :: Z      = 0
      integer :: mmtype = 0
      LIODBLE :: q1     = 0.0D0
      LIODBLE :: q2     = 0.0D0
      LIODBLE :: s1     = 0.0D0
      LIODBLE :: s2     = 0.0D0
      LIODBLE :: e1     = 0.0D0
      LIODBLE :: e2     = 0.0D0
      LIODBLE :: eps    = 0.0D0
      LIODBLE :: sig    = 0.0D0
      LIODBLE :: deps   = 0.0D0
      LIODBLE :: dsig   = 0.0D0
      integer, allocatable :: basis_id(:)

      contains
         procedure, pass :: set_eps_sig
         procedure, pass :: kill => destroy_lj
   end type lj_atom

   type mm_atom
      integer :: mmtype = 0
      LIODBLE, allocatable :: dist(:)

      contains
         procedure, pass :: kill => destroy_mm
   end type mm_atom

   LIODBLE :: k_fermi = 10.0D0
   integer :: n_lj_atoms
   type(lj_atom), allocatable :: lj_atoms(:)
   type(mm_atom), allocatable :: mm_atoms(:)

   ! AMBER type sigma/eps, these arrays are sized ntypes.
   LIODBLE, allocatable :: mmlj_eps(:)
   LIODBLE, allocatable :: mmlj_sig(:)

contains

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

