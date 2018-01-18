module fsiesta

! Support routines for siesta-as-a-subroutine in Unix/Linux.
! The routines that handle the other side of the communication are
! in module iopipes of siesta program.
! Usage:
!   call siesta_launch( label, nnodes )
!     character(len=*),intent(in) :: label  : Name of siesta process
!                                             (prefix of its .fdf file)
!     integer,optional,intent(in) :: nnodes : Number of MPI nodes
!
!   call siesta_units( length, energy )
!     character(len=*),intent(in) :: length : Physical unit of length
!     character(len=*),intent(in) :: energy : Physical unit of energy
!
!   call siesta_forces( label, na, xa, cell, energy, fa, stress )
!     character(len=*), intent(in) :: label      : Name of siesta process
!     integer,          intent(in) :: na         : Number of atoms
!     real(dp),         intent(in) :: xa(3,na)   : Cartesian coords
!     real(dp),optional,intent(in) :: cell(3,3)  : Unit cell vectors
!     real(dp),optional,intent(out):: energy     : Total energy
!     real(dp),optional,intent(out):: fa(3,na)   : Atomic forces
!     real(dp),optional,intent(out):: stress(3,3): Stress tensor
!   call siesta_quit( label )
!     character(len=*),intent(in) :: label  : Name of one siesta process,
!                                             or 'all' to stop all procs.
! Behaviour:
! - If nnodes is not present among siesta_launch arguments, or nnodes<2, 
!   a serial siesta process will be launched. Otherwise, a parallel 
!   mpirun process will be launched.
! - If siesta_units is not called, length='Ang', energy='eV' are
!   used by default. If it is called more than once, the units in the
!   last call become in effect.
! - The physical units set by siesta_units are used for all the siesta
!   processes launched
! - If siesta_forces is called without a previous call to siesta_launch
!   for that label, it assumes that the siesta process has been launched
!   (and the communication pipes created) externally in the shell.
!   In this case, siesta_forces only opens its end of the pipes and begins
!   communication through them.
! - If argument cell is not present in the call to siesta_forces, or if
!   the cell has zero volume, it is assumed that the system is a molecule,
!   and a supercell is generated automatically by siesta so that the 
!   different images do not overlap. In this case the stress returned
!   has no physical meaning.
! - The stress is defined as dE/d(strain)/Volume, with a positive sign
!   when the system tends to contract (negative pressure)
! - The following events result in a stopping error message:
!   - siesta_launch is called twice with the same label
!   - siesta_forces finds a communication error trough the pipes
!   - siesta_quit is called without a prior call to siesta_launch or
!     siesta_forces for that label
! - If siesta_quit is not called for a launched siesta process, that
!   process will stay listening indefinitedly to the pipe and will need
!   to be killed in the shell.
! - siesta_units may be called either before or after siesta_launch
! J.M.Soler and A.Garcia. Nov.2003

  implicit none

PRIVATE ! Nothing is declared public beyond this point

! Holds data on siesta processes and their communication pipes
!  type proc
!    private
!    character(len=80) :: label     ! Name of process
!    integer           :: iuc, iuf  ! I/O units for coords/forces commun.
!  end type proc

! Global module variables
!  integer, parameter :: max_procs = 100
!  integer, parameter :: dp = kind(1.d0)
!  type(proc),   save :: p(max_procs)
!  integer,      save :: np=0
!  character(len=32), save :: xunit = 'Ang'
!  character(len=32), save :: eunit = 'eV'
!  character(len=32), save :: funit = 'eV/Ang'
!  character(len=32), save :: sunit = 'eV/Ang**3'

end module fsiesta
