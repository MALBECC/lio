!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       module sysdata
!--------------------------------------------------------------------!
       implicit none
       double precision :: Energy_Total
       double precision :: Dipole_Moment(3)

       double precision :: External_Field(3)

       integer                      :: Quantum_Atoms
       double precision,allocatable :: Quantum_Coords(:,:)
       double precision,allocatable :: Quantum_Forces(:,:)

       integer                      :: Classic_Particles
       double precision,allocatable :: Classic_Coords(:,:)
       double precision,allocatable :: Classic_Forces(:,:)

       integer                      :: Basis_Functions
       double precision,allocatable :: Smat(:,:)
       double precision,allocatable :: Fmat(:,:)
       complex*16,allocatable       :: Pmat(:,:)

       contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       subroutine set_quantum(num)
!--------------------------------------------------------------------!
       implicit none
       integer,intent(in) :: num
       Quantum_Atoms=num
       if (num.gt.0) then
         allocate(Quantum_Coords(3,num))
         allocate(Quantum_Forces(3,num))
       else
         if (allocated(Quantum_Coords)) deallocate(Quantum_Coords)
         if (allocated(Quantum_Forces)) deallocate(Quantum_Forces)
       endif
       return;end subroutine
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       subroutine set_classic(num)
!--------------------------------------------------------------------!
       implicit none
       integer,intent(in) :: num
       Classic_Particles=num
       if (num.gt.0) then
         allocate(Classic_Coords(3,num))
         allocate(Classic_Forces(3,num))
       else
         if (allocated(Classic_Coords)) deallocate(Classic_Coords)
         if (allocated(Classic_Forces)) deallocate(Classic_Forces)
       endif
       return;end subroutine
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       subroutine set_basis(num)
!--------------------------------------------------------------------!
       implicit none
       integer,intent(in) :: num
       Basis_Functions=num
       if (num.gt.0) then
         allocate(Smat(num,num))
         allocate(Fmat(num,num))
         allocate(Pmat(num,num))
       else
         if (allocated(Smat)) deallocate(Smat)
         if (allocated(Fmat)) deallocate(Fmat)
         if (allocated(Pmat)) deallocate(Pmat)
       endif
       return;end subroutine
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       end module
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
