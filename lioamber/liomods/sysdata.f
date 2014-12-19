!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       module sysdata
!--------------------------------------------------------------------!
       implicit none
       real*8 :: Energy_Total
       real*8 :: Dipole_Moment(3)

       real*8 :: External_Field(3)

       integer                :: Quantum_Atoms
       real*8,allocatable     :: Quantum_Coords(:,:)
       real*8,allocatable     :: Quantum_Forces(:,:)

       integer                :: Classic_Particles
       real*8,allocatable     :: Classic_Coords(:,:)
       real*8,allocatable     :: Classic_Forces(:,:)

       integer                :: Basis_Functions
       real*8,allocatable     :: Smat(:,:)
       real*8,allocatable     :: Fmat(:,:)
       complex*16,allocatable :: Pmat(:,:)

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
