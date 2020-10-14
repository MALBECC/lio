!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% LOWDIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Performs a Löwdin Population Analysis and outputs atomic charges.            !
! Closed shell.
subroutine lowdin_cs(Pmat, Smat_sq, atom_of_func, atom_q)
   ! Pmat        : Density matrix in atomic basis. Size Nfuncs x Nfuncs.
   ! Smat_sq     : Lowdin-orthogonalised overlap matrix. Size Nfuncs x Nfuncs.
   ! atom_of_func: Atom center of a given function. Size Nfuncs.
   ! atom_q      : Atomic charges output. Size Natoms
   implicit none
   integer , intent(in)  :: atom_of_func(:)
   LIODBLE , intent(in)  :: Pmat(:,:)
   LIODBLE , intent(in)  :: Smat_sq(:,:)
   LIODBLE , intent(out) :: atom_q(:)
   
   LIODBLE :: newterm
   integer :: iatom
   integer :: ifunc, jfunc, kfunc

   do kfunc = 1, M
      newterm = 0.0D0
      
      do ifunc = 1, M
      do jfunc = 1, M
         newterm = newterm + Smat_sq(kfunc, ifunc) * &
                             Pmat(ifunc, jfunc) * Smat_sq(jfunc, kfunc)
      enddo
      enddo
      
      iatom = atom_of_func(k)
      atom_q(iatom) = atom_q(iatom) - newterm
   enddo

   return
end subroutine lowdin_cs

! Open shell case.
subroutine lowdin_os(Pmat_a, Pmat_b, Smat_sq, atom_of_func, atom_q, atom_s)
   ! Pmat_a      : Density matrix alpha in atomic basis. Size Nfuncs x Nfuncs.
   ! Pmat_b      : Density matrix beta in atomic basis. Size Nfuncs x Nfuncs.
   ! Smat        : Lowdin-orthogonalised overlap matrix. Size Nfuncs x Nfuncs.
   ! atom_of_func: Atom center of a given function. Size Nfuncs.
   ! atom_q      : Atomic charges output. Size Natoms
   ! atom_s      : Atomic spin population output. Size Natoms
   implicit none
   integer , intent(in)  :: atom_of_func(:)
   LIODBLE , intent(in)  :: Pmat_a(:,:)
   LIODBLE , intent(in)  :: Pmat_b(:,:)
   LIODBLE , intent(in)  :: Smat_sq(:,:)
   LIODBLE , intent(out) :: atom_q(:)

   LIODBLE, allocatable :: temp_q(:)
   integer :: ifunc, jfunc, iatom

   atom_q = 0.0D0
   atom_s = 0.0D0

   call lowdin_cs(Pmat_a, Smat_sq, atom_of_func, atom_q)
   atom_s = atom_q

   allocate(temp_q(size(atom_q,1)))
   temp_q = 0.0D0
   call lowdin_cs(Pmat_b, Smat_sq, atom_of_func, atom_q)

   atom_q = temp_q + atom_q
   atom_s = temp_q - atom_s

   deallocate(temp_q)
end subroutine lowdin_os
