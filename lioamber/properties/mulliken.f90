!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% MULLIKEN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Performs a Mulliken Population Analysis, outputing atomic charges.           !
! Closed shell case.
subroutine mulliken_cs(Pmat, Smat, atom_of_func, atom_q)
   ! Pmat        : Density matrix in atomic basis. Size Nfuncs x Nfuncs.
   ! Smat        : Atomic basis overlap matrix. Size Nfuncs x Nfuncs.
   ! atom_of_func: Atom center of a given function. Size Nfuncs.
   ! atom_q      : Atomic charges output. Size Natoms
   implicit none
   integer , intent(in)  :: atom_of_func(:)
   LIODBLE , intent(in)  :: Pmat(:,:)
   LIODBLE , intent(in)  :: Smat(:,:)
   LIODBLE , intent(out) :: atom_q(:)

   integer :: ifunc, jfunc, iatom

   atom_q = 0.0D0
   do ifunc =1, size(Pmat,2)
      iatom = atom_of_func(ifunc)

      do jfunc = 1, size(Pmat,2)
         atom_q(iatom) = atom_q(iatom) &
                       - Pmat(ifunc, jfunc) * Smat(ifunc, jfunc)
      enddo
   enddo
end subroutine mulliken_cs

! Open shell case.
subroutine mulliken_os(Pmat_a, Pmat_b, Smat, atom_of_func, atom_q, atom_s)
   ! Pmat_a      : Density matrix alpha in atomic basis. Size Nfuncs x Nfuncs.
   ! Pmat_b      : Density matrix beta in atomic basis. Size Nfuncs x Nfuncs.
   ! Smat        : Atomic basis overlap matrix. Size Nfuncs x Nfuncs.
   ! atom_of_func: Atom center of a given function. Size Nfuncs.
   ! atom_q      : Atomic charges output. Size Natoms
   ! atom_s      : Atomic spin population output. Size Natoms
   implicit none
   integer , intent(in)  :: atom_of_func(:)
   LIODBLE , intent(in)  :: Pmat_a(:,:)
   LIODBLE , intent(in)  :: Pmat_b(:,:)
   LIODBLE , intent(in)  :: Smat(:,:)
   LIODBLE , intent(out) :: atom_q(:)

   LIODBLE, allocatable :: temp_q(:)
   integer :: ifunc, jfunc, iatom

   atom_q = 0.0D0
   atom_s = 0.0D0

   call mulliken_cs(Pmat_a, Smat, atom_of_func, atom_q)
   atom_s = atom_q

   allocate(temp_q(size(atom_q,1)))
   temp_q = 0.0D0
   call mulliken_cs(Pmat_b, Smat, atom_of_func, atom_q)

   atom_q = temp_q + atom_q
   atom_s = temp_q - atom_s

   deallocate(temp_q)

end subroutine mulliken_os


