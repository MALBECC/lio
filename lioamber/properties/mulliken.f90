!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% MULLIKEN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Performs a Mulliken Population Analysis, outputing atomic charges.           !
! Closed shell case.
subroutine mulliken_cs(Pmat, Smat, atom_of_func, atom_z, atom_q)
   ! Pmat        : Density matrix in atomic basis. Size Nfuncs x Nfuncs.
   ! Smat        : Atomic basis overlap matrix. Size Nfuncs x Nfuncs.
   ! atom_of_func: Atom center of a given function. Size Nfuncs.
   ! atom_q      : Atomic charges output. Size Natoms
   implicit none
   integer , intent(in)  :: atom_of_func(:)
   integer , intent(in)  :: atom_z(:)
   LIODBLE , intent(in)  :: Pmat(:,:)
   LIODBLE , intent(in)  :: Smat(:,:)
   LIODBLE , intent(out) :: atom_q(:)

   integer :: ifunc, jfunc, iatom

   atom_q = dble(atom_z)
   do ifunc =1, size(Pmat,2)
      iatom = atom_of_func(ifunc)

      do jfunc = 1, size(Pmat,2)
         atom_q(iatom) = atom_q(iatom) &
                       - Pmat(ifunc, jfunc) * Smat(ifunc, jfunc)
      enddo
   enddo
end subroutine mulliken_cs

! Open shell case.
subroutine mulliken_os(Pmat_a, Pmat_b, Smat, atom_of_func, atom_z, atom_q, atom_s)
   ! Pmat_a      : Density matrix alpha in atomic basis. Size Nfuncs x Nfuncs.
   ! Pmat_b      : Density matrix beta in atomic basis. Size Nfuncs x Nfuncs.
   ! Smat        : Atomic basis overlap matrix. Size Nfuncs x Nfuncs.
   ! atom_of_func: Atom center of a given function. Size Nfuncs.
   ! atom_q      : Atomic charges output. Size Natoms
   ! atom_s      : Atomic spin population output. Size Natoms
   implicit none
   integer , intent(in)  :: atom_of_func(:)
   integer , intent(in)  :: atom_z(:)
   LIODBLE , intent(in)  :: Pmat_a(:,:)
   LIODBLE , intent(in)  :: Pmat_b(:,:)
   LIODBLE , intent(in)  :: Smat(:,:)
   LIODBLE , intent(out) :: atom_q(:)
   LIODBLE , intent(out) :: atom_s(:)

   LIODBLE, allocatable :: temp_q(:)

   atom_q = 0.0D0
   atom_s = 0.0D0

   call mulliken_cs(Pmat_a, Smat, atom_of_func, atom_z, atom_q)
   atom_s = atom_q

   allocate(temp_q(size(atom_q,1)))
   temp_q = 0.0D0
   call mulliken_cs(Pmat_b, Smat, atom_of_func, atom_z, temp_q)

   atom_q = temp_q + atom_q - dble(atom_z)
   atom_s = temp_q - atom_s

   deallocate(temp_q)

end subroutine mulliken_os

! These routines are a general interface for property calculation
! and printing.
subroutine print_mulliken_cs(Pmat, Smat, atom_of_func, atom_z, real_z)
   use properties_data, only: UIDs, fmulliken
   
   implicit none
   LIODBLE, intent(in) :: Pmat(:,:)
   LIODBLE, intent(in) :: Smat(:,:)
   integer, intent(in) :: atom_of_func(:)
   integer, intent(in) :: atom_z(:)
   integer, intent(in) :: real_z(:)

   LIODBLE, allocatable :: q(:)

   call g2g_timer_sum_start('Mulliken')
   allocate(q(size(atom_z,1)))
   
   call mulliken(Pmat, Smat, atom_of_func, atom_z, q)
   call write_population(real_z, q, 0, UIDs%mul, fmulliken)

   deallocate(q)
   call g2g_timer_sum_pause('Mulliken')
end subroutine print_mulliken_cs

subroutine print_mulliken_os(Pmat, Pmat_b, Smat, atom_of_func, atom_z, real_z)
   use properties_data, only: UIDs, fmulliken
   
   implicit none
   LIODBLE, intent(in) :: Pmat(:,:)
   LIODBLE, intent(in) :: Pmat_b(:,:)
   LIODBLE, intent(in) :: Smat(:,:)
   integer, intent(in) :: atom_of_func(:)
   integer, intent(in) :: atom_z(:)
   integer, intent(in) :: real_z(:)

   LIODBLE, allocatable :: q(:), s(:)
   character(len=100)   :: spinfile

   call g2g_timer_sum_start('Mulliken')
   allocate(q(size(atom_z,1)), s(size(atom_z,1)))
   
   call mulliken(Pmat, Pmat_b, Smat, atom_of_func, atom_z, q, s)
   call write_population(real_z, q, 0, UIDs%mul, fmulliken)

   spinfile = trim(fmulliken) // "_spin"
   call write_population(real_z, s, 1, UIDs%muls, spinfile)

   deallocate(q, s)
call g2g_timer_sum_pause('Mulliken')

end subroutine print_mulliken_os