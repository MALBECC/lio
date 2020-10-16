!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% LOWDIN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Performs a Löwdin Population Analysis and outputs atomic charges.            !
! Closed shell.
subroutine lowdin_cs(Pmat, Smat_sq, atom_of_func, atom_z, atom_q)
   ! Pmat        : Density matrix in atomic basis. Size Nfuncs x Nfuncs.
   ! Smat_sq     : Lowdin-orthogonalised overlap matrix. Size Nfuncs x Nfuncs.
   ! atom_of_func: Atom center of a given function. Size Nfuncs.
   ! atom_q      : Atomic charges output. Size Natoms
   implicit none
   integer , intent(in)  :: atom_of_func(:)
   integer , intent(in)  :: atom_z(:)
   LIODBLE , intent(in)  :: Pmat(:,:)
   LIODBLE , intent(in)  :: Smat_sq(:,:)
   LIODBLE , intent(out) :: atom_q(:)
   
   LIODBLE :: newterm
   integer :: iatom
   integer :: ifunc, jfunc, kfunc

   atom_q = dble(atom_z)
   do kfunc = 1, size(Pmat,1)
      newterm = 0.0D0
      
      do ifunc = 1, size(Pmat,1)
      do jfunc = 1, size(Pmat,1)
         newterm = newterm + Smat_sq(kfunc, ifunc) * &
                             Pmat(ifunc, jfunc) * Smat_sq(jfunc, kfunc)
      enddo
      enddo
      
      iatom = atom_of_func(kfunc)
      atom_q(iatom) = atom_q(iatom) - newterm
   enddo

   return
end subroutine lowdin_cs

! Open shell case.
subroutine lowdin_os(Pmat_a, Pmat_b, Smat_sq, atom_of_func, atom_z, atom_q, atom_s)
   ! Pmat_a      : Density matrix alpha in atomic basis. Size Nfuncs x Nfuncs.
   ! Pmat_b      : Density matrix beta in atomic basis. Size Nfuncs x Nfuncs.
   ! Smat        : Lowdin-orthogonalised overlap matrix. Size Nfuncs x Nfuncs.
   ! atom_of_func: Atom center of a given function. Size Nfuncs.
   ! atom_q      : Atomic charges output. Size Natoms
   ! atom_s      : Atomic spin population output. Size Natoms
   implicit none
   integer , intent(in)  :: atom_of_func(:)
   integer , intent(in)  :: atom_z(:)
   LIODBLE , intent(in)  :: Pmat_a(:,:)
   LIODBLE , intent(in)  :: Pmat_b(:,:)
   LIODBLE , intent(in)  :: Smat_sq(:,:)
   LIODBLE , intent(out) :: atom_q(:)
   LIODBLE , intent(out) :: atom_s(:)

   LIODBLE, allocatable :: temp_q(:)

   atom_q = 0.0D0
   atom_s = 0.0D0

   call lowdin_cs(Pmat_a, Smat_sq, atom_of_func, atom_z, atom_q)
   atom_s = atom_q

   allocate(temp_q(size(atom_q,1)))
   temp_q = 0.0D0
   call lowdin_cs(Pmat_b, Smat_sq, atom_of_func, atom_z, temp_q)

   atom_q = temp_q + atom_q - dble(atom_z)
   atom_s = temp_q - atom_s

   deallocate(temp_q)
end subroutine lowdin_os

! These routines are a general interface for property calculation
! and printing.
subroutine print_lowdin_cs(Pmat, Smat, atom_of_func, atom_z, real_z)
   use properties_data, only: UIDs, flowdin
   
   implicit none
   LIODBLE, intent(in) :: Pmat(:,:)
   LIODBLE, intent(in) :: Smat(:,:)
   integer, intent(in) :: atom_of_func(:)
   integer, intent(in) :: atom_z(:)
   integer, intent(in) :: real_z(:)

   LIODBLE, allocatable :: q(:)

   call g2g_timer_sum_start('Lowdin')
   allocate(q(size(atom_z,1)))
   
   call lowdin(Pmat, Smat, atom_of_func, atom_z, q)
   call write_population(real_z, q, 2, UIDs%mul, flowdin)

   deallocate(q)
   call g2g_timer_sum_pause('Lowdin')
end subroutine print_lowdin_cs

subroutine print_lowdin_os(Pmat, Pmat_b, Smat, atom_of_func, atom_z, real_z)
   use properties_data, only: UIDs, flowdin
   
   implicit none
   LIODBLE, intent(in) :: Pmat(:,:)
   LIODBLE, intent(in) :: Pmat_b(:,:)
   LIODBLE, intent(in) :: Smat(:,:)
   integer, intent(in) :: atom_of_func(:)
   integer, intent(in) :: atom_z(:)
   integer, intent(in) :: real_z(:)

   LIODBLE, allocatable :: q(:), s(:)
   character(len=100)   :: spinfile

   call g2g_timer_sum_start('Lowdin')
   allocate(q(size(atom_z,1)), s(size(atom_z,1)))
   
   call lowdin(Pmat, Pmat_b, Smat, atom_of_func, atom_z, q, s)
   call write_population(real_z, q, 2, UIDs%low, flowdin)

   spinfile = trim(flowdin) // "_spin"
   call write_population(real_z, s, 3, UIDs%lows, spinfile)

   deallocate(q, s)
   call g2g_timer_sum_pause('Lowdin')

end subroutine print_lowdin_os