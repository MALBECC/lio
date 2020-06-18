#include "../datatypes/datatypes.fh"
#include "../lennard_jones/ljswitch_data.f90"
#include "../lennard_jones/ljswitch.f90"

program test_ljswitch

   write(*,'(A)') "== Testing module LJ_SWITCH =="
   call test_init_end()  ! [INI]


end program test_ljswitch

! Verfies proper module initialisation and deinitialisation. [INI]
subroutine test_init_end()
   use lj_switch_data, only: n_lj_atoms, lj_atoms, mm_atoms, mmlj_eps, mmlj_sig
   use lj_switch     , only: ljs_initialise, ljs_finalise

   logical :: passed
   integer :: counter, counter2
   integer, allocatable :: lio_nuc(:), lio_iz(:), result_z(:), result_bas(:,:)
   LIODBLE, allocatable :: eps_ext(:), sig_ext(:)

   ! Variable setup
   n_lj_atoms = 3
   allocate(lj_atoms(n_lj_atoms))
   allocate(lio_iz(5))
   allocate(lio_nuc(10))
   allocate(eps_ext(10), sig_ext(10))

   lio_iz  = (/1, 4, 9, 16, 25/)
   lio_nuc = (/1, 3, 5, 5, 1, 3, 0, 3, 1, 0/)
   lj_atoms(1)%idx = 1; lj_atoms(2)%idx = 3; lj_atoms(3)%idx = 5

   do counter = 1, 10
      eps_ext(counter) = 1.0D0 / dble(counter)
      sig_ext(counter) = 2.0D0 * dble(counter)
   enddo
   write(*,'(A)') "Subroutine ljs_initialise()"
   call ljs_initialise(eps_ext, sig_ext, lio_iz, lio_nuc)
   
   ! Test if atom Z has been properly assigned.
   passed = .true.
   allocate(result_z(3))
   result_z = (/1,9,25/)
   do counter = 1, 3
      if (lj_atoms(counter)%Z /=  result_z(counter)) then
         passed = .false.
         write(*,'(A30,I1)') "ERROR: Wrong atomic Z in atom ", counter
         write(*,'(7x,A9,I2,A11,I2)') "Expected ", result_z(counter), &
                                       " and found ", lj_atoms(counter)%Z
      endif
   enddo
   deallocate(result_z)
   if (passed) write(*,'(A)') "=> Atomic Z properly assigned."

   ! Tests proper assignation of basis functions.
   passed = .true.
   allocate(result_bas(3,3))
   result_bas = reshape((/1,5,9,2,6,8,3,4,0/), shape(result_bas))
   do counter  = 1, 3
   do counter2 = 1, size(lj_atoms(counter)%basis_id,1)
      if (lj_atoms(counter)%basis_id(counter2) /=  &
          result_bas(counter2, counter) ) then
         passed = .false.
         write(*,'(A30,I1,A8,I1)') "ERROR: Wrong basis in atom ", counter, &
                                   ", index ", counter2
         write(*,'(7x,A9,I2,A11,I2)') "Expected ", &
                                       result_bas(counter2, counter), &
                                       " and found ", &
                                       lj_atoms(counter)%basis_id(counter2)
      endif
   enddo
   enddo 
   deallocate(result_bas)
   if (passed) write(*,'(A)') "=> Atomic basis properly assigned."

   write(*,'(A)') ""

   ! Test all deinitialisations.
   write(*,'(A)') "Subroutine ljs_finalise()"
   allocate(mm_atoms(10))
   call ljs_finalise()
   
   if (.not. allocated(lj_atoms)) then
      write(*,'(A)') "=> lj_atoms deallocated."
   else
      write(*,'(A)') "ERROR: lj_atoms not deallocated."
   endif

   if (.not. allocated(mm_atoms)) then
      write(*,'(A)') "=> mm_atoms deallocated."
   else
      write(*,'(A)') "ERROR: mm_atoms not deallocated."
   endif

   if (.not. allocated(mmlj_sig)) then
      write(*,'(A)') "=> mmlj_sig deallocated."
   else
      write(*,'(A)') "ERROR: mmlj_sig not deallocated."
   endif
   
   if (.not. allocated(mmlj_eps)) then
      write(*,'(A)') "=> mmlj_eps deallocated."
   else
      write(*,'(A)') "ERROR: mmlj_eps not deallocated."
   endif

   if (n_lj_atoms == 0) then
      write(*,'(A)') "=> n_lj_atoms set to zero."
   else
      write(*,'(A)') "ERROR: n_lj_atoms not set to zero."
   endif
   write(*,'(A)') ""
end subroutine test_init_end