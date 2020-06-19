#include "../datatypes/datatypes.fh"
#include "../lennard_jones/ljswitch_data.f90"
#include "../lennard_jones/ljswitch.f90"

program test_ljswitch
   write(*,'(A)') ""
   write(*,'(A)') "== Testing module LJ_SWITCH =="
   write(*,'(A)') ""

   ! The abbreviations between [] provide a way to easily
   ! navigate this file by using find (usually / or ctrl+f).
   call test_init_end()                      ! [INI]
   call test_mm_setting()                    ! [MMS]
   call test_mm_interface()                  ! [MMI]
   call test_crg_and_derivs()                ! [CRG]
   call test_fock_terms()                    ! [FOK]

end program test_ljswitch

! Initialisation and finalisation [INI]
! Verfies proper module initialisation and deinitialisation.
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
   eps_ext(5) = 0.0D0
   sig_ext(5) = 0.0D0

   write(*,'(A)') "Subroutine ljs_initialise"
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

   ! Checks input from external epsilon and sigma lists.
   passed = .true.
   do counter = 1, 10
      if (mmlj_eps(counter) /= eps_ext(counter)) then
         passed = .false.
         write(*,'(A30,I1)') "ERROR: Wrong epsilon in index ", counter
         write(*,'(7x,A9,F14.7,A11,F14.7)') "Expected ", eps_ext(counter), &
                                            " and found ", mmlj_eps(counter)
      endif
   enddo
   if (passed) write(*,'(A)') "=> Epsilon passed correctly."
   passed = .true.
   do counter = 1, 10
      if (mmlj_sig(counter) /= sig_ext(counter)) then
         passed = .false.
         write(*,'(A30,I1)') "ERROR: Wrong epsilon in index ", counter
         write(*,'(7x,A9,F14.7,A11,F14.7)') "Expected ", sig_ext(counter), &
                                            " and found ", mmlj_sig(counter)
      endif
   enddo
   if (passed) write(*,'(A)') "=> Sigma passed correctly."
   

   deallocate(lio_iz, lio_nuc, eps_ext, sig_ext)
   write(*,'(A)') ""

   ! Test all deinitialisations.
   write(*,'(A)') "Subroutine ljs_finalise"
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

! MM Setting [MMS]
! This tests the proper setup of incoming data from the MM
! software.
subroutine test_mm_setting()
   use lj_switch_data, only: lj_atoms, mm_atoms
   use lj_switch     , only: ljs_settle_mm, ljs_finalise

   implicit none
   integer, allocatable :: qm_typ(:), mm_typ(:)
   LIODBLE, allocatable :: pos(:,:), dist_res(:,:)

   integer :: nQM, nMM, counter, counter2, qm_res(2)
   LIODBLE :: sq2, sq3, diff
   logical :: passed
   
   nQM = 5
   nMM = 7
   allocate(qm_typ(nQM), mm_typ(nMM), pos(nQM+nMM,3))
   allocate(lj_atoms(2))
   lj_atoms(1)%idx = 1; lj_atoms(2)%idx = 3

   qm_typ = (/5,2,1,3,4/)
   mm_typ = (/5,4,3,2,1,1,2/)
 
   pos = reshape((/0.0D0, 1.0D0, 1.0D0, 2.0D0, 3.0D0, 1.0D0,   &
                   0.0D0, 0.0D0, 1.0D0, 1.0D0, 0.0D0, 1.0D0,   &
                   0.0D0, 1.0D0, 2.0D0, 2.0D0, 3.0D0, 0.0D0,   &
                   1.0D0, 0.0D0, 1.0D0, 0.0D0, 1.0D0, 1.0D0,   &
                   0.0D0, 1.0D0, 3.0D0, 2.0D0, 3.0D0, 0.0D0,   &
                   0.0D0, 1.0D0, 0.0D0, 1.0D0, 1.0D0, 1.0D0/), &
                   shape(pos))
  
   write(*,'(A)') "Subroutine ljs_settle_mm"         
   call ljs_settle_mm(qm_typ, mm_typ, pos)

   ! Test if atom types have been properly assigned.
   passed = .true.
   qm_res = (/5,1/)
   do counter = 1, 2
      if (lj_atoms(counter)%mmtype /=  qm_res(counter)) then
         passed = .false.
         write(*,'(A34,I1)') "ERROR: Wrong atom type in QM atom ", counter
         write(*,'(7x,A9,I2,A11,I2)') "Expected ", qm_res(counter), &
                                      " and found ", lj_atoms(counter)%mmtype
      endif
   enddo
   if (passed) write(*,'(A)') "=> QM atom types correctly assigned."
   passed = .true.
   do counter = 1, nMM
      if (mm_atoms(counter)%mmtype /=  mm_typ(counter)) then
         passed = .false.
         write(*,'(A34,I1)') "ERROR: Wrong atom type in MM atom ", counter
         write(*,'(7x,A9,I2,A11,I2)') "Expected ", mm_typ(counter), &
                                      " and found ", mm_atoms(counter)%mmtype
      endif
   enddo
   if (passed) write(*,'(A)') "=> MM atom types correctly assigned."

   ! Verifies distances between QM and MM atoms.
   passed = .true.
   allocate(dist_res(nMM,2))
   sq2 = sqrt(2.0D0)
   sq3 = sqrt(3.0D0)
   dist_res = reshape((/1.0D0, 1.0D0, 1.0D0, sq2, sq2, sq2, sq3,         &
                        sqrt(13.0D0), sqrt(11.0D0), 3.0D0, sqrt(10.0D0), &
                        2.0D0*sq2, sq2*sq3, sqrt(5.0D0)/), shape(dist_res))
   do counter = 1, nMM
   do counter2 = 1, 2
      diff = abs(mm_atoms(counter)%dist(counter2) - dist_res(counter,counter2))
      if (diff > 1.0D-10) then
         passed = .false.
         write(*,'(A39,I1,A13,I1)') "ERROR: Wrong distances between MM atom ",&
                                    counter, " and LJ atom ", counter2
         write(*,'(7x,A9,F14.7,A11,F14.7)') "Expected ", &
                                    dist_res(counter,counter2), " and found ",&
                                    mm_atoms(counter)%dist(counter2)
         
      endif
   enddo
   enddo
   deallocate(dist_res)
   if (passed) write(*,'(A)') "=> MM-QM distances calculated correctly."

   write(*,'(A)') ""
   deallocate(qm_typ, mm_typ, pos)
   call ljs_finalise()
end subroutine test_mm_setting

! MM Interface [MMI]
! Verifies the values LIO returns to the MM software as energy
! and gradients.
subroutine test_mm_interface()
   use lj_switch_data, only: lj_atoms, mmlj_eps, mmlj_sig
   use lj_switch     , only: ljs_settle_mm, ljs_finalise, &
                             ljs_substract_mm, ljs_gradients_qmmm

   implicit none
   integer, allocatable :: qm_typ(:), mm_typ(:)
   LIODBLE, allocatable :: pos(:,:), gradmm(:,:), gradqm(:,:), grad_res(:,:)

   integer :: nQM, nMM, counter, counter2
   LIODBLE :: energy, ener_res, diff
   logical :: passed
   
   nQM = 5
   nMM = 3
   allocate(qm_typ(nQM), mm_typ(nMM), pos(nQM+nMM,3))
   allocate(lj_atoms(2))
   lj_atoms(1)%idx = 1; lj_atoms(2)%idx = 3

   qm_typ = (/1,2,1,3,4/)
   mm_typ = (/1,4,3/)
 
   pos = reshape((/0.0D0, 1.0D0, 1.0D0, 2.0D0, 3.0D0, &
                   1.0D0, 0.0D0, 0.0D0, &
                   0.0D0, 1.0D0, 2.0D0, 2.0D0, 3.0D0, &
                   0.0D0, 1.0D0, 0.0D0, &
                   0.0D0, 1.0D0, 3.0D0, 2.0D0, 3.0D0, &
                   0.0D0, 0.0D0, 1.0D0/), shape(pos))
  
   call ljs_settle_mm(qm_typ, mm_typ, pos)

   allocate(mmlj_eps(5), mmlj_sig(5))
   do counter = 1, 5
      mmlj_eps(counter) = 1.0D0 / dble(counter)
      mmlj_sig(counter) = 2.0D0 + dble(counter) * 0.2D0
   enddo
   mmlj_eps(3) = 0.0D0
   mmlj_sig(3) = 0.0D0

   allocate(gradqm(3,nQM), gradmm(3,nMM))
   gradqm = 0.0D0; gradmm = 0.0D0; energy = 0.0D0

   write(*,'(A)') "Subroutine ljs_substract_mm"
   call ljs_substract_mm(energy, gradqm, gradmm, pos, nQM)

   ! Checks total energy.
   ener_res = -42421.7509682040D0
   diff = abs(energy - ener_res) / abs(energy + ener_res)
   if (diff > 1.0D-9) then
      passed = .false.
      write(*,'(A19)') "ERROR: Wrong energy"
      write(*,'(7x,A9,F15.7,A11,F15.7)') "Expected ", ener_res, &
                                         " and found ", energy
   else
      passed = .true.
      write(*,'(A)') "=> Energy calculated correctly."
   endif

   ! Checks QM region gradients.
   allocate(grad_res(3,nQM))
   passed = .true.
   grad_res = reshape((/-153579.7521485972D0, -356895.44677734373D0, 0.0D0, &
                         0.0D0, 0.0D0, 0.0D0,    &
                        -3.1673435152421874D-2, -7.4393611725054196D-2,    & 
                        -0.15910057031621411D0, &
                         0.0D0, 0.0D0, 0.0D0,    &
                         0.0D0, 0.0D0, 0.0D0/), shape(grad_res))
   do counter = 1, 3
   do counter2 = 1, nQM
      diff = abs(grad_res(counter,counter2) - gradqm(counter,counter2)) / &
             abs(grad_res(counter,counter2) + gradqm(counter,counter2)) 
      if (diff > 1.0D-13) then
         passed = .false.
         write(*,'(A34,I1,A14,I1)') "ERROR: Wrong gradient for QM atom ",&
                                    counter2, " in direction ", counter
         write(*,'(7x,A9,F15.7,A11,F15.7)') &
                                "Expected ", grad_res(counter,counter2), &
                                " and found ", gradqm(counter,counter2)
      endif
   enddo
   enddo
   deallocate(grad_res)
   if (passed) write(*,'(A)') "=> QM region gradients calculated correctly."

   ! Checks MM region gradients.
   allocate(grad_res(3,nMM))
   passed = .true.
   grad_res = reshape((/153579.7521485907D0, 4.2720176572632322D-2, &
                        6.4080264858948480D-2, &
                        3.1673435152421874D-2, 356895.47845077893D0, &
                        9.5020305457265614D-2, &
                        0.0D0, 0.0D0, 0.0D0/), shape(grad_res))
   do counter = 1, 3
   do counter2 = 1, nMM
      diff = abs(grad_res(counter,counter2) - gradmm(counter,counter2)) / &
             abs(grad_res(counter,counter2) + gradmm(counter,counter2))
             
      if (diff > 1.0D-13) then
         passed = .false.
         write(*,'(A34,I1,A14,I1)') "ERROR: Wrong gradient for MM atom ",&
                                    counter2, " in direction ", counter
         write(*,'(7x,A9,F15.7,A11,F15.7)') &
                                "Expected ", grad_res(counter,counter2), &
                                " and found ", gradmm(counter,counter2)
         
      endif
   enddo
   enddo
   deallocate(grad_res)
   if (passed) write(*,'(A)') "=> MM region gradients calculated correctly."
   write(*,'(A)') ""

   ! Makes the same calculation with the internal LJ parameters, adjusted 
   ! by charge. For the purpose of this test, the results should be the same
   ! as above.
   lj_atoms(1)%sig = 2.2D0; lj_atoms(1)%eps = 1.0D0
   lj_atoms(2)%sig = 2.2D0; lj_atoms(2)%eps = 1.0D0
   gradqm = 0.0D0; gradmm = 0.0D0

   write(*,'(A)') "Subroutine ljs_gradients_qmmm"
   call ljs_gradients_qmmm(gradqm, gradmm, pos, nQM)

   ! Checks QM region gradients.
   allocate(grad_res(3,nQM))
   passed = .true.
   grad_res = reshape((/153579.7521485972D0, 356895.44677734373D0, 0.0D0, &
                        0.0D0, 0.0D0, 0.0D0,    &
                        3.1673435152421874D-2, 7.4393611725054196D-2,    & 
                        0.15910057031621411D0, &
                        0.0D0, 0.0D0, 0.0D0,    &
                        0.0D0, 0.0D0, 0.0D0/), shape(grad_res))
   do counter = 1, 3
   do counter2 = 1, nQM
      diff = abs(grad_res(counter,counter2) - gradqm(counter,counter2)) / &
             abs(grad_res(counter,counter2) + gradqm(counter,counter2)) 
      if (diff > 1.0D-13) then
         passed = .false.
         write(*,'(A34,I1,A14,I1)') "ERROR: Wrong gradient for QM atom ",&
                                    counter2, " in direction ", counter
         write(*,'(7x,A9,F15.7,A11,F15.7)') &
                                "Expected ", grad_res(counter,counter2), &
                                " and found ", gradqm(counter,counter2)
      endif
   enddo
   enddo
   deallocate(grad_res)
   if (passed) write(*,'(A)') "=> QM region gradients calculated correctly."

   ! Checks MM region gradients.
   allocate(grad_res(3,nMM))
   passed = .true.
   grad_res = reshape((/-153579.7521485907D0, -4.2720176572632322D-2, &
                        -6.4080264858948480D-2, &
                        -3.1673435152421874D-2, -356895.47845077893D0, &
                        -9.5020305457265614D-2, &
                         0.0D0, 0.0D0, 0.0D0/), shape(grad_res))
   do counter = 1, 3
   do counter2 = 1, nMM
      diff = abs(grad_res(counter,counter2) - gradmm(counter,counter2)) / &
             abs(grad_res(counter,counter2) + gradmm(counter,counter2))
             
      if (diff > 1.0D-13) then
         passed = .false.
         write(*,'(A34,I1,A14,I1)') "ERROR: Wrong gradient for MM atom ",&
                                    counter2, " in direction ", counter
         write(*,'(7x,A9,F15.7,A11,F15.7)') &
                                "Expected ", grad_res(counter,counter2), &
                                " and found ", gradmm(counter,counter2)
         
      endif
   enddo
   enddo
   deallocate(grad_res)
   if (passed) write(*,'(A)') "=> MM region gradients calculated correctly."
   write(*,'(A)') ""

   deallocate(qm_typ, mm_typ, pos, gradqm, gradmm)
   call ljs_finalise()
end subroutine test_mm_interface


! Charges and derivatives [CRG]
! This subroutine tests the proper assignation of LJ parameters
! and their derivatives in accordance to the atomic charge.
subroutine test_crg_and_derivs()
   use lj_switch_data, only: lj_atoms, k_fermi
   use lj_switch     , only: ljs_finalise
   implicit none

   integer :: counter
   logical :: passed
   LIODBLE :: diff, res_out(5)

   k_fermi = 10.0D0
   allocate(lj_atoms(5))
   do counter = 1, 5
      lj_atoms(counter)%s1 = 2.0D0; lj_atoms(counter)%s2 = 3.0D0
      lj_atoms(counter)%e1 = 1.0D0; lj_atoms(counter)%e2 = 0.5D0
      lj_atoms(counter)%q1 = 0.0D0; lj_atoms(counter)%q2 = 1.0D0
   enddo

   write(*,'(A)') "Subroutine set_eps_sig"
   call lj_atoms(1)%set_eps_sig(-0.5D0)
   call lj_atoms(2)%set_eps_sig( 0.0D0)
   call lj_atoms(3)%set_eps_sig( 0.5D0)
   call lj_atoms(4)%set_eps_sig( 1.0D0)
   call lj_atoms(5)%set_eps_sig( 1.5D0)

   ! Verifies sigma and  epsilon
   passed = .true.
   res_out  = (/ 2.00005D0,  2.00669D0,  2.50000D0,  2.99331D0,  2.99995D0/)
   do counter = 1, 5
      diff = abs(lj_atoms(counter)%sig - res_out(counter))
      
      if (diff > 1.0D-5) then
      passed = .false.
      write(*,'(A32,I1)') "ERROR: Wrong sig value for atom ", counter
      write(*,'(7x,A9,F14.7,A11,F14.7)') "Expected ", res_out(counter), &
                                         " and found ", lj_atoms(counter)%sig
      endif
   enddo
   if (passed) write(*,'(A)') "=> Sigma calculated correctly."

   passed = .true.
   res_out  = (/ 0.99998D0,  0.99665D0,  0.75000D0,  0.50335D0,  0.50002D0/)
   do counter = 1, 5
      diff = abs(lj_atoms(counter)%eps - res_out(counter))
      
      if (diff > 1.0D-5) then
      passed = .false.
      write(*,'(A32,I1)') "ERROR: Wrong eps value for atom ", counter
      write(*,'(7x,A9,F14.7,A11,F14.7)') "Expected ", res_out(counter), &
                                         " and found ", lj_atoms(counter)%eps
      endif
   enddo
   if (passed) write(*,'(A)') "=> Epsilon calculated correctly."

   ! Verifies sigma and epsilon derivatives with respect to Q.
   passed = .true.
   res_out  = (/ 0.00045D0,  0.06648D0,  2.50000D0,  0.06648D0,  0.00045D0/)
   do counter = 1, 5
      diff = abs(lj_atoms(counter)%dsig - res_out(counter))
      
      if (diff > 1.0D-5) then
      passed = .false.
      write(*,'(A33,I1)') "ERROR: Wrong dsig value for atom ", counter
      write(*,'(7x,A9,F14.7,A11,F14.7)') "Expected ", res_out(counter), &
                                         " and found ", lj_atoms(counter)%dsig
      endif
   enddo
   if (passed) write(*,'(A)') "=> Sigma derivative calculated correctly."

   passed = .true.
   res_out  = (/-0.00023D0, -0.03324D0, -1.25000D0, -0.03324D0, -0.00023D0/)
   do counter = 1, 5
      diff = abs(lj_atoms(counter)%deps - res_out(counter))
      
      if (diff > 1.0D-5) then
      passed = .false.
      write(*,'(A3,I1)') "ERROR: Wrong deps value for atom ", counter
      write(*,'(7x,A9,F14.7,A11,F14.7)') "Expected ", res_out(counter), &
                                         " and found ", lj_atoms(counter)%deps
      endif
   enddo
   if (passed) write(*,'(A)') "=> Epsilon derivative calculated correctly."
   write(*,'(A)') ""

   call ljs_finalise()
end subroutine test_crg_and_derivs

! Fock terms [FOK]
! This subroutine tests the proper calculation of Fock matrix elements.
! It also implicitly tests the value dE/dQ, which is needed by the
! aforementioned calculation.
subroutine test_fock_terms()
   use lj_switch_data, only: lj_atoms, k_fermi, mmlj_eps, mmlj_sig, mm_atoms, &
                             n_lj_atoms
   use lj_switch     , only: ljs_finalise, ljs_add_fock_terms
   implicit none

   integer :: counter, counter2, n_bas
   logical :: passed
   LIODBLE :: diff, ener, ener_res
   LIODBLE, allocatable :: fock_dum(:,:), ovlap_dum(:,:), rho_dum(:,:), &
                           fock_res(:,:)

   k_fermi = 10.0D0
   n_lj_atoms = 2
   allocate(lj_atoms(2))
   do counter = 1, 2
      lj_atoms(counter)%s1 = 2.0D0; lj_atoms(counter)%s2 = 3.0D0
      lj_atoms(counter)%e1 = 1.0D0; lj_atoms(counter)%e2 = 0.5D0
      lj_atoms(counter)%q1 = 0.0D0; lj_atoms(counter)%q2 = 1.0D0
      lj_atoms(counter)%Z  = 5
   enddo
   allocate(lj_atoms(1)%basis_id(2))
   allocate(lj_atoms(2)%basis_id(3))
   lj_atoms(1)%basis_id = (/1,3/)
   lj_atoms(2)%basis_id = (/2,6,7/)


   n_bas = 10
   allocate(fock_dum(n_bas,n_bas), ovlap_dum(n_bas,n_bas), &
            rho_dum(n_bas,n_bas))
   fock_dum  = 0.0D0
   ovlap_dum = 0.0D0
   rho_dum   = 0.0D0

   do counter = 1, 6
      rho_dum(counter  , counter  ) = 1.00D0
      rho_dum(counter  , counter+1) = 0.25D0
      rho_dum(counter+1, counter  ) = 0.25D0
   enddo
   
   ! Charge of atom 1 should be 0.5, charge of atom 2 should be -0.5.
   do counter = 1, 4
      ovlap_dum(counter,counter) = 1.0D0
   enddo
   do counter = 5, 8
      ovlap_dum(counter,counter) = 0.5D0
   enddo
   ovlap_dum(1,2) = 3.0D0; ovlap_dum(2,1) = 3.0D0
   ovlap_dum(2,3) = 4.0D0; ovlap_dum(3,2) = 4.0D0
   ovlap_dum(3,4) = 3.0D0; ovlap_dum(4,3) = 3.0D0
   ovlap_dum(5,6) = 3.0D0; ovlap_dum(6,5) = 3.0D0
   ovlap_dum(6,7) = 3.0D0; ovlap_dum(7,6) = 3.0D0

   ! MM atoms setup
   allocate(mmlj_sig(3), mmlj_eps(3), mm_atoms(3))
   allocate(mm_atoms(1)%dist(2), mm_atoms(2)%dist(2), mm_atoms(3)%dist(2))
   mmlj_sig = (/2.0D0, 2.0D0, 3.0D0/)
   mmlj_eps = (/0.5D0, 0.5D0, 1.0D0/)
   mm_atoms(1)%mmtype = 1; mm_atoms(2)%mmtype = 2; mm_atoms(3)%mmtype = 3

   mm_atoms(1)%dist(1) = 2.0D0; mm_atoms(1)%dist(2) = 4.0D0
   mm_atoms(2)%dist(1) = 4.0D0; mm_atoms(2)%dist(2) = 2.0D0
   mm_atoms(3)%dist(1) = 3.0D0; mm_atoms(3)%dist(2) = 3.0D0

   ener = 0.0D0
   write(*,'(A)') "Subroutine ljs_add_fock_terms"
   call ljs_add_fock_terms(fock_dum, ener, rho_dum, ovlap_dum)
   
   ! Checks total energy.
   ener_res = 0.81400396217800775D0
   diff = abs(ener - ener_res) / abs(ener + ener_res)
   if (diff > 1.0D-9) then
       passed = .false.
       write(*,'(A19)') "ERROR: Wrong energy"
       write(*,'(7x,A9,F15.7,A11,F15.7)') "Expected ", ener_res, &
                                          " and found ", ener
   else
       passed = .true.
       write(*,'(A)') "=> Energy calculated correctly."
   endif
   
   ! Checks Fock matrix output
   allocate(fock_res(n_bas,n_bas))
   fock_res = 0.0D0
   fock_res(1,1) = -11.968292133056625D0
   fock_res(2,1) = -35.906198098585342D0 ; fock_res(1,2) = fock_res(2,1)
   fock_res(2,2) = -4.4056647182182118D-4
   fock_res(2,3) = -47.874930798113787D0 ; fock_res(3,2) = fock_res(2,3)
   fock_res(3,3) = -11.968292133056625D0
   fock_res(3,4) = -35.904876399169879D0 ; fock_res(4,3) = fock_res(3,4)
   fock_res(5,6) = -1.3216994154654636D-3; fock_res(6,5) = fock_res(5,6)
   fock_res(6,6) = -2.2028323591091059D-4
   fock_res(6,7) = -2.6433988309309273D-3; fock_res(7,6) = fock_res(6,7)
   fock_res(7,7) = -2.2028323591091059D-4

   passed = .true.
   do counter  = 1, n_bas
   do counter2 = 1, n_bas
      diff = abs(fock_res(counter,counter2) - fock_dum(counter,counter2)) / &
             abs(fock_res(counter,counter2) + fock_dum(counter,counter2))
             
      if (diff > 1.0D-13) then
         passed = .false.
         write(*,'(A33,I2,A1,I2)') "ERROR: Wrong Fock matrix element ",&
                                    counter, ",", counter
         write(*,'(7x,A9,F15.7,A11,F15.7)') &
                                "Expected ", fock_res(counter,counter2), &
                                " and found ", fock_dum(counter,counter2)
         
      endif
   enddo
   enddo
   if (passed) write(*,'(A)') "=> Fock matrix elements calculated correctly."
   write(*,'(A)') ""

   call ljs_finalise()
   deallocate(fock_dum, fock_res, ovlap_dum, rho_dum)
end subroutine test_fock_terms