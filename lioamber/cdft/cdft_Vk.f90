! Gets the jacobian matrix by making a small perturbation in each direction.
! This is done in order to propagate the constraint potentials Vc and Vs
! by means of Newton's method. Instead of calculating J-1, we solve an
! alternative problem resulting in ΔVi (i=c,s). 
subroutine cdft_get_deltaV(fock_a, rho_a, fock_b, rho_b)
   use typedef_operator, only: operator
   use cdft_data       , only: cdft_c

   implicit none
   type(operator), intent(inout) :: fock_a, rho_a, fock_b, rho_b

   if (cdft_c%dual) then
      ! We consider the special case of only two regions.
      call cdft_get_deltaV_dual(fock_a, rho_a, fock_b, rho_b)
   else
      call cdft_get_deltaV_regular(fock_a, rho_a, fock_b, rho_b)
   endif
end subroutine cdft_get_deltaV

! We will always have n regions encomprising all of the atoms.
! As such, one of the Jacobian terms will be linearly dependent
! on the others, since Sum(Vk) = 0 for a given type of constraint
! (charge or spin).
subroutine cdft_get_deltaV_regular(fock_a, rho_a, fock_b, rho_b)
   use typedef_operator, only: operator
   use cdft_data       , only: cdft_c, cdft_reg
   
   implicit none
   type(operator), intent(inout) :: fock_a, rho_a, fock_b, rho_b

   integer :: ii, jj
   LIODBLE :: energ, dV

   ! Variables for LAPACK
   integer              :: LWORK, INFO
   LIODBLE, allocatable :: WORK(:)

   cdft_c%jacob = 0.0D0
   call cdft_get_constraints()
   cdft_reg%cst_old = cdft_reg%cst

   if (cdft_c%do_chrg) then
      cdft_reg%Vc_old = cdft_reg%Vc

      do ii = 1, cdft_c%n_regions -1
         dV = cdft_reg%cst(ii)
         if (abs(dV) > 0.01D0) dV = 0.01D0
         cdft_reg%Vc(ii) = cdft_reg%Vc(ii) + dV
         cdft_reg%Vc(cdft_c%n_regions) = 0.0D0
         do jj = 1, cdft_c%n_regions -1
            cdft_reg%Vc(cdft_c%n_regions) = &
               cdft_reg%Vc(cdft_c%n_regions) - cdft_reg%Vc(jj)
         enddo

         call g2g_cdft_set_v(cdft_reg%Vc, cdft_reg%Vs)
         call SCF(energ, fock_a, rho_a, fock_b, rho_b)

         call cdft_get_constraints()
         do jj = 1, cdft_c%n_regions -1
            cdft_c%jacob(jj,ii) = (cdft_reg%cst(jj) - cdft_reg%cst_old(jj)) &
                                  / dV
         enddo
         if (cdft_c%do_spin) then
            do jj = cdft_c%sp_idx +1, cdft_c%sp_idx + cdft_c%n_regions -1
               cdft_c%jacob(jj,ii) = (cdft_reg%cst(jj) - cdft_reg%cst_old(jj)) &
                                     / dV
            enddo
         endif
         cdft_reg%Vc = cdft_reg%Vc_old    
      enddo
   endif

   if (cdft_c%do_spin) then
      cdft_reg%Vs_old = cdft_reg%Vs

      do ii = 1, cdft_c%n_regions
         dV = cdft_reg%cst(ii+cdft_c%sp_idx)
         if (abs(dV) > 0.01D0) dV = 0.01D0
         cdft_reg%Vs(ii) = cdft_reg%Vs(ii) + dV
         cdft_reg%Vs(cdft_c%n_regions) = 0.0D0
         do jj = 1, cdft_c%n_regions -1
            cdft_reg%Vs(cdft_c%n_regions) = &
               cdft_reg%Vs(cdft_c%n_regions) - cdft_reg%Vs(jj)
         enddo

         call g2g_cdft_set_V(cdft_reg%Vc, cdft_reg%Vs)
         call SCF(energ, fock_a, rho_a, fock_b, rho_b)
         
         call cdft_get_constraints()
         if (cdft_c%do_chrg) then
            do jj = 1, cdft_c%n_regions -1
               cdft_c%jacob(jj,ii+cdft_c%sp_idx) = (cdft_reg%cst(jj) - &
                                                   cdft_reg%cst_old(jj)) / dV
            enddo
         endif
         
         do jj = cdft_c%sp_idx +1, cdft_c%sp_idx + cdft_c%n_regions -1
            cdft_c%jacob(jj,ii+cdft_c%sp_idx) = (cdft_reg%cst(jj) - &
                                   cdft_reg%cst_old(jj)) / dV
         enddo
         cdft_reg%Vs = cdft_reg%Vs_old
      enddo
   endif
   cdft_reg%cst = cdft_reg%cst_old

   ! Alternative to invert J-1: Solving A*x = B where A is the jacobian
   ! and B is the negative constraints array. The result is an array containing
   ! Xn+1 - Xn for each constraint. This is totally equivalent to solve
   ! Xn+1 = Xn - J^(-1) * Cst, with A = J and B = -Cst, but much less costly.
   if (size(cdft_c%jacob,1) > 1) then
      cdft_reg%Vmix = -cdft_reg%cst 
      allocate(WORK(1))
      call dgels('N', size(cdft_c%jacob,1), size(cdft_c%jacob,1), 1, &
                 cdft_c%jacob, size(cdft_c%jacob,1), cdft_reg%Vmix,  &
                 size(cdft_c%jacob,1), WORK, -1, INFO)
      LWORK = int(WORK(1))
      deallocate(WORK)
      allocate(WORK(LWORK))
      call dgels('N', size(cdft_c%jacob,1), size(cdft_c%jacob,1), 1, &
                 cdft_c%jacob, size(cdft_c%jacob,1), cdft_reg%Vmix,  &
                 size(cdft_c%jacob,1), WORK, LWORK, INFO)
      deallocate(WORK)
   endif
end subroutine cdft_get_deltaV_regular

subroutine cdft_get_deltaV_dual(fock_a, rho_a, fock_b, rho_b)
   use typedef_operator, only: operator
   use cdft_data       , only: cdft_c, cdft_reg
   
   implicit none
   type(operator), intent(inout) :: fock_a, rho_a, fock_b, rho_b
   LIODBLE :: energ, dV

   cdft_c%jacob = 0.0D0
   ! Even though we have a 2x2 jacobian, we only need one element
   ! since there is only one (symmetrical) field, and thus only
   ! one real constraint (the other one has linear dependency
   ! with the first).
   call cdft_get_constraints()
   cdft_reg%cst_old = cdft_reg%cst

   if (cdft_c%do_chrg) then
      cdft_reg%Vc_old = cdft_reg%Vc

      dV = cdft_reg%cst(1)
      if (abs(dV) > 0.01D0) dV = 0.01D0
      cdft_reg%Vc(1) = cdft_reg%Vc(1) + dV
      cdft_reg%Vc(2) = cdft_reg%Vc(2) - dV

      call g2g_cdft_set_v(cdft_reg%Vc, cdft_reg%Vs)
      call SCF(energ, fock_a, rho_a, fock_b, rho_b)

      call cdft_get_constraints()
      cdft_c%jacob(1,1) = (cdft_reg%cst(1) - cdft_reg%cst_old(1))  / dV
      if (cdft_c%do_spin) &
         cdft_c%jacob(3,1) = (cdft_reg%cst(3) - cdft_reg%cst_old(3)) / dV
      cdft_reg%Vc = cdft_reg%Vc_old
   endif

   if (cdft_c%do_spin) then
      cdft_reg%Vs_old = cdft_reg%Vs

      dV = cdft_reg%cst(3)
      if (abs(dV) > 0.01D0) dV = 0.01D0
      cdft_reg%Vs(1) = cdft_reg%Vs(1) + dV
      cdft_reg%Vs(2) = cdft_reg%Vs(2) - dV

      call g2g_cdft_set_V(cdft_reg%Vc, cdft_reg%Vs)
      call SCF(energ, fock_a, rho_a, fock_b, rho_b)
      
      call cdft_get_constraints()
      if (cdft_c%do_chrg) &
         cdft_c%jacob(1,1+cdft_c%sp_idx) = (cdft_reg%cst(1) - &
                                            cdft_reg%cst_old(1)) / dV
       
      cdft_c%jacob(1+cdft_c%sp_idx,1+cdft_c%sp_idx) = &
            (cdft_reg%cst(1+cdft_c%sp_idx) - &
             cdft_reg%cst_old(1+cdft_c%sp_idx)) / dV
      cdft_reg%Vs = cdft_reg%Vs_old
   endif
   cdft_reg%cst = cdft_reg%cst_old

end subroutine cdft_get_deltaV_dual

! Propagates the constraint potentials by means of Newton's method. Vmix,
! which contains ΔVc and ΔVs, is obtained in the previous routine.
subroutine cdft_set_potential()
   use cdft_data, only: cdft_c, cdft_reg
   implicit none
   LIODBLE :: Vdiff
   integer :: ii

   if (cdft_c%dual) then
      if (cdft_c%do_chrg) then
         Vdiff = cdft_reg%cst(1) / cdft_c%jacob(1,1)
         cdft_reg%Vc(1) = cdft_reg%Vc_old(1) - Vdiff
         cdft_reg%Vc(2) = cdft_reg%Vc_old(2) + Vdiff
      endif
      if (cdft_c%do_spin) then
         Vdiff = cdft_reg%cst(1+cdft_c%sp_idx) &
                 / cdft_c%jacob(1+cdft_c%sp_idx,1+cdft_c%sp_idx)
         cdft_reg%Vs(1) = cdft_reg%Vs_old(1) - Vdiff
         cdft_reg%Vs(2) = cdft_reg%Vs_old(2) + Vdiff
      endif
   else
      if (cdft_c%do_chrg) cdft_reg%Vm_old(1:cdft_c%n_regions-1) = &
                          cdft_reg%Vc_old(1:cdft_c%n_regions-1)
      if (cdft_c%do_spin) &
               cdft_reg%Vm_old((cdft_c%sp_idx+1):(cdft_c%sp_idx+cdft_c%n_regions-1)) = &
               cdft_reg%Vs_old(1:cdft_c%n_regions-1)

      cdft_reg%Vmix = cdft_reg%Vmix + cdft_reg%Vm_old

      if (cdft_c%do_chrg) then
         cdft_reg%Vc(1:cdft_c%n_regions-1) = cdft_reg%Vmix(1:cdft_c%n_regions-1)
         cdft_reg%Vc(cdft_c%n_regions) = 0.0D0
         do ii = 1, cdft_c%n_regions -1
            cdft_reg%Vc(cdft_c%n_regions) = cdft_reg%Vc(cdft_c%n_regions) - &
                                            cdft_reg%Vc(ii) 
         enddo
      endif
      if (cdft_c%do_spin) then
         cdft_reg%Vs(1:cdft_c%n_regions-1) = cdft_reg%Vmix((cdft_c%sp_idx+1):&
                                                (cdft_c%sp_idx+cdft_c%n_regions-1))    
         cdft_reg%Vs(cdft_c%n_regions) = 0.0D0
         do ii = 1, cdft_c%n_regions -1
            cdft_reg%Vs(cdft_c%n_regions) = cdft_reg%Vs(cdft_c%n_regions) - &
                                             cdft_reg%Vs(ii) 
         enddo
      endif
   endif

   call g2g_cdft_set_V(cdft_reg%Vc, cdft_reg%Vs)
end subroutine cdft_set_potential