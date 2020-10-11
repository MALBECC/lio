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
      call cdft_get_deltaV_dual(fock_a, rho_a, fock_b, rho_b)
   else
      call cdft_get_deltaV_regular(fock_a, rho_a, fock_b, rho_b)
   endif


end subroutine cdft_get_deltaV

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

      do ii = 1, cdft_c%n_regions
         dV = cdft_reg%cst(ii)
         if (abs(dV) > 0.01D0) dV = 0.01D0
         cdft_reg%Vc(ii) = cdft_reg%Vc(ii) + dV

         call g2g_cdft_set_v(cdft_reg%Vc, cdft_reg%Vs)
         call SCF(energ, fock_a, rho_a, fock_b, rho_b)

         call cdft_get_constraints()
         do jj = 1, cdft_c%n_regions
            cdft_c%jacob(jj,ii) = (cdft_reg%cst(jj) - cdft_reg%cst_old(jj)) &
                                  / dV
         enddo
         if (cdft_c%do_spin) then
            do jj = 1+cdft_c%sp_idx, cdft_c%n_regions+cdft_c%sp_idx
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

         call g2g_cdft_set_V(cdft_reg%Vc, cdft_reg%Vs)
         call SCF(energ, fock_a, rho_a, fock_b, rho_b)
         
         call cdft_get_constraints()
         if (cdft_c%do_chrg) then
            do jj = 1, cdft_c%n_regions
               cdft_c%jacob(jj,ii+cdft_c%sp_idx) = (cdft_reg%cst(jj) - &
                                      cdft_reg%cst_old(jj)) / dV
            enddo
         endif
         
         do jj = 1+cdft_c%sp_idx, cdft_c%n_regions+cdft_c%sp_idx
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
      if (size(cdft_c%jacob,1) == 1) then
         if (cdft_c%do_chrg) then
          cdft_reg%Vc(1) = cdft_reg%Vc_old(1) - &
                             cdft_reg%cst(1) / cdft_c%jacob(1,1)
         else if (cdft_c%do_spin) then
            cdft_reg%Vs(1) = cdft_reg%Vs_old(1) - &
                             cdft_reg%cst(1) / cdft_c%jacob(1,1)
       endif
      else
         if (cdft_c%do_chrg) cdft_reg%Vm_old(1:size(cdft_reg%Vc,1)) = &
                             cdft_reg%Vc_old(:)
         if (cdft_c%do_spin) &
                  cdft_reg%Vm_old((cdft_c%sp_idx+1):(cdft_c%sp_idx+size(cdft_reg%Vs,1))) = &
                  cdft_reg%Vs_old(:)

         cdft_reg%Vmix = cdft_reg%Vmix + cdft_reg%Vm_old

         if (cdft_c%do_chrg) cdft_reg%Vc = cdft_reg%Vmix(1:size(cdft_reg%Vc,1))
         if (cdft_c%do_spin) cdft_reg%Vs = &
                        cdft_reg%Vmix((cdft_c%sp_idx+1):(cdft_c%sp_idx+size(cdft_reg%Vs,1)))
      endif
   endif

   call g2g_cdft_set_V(cdft_reg%Vc, cdft_reg%Vs)
end subroutine cdft_set_potential

! Checks if CDFT converged.
subroutine cdft_check_conver(rho_new, rho_old, converged, cdft_iter, ener, &
                             rho_crit)
   use cdft_data, only: cdft_reg
   implicit none
   LIODBLE, intent(in)    :: rho_new(:), rho_old(:), rho_crit
   integer, intent(in)    :: cdft_iter
   logical, intent(out)   :: converged
   LIODBLE, intent(inout) :: ener

   LIODBLE :: rho_diff, c_max
   integer      :: jj
   
   rho_diff = 0.0D0
   do jj = 1 , size(rho_new,1)
      rho_diff  = rho_diff + (rho_new(jj) - rho_old(jj)) * &
                             (rho_new(jj) - rho_old(jj))
   enddo
   rho_diff = sqrt(rho_diff) / dble(size(rho_new,1))

   call cdft_get_constraints()
   c_max = maxval(abs(cdft_reg%cst))

   call cdft_add_energy(ener)
   write(*,'(A)') "CDFT Convergence status:" 
   write(*,*) "Iteration n°:      ", cdft_iter
   write(*,*) "Energy:            ", ener
   write(*,*) "ΔRho:              ", rho_diff
   write(*,*) "Constraint values: ", cdft_reg%cst
   write(*,*) "Charge potential:  ", cdft_reg%Vc
   write(*,*) "Spin potential:    ", cdft_reg%Vs
   converged = .false.
   if ((rho_diff < rho_crit) .and. (c_max < 1D-5)) converged = .true.
end subroutine cdft_check_conver

! Adds CDFT terms to total energy.
subroutine cdft_add_energy(energ)
   use cdft_data, only: cdft_c, cdft_reg
   implicit none
   LIODBLE, intent(inout) :: energ
   integer :: ii

   call cdft_get_constraints()
   if (cdft_c%dual) then
      if (cdft_c%do_chrg) energ = energ + abs(cdft_reg%Vc(1)) * &
                                  (cdft_reg%chrg(2) - cdft_reg%chrg(1))
      if (cdft_c%do_chrg) energ = energ + abs(cdft_reg%Vs(1)) * &
                                  (cdft_reg%spin(2) - cdft_reg%spin(1))
   else
      if (cdft_c%do_chrg) then
         do ii = 1, cdft_c%n_regions
            energ = energ + cdft_reg%Vc(ii) * &
                            (dble(cdft_reg%nelecs(ii)) - cdft_reg%chrg(ii))
         enddo
      end if

      if (cdft_c%do_spin) then
         do ii = 1, cdft_c%n_regions
            energ = energ + cdft_reg%spin(ii) * cdft_reg%Vs(ii)
         enddo
      endif
   endif
end subroutine cdft_add_energy

! Gets Becke atomic spin and calculates its difference
! with the constrained total value.
subroutine cdft_get_constraints()
   use cdft_data, only: cdft_reg, cdft_c
   implicit none
   integer :: c_index, region

   cdft_reg%cst = 0.0D0
   if (cdft_c%do_chrg) then
      call g2g_get_becke_dens(cdft_c%at_chrg)

      do region = 1, cdft_c%n_regions
         do c_index = 1, cdft_reg%natom(region)
            cdft_reg%cst(region) = cdft_reg%cst(region) + &
                                   cdft_c%at_chrg(cdft_reg%atoms(region,c_index))
         enddo
         cdft_reg%cst(region) = cdft_reg%chrg(region) - cdft_reg%cst(region)
      enddo
   endif

   if (cdft_c%do_spin) then
      call g2g_get_becke_spin(cdft_c%at_spin)

      do region = 1, cdft_c%n_regions
         do c_index = 1, cdft_reg%natom(region)
            cdft_reg%cst(region+cdft_c%sp_idx)= cdft_reg%cst(region+cdft_c%sp_idx) + &
                                         cdft_c%at_spin(cdft_reg%atoms(region,c_index))
         enddo
         cdft_reg%cst(region+cdft_c%sp_idx) = cdft_reg%cst(region+cdft_c%sp_idx) - &
                                       cdft_reg%spin(region)
                                       
      enddo
   endif
end subroutine cdft_get_constraints

! Reorders data to better treat the cases of only 1 region or 2 regions
! that encompass the whole system.
subroutine cdft_rearrange_regions(n_atoms, charge, nunp)
   use cdft_data, only: cdft_c, cdft_reg
   implicit none
   integer, intent(in) :: charge, n_atoms, nunp

   LIODBLE, allocatable :: chrg(:), spin(:), &
                           chrg2(:), spin2(:)
   integer, allocatable :: natom(:), atoms(:,:)
   integer :: ii, jj

   if (cdft_c%n_regions > 2) return
   if (cdft_c%n_regions == 1) then
      allocate(chrg(2), spin(2), chrg2(2), spin2(2))
      allocate(natom(2),atoms(2,n_atoms))
      atoms = 0

      chrg(1)  = cdft_reg%chrg(1)
      spin(1)  = cdft_reg%spin(1)
      natom(1) = cdft_reg%natom(1) 
      atoms(1,1:natom(1)) = cdft_reg%atoms(1,1:natom(1))
      chrg(2)  = dble(charge) - cdft_reg%chrg(1)
      spin(2)  = dble(nunp)   - cdft_reg%spin(1)

      if (cdft_c%mixed) then
         chrg2(1) = cdft_reg%chrg2(1)
         spin2(1) = cdft_reg%spin2(1)
         chrg2(2) = dble(charge) - cdft_reg%chrg2(1)
         spin2(2) = dble(nunp)   - cdft_reg%spin2(1)
      endif
      natom(2) = n_atoms - cdft_reg%natom(1) 
      
      ! Adds the atoms that are not in region 1 to region 2.
      jj = 1
      do ii = 1, n_atoms
         if (.not. (any(atoms(1,:) == ii)) ) then
            atoms(2,jj) = ii
            jj = jj +1
         endif
      enddo

      ! Reallocates everything
      cdft_c%dual      = .true.
      cdft_c%n_regions = 2
      cdft_c%max_nat   = maxval(natom,1) 
      deallocate(cdft_reg%chrg, cdft_reg%spin)
      deallocate(cdft_reg%natom, cdft_reg%atoms)
      allocate(cdft_reg%chrg(2), cdft_reg%spin(2))
      allocate(cdft_reg%natom(2), cdft_reg%atoms(2,cdft_c%max_nat))
      cdft_reg%chrg  = chrg
      cdft_reg%spin  = spin
      cdft_reg%natom = natom
      cdft_reg%atoms = atoms(:,1:cdft_c%max_nat)

      if (cdft_c%mixed) then
         deallocate(cdft_reg%chrg2, cdft_reg%spin2)
         allocate(cdft_reg%chrg2(2), cdft_reg%spin2(2))
         cdft_reg%chrg2  = chrg2
         cdft_reg%spin2  = spin2 
      endif
      deallocate(chrg, spin, chrg2, spin2, natom, atoms)
   elseif ( n_atoms == (cdft_reg%natom(1) + cdft_reg%natom(2)) ) then
      cdft_c%dual = .true.
   endif
end subroutine cdft_rearrange_regions