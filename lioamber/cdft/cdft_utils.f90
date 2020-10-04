! Gets the jacobian matrix by making a small perturbation in each direction.
! This is done in order to propagate the constraint potentials Vc and Vs
! by means of Newton's method. Instead of calculating J-1, we solve an
! alternative problem resulting in ΔVi (i=c,s). 
subroutine cdft_get_deltaV(fock_a, rho_a, fock_b, rho_b)
   use typedef_operator, only: operator
   use cdft_data       , only: cdft_chrg, cdft_spin, cdft_reg, &
                               jacob, sp_idx
   
   implicit none
   type(operator), intent(inout) :: fock_a, rho_a, fock_b, rho_b
   integer      :: ii, jj
   LIODBLE :: energ, dV

   ! Variables for LAPACK
   integer                   :: LWORK, INFO
   LIODBLE, allocatable :: WORK(:)

   jacob = 0.0D0
   call cdft_get_constraints()
   cdft_reg%cst_old = cdft_reg%cst

   if (cdft_chrg) then
      cdft_reg%Vc_old = cdft_reg%Vc

      do ii = 1, cdft_reg%n_regions
         dV = cdft_reg%cst(ii)
         if (dV > 0.01D0) dV = 0.01D0
         cdft_reg%Vc(ii) = cdft_reg%Vc(ii) + dV

         call g2g_cdft_set_v(cdft_reg%Vc, cdft_reg%Vs)
         call SCF(energ, fock_a, rho_a, fock_b, rho_b)

         call cdft_get_constraints()
         do jj = 1, cdft_reg%n_regions
            jacob(jj,ii) = (cdft_reg%cst(jj) - cdft_reg%cst_old(jj)) / dV
         enddo
         if (cdft_spin) then
            do jj = 1+sp_idx, cdft_reg%n_regions+sp_idx
               jacob(jj,ii) = (cdft_reg%cst(jj) - cdft_reg%cst_old(jj)) / dV
            enddo
         endif
         cdft_reg%Vc = cdft_reg%Vc_old
      enddo
   endif

   if (cdft_spin) then
      cdft_reg%Vs_old = cdft_reg%Vs

      do ii = 1, cdft_reg%n_regions
         dV = cdft_reg%cst(ii+sp_idx)
         if (dV > 0.01D0) dV = 0.01D0
         cdft_reg%Vs(ii) = cdft_reg%Vs(ii) + dV

         call g2g_cdft_set_V(cdft_reg%Vc, cdft_reg%Vs)
         call SCF(energ, fock_a, rho_a, fock_b, rho_b)
         
         call cdft_get_constraints()
         if (cdft_chrg) then
            do jj = 1, cdft_reg%n_regions
               jacob(jj,ii+sp_idx) = (cdft_reg%cst(jj) - &
                                      cdft_reg%cst_old(jj)) / dV
            enddo
         endif
         
         do jj = 1+sp_idx, cdft_reg%n_regions+sp_idx
            jacob(jj,ii+sp_idx) = (cdft_reg%cst(jj) - &
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
   if (size(jacob,1) > 1) then
      cdft_reg%Vmix = -cdft_reg%cst 
      allocate(WORK(1))
      call dgels('N', size(jacob,1), size(jacob,1), 1, jacob, size(jacob,1), &
                 cdft_reg%Vmix, size(jacob,1), WORK, -1, INFO)
      LWORK = int(WORK(1))
      deallocate(WORK)
      allocate(WORK(LWORK))
      call dgels('N', size(jacob,1), size(jacob,1), 1, jacob, size(jacob,1), &
                 cdft_reg%Vmix, size(jacob,1), WORK, LWORK, INFO)
      deallocate(WORK)
   endif
end subroutine cdft_get_deltaV

! Propagates the constraint potentials by means of Newton's method. Vmix,
! which contains ΔVc and ΔVs, is obtained in the previous routine.
subroutine cdft_set_potential()
   use cdft_data, only: cdft_chrg, cdft_spin, jacob, cdft_reg, sp_idx
   implicit none

   if (size(jacob,1) == 1) then
      if (cdft_chrg) then
         cdft_reg%Vc(1) = cdft_reg%Vc_old(1) - cdft_reg%cst(1) / jacob(1,1)
      else if (cdft_spin) then
         cdft_reg%Vs(1) = cdft_reg%Vs_old(1) - cdft_reg%cst(1) / jacob(1,1)
      endif
   else
      if (cdft_chrg) cdft_reg%Vm_old(1:size(cdft_reg%Vc,1)) = &
                     cdft_reg%Vc_old(:)
      if (cdft_spin) &
               cdft_reg%Vm_old((sp_idx+1):(sp_idx+size(cdft_reg%Vs,1))) = &
               cdft_reg%Vs_old(:)

      cdft_reg%Vmix = cdft_reg%Vmix + cdft_reg%Vm_old

      if (cdft_chrg) cdft_reg%Vc = cdft_reg%Vmix(1:size(cdft_reg%Vc,1))
      if (cdft_spin) cdft_reg%Vs = &
                     cdft_reg%Vmix((sp_idx+1):(sp_idx+size(cdft_reg%Vs,1)))
   endif

   call g2g_cdft_set_V(cdft_reg%Vc, cdft_reg%Vs)
end subroutine cdft_set_potential

! Checks if CDFT converged.
subroutine cdft_check_conver(rho_new, rho_old, converged, cdft_iter, ener, &
                             rho_crit)
   use cdft_data, only: cdft_reg
   implicit none
   LIODBLE, intent(in)    :: rho_new(:), rho_old(:), rho_crit
   integer     , intent(in)    :: cdft_iter
   logical     , intent(out)   :: converged
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
   use cdft_data, only: cdft_chrg, cdft_spin, cdft_reg, sp_idx
   implicit none
   LIODBLE, intent(inout) :: energ
   integer :: ii

   call cdft_get_constraints()
   if (cdft_chrg) then
      do ii = 1, cdft_reg%n_regions
         energ = energ + cdft_reg%cst(ii) * cdft_reg%Vc(ii)
      enddo
   endif

   if (cdft_spin) then
      do ii = 1, cdft_reg%n_regions
         energ = energ + cdft_reg%cst(ii+sp_idx) * cdft_reg%Vs(ii)
      enddo
   endif
end subroutine cdft_add_energy

! Gets Becke atomic spin and calculates its difference
! with the constrained total value.
subroutine cdft_get_constraints()
   use cdft_data, only: cdft_reg, at_spin, at_chrg, &
                        cdft_spin, cdft_chrg, sp_idx
   implicit none
   integer :: c_index, region

   cdft_reg%cst = 0.0D0
   if (cdft_chrg) then
      call g2g_get_becke_dens(at_chrg)

      do region = 1, cdft_reg%n_regions
         do c_index = 1, cdft_reg%natom(region)
            cdft_reg%cst(region) = cdft_reg%cst(region) + &
                                   at_chrg(cdft_reg%atoms(region,c_index))
         enddo
         cdft_reg%cst(region) = cdft_reg%chrg(region) - cdft_reg%cst(region)
      enddo
   endif

   if (cdft_spin) then
      call g2g_get_becke_spin(at_spin)

      do region = 1, cdft_reg%n_regions
         do c_index = 1, cdft_reg%natom(region)
            cdft_reg%cst(region+sp_idx)= cdft_reg%cst(region+sp_idx) + &
                                         at_spin(cdft_reg%atoms(region,c_index))
         enddo
         cdft_reg%cst(region+sp_idx) = cdft_reg%cst(region+sp_idx) - &
                                       cdft_reg%spin(region)
                                       
      enddo
   endif
end subroutine cdft_get_constraints