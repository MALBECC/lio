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
   if ((rho_diff < rho_crit) .and. (c_max < 1D-4)) converged = .true.
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

      do region = 1, cdft_c%n_regions -1
         do c_index = 1, cdft_reg%natom(region)
            cdft_reg%cst(region) = cdft_reg%cst(region) + &
                                   cdft_c%at_chrg(cdft_reg%atoms(region,c_index))
         enddo
         cdft_reg%cst(region) = cdft_reg%chrg(region) - cdft_reg%cst(region)
      enddo
   endif

   if (cdft_c%do_spin) then
      call g2g_get_becke_spin(cdft_c%at_spin)

      do region = 1, cdft_c%n_regions -1
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
   integer :: ii, jj, n_old, t_atoms

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
      return
   elseif ( n_atoms == (cdft_reg%natom(1) + cdft_reg%natom(2)) ) then
      cdft_c%dual = .true.
      return
   endif

   ! Checks if all atoms are included in region.
   t_atoms = 0
   do ii = 1, cdft_c%n_regions
      t_atoms = t_atoms + cdft_reg%natom(ii)
   enddo

   if (.not. (t_atoms == n_atoms)) then
      n_old = cdft_c%n_regions
      

      allocate(chrg(n_old +1) , spin(n_old +1) )
      allocate(chrg2(n_old +1), spin2(n_old +1))
      allocate(natom(n_old +1), atoms(n_old +1,n_atoms))
      atoms = 0

      chrg(1:n_old)  = cdft_reg%chrg
      spin(1:n_old)  = cdft_reg%spin
      natom(1:n_old) = cdft_reg%natom
      atoms(1:n_old,1:natom(1)) = cdft_reg%atoms(1:n_old,1:natom(1))
      chrg(n_old+1)  = dble(charge)
      spin(n_old+1)  = dble(nunp)
      do ii = 1, n_old
         chrg(n_old+1)  = - cdft_reg%chrg(ii)
         spin(n_old+1)  = - cdft_reg%spin(ii)
      enddo

      if (cdft_c%mixed) then
         chrg2(1:n_old) = cdft_reg%chrg2(1:n_old)
         spin2(1:n_old) = cdft_reg%spin2(1:n_old)
         chrg2(n_old+1) = dble(charge)
         spin2(n_old+1) = dble(nunp)
         do ii = 1, n_old
            chrg2(n_old+1)  = - cdft_reg%chrg2(ii)
            spin2(n_old+1)  = - cdft_reg%spin2(ii)
         enddo
      endif

      natom(n_old+1) = n_atoms - t_atoms 
      
      ! Adds the atoms that are not in region 1 to region 2.
      jj = 1
      do ii = 1, n_atoms
         if (.not. (any(cdft_reg%atoms == ii)) ) then
            atoms(2,jj) = ii
            jj = jj +1
         endif
      enddo

      ! Reallocates everything
      cdft_c%max_nat   = maxval(natom,1) 
      cdft_c%n_regions = n_old +1
      deallocate(cdft_reg%chrg, cdft_reg%spin)
      deallocate(cdft_reg%natom, cdft_reg%atoms)
      allocate(cdft_reg%chrg(cdft_c%n_regions))
      allocate(cdft_reg%spin(cdft_c%n_regions))
      allocate(cdft_reg%natom(cdft_c%n_regions))
      allocate(cdft_reg%atoms(cdft_c%n_regions,cdft_c%max_nat))
      cdft_reg%chrg  = chrg
      cdft_reg%spin  = spin
      cdft_reg%natom = natom
      cdft_reg%atoms = atoms(:,1:cdft_c%max_nat)

      if (cdft_c%mixed) then
         deallocate(cdft_reg%chrg2, cdft_reg%spin2)
         allocate(cdft_reg%chrg2(cdft_c%n_regions))
         allocate(cdft_reg%spin2(cdft_c%n_regions))
         cdft_reg%chrg2  = chrg2
         cdft_reg%spin2  = spin2 
      endif
      deallocate(chrg, spin, chrg2, spin2, natom, atoms)
   endif
end subroutine cdft_rearrange_regions