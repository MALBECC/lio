
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine converger_options_check(energ_all_iter)
   use converger_data, only: damping_factor, gOld, DIIS, hybrid_converg, nDIIS,&
                             conver_criter, rho_LS
   
   logical, intent(inout) :: energ_all_iter

   if (abs(gOld - 10.0D0) > 1e-12) then
      damping_factor = gOld
   else if (abs(damping_factor - 10.0D0) > 1e-12) then
      gOld = damping_factor
   endif

   if (conver_criter == 2) then
      if (hybrid_converg) then
         diis          = .true.
         conver_criter = 3
      else if (diis) then
         conver_criter = 2
      else
         conver_criter = 1
      endif
   else if (conver_criter /= 1) then
      diis = .true.
   endif

   if ((Rho_LS > 1) .and. (conver_criter /=1)) then
      hybrid_converg = .false.
      DIIS           = .false.
      conver_criter  = 1
      write(*,'(A)') &
      '  WARNING - Turning to damping-only convergence for linear search.'
    end if
   
   ! For biased DIIS, full energy is needed.
   if ((conver_criter == 4) .or. (conver_criter == 5)) energ_all_iter = .true.
end subroutine converger_options_check

subroutine converger_init( M_in, OPshell )
   use converger_data, only: fock_damped, told, etold
   implicit none
   integer         , intent(in) :: M_in
   logical         , intent(in) :: OPshell

   write(*,'(2x,A)') "Convergence criteria are: "
   write(*,'(2x,ES8.2,A33,ES8.2,A26)')                 &
            told, " in Rho mean squared difference, ", &
            Etold, " Eh in energy differences."

   if(OPshell) then
      if (.not. allocated(fock_damped) ) allocate(fock_damped(M_in, M_in, 2))
   else
      if (.not. allocated(fock_damped) ) allocate(fock_damped(M_in, M_in, 1))
   end if
   fock_damped(:,:,:) = 0.0D0

   call diis_init(M_in, OPshell)
   call ediis_init(M_in, OPshell)
end subroutine converger_init

subroutine converger_finalise()
   use converger_data, only: fockm, FP_PFm, conver_criter, fock_damped, &
                             bcoef, ndiis, EMAT2, energy_list
   implicit none
   
   if (allocated(fock_damped)) deallocate(fock_damped)

   call diis_finalise()
   call ediis_finalise()
   call rho_ls_finalise()
end subroutine converger_finalise

subroutine conver_fock(niter, M_in, dens_op, fock_op, spin, energy, n_orbs, &
#ifdef CUBLAS
                      devPtrX, devPtrY)
#else
                      Xmat, Ymat)
#endif
   use converger_data  , only: damping_factor, fockm, FP_PFm, ndiis,     &
                               fock_damped, bcoef, EMAT2, conver_criter, &
                               good_cut, level_shift, lvl_shift_en,      &
                               lvl_shift_cut, rho_ls, rho_diff, nediis,  &
                               EDIIS_start, bDIIS_start, DIIS_start
   use typedef_operator, only: operator
   use fileio_data     , only: verbose
   
   implicit none
   ! Spin allows to store correctly alpha or beta information. - Carlos
   integer       , intent(in)    :: niter, M_in, spin, n_orbs   
#ifdef  CUBLAS
   integer*8     , intent(in)    :: devPtrX, devPtrY
#else
   real(kind=8)  , intent(in)    :: Xmat(M_in,M_in), Ymat(M_in,M_in)
#endif
   real(kind=8)  , intent(in)    :: energy
   type(operator), intent(inout) :: dens_op, fock_op

   logical :: diis_on, ediis_on, bdiis_on
   integer :: ndiist, nediist, ii
   real(kind=8), allocatable :: fock00(:,:), EMAT(:,:), suma(:,:), &
                                fock(:,:), rho(:,:), BMAT(:,:)


! INITIALIZATION
! If DIIS is turned on, update fockm with the current transformed F' (into ON
! basis) and update FP_PFm with the current transformed [F',P']
!
! (1)     Calculate F' and [F',P']
!       update fockm with F'
! now, scratch1 = A = F' * P'; scratch2 = A^T
! [F',P'] = A - A^T
! BASE CHANGE HAPPENS INSIDE OF FOCK_COMMUTS

   allocate(fock00(M_in,M_in), fock(M_in,M_in), rho(M_in,M_in))
   fock00 = 0.0D0
   fock   = 0.0D0
   rho    = 0.0D0

   ! Saving rho and the first fock AO
   call dens_op%Gets_data_AO(rho)
   call fock_op%Gets_data_AO(fock00)

   diis_on  = .false.
   ediis_on = .false.
   bdiis_on = .false.
   ! Checks convergence criteria.
   select case (conver_criter)
      ! Always do damping
      case (1)

      ! Damping the first two steps, diis afterwards
      case (2)
         if (niter > 2) diis_on = .true.
      case (4)
         bdiis_on = .true.
         if (niter > 2) diis_on  = .true.

      ! Damping until good enough, diis afterwards
      case(3)
         if ((rho_diff < good_cut) .and. (niter > 2) .and. (rho_LS < 2)) &
            diis_on = .true.
      case(5)
         bdiis_on = .true.
         if ((rho_diff < good_cut) .and. (niter > 2) .and. (rho_LS < 2)) &
            diis_on = .true.
      ! Mix of all criteria, according to rho_diff value.
      case(6:)
         if (rho_diff  < EDIIS_start) ediis_on = .true.
         if ((rho_diff < bDIIS_start) .or. (rho_diff < DIIS_start)) then
            ediis_on = .false.
            diis_on  = .true.
         endif
         if ((rho_diff < bDIIS_start) .and. (rho_diff > DIIS_start)) &
            bdiis_on = .true.
      case default
         write(*,'(A,I4)') 'ERROR - Wrong conver_criter = ', conver_criter
         stop
   endselect

   ! Turn off diis is calculation when change to lineal search, Nick
   if (rho_ls > 1) then 
      diis_on  = .false.
      ediis_on = .false.
   endif

   ndiist  = min(niter, ndiis )
   nediist = min(niter, nediis)
   if (conver_criter /= 1) then
#ifdef CUBLAS
      call dens_op%BChange_AOtoON(devPtrY, M_in, 'r')
      call fock_op%BChange_AOtoON(devPtrX, M_in, 'r')
#else
      call dens_op%BChange_AOtoON(Ymat, M_in, 'r')
      call fock_op%BChange_AOtoON(Xmat, M_in, 'r')
#endif
      call diis_fock_commut(dens_op, fock_op, rho, M_in, spin, ndiist)
   endif

   ! THIS IS DAMPING 
   ! THIS IS SPARTA!
   ! If we are not doing diis this iteration, apply damping to F, save this
   ! F in fock_damped for next iteration's damping and put F' = X^T * F * X in
   ! fock the newly constructed damped matrix is stored, for next iteration in
   ! fock_damped
   if ((.not. diis_on) .and. (.not. ediis_on)) then
      fock = fock00

      if (niter > 1) &
         fock = (fock  + damping_factor * fock_damped(:,:,spin)) / &
                (1.0D0 + damping_factor)
      fock_damped(:,:,spin) = fock
      call fock_op%Sets_data_AO(fock)

#ifdef CUBLAS
      call fock_op%BChange_AOtoON(devPtrX,M_in,'r')
#else
      call fock_op%BChange_AOtoON(Xmat,M_in,'r')
#endif
   endif

   ! DIIS and EDIIS
   if (conver_criter > 1) then

      ! DIIS
      allocate(EMAT(ndiist+1,ndiist+1))
      call diis_update_emat(EMAT, niter, ndiist, spin, M_in)
      ! Stores energy for bDIIS
      if ((conver_criter > 3) .and. (niter > 1)) call diis_update_energy(energy, spin)
   
      if (diis_on) then
         if (bdiis_on) call diis_emat_bias(EMAT, ndiist)

         call diis_get_new_fock(fock, EMAT, ndiist, M_in, spin)
         call fock_op%Sets_data_ON(fock)
      endif
      deallocate(EMAT)

      ! EDIIS
      if (conver_criter > 5) then

         allocate(BMAT(nediist, nediist))
         call ediis_update_energy_fock_rho(energy, fock_op, dens_op, nediist, &
                                           spin, niter)
         call ediis_update_bmat(BMAT, nediist, niter, M_in, spin)


         if (ediis_on) then
            call ediis_get_new_fock(fock, BMAT, nediist, spin)
            call fock_op%Sets_data_ON(fock)
         endif

         deallocate(BMAT)
      endif
   endif

   if ((rho_diff > lvl_shift_cut) .and. (level_shift)) then
      call fock_op%Shift_diag_ON(lvl_shift_en, n_orbs+1)
   endif
   
end subroutine conver_fock

! The following subs are only internal.
subroutine check_convergence(rho_old, rho_new, energy_old, energy_new,&
                             n_iterations, is_converged)
   use fileio        , only: write_energy_convergence
   use converger_data, only: told, Etold, nMax, rho_ls, may_conv, &
                             rho_diff
   ! Calculates convergence criteria in density matrix, and
   ! store new density matrix in Pmat_vec.
   implicit none
   integer     , intent(in)  :: n_iterations
   real(kind=8), intent(in)  :: energy_old, energy_new
   real(kind=8), intent(in)  :: rho_new(:,:), rho_old(:)
   logical     , intent(out) :: is_converged

   integer      :: jj, kk, Rposition, M2
   real(kind=8) :: del, e_diff
   
   M2       = 2 * size(rho_new,1)
   e_diff   = abs(energy_new - energy_old)

   rho_diff = 0.0D0
   do jj = 1 , size(rho_new,1)
   do kk = jj, size(rho_new,1)
         Rposition = kk + (M2 - jj) * (jj -1) / 2
         del       = (rho_new(jj,kk) - rho_old(Rposition)) * sqrt(2.0D0)
         rho_diff  = rho_diff + del * del
      enddo
   enddo
   rho_diff = sqrt(rho_diff) / dble(size(rho_new,1))

   if (.not. may_conv) rho_diff = -1.0D0

   is_converged = .false.
   if ((rho_diff < told) .and. (e_diff < Etold)) is_converged = .true.
   call write_energy_convergence(n_iterations, energy_new, rho_diff, told, &
                                 e_diff, etold)

   if (.not. may_conv) rho_diff = 100.0D0
end subroutine check_convergence
   