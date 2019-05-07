
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
   use converger_data, only: fockm, FP_PFm, conver_criter, fock_damped,  &
                             hagodiis, bcoef, ndiis, EMAT2, energy_list, &
                             told, etold
   implicit none
   integer         , intent(in) :: M_in
   logical         , intent(in) :: OPshell

   hagodiis = .false.

   write(*,'(2x,A)') "Convergence criteria are: "
   write(*,'(2x,ES8.2,A33,ES8.2,A26)')                 &
            told, " in Rho mean squared difference, ", &
            Etold, " Eh in energy differences."

   ! Added to change from damping to DIIS. - Nick
   if (conver_criter /= 1) then
      if(OPshell) then
         if (.not. allocated(fockm)  ) allocate(fockm (M_in, M_in, ndiis, 2))
         if (.not. allocated(FP_PFm) ) allocate(FP_PFm(M_in, M_in, ndiis, 2))
         if (.not. allocated(bcoef)  ) allocate(bcoef(ndiis+1, 2) )
         if (.not. allocated(EMAT2)  ) allocate(EMAT2(ndiis+1,ndiis+1,2))
      else
         if (.not. allocated(fockm)  ) allocate(fockm (M_in, M_in, ndiis, 1))
         if (.not. allocated(FP_PFm) ) allocate(FP_PFm(M_in, M_in, ndiis, 1))
         if (.not. allocated(bcoef)  ) allocate(bcoef (ndiis+1, 1))
         if (.not. allocated(EMAT2)  ) allocate(EMAT2(ndiis+1,ndiis+1,1))
      end if
      fockm   = 0.0D0
      FP_PFm  = 0.0D0
      bcoef   = 0.0D0
      EMAT2   = 0.0D0
      if ((conver_criter == 4) .or. (conver_criter == 5)) then
         if (.not. allocated(energy_list)) allocate(energy_list(ndiis))
         energy_list = 0.0D0
      endif
   endif

   if(OPshell) then
      if (.not. allocated(fock_damped) ) allocate(fock_damped(M_in, M_in, 2))
   else
      if (.not. allocated(fock_damped) ) allocate(fock_damped(M_in, M_in, 1))
   end if
   fock_damped(:,:,:) = 0.0D0
end subroutine converger_init

subroutine converger_finalise()
   use converger_data, only: fockm, FP_PFm, conver_criter, fock_damped, &
                             hagodiis, bcoef, ndiis, EMAT2, energy_list
   implicit none
   
   if (conver_criter /= 1) then
      if (allocated(fockm) ) deallocate(fockm )
      if (allocated(FP_PFm)) deallocate(FP_PFm)
      if (allocated(bcoef) ) deallocate(bcoef )
      if (allocated(EMAT2) ) deallocate(EMAT2 )
      if ((conver_criter == 4) .or. (conver_criter == 5)) then
         if (allocated(energy_list)) deallocate(energy_list)
      endif
   endif

   if (allocated(fock_damped)) deallocate(fock_damped)
   call rho_ls_finalise()
end subroutine converger_finalise

subroutine conver_fock(niter, M_in, rho_op, fock_op, spin, energy, n_orbs, &
#ifdef CUBLAS
                      devPtrX, devPtrY)
#else
                      Xmat, Ymat)
#endif
   use converger_data  , only: damping_factor, hagodiis, fockm, FP_PFm, ndiis, &
                               fock_damped, bcoef, EMAT2, conver_criter,       &
                               good_cut, energy_list, DIIS_bias, level_shift,  &
                               lvl_shift_en, lvl_shift_cut, rho_ls, rho_diff
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
   type(operator), intent(inout) :: rho_op, fock_op

   integer      :: ndiist, ii
   real(kind=8), allocatable :: fock00(:,:), EMAT(:,:), suma(:,:), &
                                fock(:,:), rho(:,:)


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
   call rho_op%Gets_data_AO(rho)
   call fock_op%Gets_data_AO(fock00)

   ! Checks convergence criteria.
   select case (conver_criter)
      ! Always do damping
      case (1)
         hagodiis = .false.

      ! Damping the first two steps, diis afterwards
      case (2,4)
         if (niter > 2) then
            hagodiis = .true.
         else
            hagodiis = .false.
         endif

      ! Damping until good enough, diis afterwards
      case(3,5)
         if (.not. hagodiis) then
            if ((rho_diff < good_cut) .and. (niter > 2) .and. (rho_LS < 2)) then
               if (verbose > 3) &
                  write(6,'(A,I4)') "  Changing to DIIS at step: ", niter
               hagodiis = .true.
            endif
         endif

      case default
         write(*,'(A,I4)') 'ERROR - Wrong conver_criter = ', conver_criter
         stop
   endselect

   ! Turn off diis is calculation when change to lineal search, Nick
   if (rho_ls > 1) hagodiis = .false.

   ndiist = min( niter, ndiis )
   if (conver_criter /= 1) call diis_fock_commut(rho_op, fock_op, rho, M_in, &
#ifdef CUBLAS
      spin, ndiist, devPtrX, devPtrY)
#else
      spin, ndiist, Xmat, Ymat)
#endif

   ! THIS IS DAMPING 
   ! THIS IS SPARTA!
   ! If we are not doing diis this iteration, apply damping to F, save this
   ! F in fock_damped for next iteration's damping and put F' = X^T * F * X in
   ! fock the newly constructed damped matrix is stored, for next iteration in
   ! fock_damped
   if (.not. hagodiis) then
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

   ! DIIS
   if (conver_criter /= 1) then
      if (((conver_criter == 4) .or. (conver_criter == 5)) .and. (spin == 1) &
          .and. (niter > 1)) then
         do ii = 1, ndiis -1
            energy_list(ii) = energy_list(ii+1)
         enddo
         energy_list(ndiis) = energy
      endif

      allocate(EMAT(ndiist+1,ndiist+1))
      call diis_update_emat(EMAT, niter, ndiist, spin, M_in)
   

      if (hagodiis) then
         if ((conver_criter == 4) .or. (conver_criter == 5)) &
            call diis_emat_bias(EMAT, ndiist)

         call diis_get_new_fock(fock, EMAT, ndiist, M_in, spin)
         call fock_op%Sets_data_ON(fock)
      endif
   endif

   if ((rho_diff > lvl_shift_cut) .and. (level_shift)) then
   if (hagodiis) &
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
   