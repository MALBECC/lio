module converger_subs

   implicit none
contains

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

   if (Rho_LS > 0 .and. (conver_criter /=1)) then
      hybrid_converg = .false.
      DIIS           = .false.
      conver_criter  = 1
      write(*,'(A)') &
      '  WARNING - Turning to damping-only convergence (linear search).'
    end if
   

end subroutine converger_options_check

subroutine converger_init( M_in, OPshell )
   use converger_data, only: fockm, FP_PFm, conver_criter, fock_damped, &
                             hagodiis, bcoef, ndiis, EMAT2, energy_list
   implicit none
   integer         , intent(in) :: M_in
   logical         , intent(in) :: OPshell

   hagodiis = .false.
   
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
   call P_linearsearch_fin()
end subroutine converger_finalise

subroutine conver (niter, good, M_in, rho_op, fock_op, spin, energy, n_orbs,&
#ifdef CUBLAS
                      devPtrX, devPtrY)
#else
                      Xmat, Ymat)
#endif
   use converger_data  , only: damping_factor, hagodiis, fockm, FP_PFm, ndiis, &
                               fock_damped, bcoef, EMAT2, conver_criter,       &
                               good_cut, energy_list, DIIS_bias, level_shift,  &
                               lvl_shift_en, lvl_shift_cut, changed_to_LS
   use typedef_operator, only: operator
   use fileio_data     , only: verbose
   use linear_algebra  , only: matmuldiag
   implicit none
   ! Spin allows to store correctly alpha or beta information. - Carlos
   integer       , intent(in)    :: niter, M_in, spin, n_orbs   
#ifdef  CUBLAS
   integer*8     , intent(in)    :: devPtrX, devPtrY
#else
   real(kind=8)  , intent(in)    :: Xmat(M_in,M_in), Ymat(M_in,M_in)
#endif
   real(kind=8)  , intent(in)    :: good, energy
   type(operator), intent(inout) :: rho_op, fock_op

   integer      :: ndiist, ii, jj, kk, kknew, lwork, info, Emin_index
   real(kind=8), allocatable :: fock00(:,:), EMAT(:,:), diag1(:,:),      &
                                suma(:,:), scratch1(:,:), scratch2(:,:), &
                                fock(:,:), rho(:,:), work(:)


! INITIALIZATION
! If DIIS is turned on, update fockm with the current transformed F' (into ON
! basis) and update FP_PFm with the current transformed [F',P']
!
! (1)     Calculate F' and [F',P']
!       update fockm with F'
! now, scratch1 = A = F' * P'; scratch2 = A^T
! [F',P'] = A - A^T
! BASE CHANGE HAPPENS INSIDE OF FOCK_COMMUTS

   allocate(fock00(M_in,M_in), fock(M_in,M_in), rho(M_in,M_in), work(1000))
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
            if ((good < good_cut) .and. (niter > 2) .and. &
               (.not. changed_to_LS)) then
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
   if (changed_to_LS) hagodiis = .false.

   ndiist = min( niter, ndiis )
   if (conver_criter /= 1) then      
      allocate( suma(M_in, M_in), diag1(M_in, M_in) )
      allocate( scratch1(M_in, M_in), scratch2(M_in, M_in) )
      suma = 0.0D0
      diag1 = 0.0D0
      scratch1 = 0.0D0
      scratch2 = 0.0D0
      
      ! If DIIS is turned on, update fockm with the current transformed F' (into ON
      ! basis) and update FP_PFm with the current transformed [F',P']

      do jj = ndiis-(ndiist-1), ndiis-1
         fockm(:,:,jj,spin)  = fockm(:,:,jj+1,spin)
         FP_PFm(:,:,jj,spin) = FP_PFm(:,:,jj+1,spin)
      enddo
      
#ifdef CUBLAS
      call rho_op%BChange_AOtoON(devPtrY, M_in, 'r')
      call fock_op%BChange_AOtoON(devPtrX, M_in, 'r')
#else
      call rho_op%BChange_AOtoON(Ymat, M_in, 'r')
      call fock_op%BChange_AOtoON(Xmat,M_in, 'r')
#endif
      call rho_op%Gets_data_ON(rho)
      call fock_op%Commut_data_r(rho, scratch1, M_in)

      FP_PFm(:,:,ndiis,spin) = scratch1(:,:)
      call fock_op%Gets_data_ON( fockm(:,:,ndiis,spin) )
   endif

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

      ! Before ndiis iterations, we just start from the old EMAT
      if ((niter > 1) .and. (niter <= ndiis)) then
         EMAT = 0.0D0
         do jj = 1, ndiist-1
         do ii = 1, ndiist-1
            EMAT(ii,jj) = EMAT2(ii,jj,spin)
         enddo
         enddo
      ! After ndiis iterations, we start shifting the oldest iteration stored
      else if (niter > ndiis) then
         EMAT = 0.0D0
         do jj = 1, ndiist-1
         do ii = 1, ndiist-1
            EMAT(ii,jj) = EMAT2(ii+1,jj+1,spin)
         enddo
         enddo
      endif

      ! scratch1 and scratch2 store the commutations from different iterations.
      do kk = 1, ndiist
         kknew = kk + (ndiis - ndiist)
         scratch1(:,:) = FP_PFm(:,:,ndiis,spin)
         scratch2(:,:) = FP_PFm(:,:,kknew,spin)

         call matmuldiag( scratch1, scratch2, diag1, M_in )
         EMAT(ndiist,kk) = 0.0d0

         if (kk /= ndiist) EMAT(kk,ndiist) = 0.0d0
         do ii = 1, M_in
            EMAT(ndiist,kk) = EMAT(ndiist,kk) + diag1(ii,ii)
            if (kk /= ndiist) then
               EMAT(kk,ndiist) = EMAT(ndiist,kk)
            endif
         enddo
      enddo

      do kk = 1, ndiist
         EMAT(kk,ndiist+1) = -1.0d0
         EMAT(ndiist+1,kk) = -1.0d0
      enddo
      EMAT(ndiist+1, ndiist+1)= 0.0d0
      EMAT2(1:ndiist+1,1:ndiist+1,spin) = EMAT

      !   THE MATRIX EMAT SHOULD HAVE THE FOLLOWING SHAPE:
      !      |<E(1)*E(1)>  <E(1)*E(2)> ...   -1.0|
      !      |<E(2)*E(1)>  <E(2)*E(2)> ...   -1.0|
      !      |<E(3)*E(1)>  <E(3)*E(2)> ...   -1.0|
      !      |<E(4)*E(1)>  <E(4)*E(2)> ...   -1.0|
      !      |     .            .      ...     . |
      !      |   -1.0         -1.0     ...    0. |
      !   WHERE <E(I)*E(J)> IS THE SCALAR PRODUCT OF [F*P] FOR ITERATION I
      !   TIMES [F*P] FOR ITERATION J.   

      if (hagodiis) then
         if ((conver_criter == 4) .or. (conver_criter == 5)) then
            if (niter > ndiis) then
               Emin_index = minloc(energy_list,1)
            else
               Emin_index = minloc(energy_list((ndiis - ndiist +1):ndiis),1)
            endif

            do ii = 1, Emin_index -1
               Emat(ii,ii) = Emat(ii,ii) * DIIS_bias
            enddo
            do ii = Emin_index +1, ndiist
               Emat(ii,ii) = Emat(ii,ii) * DIIS_bias
            enddo
         endif

         do ii = 1, ndiist
            bcoef(ii,spin) = 0.0d0
         enddo
         bcoef(ndiist+1,spin) = -1.0d0

         ! First call to DGELS sets optimal WORK size. Second call solves the
         ! A*X = B (EMAT * Ci = bCoef) problem, with bCoef also storing the
         ! result.
         LWORK = -1
         CALL DGELS( 'No transpose',ndiist+1, ndiist+1, 1, EMAT, &
                     ndiist+1, bcoef(:,spin), ndiist+1, WORK, LWORK, INFO )

         LWORK = MIN( 1000, INT( WORK( 1 ) ) )
         CALL DGELS( 'No transpose',ndiist+1, ndiist+1, 1, EMAT, &
                     ndiist+1, bcoef(:,spin), ndiist+1, WORK, LWORK, INFO )

         ! Build new Fock as a linear combination of previous steps.
         suma = 0.0D0
         do kk = 1, ndiist
            kknew = kk + (ndiis - ndiist)
            do ii = 1, M_in
            do jj = 1, M_in
               suma(ii,jj) = suma(ii,jj) + bcoef(kk,spin) * &
                                           fockm(ii,jj,kknew,spin)
            enddo
            enddo
         enddo
         fock = suma
         call fock_op%Sets_data_ON(fock)

      endif
   endif

   if ((good > lvl_shift_cut) .and. (level_shift)) then
   if (hagodiis) &
       call fock_op%Shift_diag_ON(lvl_shift_en, n_orbs+1)
   endif
   deallocate (work)
   
end subroutine conver

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Lineal search subroutines
!
! This procedure improve convergency in systems finding best lineal
! combinations of density matrix in steps n and n-1. To do so, it
! evaluates different combinations of rho(n) and rho(n-1), using the
! variable lambda as a weight. The possibilities it evaluates also
! depend on the performance of the algorithm in previous steps,
! which is regulated by the parameters Elast and Pstepsize.
!
! FFR comments: Rho_LS, changed_to_LS and P_oscilation_analysis are
! communication variables that appear here and in some external subs.
! P_oscilation_analysis and changed_to_LS are not even used inside here.
! These should be dealt with differently.
!
!------------------------------------------------------------------------------!
! LOG:
! V 1.00 September 2018 Final version - Nicolas Foglia
! V 1.01 September 2018 adaptation - FFR
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine P_linearsearch_init(MM, open_shell, rho_alp, rho_bet)
   use converger_data, only: rho_lambda0, rho_lambda1, rho_lambda0_alpha, &
                             rho_lambda1_alpha, rho_lambda0_betha,        &
                             rho_lambda1_betha, rho_LS
   implicit none
   integer     , intent(in) :: MM
   logical     , intent(in) :: open_shell
   real(kind=8), intent(in) :: rho_alp(MM), rho_bet(MM)

   if (rho_LS < 1) return

   if (.not. allocated(rho_lambda1)) allocate(rho_lambda1(MM))
   if (.not. allocated(rho_lambda0)) allocate(rho_lambda0(MM))
   if (open_shell) then
      if (.not. allocated(rho_lambda1_alpha)) allocate(rho_lambda1_alpha(MM))
      if (.not. allocated(rho_lambda0_alpha)) allocate(rho_lambda0_alpha(MM))
      if (.not. allocated(rho_lambda1_betha)) allocate(rho_lambda1_betha(MM))
      if (.not. allocated(rho_lambda0_betha)) allocate(rho_lambda0_betha(MM))
      rho_lambda0_alpha = rho_alp
      rho_lambda0_betha = rho_bet
   end if
end subroutine P_linearsearch_init

subroutine P_linearsearch_fin()
   use converger_data, only: rho_lambda0, rho_lambda1, rho_lambda0_alpha, &
                             rho_lambda1_alpha, rho_lambda0_betha,        &
                             rho_lambda1_betha, rho_LS
   
   if (rho_LS < 1) return
   if (allocated(rho_lambda1))       deallocate(rho_lambda1)
   if (allocated(rho_lambda0))       deallocate(rho_lambda0)
   if (allocated(rho_lambda1_alpha)) deallocate(rho_lambda1_alpha)
   if (allocated(rho_lambda0_alpha)) deallocate(rho_lambda0_alpha)
   if (allocated(rho_lambda1_betha)) deallocate(rho_lambda1_betha)
   if (allocated(rho_lambda0_betha)) deallocate(rho_lambda0_betha)
end subroutine P_linearsearch_fin

subroutine P_conver(niter, En, E1, E2, Ex, good, rho_new, rho_old, Hmat_vec, &
                    Fmat_vec, Fmat_vec2, Gmat_vec, Ginv_vec, int_memo,       &
                    rhoa_new, rhob_new, rhoa_old, rhob_old)
   ! rho_LS values:
   !   = 0 calculate convergence criteria for actual density matrix.
   !   = 1 do linear search for density matrix if energy > previous energy.
   !   = 2 do linear search for density matrix in all steps.
   use converger_data, only: rho_LS

   implicit none
   ! Step number and energies.
   integer     , intent(in)    :: niter
   real(kind=8), intent(inout) :: E1, E2, Ex, En
   
   ! Variables for electron integrals.
   logical     , intent(in)    :: int_memo
   real(kind=8), intent(inout) :: Hmat_vec(:), Fmat_vec(:), Fmat_vec2(:), &
                                  Gmat_vec(:), Ginv_vec(:)

   ! Convergence criteria
   real(kind=8), intent(out)   :: good

   ! Density matrices of current and previous step. Due to how SCF works,
   ! the current step is in matrix form and the previous in vector form.
   real(kind=8), intent(inout) :: rho_old(:)
   real(kind=8), intent(inout) :: rho_new(:,:)
   real(kind=8), optional, intent(inout) :: rhoa_old(:)  , rhob_old(:)
   real(kind=8), optional, intent(inout) :: rhoa_new(:,:), rhob_new(:,:)

   ! True if predicted density != density of previous steep
   logical :: may_conv = .true.

   select case (rho_LS)
      case (0,-1)
      case (1, 2)

         ! Makes separate calls for open and closed shell.
         if (present(rhoa_old)) then
            call P_linear_calc(niter, En, E1, E2, Ex, may_conv, rho_new, &
                               rho_old, Hmat_vec, Fmat_vec, Fmat_vec2,   &
                               Gmat_vec, Ginv_vec, int_memo)
         else
            call P_linear_calc(niter, En, E1, E2, Ex, may_conv, rho_new, &
                               rho_old, Hmat_vec, Fmat_vec, Fmat_vec2,   &
                               Gmat_vec, Ginv_vec, int_memo,             &
                               rhoa_new, rhob_new, rhoa_old, rhob_old)
         endif
      case default
         write(*,'(A,I2)') &
            "ERROR - P_conver: Wrong Rho_LS value, current value is ", Rho_LS
         stop
   end select
      
   call P_calc_fluctuation(rho_old, rho_new, good)
   if (.not. may_conv) good = -1.0D0
end subroutine P_conver

! The following subs are only internal.
subroutine P_calc_fluctuation(rho_old, rho_new, good)
   ! Calculates convergence criteria in density matrix, and
   ! store new density matrix in Pmat_vec.
   real(kind=8), intent(in)    :: rho_new(:,:), rho_old(:)
   real(kind=8), intent(out)   :: good

   integer      :: jj, kk, Rposition, M2
   real(kind=8) :: del
   
   M2   = 2 * size(rho_new,1)
   good = 0.0D0

   do jj = 1 , size(rho_new,1)
   do kk = jj, size(rho_new,1)
         Rposition = kk + (M2 - jj) * (jj -1) / 2
         del  = (rho_new(jj,kk) - rho_old(Rposition)) * sqrt(2.0D0)
         good = good + del * del
      enddo
   enddo
   good = sqrt(good) / dble(size(rho_new,1))
end subroutine P_calc_fluctuation

subroutine P_linear_calc(niter, En, E1, E2, Ex, may_conv, rho_new, rho_old, &
                         Hmat_vec, Fmat_vec, Fmat_vec2, Gmat_vec, Ginv_vec, &
                         int_memo, rhoa_new, rhob_new, rhoa_old, rhob_old)
   use liosubs       , only: line_search
   use converger_data, only: rho_lambda0, rho_lambda1, rho_lambda0_alpha,   &
                             rho_lambda1_alpha, rho_lambda0_betha,          &
                             rho_lambda1_betha, Elast, pstepsize, rho_LS
   implicit none
   integer     , intent(in)    :: niter
   logical     , intent(in)    :: int_memo
   logical     , intent(inout) :: may_conv
   real(kind=8), intent(inout) :: En, E1, E2, Ex
   real(kind=8), intent(inout) :: Hmat_vec(:), Fmat_vec(:), Fmat_vec2(:), &
                                  Gmat_vec(:), Ginv_vec(:)
   ! Density matrices of the current and previous steps.
   real(kind=8), intent(inout) :: rho_old(:)
   real(kind=8), intent(inout) :: rho_new(:,:)
   real(kind=8), optional, intent(inout) :: rhoa_old(:)  , rhob_old(:)
   real(kind=8), optional, intent(inout) :: rhoa_new(:,:), rhob_new(:,:)

   ! Values for combination of density matrices 
   real(kind=8) :: dlambda, Blambda 
   
   ! Auxiliars
   integer :: M, MM, M2, jj, kk, Rposition, ilambda
   logical :: open_shell = .false.
   real(kind=8), allocatable :: RMM_temp(:),E_lambda(:)

   M  = size(rho_new,1)
   MM = size(rho_old,1)
   M2 = 2 * M
   if (present(rhoa_old)) open_shell = .true.

   allocate(E_lambda(0:10))
   allocate(RMM_temp(1:MM))

   if (niter == 1)  then
      Pstepsize   = 1.d0
      rho_lambda0 = rho_old
      if (open_shell) rho_lambda0_alpha = rhoa_old
      if (open_shell) rho_lambda0_betha = rhob_old
   end if
      
   do jj=1,M
   do kk=jj,M
      Rposition = kk + (M2 - jj) * (jj - 1) / 2
      rho_lambda1(Rposition) = rho_new(jj,kk)
      if (open_shell) rho_lambda1_alpha(Rposition) = rhoa_new(jj,kk)
      if (open_shell) rho_lambda1_betha(Rposition) = rhob_new(jj,kk)
   enddo
   enddo
   
   rho_old = rho_lambda1
   if (open_shell) rhoa_old = rho_lambda1_alpha
   if (open_shell) rhob_old = rho_lambda1_betha
   
   call give_me_energy(E_lambda(10), En, E1, E2, Ex, rho_old, Hmat_vec,    &
                       Fmat_vec, Fmat_vec2, Gmat_vec, Ginv_vec, open_shell,&
                       int_memo)
   
   if ((Elast < E_lambda(10)) .or. (Rho_LS == 2)) then
      write(*,*) "Lambda this step: ", E_lambda(10), ", last step: ", Elast
      write(*,*) "Doing lineal interpolation in Rho."
      do ilambda = 0, 10
         dlambda = Pstepsize * dble(ilambda) / 10.d0
         if (dlambda > 1.d0) STOP "dlambda > 1.d0"

         rho_old = rho_lambda0 * (1.d0 - dlambda) + rho_lambda1 * dlambda 
         if (open_shell) rhoa_old = rho_lambda0_alpha * (1.d0 - dlambda) + &
                                    rho_lambda1_alpha * dlambda
         if (open_shell) rhob_old = rho_lambda0_betha * (1.d0 - dlambda) + &
                                    rho_lambda1_betha * dlambda
         call give_me_energy(E_lambda(ilambda), En, E1, E2, Ex, rho_old,       &
                             Hmat_vec, Fmat_vec, Fmat_vec2, Gmat_vec, Ginv_vec,&
                             open_shell, int_memo)
         write(*,*) "Step n°", ilambda, ", energy: ", E_lambda(ilambda)
      end do
   
      call line_search(11, E_lambda, 1d0, Blambda)
      if (Blambda >= 1.d0) Blambda = Blambda - 1.0d0
      write(*,*) "Best lambda: ", Blambda
      Blambda = Blambda * Pstepsize / 10.d0
      write(*,*) "Fluctuation: ", Blambda
   else
      Blambda = Pstepsize
   end if
   
   rho_old = rho_lambda0 * (1.d0 - Blambda) + rho_lambda1 * Blambda
   if (open_shell) rhoa_old = rho_lambda0_alpha * (1.d0 - Blambda) + &
                              rho_lambda1_alpha * Blambda
   if (open_shell) rhob_old = rho_lambda0_betha * (1.d0 - Blambda) + &
                              rho_lambda1_betha * Blambda
   
   do jj = 1, M
      do kk = jj, M
         Rposition = kk + (M2 - jj) * (jj -1) / 2
         rho_new(jj,kk) = rho_old(Rposition)
         if (open_shell) rhoa_new(jj,kk) = rhoa_old(Rposition)
         if (open_shell) rhob_new(jj,kk) = rhob_old(Rposition)
      enddo
   enddo

   call give_me_energy(Elast, En, E1, E2, Ex, rho_old, Hmat_vec, &
                       Fmat_vec, Fmat_vec2, Gmat_vec, Ginv_vec,  &
                       open_shell, int_memo)
      
   RMM_temp    = rho_old
   rho_old     = rho_lambda0
   rho_lambda0 = RMM_temp
   
   if (open_shell) then
      RMM_temp          = rhoa_old
      rhoa_old          = rho_lambda0_alpha
      rho_lambda0_alpha = RMM_temp
   
      RMM_temp          = rhob_old
      rhob_old          = rho_lambda0_betha
      rho_lambda0_betha = RMM_temp
   end if
   
   if (Blambda <= 4.d-1 * Pstepsize) Pstepsize = Pstepsize * 0.5d0 
   if (Blambda >= 8.d-1 * Pstepsize) Pstepsize = Pstepsize * 1.2d0
   if (Pstepsize > 1.d0) Pstepsize = 1.d0
   if ((Blambda <= 2.d-1*Pstepsize) .and. (Pstepsize > 1d-4)) may_conv = .false.

   deallocate(E_lambda, RMM_temp)
end subroutine P_linear_calc

subroutine give_me_energy(E, En, E1, E2, Ex, Pmat_vec, Hmat_vec, Fmat_vec, &
                          Fmat_vec2, Gmat_vec, Ginv_vec, open_shell, int_memo)
   !  return Energy components for a density matrix stored in Pmat_vec
   use faint_cpu, only: int3lu
   implicit none
   logical     , intent(in)    :: open_shell, int_memo
   real(kind=8), intent(in)    :: En
   real(kind=8), intent(out)   :: E, E1, E2, Ex
   real(kind=8), intent(inout) :: Pmat_vec(:), Hmat_vec(:), Fmat_vec(:), &
                                  Fmat_vec2(:), Gmat_vec(:), Ginv_vec(:)
   integer :: kk
      
   E  = 0.0D0; E1 = 0.0D0
   E2 = 0.0D0; Ex = 0.0D0
   
   do kk = 1, size(Pmat_vec,1)
      E1 = E1 + Pmat_vec(kk) * Hmat_vec(kk) !Computes 1e energy
   enddo

   ! Computes Coulomb part of Fock, and energy on E2.
   call int3lu(E2, Pmat_vec, Fmat_vec2, Fmat_vec, Gmat_vec, Ginv_vec, &
               Hmat_vec, open_shell, int_memo)

   ! Computes XC integration / Fock elements.
   call g2g_solve_groups(0,Ex,0)

   ! Adds all energy components.
   E = E1 + E2 + En + Ex

end subroutine give_me_energy
   
end module converger_subs
