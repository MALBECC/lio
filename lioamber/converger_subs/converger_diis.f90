subroutine diis_init(M_in, OPshell)
   use converger_data, only: fockm, FP_PFm, conver_method, bcoef, ndiis, &
                             EMAT, energy_list
   implicit none
   integer, intent(in) :: M_in
   logical, intent(in) :: OPshell

   if (conver_method < 1) return

   if (.not. allocated(bcoef)) allocate(bcoef(ndiis+1) )
   if (.not. allocated(EMAT )) allocate(EMAT (ndiis,ndiis))

   if (OPshell) then
      if (.not. allocated(fockm) ) allocate(fockm (M_in, M_in, ndiis, 2))
      if (.not. allocated(FP_PFm)) allocate(FP_PFm(M_in, M_in, ndiis, 2))
   else
      if (.not. allocated(fockm) ) allocate(fockm (M_in, M_in, ndiis, 1))
      if (.not. allocated(FP_PFm)) allocate(FP_PFm(M_in, M_in, ndiis, 1))
   end if
   fockm   = 0.0D0
   FP_PFm  = 0.0D0
   bcoef   = 0.0D0
   EMAT    = 0.0D0

   if (conver_method > 3) then
      if (.not. allocated(energy_list)) allocate(energy_list(ndiis))
      energy_list = 0.0D0
   endif
end subroutine diis_init

subroutine diis_finalise()
   use converger_data, only: fockm, FP_PFm, conver_method, bcoef, EMAT, &
                             energy_list

   if (conver_method < 1) return

   if (allocated(fockm) ) deallocate(fockm )
   if (allocated(FP_PFm)) deallocate(FP_PFm)
   if (allocated(bcoef) ) deallocate(bcoef )
   if (allocated(EMAT ) ) deallocate(EMAT  )
   if ((conver_method > 3) .and. allocated(energy_list)) deallocate(energy_list)

end subroutine diis_finalise

subroutine diis_fock_commut(dens_op, fock_op, dens, M_in, spin, ndiist)
   use converger_data  , only: fockm, FP_PFm, ndiis
   use typedef_operator, only: operator
   
   implicit none
   integer       , intent(in)    :: M_in, ndiist, spin
   real(kind=8)  , intent(inout) :: dens(:,:)
   type(operator), intent(inout) :: dens_op, fock_op

   integer :: jj

   ! If DIIS is turned on, update fockm with the current transformed F' (into ON
   ! basis) and update FP_PFm with the current transformed [F',P']
   do jj = ndiis-(ndiist-1), ndiis-1
      fockm(:,:,jj,spin)  = fockm(:,:,jj+1,spin)
      FP_PFm(:,:,jj,spin) = FP_PFm(:,:,jj+1,spin)
   enddo
   
   call dens_op%Gets_data_ON(dens)
   call fock_op%Commut_data_r(dens, FP_PFm(:,:,ndiis,spin), M_in)
   call fock_op%Gets_data_ON( fockm(:,:,ndiis,spin) )
   
end subroutine diis_fock_commut

subroutine diis_get_error(M_in, spin, verbose)
   use converger_data, only: ndiis, FP_PFm, diis_error
   integer, intent(in)  :: M_in, spin, verbose

   integer      :: ii, jj
   real(kind=8) :: max_error, avg_error
   
   max_error = maxval(abs(FP_PFm(:,:,ndiis,spin)))

   if (verbose > 3) then
      avg_error = 0.0D0
      do ii=1, M_in
      do jj=1, M_in
         avg_error = avg_error + &
                     FP_PFm(ii,jj,ndiis,spin) * FP_PFm(ii,jj,ndiis,spin)
      enddo
      enddo
      avg_error = sqrt(avg_error / M_in)

      if (spin == 1) write(*,'(2x,A22,ES14.7,A18,ES14.7)')        &
                           "Max Alpha DIIS error: ", max_error, &
                           " | Average error: ", avg_error
      if (spin == 2) write(*,'(2x,A22,ES14.7,A18,ES14.7)')        &
                           "Max Beta  DIIS error: ", max_error, &
                           " | Average error: ", avg_error
   endif

   ! Updates DIIS error with the greatest value for alpha or beta.
   if (spin == 1) then
      diis_error = max_error
   else if (max_error > diis_error) then
      diis_error = max_error
   endif

end subroutine diis_get_error

subroutine diis_update_energy(energy)
   use converger_data, only: energy_list, ndiis
   implicit none
   real(kind=8), intent(in) :: energy
   integer :: ii

   do ii = 1, ndiis -1
      energy_list(ii) = energy_list(ii+1)
   enddo

   energy_list(ndiis) = energy
end subroutine diis_update_energy

subroutine diis_update_emat(niter, ndiist, M_in, open_shell)
   use converger_data, only: EMAT, ndiis, FP_PFm
   use linear_algebra, only: matmuldiag

   implicit none
   integer     , intent(in)  :: niter, ndiist, M_in
   logical     , intent(in)  :: open_shell

   integer                   :: ii, jj, kind
   real(kind=8), allocatable :: diag1(:,:)

   ! Before ndiis iterations, we just start from the old EMAT
   ! After ndiis iterations, we start shifting the oldest iteration stored
   if (niter > ndiis) then
      do jj = 1, ndiist-1
      do ii = 1, ndiist-1
         EMAT(ii,jj) = EMAT(ii+1,jj+1)
      enddo
      enddo
   endif

   allocate( diag1(M_in, M_in) )
   diag1 = 0.0D0
   do ii = 1, ndiist
      kind = ii + (ndiis - ndiist)
      EMAT(ndiist,ii) = 0.0d0

      ! Make diagonal-only multiplication for the commutations of different
      ! iterations.
      call matmuldiag(FP_PFm(:,:, ndiis, 1), FP_PFm(:,:, kind, 1), &
                      diag1, M_in)
      do jj = 1, M_in
         EMAT(ndiist,ii) = EMAT(ndiist,ii) + diag1(jj,jj)
      enddo
      if (open_shell) then
         call matmuldiag(FP_PFm(:,:, ndiis, 2), FP_PFm(:,:, kind, 2), &
                         diag1, M_in)
         do jj = 1, M_in
            EMAT(ndiist,ii) = EMAT(ndiist,ii) + diag1(jj,jj)
         enddo
      endif
      EMAT(ii,ndiist) = EMAT(ndiist,ii)
   enddo
   deallocate(diag1)

end subroutine diis_update_emat

subroutine diis_emat_bias(ndiist)
   use converger_data, only: EMAT, ndiis, energy_list, DIIS_bias
   
   implicit none
   integer, intent(in)    :: ndiist
   integer :: Emin_index, ii
   
   Emin_index = minloc(energy_list((ndiis - ndiist +1):ndiis),1)   
   do ii = 1, Emin_index -1
      Emat(ii,ii) = Emat(ii,ii) * DIIS_bias
   enddo
   do ii = Emin_index +1, ndiist
      Emat(ii,ii) = Emat(ii,ii) * DIIS_bias
   enddo

end subroutine diis_emat_bias

subroutine diis_get_new_fock(fock, ndiist, M_in, spin)
   use converger_data, only: EMAT, ndiis, bcoef, fockm
   implicit none
   integer     , intent(in)  :: ndiist, M_in, spin
   real(kind=8), intent(out) :: fock(:,:)

   integer :: ii, jj, kk, kknew, LWORK, INFO
   real(kind=8), allocatable :: work(:), EMAT_aux(:,:)

   
   ! First call to this routines gets DIIS coefficients.
   if (spin == 1) then
      allocate(EMAT_aux(ndiist+1,ndiist+1))
      do ii = 1, ndiist
         do jj = 1, ndiist
            EMAT_aux(ii,jj) = EMAT(ii,jj)
         enddo
         EMAT_aux(ii,ndiist+1) = -1.0D0
         EMAT_aux(ndiist+1,ii) = -1.0D0
      enddo
      EMAT_aux(ndiist+1, ndiist+1) = 0.0D0
      
      !   THE MATRIX EMAT SHOULD HAVE THE FOLLOWING SHAPE:
      !      |<E(1)*E(1)>  <E(1)*E(2)> ...   -1.0|
      !      |<E(2)*E(1)>  <E(2)*E(2)> ...   -1.0|
      !      |<E(3)*E(1)>  <E(3)*E(2)> ...   -1.0|
      !      |<E(4)*E(1)>  <E(4)*E(2)> ...   -1.0|
      !      |     .            .      ...     . |
      !      |   -1.0         -1.0     ...    0. |
      !   WHERE <E(I)*E(J)> IS THE SCALAR PRODUCT OF [F*P] FOR ITERATION I
      !   TIMES [F*P] FOR ITERATION J.

      do ii = 1, ndiist
         bcoef(ii) = 0.0d0
      enddo
      bcoef(ndiist+1) = -1.0d0

      ! First call to DGELS sets optimal WORK size. Second call solves the
      ! A*X = B (EMAT * Ci = bCoef) problem, with bCoef also storing the
      ! result.
      allocate(work(1000))
      LWORK = -1
      CALL DGELS('No transpose',ndiist+1, ndiist+1, 1, EMAT_aux, ndiist+1, &
                 bcoef, ndiist+1, WORK, LWORK, INFO )

      LWORK = MIN( 1000, INT( WORK( 1 ) ) )
      CALL DGELS('No transpose',ndiist+1, ndiist+1, 1, EMAT_aux, ndiist+1, &
                 bcoef, ndiist+1, WORK, LWORK, INFO )
      deallocate (work, EMAT_aux)
   endif

   ! Build new Fock as an extrapolation of previous steps.
   fock = 0.0D0
   do kk = 1, ndiist
      kknew = kk + (ndiis - ndiist)
      do ii = 1, M_in
      do jj = 1, M_in
         fock(ii,jj) = fock(ii,jj) + bcoef(kk) * fockm(ii,jj,kknew,spin)
      enddo
      enddo
   enddo

end subroutine diis_get_new_fock