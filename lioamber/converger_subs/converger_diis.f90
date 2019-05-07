subroutine diis_fock_commut(dens_op, fock_op, dens, M_in, spin, ndiist, &
#ifdef CUBLAS   
                            devPtrX, devPtrY)
#else
                            Xmat, Ymat)
#endif
   use converger_data  , only: fockm, FP_PFm, ndiis
   use typedef_operator, only: operator
   
   implicit none
#ifdef  CUBLAS
   integer*8     , intent(in)    :: devPtrX, devPtrY
#else
   real(kind=8)  , intent(in)    :: Xmat(M_in,M_in), Ymat(M_in,M_in)
#endif
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
   
#ifdef CUBLAS
   call dens_op%BChange_AOtoON(devPtrY, M_in, 'r')
   call fock_op%BChange_AOtoON(devPtrX, M_in, 'r')
#else
   call dens_op%BChange_AOtoON(Ymat, M_in, 'r')
   call fock_op%BChange_AOtoON(Xmat, M_in, 'r')
#endif
   call dens_op%Gets_data_ON(dens)
   call fock_op%Commut_data_r(dens, FP_PFm(:,:,ndiis,spin), M_in)
   call fock_op%Gets_data_ON( fockm(:,:,ndiis,spin) )

end subroutine diis_fock_commut

subroutine diis_update_emat(EMAT, niter, ndiist, spin, M_in)
   use converger_data, only: EMAT2, ndiis, FP_PFm
   use linear_algebra, only: matmuldiag

   implicit none
   integer     , intent(in)    :: niter, ndiist, spin, M_in
   real(kind=8), intent(inout) :: EMAT(:,:)

   integer                   :: ii, jj 
   real(kind=8), allocatable :: diag1(:,:)

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

   allocate( diag1(M_in, M_in) )
   diag1 = 0.0D0
   do ii = 1, ndiist
      jj = ii + (ndiis - ndiist)

      ! Make diagonal-only multiplication for the commutations of different
      ! iterations.
      call matmuldiag( FP_PFm(:,:, ndiis, spin), FP_PFm(:,:, jj   , spin), &
                       diag1, M_in )
      EMAT(ndiist,ii) = 0.0d0

      if (ii /= ndiist) EMAT(ii,ndiist) = 0.0d0
      do jj = 1, M_in
         EMAT(ndiist,ii) = EMAT(ndiist,ii) + diag1(jj,jj)
         if (ii /= ndiist) then
            EMAT(ii,ndiist) = EMAT(ndiist,ii)
         endif
      enddo
   enddo

   do ii = 1, ndiist
      EMAT(ii,ndiist+1) = -1.0d0
      EMAT(ndiist+1,ii) = -1.0d0
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
   deallocate(diag1)

end subroutine diis_update_emat

subroutine diis_emat_bias(EMAT, ndiist)
   use converger_data, only: ndiis, energy_list, DIIS_bias
   
   implicit none
   integer     , intent(in)    :: ndiist
   real(kind=8), intent(inout) :: EMAT(:,:)

   integer :: Emin_index, ii
   
   Emin_index = minloc(energy_list((ndiis - ndiist +1):ndiis),1)   
   do ii = 1, Emin_index -1
      Emat(ii,ii) = Emat(ii,ii) * DIIS_bias
   enddo
   do ii = Emin_index +1, ndiist
      Emat(ii,ii) = Emat(ii,ii) * DIIS_bias
   enddo

end subroutine diis_emat_bias


subroutine diis_get_new_fock(fock, EMAT, ndiist, M_in, spin)
   use converger_data, only: ndiis, bcoef, fockm
   implicit none
   integer     , intent(in)  :: ndiist, M_in, spin
   real(kind=8), intent(in)  :: EMAT(:,:)
   real(kind=8), intent(out) :: fock(:,:)

   integer :: ii, jj, kk, kknew, LWORK, INFO
   real(kind=8), allocatable :: work(:)


   do ii = 1, ndiist
      bcoef(ii,spin) = 0.0d0
   enddo
   bcoef(ndiist+1,spin) = -1.0d0

! First call to DGELS sets optimal WORK size. Second call solves the
! A*X = B (EMAT * Ci = bCoef) problem, with bCoef also storing the
! result.
   allocate(work(1000))
   LWORK = -1
   CALL DGELS( 'No transpose',ndiist+1, ndiist+1, 1, EMAT, &
               ndiist+1, bcoef(:,spin), ndiist+1, WORK, LWORK, INFO )

   LWORK = MIN( 1000, INT( WORK( 1 ) ) )
   CALL DGELS( 'No transpose',ndiist+1, ndiist+1, 1, EMAT, &
               ndiist+1, bcoef(:,spin), ndiist+1, WORK, LWORK, INFO )
   deallocate (work)

   ! Build new Fock as a linear combination of previous steps.
   fock = 0.0D0
   do kk = 1, ndiist
      kknew = kk + (ndiis - ndiist)
      do ii = 1, M_in
      do jj = 1, M_in
         fock(ii,jj) = fock(ii,jj) + bcoef(kk,spin) * &
                                    fockm(ii,jj,kknew,spin)
      enddo
      enddo
   enddo

end subroutine diis_get_new_fock