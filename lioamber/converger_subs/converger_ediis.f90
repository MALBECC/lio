subroutine EDIIS_init(M_in, OP_shell)

   use converger_data, only: nediis, damp_ediis, fock_ediis_mat, rho_ediis_mat,&
                             BMAT, EDIIS_E, fock_damp

   implicit none
   logical, intent(in) :: OP_shell
   integer, intent(in) :: M_in

   if (OP_shell) then
      allocate(fock_ediis_mat(M_in,M_in,nediis,2),                             &
               rho_ediis_mat(M_in,M_in,nediis,2), BMAT(nediis,nediis,2),       &
               EDIIS_E(nediis), fock_damp(M_in,M_in,2))
   else
      allocate(fock_ediis_mat(M_in,M_in,nediis,1),                             &
               rho_ediis_mat(M_in,M_in,nediis,1), BMAT(nediis,nediis,1),       &
               EDIIS_E(nediis), fock_damp(M_in,M_in,1))
   endif

   BMAT           = 0.0d0
   fock_ediis_mat = 0.0d0
   rho_ediis_mat  = 0.0d0
   EDIIS_E        = 0.0d0
   fock_damp      = 0.0d0

end subroutine EDIIS_init

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#ifdef CUBLAS
subroutine ediis_conver(niter, M_in, energy, devPtrX, devPTrY, open_shell, &
                        fock_aop, dens_aop, fock_bop, dens_bop)
#else
subroutine ediis_conver(niter, M_in, energy ,Xmat, Ymat, open_shell, &
                        fock_aop, dens_aop, fock_bop, dens_bop)
#endif

   use converger_data  , only: nediis, damp_ediis, fock_ediis_mat,        &
                               rho_ediis_mat, BMAT, step_nediis, EDIIS_E, &
                               fock_damp
   use typedef_operator, only: operator

   implicit none

   type(operator), intent(inout)           :: fock_aop, dens_aop
   type(operator), intent(inout), optional :: fock_bop, dens_bop

   logical     , intent(in) :: open_shell
   integer     , intent(in) :: niter
   integer     , intent(in) :: M_in
   real(kind=8), intent(in)  :: energy

   logical                   :: acceptable
   integer                   :: position, ii, jj, kk
   real(kind=8)              :: trace
   real(kind=8)              :: norm_factor
   real(kind=8), allocatable :: mat_aux1(:,:,:), mat_aux2(:,:,:)
   real(kind=8), allocatable :: fock_new(:,:,:)
   real(kind=8), allocatable :: BMAT_aux(:,:,:)
   real(kind=8), allocatable :: EDIIS_coef(:,:)

#ifdef  CUBLAS
   integer*8   , intent(in) :: devPtrX
   integer*8   , intent(in) :: devPtrY
#else
   real(kind=8), intent(in) :: Xmat(M_in,M_in)
   real(kind=8), intent(in) :: Ymat(M_in,M_in)
#endif

   step_nediis = niter
   position    = min(step_nediis, nediis)

   if (allocated(mat_aux1)  ) deallocate(mat_aux1)
   if (allocated(mat_aux2)  ) deallocate(mat_aux2)
   if (allocated(EDIIS_coef)) deallocate(EDIIS_coef)
   if (allocated(fock_new)  ) deallocate(fock_new)


   if (open_shell) then
      allocate(mat_aux1(M_in,M_in,2),mat_aux2(M_in,M_in,2), &
               EDIIS_coef(position,2), fock_new(M_in,M_in,2))
   else
      allocate(mat_aux1(M_in,M_in,1),mat_aux2(M_in,M_in,1), &
               EDIIS_coef(position,1),fock_new(M_in,M_in,1))
   endif

   mat_aux1   = 0.0d0
   mat_aux2   = 0.0d0
   EDIIS_coef = 0.0d0
   fock_new   = 0.0d0
   acceptable = .true.


   ! Updating Fock, rho and the energy:
   if (step_nediis > nediis) then
      do ii = 1, nediis
         fock_ediis_mat(:,:,ii,1) = fock_ediis_mat(:,:,ii+1,1)
         rho_ediis_mat(:,:,ii,1)  = rho_ediis_mat(:,:,ii+1,1)
         if (open_shell) then
            fock_ediis_mat(:,:,ii,2) = fock_ediis_mat(:,:,ii+1,2)
            rho_ediis_mat(:,:,ii,2)  = rho_ediis_mat(:,:,ii+1,2)
         endif
         EDIIS_E(ii) = EDIIS_E(ii+1)
      enddo
   endif

#ifdef CUBLAS
    call dens_aop%BChange_AOtoON(devPtrY, M_in, 'r')
    call fock_aop%BChange_AOtoON(devPtrX, M_in, 'r')
    if(open_shell)then
      call dens_bop%BChange_AOtoON(devPtrY, M_in, 'r')
      call fock_bop%BChange_AOtoON(devPtrX, M_in, 'r')
    endif
#else
    call dens_aop%BChange_AOtoON(Ymat, M_in, 'r')
    call fock_aop%BChange_AOtoON(Xmat, M_in, 'r')
    if (open_shell) then
       call dens_bop%BChange_AOtoON(Ymat, M_in, 'r')
       call fock_bop%BChange_AOtoON(Xmat, M_in, 'r')
    endif
#endif

   call fock_aop%Gets_data_ON(mat_aux1(:,:,1))
   call dens_aop%Gets_data_ON(mat_aux2(:,:,1))

   fock_ediis_mat(:,:,position,1) = mat_aux1(:,:,1)
   rho_ediis_mat(:,:,position,1)  = mat_aux2(:,:,1)

   if (open_shell) then
      call fock_bop%Gets_data_ON(mat_aux1(:,:,2))
      call dens_bop%Gets_data_ON(mat_aux2(:,:,2))

      fock_ediis_mat(:,:,position,1) = mat_aux1(:,:,2)
      rho_ediis_mat(:,:,position,1)  = mat_aux2(:,:,2)
   endif


   EDIIS_E(position) = energy


   ! First and second steps
   if (step_nediis == 1) then
      fock_new  = mat_aux1
      fock_damp = fock_new
   else if (step_nediis == 2) then
      fock_new  = (mat_aux1 + 80.0d0 * fock_damp) / 81.0d0
      fock_damp = fock_new
   endif

   ! Updating BMAT
   if (allocated(BMAT_aux)) deallocate(BMAT_aux)
   if (open_shell) then
      allocate(BMAT_aux(position,position,2))
   else
      allocate(BMAT_aux(position,position,1))
   endif

   BMAT_aux = 0.0d0
   if ((step_nediis > 1) .and. (step_nediis <= nediis)) then
      do jj = 1, position-1
      do ii = 1, position-1
         BMAT_aux(ii,jj,1) = BMAT(ii,jj,1)
         if (open_shell) BMAT_aux(ii,jj,2) = BMAT(ii,jj,2)
      enddo
      enddo

   else if (step_nediis > nediis) then
      do jj = 1, position-1
      do ii = 1, position-1
         BMAT_aux(ii,jj,1) = BMAT(ii+1,jj+1,1)
         if (open_shell) BMAT_aux(ii,jj,2) = BMAT(ii+1,jj+1,2)
      enddo
      enddo
   endif

   mat_aux1 = 0.0d0
   mat_aux2 = 0.0d0

   do jj = 1, position-1
      mat_aux1(:,:,1) = fock_ediis_mat(:,:,jj,1) - &
                        fock_ediis_mat(:,:,position,1)
      mat_aux2(:,:,1) = rho_ediis_mat(:,:,jj,1)  - &
                        rho_ediis_mat(:,:,position,1)

      call matmul_trace(mat_aux1(:,:,1), mat_aux2(:,:,1), M_in, trace)
      BMAT_aux(jj,position,1) = trace

      if (open_shell) then
         mat_aux1(:,:,2) = fock_ediis_mat(:,:,jj,2) - &
                           fock_ediis_mat(:,:,position,2)
         mat_aux2(:,:,2) = rho_ediis_mat(:,:,jj,2)  - &
                           rho_ediis_mat(:,:,position, 2)
         call matmul_trace(mat_aux1(:,:,2), mat_aux2(:,:,2), M_in, trace)
         BMAT_aux(jj,position,2) = trace
      endif
      BMAT_aux(position, jj,:) = BMAT_aux(jj, position,:)
   enddo


   BMAT(1:position,1:position,:) = BMAT_aux(1:position,1:position,:)
   do ii = 1, position
      EDIIS_coef(ii,:) = EDIIS_E(ii)
   enddo


   ! Solving linear equation and getting new Fock matrix:
   if (step_nediis > 2) then
      call solve_linear_constraints(EDIIS_coef(:,1), EDIIS_E(1:position), &
                                    BMAT_aux(:,:,1), position)
      if (open_shell) call solve_linear_constraints(EDIIS_coef(:,2),      &
                                              EDIIS_E(1:position),        &
                                              BMAT_aux(:,:,2), position)

      do ii = 1, position
         fock_new(:,:,1) = fock_new(:,:,1) + EDIIS_coef(ii,1) * &
                           fock_ediis_mat(:,:,ii,1)
         if (open_shell) fock_new(:,:,2) = fock_new(:,:,2) + EDIIS_coef(ii,2) *&
                                           fock_ediis_mat(:,:,ii,2)
      enddo
   endif

   call fock_aop%Sets_data_ON(fock_new(:,:,1))
   if (open_shell) call fock_bop%Sets_data_ON(fock_new(:,:,2))
end subroutine ediis_conver

subroutine matmul_trace(mat1,mat2, M_in, trace)
   implicit none
   integer     , intent(in)  :: M_in
   real(kind=8), intent(in)  :: mat1(M_in,M_in), mat2(M_in, M_in)
   real(kind=8), intent(out) :: trace
   integer      :: ii, jj
   real(kind=8) :: mat3(M_in)

   mat3  = 0.0d0
   trace = 0.0d0

   do ii=1, M_in
   do jj=1, M_in
      mat3(ii) = mat1(ii,jj) * mat2(jj,ii) + mat3(ii)
   enddo
   enddo

   do ii=1, M_in
      trace = trace + mat3(ii)
   enddo
end subroutine matmul_trace

subroutine solve_linear_constraints(coef, Ener, BMAT, ndim)

   implicit none
   integer     , intent(in)  :: ndim
   real(kind=8), intent(in)  :: Ener(ndim), BMAT(ndim,ndim)
   real(kind=8), intent(out) :: coef(ndim)
   logical      :: converged, big_alpha1, big_alpha2, update
   integer      :: ii, jj, conv_steps, yind, zind(ndim-1), lastindex, newindex
   real(kind=8) :: aux_coef(ndim), new_coef(ndim), grad(ndim), delta(ndim), &
                   r_grad(ndim-1), alpha1, alpha2, alpha3, alpha_aux,       &
                   vec_alpha2(ndim-1), result1, result2, result3

   coef       = 1.0d0 / dble(ndim)
   new_coef   = 0.0d0
   delta      = 0.0d0
   aux_coef   = 0.0d0
   r_grad     = 0.0d0
   lastindex  = 0
   newindex   = 1
   yind       = 1
   alpha1     = 0.0d0
   alpha2     = 0.0d0
   alpha3     = 0.0d0
   alpha_aux  = 0.0d0
   converged  = .true.
   conv_steps = 0

   do ii = 2, ndim
      zind(ii-1) = ii
   enddo

   do while (converged .and. (conv_steps <= 100000))
      conv_steps = conv_steps +1

      big_alpha1 = .false.
      big_alpha2 = .false.
      update     = .false.
      converged  = .false.

      call gradient(coef, grad, Ener, BMAT ,ndim)
      call displacement(grad, zind, yind, delta, coef, ndim)

      do ii = 1, ndim-1
         if (abs(delta(zind(ii))) > 1.0D-8) then
            converged = .true.
            exit
         endif
      enddo

      if (converged .eqv. .false.) exit

      if (delta(yind) < 0.0d0) then
         alpha1 = - coef(yind) / delta(yind)
      else
         big_alpha1 = .true.
      endif

      do ii = 1, ndim-1
         vec_alpha2(ii) = -delta(zind(ii)) / coef(zind(ii))
      enddo

      alpha2 = maxval(vec_alpha2)
      alpha2 = 1.0d0 / alpha2
      if (alpha2 <= 0.0d0) big_alpha2 = .true.

      call min_alpha(Ener,BMAT, coef,delta, alpha3, ndim)

      if (big_alpha1 .and. big_alpha2) then
         call f_coef(Ener,BMAT, coef+alpha3*delta, result1, ndim)
         call f_coef(Ener,BMAT, coef, result2, ndim)

      else if (big_alpha1) then
         if (alpha3 > alpha2) alpha3 = alpha2
         call f_coef(Ener, BMAT, coef + alpha2 * delta, result1, ndim)
         call f_coef(Ener, BMAT, coef + alpha3 * delta, result2, ndim)
         call f_coef(Ener, BMAT, coef                 , result3, ndim)
         if (result1 < result2) alpha3 = alpha2

      else if (big_alpha2) then
         if (alpha3 > alpha1) then
            alpha3 = alpha1
            update = .true.
         endif
         call f_coef(Ener, BMAT, coef + alpha1 * delta, result1, ndim)
         call f_coef(Ener, BMAT, coef + alpha3 * delta, result2, ndim)
         call f_coef(Ener, BMAT, coef                 , result3, ndim)
         if (result1 < result2) then
            alpha3 = alpha1
            update = .true.
         endif
      else
         alpha_aux = alpha2
         if (alpha1 < alpha2)    alpha_aux = alpha1
         if (alpha3 > alpha_aux) alpha3    = alpha_aux

         call f_coef(Ener,BMAT, coef + alpha_aux * delta, result1, ndim)
         call f_coef(Ener,BMAT, coef + alpha3    * delta, result2, ndim)
         call f_coef(Ener,BMAT, coef                    , result3, ndim)
         if (result1 < result2) alpha3 = alpha_aux
         if (alpha3  >= alpha1) update = .true.
      endif

      new_coef = coef + alpha3 * delta
      coef     = new_coef

      if (update) then
         newindex = newindex+1
         if (newindex == ndim+1) newindex = 2

         lastindex        = yind
         yind             = zind(newindex-1)
         zind(newindex-1) = lastindex
      endif

      ! Checking the restriction.
      result1 = 0.0d0
      do ii = 1, ndim
         result1 = coef(ii) + result1
      enddo
      if (abs(result1-1.0d0) > 1.0d-8) then
         print*,"EDIIS restriction not complied."
         stop
      endif
   enddo
end subroutine solve_linear_constraints

subroutine gradient(coef, grad, Ener, BMAT ,ndim)

   implicit none
   integer     , intent(in)  :: ndim
   real(kind=8), intent(in)  :: coef(ndim), Ener(ndim), BMAT(ndim,ndim)
   real(kind=8), intent(out) :: grad(ndim)
   integer :: ii, jj

   grad = 0.0d0
   do ii = 1, ndim
      do jj = 1, ndim
         grad(ii) = - BMAT(ii,jj) * coef(jj) + grad(ii)
      enddo
      grad(ii) = grad(ii) + Ener(ii)
   enddo

end subroutine gradient

subroutine displacement(grad, zind, yind, delta, coef, ndim)
  implicit none
   integer     , intent(in)  :: ndim, yind, zind(ndim-1)
   real(kind=8), intent(in)  :: coef(ndim), grad(ndim)
   real(kind=8), intent(out) :: delta(ndim)
   integer      :: ii
   real(kind=8) :: r_grad(ndim-1)

   delta = 0.0d0
   do ii = 1, ndim-1
      r_grad(ii) = grad(zind(ii)) - grad(yind)
      delta(zind(ii)) = 0.0d0
      if ((r_grad(ii) < 0.0d0) .or. (coef(zind(ii)) > 0.0d0)) &
                                  delta(zind(ii)) = - r_grad(ii)
   enddo

   do ii = 1, ndim-1
       delta(yind) = - delta(zind(ii)) + delta(yind)
   enddo
end subroutine displacement

subroutine min_alpha(Ener, BMAT, coef,delta, alpha, ndim)
   implicit none
   integer     , intent(in)  :: ndim
   real(kind=8), intent(in)  :: Ener(ndim), BMAT(ndim,ndim), coef(ndim), &
                                delta(ndim)
   real(kind=8), intent(out) :: alpha
   integer      :: ii, jj
   real(kind=8) :: num1, num2, den1

   num1  = 0.0d0
   num2  = 0.0d0
   den1  = 0.0d0
   alpha = 0.0d0

   do ii = 1, ndim
      num1 = num1 + Ener(ii) * delta(ii)
   enddo

   do jj = 1, ndim
   do ii = 1, ndim
      num2 = num2 + BMAT(ii,jj) * delta(ii) * coef(jj)
      den1 = den1 + BMAT(ii,jj) * delta(ii) * delta(jj)
   enddo
   enddo
   alpha = (num1 - num2) / den1
end subroutine min_alpha

subroutine f_coef(Ener, BMAT, coef, result, ndim)

   implicit none
   integer     , intent(in)  :: ndim
   real(kind=8), intent(in)  :: Ener(ndim), BMAT(ndim,ndim), coef(ndim)
   real(kind=8), intent(out) :: result
   integer      :: ii, jj
   real(kind=8) :: sum1, sum2

   sum1   = 0.0d0
   sum2   = 0.0d0
   result = 0.0d0

   do ii = 1, ndim
      sum1 = sum1 + Ener(ii) * coef(ii)
   enddo

   do jj = 1, ndim
   do ii = 1, ndim
      sum2 = sum2 + BMAT(ii,jj) * coef(ii) * coef(jj)
   enddo
   enddo

   result = sum1 - 0.5d0 * sum2
end subroutine