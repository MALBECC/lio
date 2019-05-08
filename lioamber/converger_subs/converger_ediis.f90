subroutine ediis_init(M_in, open_shell)
   use converger_data, only: nediis, ediis_fock, ediis_dens, BMAT, EDIIS_E,&
                             conver_criter

   implicit none
   logical, intent(in) :: open_shell
   integer, intent(in) :: M_in

   if (conver_criter < 6) return
   if (.not. allocated(EDIIS_E)) allocate(EDIIS_E(nediis))
   if (.not. open_shell) then
      if (.not. allocated(ediis_fock)) allocate(ediis_fock(M_in,M_in,nediis,1))
      if (.not. allocated(ediis_dens)) allocate(ediis_dens(M_in,M_in,nediis,1))
      if (.not. allocated(BMAT)      ) allocate(BMAT(nediis,nediis,1))
   else
      if (.not. allocated(ediis_fock)) allocate(ediis_fock(M_in,M_in,nediis,2))
      if (.not. allocated(ediis_dens)) allocate(ediis_dens(M_in,M_in,nediis,2))
      if (.not. allocated(BMAT)      ) allocate(BMAT(nediis,nediis,2))
   endif

   BMAT       = 0.0d0
   ediis_fock = 0.0d0
   ediis_dens = 0.0d0
   EDIIS_E    = 0.0d0
end subroutine ediis_init

subroutine ediis_finalise()
   use converger_data, only: ediis_fock, ediis_dens, BMAT, EDIIS_E, &
                             conver_criter
   implicit none

   if (conver_criter < 6) return
   if (allocated(EDIIS_E)   ) deallocate(EDIIS_E)
   if (allocated(ediis_fock)) deallocate(ediis_fock)
   if (allocated(ediis_dens)) deallocate(ediis_dens)
   if (allocated(BMAT)      ) deallocate(BMAT)

end subroutine ediis_finalise

subroutine ediis_update_energy_fock_rho(energy, fock_op, dens_op, position, &
                                        spin, niter)
   use converger_data  , only: EDIIS_E, ediis_fock, ediis_dens, nediis
   use typedef_operator, only: operator
   
   integer       , intent(in) :: position, spin, niter
   real(kind=8)  , intent(in) :: energy
   type(operator), intent(in) :: fock_op, dens_op
   integer :: ii

   if (niter > nediis) then
      do ii = 1, nediis-1
         ediis_fock(:,:,ii,spin) = ediis_fock(:,:,ii+1,spin)
         ediis_dens(:,:,ii,spin) = ediis_dens(:,:,ii+1,spin)
         EDIIS_E(ii) = EDIIS_E(ii+1)
      enddo
   endif

   call fock_op%Gets_data_ON(ediis_fock(:,:,position,spin))
   call dens_op%Gets_data_ON(ediis_dens(:,:,position,spin))
   EDIIS_E(position) = energy
endsubroutine ediis_update_energy_fock_rho

subroutine ediis_update_bmat(BMAT_aux, position, niter, M_in, spin)
   use converger_data  , only: nediis, ediis_fock, ediis_dens, BMAT
   use typedef_operator, only: operator

   implicit none
   integer     , intent(in)    :: niter, position, spin, M_in
   real(kind=8), intent(inout) :: BMAT_aux(:,:)

   integer                   :: ii, jj
   real(kind=8)              :: trace
   real(kind=8), allocatable :: mat_aux1(:,:), mat_aux2(:,:)

   allocate(mat_aux1(M_in,M_in),mat_aux2(M_in,M_in))
   mat_aux1   = 0.0d0
   mat_aux2   = 0.0d0

   ! Updating BMAT
   BMAT_aux = 0.0D0
   if (niter > nediis) then
      do jj = 1, position-1
      do ii = 1, position-1
         BMAT_aux(ii,jj) = BMAT(ii+1,jj+1,spin)
      enddo
      enddo
   else if (niter > 1) then
      do jj = 1, position-1
      do ii = 1, position-1
         BMAT_aux(ii,jj) = BMAT(ii,jj,spin)
      enddo
      enddo
   endif

   do jj = 1, position-1
      mat_aux1(:,:) = ediis_fock(:,:,jj,spin) - &
                      ediis_fock(:,:,position,spin)
      mat_aux2(:,:) = ediis_dens(:,:,jj,spin)  - &
                      ediis_dens(:,:,position,spin)

      call matmul_trace(mat_aux1(:,:), mat_aux2(:,:), M_in, trace)
      BMAT_aux(jj,position) = trace
      BMAT_aux(position,jj) = BMAT_aux(jj,position)
   enddo

   BMAT(1:position,1:position,spin) = BMAT_aux(1:position,1:position)

   deallocate(mat_aux1, mat_aux2)

end subroutine ediis_update_bmat

subroutine ediis_get_new_fock(fock, BMAT_aux, position, spin)
   use converger_data, only: EDIIS_E, ediis_fock

   implicit none
   integer     , intent(in)  :: position, spin
   real(kind=8), intent(in)  :: BMAT_aux(:,:)
   real(kind=8), intent(out) :: fock(:,:)

   integer :: ii
   real(kind=8), allocatable :: EDIIS_coef(:)


   allocate(EDIIS_coef(position))
   do ii = 1, position
      EDIIS_coef(ii) = EDIIS_E(ii)
   enddo

   ! Solving linear equation and getting new Fock matrix:
   fock = 0.0D0
   call solve_linear_constraints(EDIIS_coef(:), EDIIS_E(1:position), &
                                 BMAT_aux(:,:), position)
   do ii = 1, position
      fock(:,:) = fock(:,:) + EDIIS_coef(ii) * ediis_fock(:,:,ii,spin)
   enddo

   deallocate(EDIIS_coef)
end subroutine ediis_get_new_fock

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
   integer      :: ii, conv_steps, yind, zind(ndim-1), lastindex, newindex
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

   ! Initial guess
   coef = 0.0D0
   coef(ndim) = 1.0D0
   do ii = ndim, 2, -1
      coef(ii)   = coef(ii) / 2.0D0
      coef(ii-1) = coef(ii)
   enddo

   ! Main algorithm.
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
      if (abs(sum(coef) - 1.0D0) > 1.0D-8) then
         print*,"EDIIS restriction not complied."
         stop
      endif
   enddo

   if (coef(ndim) < 1e-37) then
      coef = 0.0D0
      coef(ndim) = 1.0D0
      do ii = ndim, 2, -1
         coef(ii)   = coef(ii) / 2.0D0
         coef(ii-1) = coef(ii)
      enddo
   endif
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
         grad(ii) = grad(ii) - BMAT(ii,jj) * coef(jj)
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