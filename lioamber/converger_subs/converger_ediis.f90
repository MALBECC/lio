!#############################################################################!
! These routines are intended for EDIIS and ADIIS convergence. However,       !
! they do not currently work as intended. Lose time here at your own risk.    !
!                                                                             !
! Externally called subroutines (from converger_commons):                     !
!   * ediis_init                                                              !
!   * ediis_finalise                                                          !
!   * ediis_update_energy_fock_rho                                            !
!   * ediis_update_bmat                                                       !
!   * ediis_get_new_fock                                                      !
!                                                                             !
! Internal subroutines:                                                       !
!   * matmul_trace                                                            !
!   * get_coefs_dfp                                                           !
!   * get_coefs_charly                                                        !
!   * f_coef                                                                  !
!   * min_alpha                                                               !
!   * gradient                                                                !
!   * displacement                                                            !
!                                                                             !
! Useful references are:                                                      !
!   * Kudin et al., JCP 116 (2002), DOI:10.1063/1.1470195 (for EDIIS)         !
!   * Hu et al.   , JCP 132 (2010), DOI:10.1063/1.3304922 (for ADIIS)         !
!                                                                             ! 
! Last one losing time here: Federico Pedron, May/2019                        !
!#############################################################################!
subroutine ediis_init(M_in, open_shell)
   use converger_data, only: nediis, ediis_fock, ediis_dens, BMAT, EDIIS_E,&
                             conver_method, EDIIS_coef

   implicit none
   logical, intent(in) :: open_shell
   integer, intent(in) :: M_in

   if (conver_method < 6) return
   if (.not. allocated(BMAT))       allocate(BMAT(nediis,nediis))
   if (.not. allocated(EDIIS_E))    allocate(EDIIS_E(nediis))
   if (.not. allocated(EDIIS_coef)) allocate(EDIIS_coef(nediis))
   if (.not. open_shell) then
      if (.not. allocated(ediis_fock)) allocate(ediis_fock(M_in,M_in,nediis,1))
      if (.not. allocated(ediis_dens)) allocate(ediis_dens(M_in,M_in,nediis,1))
   else
      if (.not. allocated(ediis_fock)) allocate(ediis_fock(M_in,M_in,nediis,2))
      if (.not. allocated(ediis_dens)) allocate(ediis_dens(M_in,M_in,nediis,2))
   endif

   BMAT       = 0.0d0
   ediis_fock = 0.0d0
   ediis_dens = 0.0d0
   EDIIS_E    = 0.0d0
   EDIIS_coef = 0.0d0
end subroutine ediis_init

subroutine ediis_finalise()
   use converger_data, only: ediis_fock, ediis_dens, BMAT, EDIIS_E, &
                             conver_method, EDIIS_coef
   implicit none

   if (conver_method < 6) return
   if (allocated(EDIIS_E)   ) deallocate(EDIIS_E)
   if (allocated(EDIIS_coef)) deallocate(EDIIS_coef)
   if (allocated(ediis_fock)) deallocate(ediis_fock)
   if (allocated(ediis_dens)) deallocate(ediis_dens)
   if (allocated(BMAT)      ) deallocate(BMAT)

end subroutine ediis_finalise

subroutine ediis_update_energy_fock_rho(fock_op, dens_op, position, spin, &
                                        niter)
   use converger_data  , only: ediis_fock, ediis_dens, nediis
   use typedef_operator, only: operator
   
   integer       , intent(in) :: position, spin, niter
   
   type(operator), intent(in) :: fock_op, dens_op
   integer :: ii

   if (niter > nediis) then
      do ii = 1, nediis-1
         ediis_fock(:,:,ii,spin) = ediis_fock(:,:,ii+1,spin)
         ediis_dens(:,:,ii,spin) = ediis_dens(:,:,ii+1,spin)
      enddo
   endif

   call fock_op%Gets_data_ON(ediis_fock(:,:,position,spin))
   call dens_op%Gets_data_ON(ediis_dens(:,:,position,spin))
endsubroutine ediis_update_energy_fock_rho

subroutine ediis_update_bmat(energy, nlast, niter, M_in, open_shell)
   use converger_data  , only: nediis, ediis_fock, ediis_dens, BMAT, EDIIS_E, &
                               EDIIS_not_ADIIS
   use typedef_operator, only: operator

   implicit none
   integer     , intent(in)  :: niter, nlast, M_in
   logical     , intent(in)  :: open_shell
   LIODBLE, intent(in)  :: energy

   integer                   :: ii, jj
   LIODBLE              :: trace1, trace2, trace3, trace4
   trace1 = 0.0D0 ; trace2 = 0.0D0
   trace3 = 0.0D0 ; trace4 = 0.0D0 

   if (EDIIS_not_ADIIS) then
      ! Calculates EDIIS components.

      ! Updating BMAT if there are enough iterations.
      if (niter > nediis) then
         do jj = 1, nlast-1
            EDIIS_E(jj) = EDIIS_E(jj+1)
            do ii = 1, nlast-1
               BMAT(ii,jj) = BMAT(ii+1,jj+1)
            enddo
         enddo
      endif

      EDIIS_E(nlast) = energy
      do ii = 1, nlast-1
         call matmul_trace(ediis_dens(:,:,ii,1) - ediis_dens(:,:,nlast,1), &
                           ediis_fock(:,:,ii,1) - ediis_fock(:,:,nlast,1), &
                           M_in, trace1)
                           
         if (open_shell) then
            call matmul_trace(ediis_dens(:,:,ii,2) - ediis_dens(:,:,nlast,2),&
                              ediis_fock(:,:,ii,2) - ediis_fock(:,:,nlast,2),&
                              M_in, trace2)
         endif

         BMAT(ii,nlast) = trace1 + trace2
         BMAT(nlast,ii) = BMAT(ii,nlast)
      enddo   
   else
      ! Calculates ADIIS components.
      EDIIS_E  = 0.0D0
      do ii = 1, nlast-1
         call matmul_trace(ediis_dens(:,:,ii,1) - ediis_dens(:,:,nlast,1), &
                           ediis_fock(:,:,nlast,1), M_in, trace1)

         if (open_shell) then
            call matmul_trace(ediis_dens(:,:,ii,2) - ediis_dens(:,:,nlast,2), &
                              ediis_fock(:,:,nlast,2), M_in, trace2)
         endif
         EDIIS_E(ii) = trace2 + trace1
         do jj = 1, nlast-1
            call matmul_trace(ediis_dens(:,:,ii,1) - ediis_dens(:,:,nlast,1), &
                              ediis_fock(:,:,jj,1) - ediis_fock(:,:,nlast,1), &
                              M_in, trace3)
            if (open_shell) then
               call matmul_trace(ediis_dens(:,:,ii,2)-ediis_dens(:,:,nlast,2), &
                                 ediis_fock(:,:,jj,2)-ediis_fock(:,:,nlast,2), &
                                 M_in, trace4)
            endif
            BMAT(ii,jj) = trace1 + trace2 + trace3 + trace4
            BMAT(jj,ii) = BMAT(ii,jj)
         enddo
      enddo
   endif

end subroutine ediis_update_bmat

subroutine ediis_get_new_fock(fock, position, spin)
   use converger_data, only: EDIIS_E, ediis_fock, BMAT, EDIIS_coef

   implicit none
   integer     , intent(in)  :: position, spin
   LIODBLE, intent(out) :: fock(:,:)

   integer :: ii, solver

   ! Solves linear equation and gets new Fock matrix.
   if (spin == 1) then
      solver = 0
      if (solver == 1) then
         call get_coefs_dfp(EDIIS_coef(1:position), EDIIS_E(1:position), &
                           BMAT(1:position,1:position), position)
      else
         call get_coefs_charly(EDIIS_coef(1:position), EDIIS_E(1:position), &
                              BMAT(1:position,1:position), position)
      endif
   endif

   fock = 0.0D0
   do ii = 1, position
      fock(:,:) = fock(:,:) + EDIIS_coef(ii) * ediis_fock(:,:,ii,spin)
   enddo
end subroutine ediis_get_new_fock

!#############################################################################!
subroutine matmul_trace(mat1, mat2, M_in, trace)
   implicit none
   integer     , intent(in)  :: M_in
   LIODBLE, intent(in)  :: mat1(M_in,M_in), mat2(M_in, M_in)
   LIODBLE, intent(out) :: trace
   integer      :: ii, jj
   LIODBLE :: mat3(M_in)

   mat3  = 0.0D0
   trace = 0.0D0

   do ii=1, M_in
   do jj=1, M_in
      mat3(ii) = mat1(ii,jj) * mat2(jj,ii) + mat3(ii)
   enddo
   enddo

   do ii=1, M_in
      trace = trace + mat3(ii)
   enddo
end subroutine matmul_trace

subroutine get_coefs_dfp(coef, Ener, BMAT, ndim)
   implicit none
   integer     , intent(in)    :: ndim
   LIODBLE, intent(in)    :: Ener(ndim), BMAT(ndim,ndim)
   LIODBLE, intent(inout) :: coef(ndim)


   integer      :: n_iter, icount, jcount
   LIODBLE :: funct_val, sum_p, step_max, temp_val, funct_ls, check_val, &
                   fac, fae, fad, sumdg, sumxi
   LIODBLE, allocatable :: hessian_inv(:,:), grad(:), coef_new(:), &
                                dgrad(:), hdgrad(:), xi(:)
     
   allocate(hessian_inv(ndim, ndim), xi(ndim), grad(ndim), coef_new(ndim), &
            dgrad(ndim), hdgrad(ndim))

   ! Initial guess
   coef = 0.0d0
   jcount = 1
   do icount = 2, ndim
      sumdg = Ener(icount) - BMAT(icount,icount)
      sumxi = Ener(jcount) - BMAT(jcount,jcount)
      if (sumdg < sumxi) jcount = icount
   enddo
   coef(jcount) = 0.9D0
   do icount = 1, jcount-1
      coef(icount) = (1.0D0 - coef(jcount)) / (ndim - 1)
   enddo
   do icount = jcount+1, ndim
      coef(icount) = (1.0D0 - coef(jcount)) / (ndim - 1)
   enddo
   coef = coef / sum(coef)
   
   call f_coef(Ener, BMAT, coef, funct_val, ndim)
   call gradient(coef, grad, Ener, BMAT, ndim)

   hessian_inv = 0.0D0
   xi = -grad
   sum_p = 0.0D0
   do icount = 1, ndim
      hessian_inv(icount,icount) = 1.0D0
      sum_p = sum_p + coef(icount) * coef(icount)
   enddo

   step_max = 100.0D0 * max(sqrt(sum_p), dble(ndim))
   do n_iter = 1, 500
      call line_search(ndim, coef, funct_val, grad, xi, coef_new, funct_ls, &
                      step_max, Ener, BMAT)

      ! Updates coefficients
      xi   = coef_new - coef
      coef = coef_new

      ! Checks convergence in coefficients.
      check_val = 0.0D0
      do icount = 1, ndim
         temp_val = abs(xi(icount)) / max(abs(coef(icount)), 1.0D0)
         if (check_val < temp_val) check_val = temp_val
      enddo
      if (check_val > 1.0D-7) exit
      
      dgrad = grad
      call gradient(coef, grad, Ener, BMAT, ndim)

      ! Checks convergence in gradient.
      check_val = 0.0D0
      do icount = 1, ndim
         temp_val = abs(grad(icount)) * max(abs(coef(icount)), 1.0D0) &
                    / max(abs(funct_ls), 1.0D0)
         if (check_val < temp_val) check_val = temp_val
      enddo
      if (check_val > 1.0D-7) exit

      dgrad  = grad - dgrad
      hdgrad = 0.0D0
      do icount = 1, ndim
      do jcount = 1, ndim
         hdgrad(icount) = hdgrad(icount) + hessian_inv(icount, jcount) &
                                         * dgrad(jcount)
      enddo
      enddo

      fac   = 0.0D0
      fae   = 0.0D0
      sumdg = 0.0D0
      sumxi = 0.0D0
      do icount = 1, ndim
         fac   = fac + dgrad(icount) * xi(icount)
         fae   = fae + dgrad(icount) * hdgrad(icount)
         sumdg = sumdg + dgrad(icount) * dgrad(icount)
         sumxi = sumxi + xi(icount)    * xi(icount)
      enddo

      ! Skips update if fac is not big enough
      if ((fac*fac) > (1D-7 * sumdg * sumxi)) then
         fac = 1.0D0 / fac
         fad = 1.0D0 / fae
         dgrad = fac * xi - fad * hdgrad
         
         do icount = 1, ndim
         do jcount = 1, ndim
            hessian_inv(icount,jcount) = hessian_inv(icount,jcount)            &
                                       + fac * xi(icount)     * xi(jcount)     &
                                       - fad * hdgrad(icount) * hdgrad(jcount) &
                                       + fae * dgrad(icount)  * dgrad(jcount)
         enddo
         enddo
      endif

      xi = 0.0D0
      do icount = 1, ndim
      do jcount = 1, ndim
         xi(icount) = xi(icount) - hessian_inv(icount, jcount) * grad(jcount)
      enddo
      enddo
   enddo
  
   if (n_iter > 500) stop "NO CONVERGENCE IN EDIIS"

   ! Standardizes coefficients so that the comply with the constrains.

   if (coef(ndim) < 1e-6) coef(ndim) = 1.0D0
   sumdg = 0.0D0
   do icount = 1, ndim
      sumdg = sumdg + coef(icount) * coef(icount)
   enddo
   do icount = 1, ndim
      coef(icount) = coef(icount) * coef(icount) / sumdg
   enddo

   deallocate(hessian_inv, xi, grad, coef_new, dgrad, hdgrad)
end subroutine get_coefs_dfp

subroutine get_coefs_charly(coef, Ener, BMAT, ndim)
   implicit none
   
   integer     , intent(in)  :: ndim
   LIODBLE, intent(in)  :: BMAT(ndim,ndim), Ener(ndim)
   LIODBLE, intent(out) :: coef(ndim)

   logical      :: converged, big_alpha1, big_alpha2, update
   integer      :: ii, conv_steps, yind, zind(ndim-1), lastindex, newindex
   LIODBLE :: aux_coef(ndim), new_coef(ndim), grad(ndim), &
                   delta(ndim), r_grad(ndim-1), alpha1, alpha2,&
                   alpha3, alpha_aux, vec_alpha2(ndim-1)
   LIODBLE :: result1, result2

   coef         = 1.0d0/dble(ndim)
   new_coef     = 0.0d0
   delta        = 0.0d0
   aux_coef     = 0.0d0
   r_grad       = 0.0d0
   lastindex    = 0
   newindex     = 1
   yind         = 1
   alpha1       = 0.0d0
   alpha2       = 0.0d0
   alpha3       = 0.0d0
   alpha_aux    = 0.0d0
   converged    = .true.
   conv_steps   = 0

   do ii=2, ndim
      zind(ii-1) = ii
   enddo

   do while (converged .and. (conv_steps <= 100000))
      conv_steps = conv_steps +1
      big_alpha1 = .false.
      big_alpha2 = .false.
      update     = .false.
      converged  = .false.

      call gradient(coef, grad, Ener, BMAT ,ndim)
      call displacement (grad, zind, yind, delta, coef, ndim)

      do ii=1, ndim-1
         if (abs(delta(zind(ii))) > 1.0D-8) then
            converged = .true.
            exit
         endif
      enddo
      if (converged .eqv. .false.) then
        exit
      endif

      if (delta(yind) < 0.0d0) then
         alpha1 = -coef(yind) / delta(yind)
      else
         big_alpha1 = .true.
      endif

      do ii=1, ndim-1
         vec_alpha2(ii) = -delta(zind(ii)) / coef(zind(ii))
      enddo
      alpha2 = maxval(vec_alpha2)
      alpha2 = 1.0d0 / alpha2

      if (alpha2 <= 0.0d0) then
         big_alpha2 = .true.
      endif

      call min_alpha(Ener, BMAT, coef,delta, alpha3, ndim)

      if (big_alpha1 .and. big_alpha2) then
      elseif (big_alpha1) then
         if (alpha3 > alpha2) alpha3 = alpha2
         
         call f_coef(Ener,BMAT, coef + alpha2 * delta, result1, ndim)
         call f_coef(Ener,BMAT, coef + alpha3 * delta, result2, ndim)

         if (result1 < result2) alpha3 = alpha2

      elseif (big_alpha2) then
         if (alpha3 > alpha1) then
            alpha3  = alpha1
            update  = .true.
         endif

         call f_coef(Ener,BMAT, coef + alpha1 * delta, result1, ndim)
         call f_coef(Ener,BMAT, coef + alpha3 * delta, result2, ndim)
         
         if (result1 < result2) then
            alpha3 = alpha1
            update = .true.
         endif

      else
         if (alpha1 < alpha2) then
            alpha_aux = alpha1
         else
            alpha_aux = alpha2
         endif
         if (alpha3 > alpha_aux) alpha3 = alpha_aux
         call f_coef(Ener,BMAT, coef + alpha_aux * delta, result1, ndim)
         call f_coef(Ener,BMAT, coef + alpha3    * delta, result2, ndim)
         if (result1 <  result2) alpha3 = alpha_aux
         if (alpha3  >= alpha1 ) update = .true.
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
   enddo

   ! Normalizes coefficients
   if (coef(ndim) < 1e-6) coef(ndim) = 1.0D0
   result1 = 0.0D0
   do ii = 1, ndim
      result1 = result1 + coef(ii) * coef(ii)
   enddo
   do ii = 1, ndim
      coef(ii) = coef(ii) * coef(ii) / result1
   enddo
end subroutine get_coefs_charly

subroutine line_search(ndim, x_old, f_old, grad, xi, x_new, f_new, step_max, &
                       Ener, BMAT)
   implicit none
   integer     , intent(in)    :: ndim
   LIODBLE, intent(in)    :: x_old(ndim), grad(ndim), f_old, step_max, &
                                  Ener(ndim), BMAT(ndim, ndim)
   LIODBLE, intent(out)   :: x_new(ndim), f_new
   LIODBLE, intent(inout) :: xi(ndim)

   integer      :: icount
   logical      :: first_step
   LIODBLE :: sum_x, slope, lambda_min, temp, lambda, lambda2, f_new2, &
                   f_old2, lambda_temp, rhs1, rhs2, temp2


   sum_x = 0.0D0
   do icount = 1, ndim
      sum_x = sum_x + xi(icount) * xi(icount)
   enddo
   sum_x = sqrt(sum_x)
   if (sum_x > step_max) xi = xi * step_max / sum_x

   slope = 0.0D0
   do icount = 1, ndim
      slope = slope + grad(icount) * xi(icount)
   enddo

   lambda_min = 0.0D0
   do icount = 1, ndim
      temp = abs(xi(icount)) / max(abs(x_old(icount)), 1.0D0)
      if (temp > lambda_min) lambda_min = temp
   enddo
   lambda_min = 1.0D-7 / lambda_min

   lambda  = 1.0D0
   lambda2 = 0.0D0
   f_new2  = 0.0D0
   f_old2  = 0.0D0
   first_step = .true.
   do
      x_new = x_old + lambda * xi
      call f_coef(Ener, BMAT, x_new, f_new, ndim)

      if (lambda < lambda_min) then
         ! Exits cycle if step is too small.
         x_new = x_old
         exit
      else if (f_new < (f_old + 1.0D-4 * lambda * slope)) then
         ! Exits cycle if function f is decreased sufficiently.
         exit
      else
         ! Runs normally
         if (first_step) then
            ! First step of the cycle
            lambda_temp = - slope / (2.0D0 * (f_new - f_old - slope))
            first_step = .false.
         else
            rhs1 = f_new  - f_old  - lambda  * slope
            rhs2 = f_new2 - f_old2 - lambda2 * slope

            temp  = (rhs1 / (lambda * lambda) - rhs2 / (lambda2 * lambda2)) / &
                    (lambda - lambda2)
            temp2 = (- lambda2 * rhs1 / (lambda  * lambda)   &
                     + lambda  * rhs2 / (lambda2 * lambda2)) &
                     / (lambda - lambda2)
            if (abs(temp) < 1E-37) then
               lambda_temp = - slope / (2.0D0 * temp2)
            else
               lambda_temp = (-temp2 + &
                              sqrt(temp2 * temp2 - 3.0D0 * temp * slope)) &
                              / (3.0D0 * temp)
            endif
            if (lambda_temp > 0.5D0 * lambda) lambda_temp = 0.5D0 * lambda
         endif
      endif

      lambda2 = lambda
      f_new2  = f_new
      f_old2  = f_old
      lambda  = max(lambda_temp, 0.1D0 * lambda)
   enddo
end subroutine line_search

subroutine gradient(coef, grad, Ener, BMAT ,ndim)
   use converger_data, only: EDIIS_not_ADIIS

   implicit none
   integer     , intent(in)  :: ndim
   LIODBLE, intent(in)  :: coef(ndim), Ener(ndim), BMAT(ndim,ndim)
   LIODBLE, intent(out) :: grad(ndim)
   integer :: ii, jj

   grad = 0.0d0
   if (EDIIS_not_ADIIS) then
      do ii = 1, ndim
         do jj = 1, ndim
            grad(ii) = grad(ii) - BMAT(ii,jj) * coef(jj)
         enddo
         grad(ii) = grad(ii) + Ener(ii)
      enddo
   else
      do ii = 1, ndim
         do jj = 1, ndim
            grad(ii) = grad(ii) + BMAT(ii,jj) * coef(jj)
         enddo
         grad(ii) = grad(ii) + Ener(ii)
      enddo
      grad = 2.0D0 * grad
   endif

end subroutine gradient

subroutine f_coef(Ener, BMAT, coef, result, ndim)
   use converger_data, only: EDIIS_not_ADIIS

   implicit none
   integer     , intent(in)  :: ndim
   LIODBLE, intent(in)  :: Ener(ndim), BMAT(ndim,ndim), coef(ndim)
   LIODBLE, intent(out) :: result
   integer      :: ii, jj
   LIODBLE :: sum1, sum2

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

   if (EDIIS_not_ADIIS) then
      result = sum1 - sum2
   else
      result = 2.0D0 * sum1 + sum2
   endif
end subroutine f_coef

subroutine displacement(grad, zind, yind, delta, coef, ndim)
   implicit none
   integer     , intent(in)  :: ndim, yind, zind(ndim-1)
   LIODBLE, intent(in)  :: coef(ndim), grad(ndim)
   LIODBLE, intent(out) :: delta(ndim)

   integer      :: ii
   LIODBLE :: r_grad(ndim-1)

   delta = 0.0D0
   do ii = 1, ndim-1
      r_grad(ii) = grad(zind(ii)) - grad(yind)
      if ((r_grad(ii) < 0.0D0) .or. (coef(zind(ii)) > 0.0D0)) then
         delta(zind(ii)) = -r_grad(ii)
      else
         delta(zind(ii)) = 0.0D0
      endif
   enddo

   do ii = 1, ndim-1
      delta(yind) = - delta(zind(ii)) + delta(yind)
   enddo
end subroutine displacement

subroutine min_alpha(Ener, BMAT, coef, delta, alpha, ndim)
   implicit none
   integer     , intent(in)   :: ndim
   LIODBLE, intent(in)   :: Ener(ndim), BMAT(ndim,ndim), &
                                 coef(ndim), delta(ndim)
   LIODBLE,  intent(out) :: alpha

   integer      :: ii, jj
   LIODBLE :: num1, num2, den1

   num1  = 0.0d0
   num2  = 0.0d0
   den1  = 0.0d0
   alpha = 0.0d0

   do ii=1, ndim
      num1 = num1 + Ener(ii) * delta(ii)
   enddo

   do jj=1, ndim
   do ii=1, ndim
      num2 = num2 + BMAT(ii,jj) * delta(ii) * coef(jj)
      den1 = den1 + BMAT(ii,jj) * delta(ii) * delta(jj)
   enddo
   enddo

   alpha = (num1 - num2)/den1
end subroutine min_alpha
!#############################################################################!