! Performs the CDFT iterative procedure.
! Variables Pmat_v, coefs, and coefs_b get allocated elsewhere and overwritten
! by SCF in each iteration.
subroutine CDFT(fock_a, rho_a, fock_b, rho_b, Pmat_v, coefs, coefs_b, natom, &
                op_shell)
   use typedef_operator, only: operator
   use converger_data  , only: told
   use cdft_data       , only: cdft_c

   implicit none
   integer, intent(in)                 :: natom
   logical, intent(in)                 :: op_shell
   LIODBLE, intent(inout)              :: Pmat_v(:), coefs(:,:), coefs_b(:,:)
   type(operator), intent(inout)       :: fock_a, rho_a, fock_b, rho_b

   integer :: cdft_iter, max_cdft_iter
   logical :: cdft_converged = .false.
   LIODBLE :: energ
   LIODBLE, allocatable :: Pmat_old(:)


   max_cdft_iter = 100
   cdft_iter     = 0
   allocate(Pmat_old(size(Pmat_v,1)))

   call cdft_initialise(natom)
   if (cdft_c%mixed) call cdft_mixed_initialise(size(Pmat_v,1), op_shell)

   do while ((.not. cdft_converged) .and. (cdft_iter < max_cdft_iter))
      cdft_iter = cdft_iter +1
      Pmat_old  = Pmat_v
      call SCF(energ, fock_a, rho_a, fock_b, rho_b)
      call cdft_check_conver(Pmat_v, Pmat_old, cdft_converged, &
                             cdft_iter, energ, told)

      if (.not. cdft_converged) then
         ! Calculates perturbations and Jacobian.
         call cdft_get_deltaV(fock_a, rho_a, fock_b, rho_b)
         call cdft_set_potential()
      elseif (cdft_c%mixed) then
         ! Gets Wmat and retrieves MO
         call cdft_mixed_set_coefs(coefs, .true., 1)
         if (op_shell) call cdft_mixed_set_coefs(coefs_b, .false., 1)
      endif
   enddo

   ! If it is a mixed calculation, repeats everything for the second state.
   if (cdft_c%mixed) then
      cdft_converged = .false.
      cdft_iter      = 0
      call cdft_mixed_switch()

      do while ((.not. cdft_converged) .and. (cdft_iter < max_cdft_iter))
         cdft_iter = cdft_iter +1
         Pmat_old  = Pmat_v
         call SCF(energ, fock_a, rho_a, fock_b, rho_b)
         call cdft_check_conver(Pmat_v, Pmat_old, cdft_converged, &
                                cdft_iter, energ, told)
   
         if (.not. cdft_converged) then
            ! Calculates perturbations and Jacobian.
            call cdft_get_deltaV(fock_a, rho_a, fock_b, rho_b)
            call cdft_set_potential()
         else
            ! Gets Wmat and retrieves MO
            call cdft_mixed_set_coefs(coefs, .true., 2)
            if (op_shell) call cdft_mixed_set_coefs(coefs_b, .false., 2)
         endif
      enddo
   endif

   call cdft_finalise()
   call cdft_mixed_finalise()
   deallocate(Pmat_old)
end subroutine CDFT