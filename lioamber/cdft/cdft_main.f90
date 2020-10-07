! Performs the CDFT iterative procedure.
! Variables Pmat_v, coefs, and coefs_b get allocated elsewhere and overwritten
! by SCF in each iteration.
subroutine CDFT(fock_a, rho_a, fock_b, rho_b, Pmat_v, coefs, coefs_b, overlap, &
                natom, nbasis, nOcc, nOcc_b, op_shell)
   use typedef_operator, only: operator
   use converger_data  , only: told
   use cdft_data       , only: cdft_c

   implicit none
   integer, intent(in)                 :: natom, nbasis, nOcc, nOcc_b
   logical, intent(in)                 :: op_shell
   LIODBLE, intent(inout)              :: Pmat_v(:), coefs(:,:), coefs_b(:,:),&
                                          overlap(:,:)
   type(operator), intent(inout)       :: fock_a, rho_a, fock_b, rho_b

   integer :: cdft_iter, max_cdft_iter
   logical :: cdft_converged = .false.
   LIODBLE :: energ, energ2
   LIODBLE, allocatable :: Pmat_old(:), Wmat_vec(:), Wmat(:,:)


   max_cdft_iter = 100
   cdft_iter     = 0
   allocate(Pmat_old(size(Pmat_v,1)))

   call cdft_initialise(natom)
   if (cdft_c%mixed) call cdft_mixed_initialise(nbasis, nOcc, nOcc_b, op_shell)
   
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
         ! Gets W for state 1, retrieves MO
         call cdft_mixed_set_coefs(coefs, .true., 1)
         if (op_shell) call cdft_mixed_set_coefs(coefs_b, .false., 1)

         allocate(Wmat_vec(size(Pmat_v,1)))
         Wmat_vec = 0.0D0
         call g2g_cdft_w(Wmat_vec)
      endif
   enddo

   ! If it is a mixed calculation, repeats everything for the second state.
   if (cdft_c%mixed) then
      cdft_converged = .false.
      cdft_iter      = 0
      call cdft_mixed_switch()

      allocate(Wmat(1,1)) ! Compilator warnings...
      do while ((.not. cdft_converged) .and. (cdft_iter < max_cdft_iter))
         cdft_iter = cdft_iter +1
         Pmat_old  = Pmat_v
         call SCF(energ2, fock_a, rho_a, fock_b, rho_b)
         call cdft_check_conver(Pmat_v, Pmat_old, cdft_converged, &
                                cdft_iter, energ2, told)
   
         if (.not. cdft_converged) then
            ! Calculates perturbations and Jacobian.
            call cdft_get_deltaV(fock_a, rho_a, fock_b, rho_b)
            call cdft_set_potential()
         else
            ! Retrieves MO and W for state 2
            call cdft_mixed_set_coefs(coefs, .true., 2)
            if (op_shell) call cdft_mixed_set_coefs(coefs_b, .false., 2)

            deallocate(Wmat)
            allocate(Wmat(nbasis,nbasis))
            ! We accumulate 1+2 over Wmat_vec and then extract it.
            call g2g_cdft_w(Wmat_vec)

            Wmat = 0.0D0
            call spunpack('L', nbasis, Wmat_vec, Wmat)
            deallocate(Wmat_vec)
         endif
      enddo

      call cdft_mixed_hab(energ, energ2, Wmat, overlap, op_shell)
   endif

   call cdft_finalise()
   call cdft_mixed_finalise()
   deallocate(Pmat_old)
   if (cdft_c%mixed) deallocate(Wmat)
end subroutine CDFT