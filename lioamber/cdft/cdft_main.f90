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
   LIODBLE, allocatable :: Pmat_old(:), Hmat(:,:)
   LIODBLE, allocatable :: Wmat_vec(:), Wmat(:,:)
   LIODBLE, allocatable :: Wmat_vec_b(:), Wmat_b(:,:) ! For open shell...


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
      endif
   enddo
   if (cdft_c%mixed .and. cdft_converged) then
      ! Gets W for state 1, retrieves MO
      call cdft_mixed_set_coefs(coefs, .true., 1)

      allocate(Wmat_vec(size(Pmat_v,1)))
      Wmat_vec = 0.0D0
      call g2g_cdft_w(Wmat_vec)

      if (op_shell) then
         call cdft_mixed_set_coefs(coefs_b, .false., 1)
         call cdft_mixed_invert_spin()
         
         allocate(Wmat_vec_b(size(Pmat_v,1)))
         Wmat_vec_b = 0.0D0
         call g2g_cdft_w(Wmat_vec_b)
      endif
   endif

   ! If it is a mixed calculation, repeats everything for the second state.
   if (cdft_c%mixed) then
      cdft_converged = .false.
      cdft_iter      = 0
      call cdft_mixed_switch()

      allocate(Wmat(1,1), Wmat_b(1,1)) ! Compilator warnings...
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
            
            ! We accumulate 1+2 over Wmat_vec and then extract it.
            call g2g_cdft_w(Wmat_vec)

            deallocate(Wmat); allocate(Wmat(nbasis,nbasis))
            Wmat = 0.0D0
            call spunpack('L', nbasis, Wmat_vec, Wmat)
            deallocate(Wmat_vec)


            if (op_shell) then
               call cdft_mixed_set_coefs(coefs_b, .false., 2)
               call cdft_mixed_invert_spin()
                 
               call g2g_cdft_w(Wmat_vec_b)

               deallocate(Wmat_b); allocate(Wmat_b(nbasis,nbasis))
               Wmat_b = 0.0D0
               call spunpack('L', nbasis, Wmat_vec_b, Wmat_b)
               deallocate(Wmat_vec_b)
            endif
         endif
      enddo

      if (cdft_converged) then
         ! Retrieves MO and W for state 2
         call cdft_mixed_set_coefs(coefs, .true., 2)
            
         ! We accumulate 1+2 over Wmat_vec and then extract it.
         call g2g_cdft_w(Wmat_vec)

         deallocate(Wmat); allocate(Wmat(nbasis,nbasis))
         Wmat = 0.0D0
         call spunpack('L', nbasis, Wmat_vec, Wmat)
         deallocate(Wmat_vec)


         if (op_shell) then
            call cdft_mixed_set_coefs(coefs_b, .false., 2)
            call cdft_mixed_invert_spin()
              
            call g2g_cdft_w(Wmat_vec_b)

            deallocate(Wmat_b); allocate(Wmat_b(nbasis,nbasis))
            Wmat_b = 0.0D0
            call spunpack('L', nbasis, Wmat_vec_b, Wmat_b)
            deallocate(Wmat_vec_b)
         endif
      endif

      allocate(Hmat(2,2))
      call cdft_mixed_hab(energ, energ2, Wmat, Wmat_b, overlap, op_shell, Hmat)
      call cdft_mixed_print(Hmat)
      deallocate(Hmat)
   endif


   call cdft_finalise()
   call cdft_mixed_finalise()
   deallocate(Pmat_old)
   if (allocated(Wmat))   deallocate(Wmat)
   if (allocated(Wmat_b)) deallocate(Wmat_b)
end subroutine CDFT