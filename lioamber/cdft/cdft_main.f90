! Performs the CDFT iterative procedure.
subroutine CDFT(fock_a, rho_a, fock_b, rho_b, Pmat_vec, natom)
   use typedef_operator, only: operator
   use converger_data  , only: told

   implicit none
   integer       , intent(in)    :: natom
   LIODBLE  , intent(inout) :: Pmat_vec(:)
   type(operator), intent(inout) :: fock_a, rho_a, fock_b, rho_b

   integer      :: cdft_iter, max_cdft_iter
   logical      :: cdft_converged = .false.
   LIODBLE :: energ
   LIODBLE, allocatable :: Pmat_old(:)

   max_cdft_iter = 100
   cdft_iter     = 0
   allocate(Pmat_old(size(Pmat_vec,1)))
   call cdft_initialise(natom)
   do while ((.not. cdft_converged) .and. (cdft_iter < max_cdft_iter))
      cdft_iter = cdft_iter +1
      Pmat_old  = Pmat_vec
      call SCF(energ, fock_a, rho_a, fock_b, rho_b)
      call cdft_check_conver(Pmat_vec, Pmat_old, cdft_converged, &
                             cdft_iter, energ, told)

      if (.not. cdft_converged) then
         ! Calculates perturbations and Jacobian.
         call cdft_get_deltaV(fock_a, rho_a, fock_b, rho_b)
         call cdft_set_potential()
      endif
   enddo

   call cdft_finalise()
   deallocate(Pmat_old)
end subroutine CDFT