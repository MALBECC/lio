subroutine cdft_mixed_initialise(matrix_size, op_shell)
   use cdft_data, only: cdft_mc

   implicit none
   integer, intent(in) :: matrix_size
   logical, intent(in) :: op_shell

   if (allocated(cdft_mc%coefs_a1)) deallocate(cdft_mc%coefs_a1)
   if (allocated(cdft_mc%coefs_a2)) deallocate(cdft_mc%coefs_a2)
   allocate(cdft_mc%coefs_a1(matrix_size, matrix_size))
   allocate(cdft_mc%coefs_a2(matrix_size, matrix_size))

   if (op_shell) then
      if (allocated(cdft_mc%coefs_b1)) deallocate(cdft_mc%coefs_b1)
      if (allocated(cdft_mc%coefs_b2)) deallocate(cdft_mc%coefs_b2)
      allocate(cdft_mc%coefs_b1(matrix_size, matrix_size))
      allocate(cdft_mc%coefs_b2(matrix_size, matrix_size))
   endif

   if (allocated(cdft_mc%Wmat1)) deallocate(cdft_mc%Wmat1)
   if (allocated(cdft_mc%Wmat2)) deallocate(cdft_mc%Wmat2)
   allocate(cdft_mc%Wmat1(matrix_size, matrix_size))
   allocate(cdft_mc%Wmat2(matrix_size, matrix_size))
end subroutine cdft_mixed_initialise

subroutine cdft_mixed_finalise()
   use cdft_data, only: cdft_mc

   implicit none

   if (allocated(cdft_mc%coefs_a1)) deallocate(cdft_mc%coefs_a1)
   if (allocated(cdft_mc%coefs_a2)) deallocate(cdft_mc%coefs_a2)
   if (allocated(cdft_mc%coefs_b1)) deallocate(cdft_mc%coefs_b1)
   if (allocated(cdft_mc%coefs_b2)) deallocate(cdft_mc%coefs_b2)
   if (allocated(cdft_mc%Wmat1)   ) deallocate(cdft_mc%Wmat1)
   if (allocated(cdft_mc%Wmat2)   ) deallocate(cdft_mc%Wmat2)
end subroutine cdft_mixed_finalise

subroutine cdft_mixed_switch()
   use cdft_data, only: cdft_c, cdft_reg

   implicit none
   LIODBLE :: tmpchrg, tmpspin
   integer :: ireg

   do ireg = 1, cdft_c%n_regions
      tmpchrg = cdft_reg%chrg2(ireg)
      tmpspin = cdft_reg%spin2(ireg)

      cdft_reg%chrg2(ireg) = cdft_reg%chrg(ireg)
      cdft_reg%spin2(ireg) = cdft_reg%spin(ireg)

      cdft_reg%chrg(ireg) = tmpchrg
      cdft_reg%spin(ireg) = tmpspin
   enddo
end subroutine cdft_mixed_switch

subroutine cdft_mixed_set_coefs(coefs, alpha, state)
   use cdft_data, only: cdft_mc

   implicit none
   LIODBLE, intent(in) :: coefs(:,:)
   logical, intent(in) :: alpha
   integer, intent(in) :: state

   if (alpha) then
      if (state == 1) cdft_mc%coefs_a1 = coefs
      if (state == 2) cdft_mc%coefs_a2 = coefs
   else
      if (state == 1) cdft_mc%coefs_b1 = coefs
      if (state == 2) cdft_mc%coefs_b2 = coefs
   endif
end subroutine cdft_mixed_set_coefs
