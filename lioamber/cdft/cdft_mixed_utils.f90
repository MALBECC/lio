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
end subroutine cdft_mixed_initialise

subroutine cdft_mixed_finalise()
   use cdft_data, only: cdft_mc

   implicit none

   if (allocated(cdft_mc%coefs_a1)) deallocate(cdft_mc%coefs_a1)
   if (allocated(cdft_mc%coefs_a2)) deallocate(cdft_mc%coefs_a2)
   if (allocated(cdft_mc%coefs_b1)) deallocate(cdft_mc%coefs_b1)
   if (allocated(cdft_mc%coefs_b2)) deallocate(cdft_mc%coefs_b2)
end subroutine cdft_mixed_finalise

subroutine cdft_mixed_switch()
   ! This sub exchanges data between the two states, storing
   ! the data of the first run in the second state.
   use cdft_data, only: cdft_c, cdft_reg

   implicit none
   LIODBLE :: tmpchrg, tmpspin, tmpVc, tmpVs
   integer :: ireg

   do ireg = 1, cdft_c%n_regions
      tmpchrg = cdft_reg%chrg2(ireg)
      tmpspin = cdft_reg%spin2(ireg)

      cdft_reg%chrg2(ireg) = cdft_reg%chrg(ireg)
      cdft_reg%spin2(ireg) = cdft_reg%spin(ireg)
 
      cdft_reg%chrg(ireg) = tmpchrg
      cdft_reg%spin(ireg) = tmpspin
   enddo

   ! This might be a little confusing since the exchange is inverse
   ! to the previous (1->2 instead of 2->1). This is so that, 
   ! while we prepare the new targets for CDFT, since it is always
   ! done in state 1, we also store the previous Vcs in state 2.
   do ireg = 1, cdft_c%n_regions
      tmpVc = cdft_reg%Vc(ireg)
      tmpVs = cdft_reg%Vs(ireg)

      cdft_reg%Vc(ireg) = cdft_reg%Vc2(ireg)
      cdft_reg%Vs(ireg) = cdft_reg%Vc2(ireg)

      cdft_reg%Vc2(ireg) = tmpVc
      cdft_reg%Vs2(ireg) = tmpVs
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

subroutine cdft_mixed_hab(e1, e2, Wat, Sat)
   use cdft_data, only: cdft_mc
   implicit none
   LIODBLE, intent(in)    :: e1, e2
   LIODBLE, intent(inout) :: Wat(:,:), Sat(:,:)
end subroutine cdft_mixed_hab