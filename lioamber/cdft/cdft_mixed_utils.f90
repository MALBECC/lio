subroutine cdft_mixed_initialise(n_basis, n_occ, n_occ2, op_shell)
   use cdft_data, only: cdft_mc

   implicit none
   integer, intent(in) :: n_basis, n_occ, n_occ2
   logical, intent(in) :: op_shell

   if (allocated(cdft_mc%coefs_a1)) deallocate(cdft_mc%coefs_a1)
   if (allocated(cdft_mc%coefs_a2)) deallocate(cdft_mc%coefs_a2)
   allocate(cdft_mc%coefs_a1(n_basis, n_occ))
   allocate(cdft_mc%coefs_a2(n_basis, n_occ))

   if (op_shell) then
      if (allocated(cdft_mc%coefs_b1)) deallocate(cdft_mc%coefs_b1)
      if (allocated(cdft_mc%coefs_b2)) deallocate(cdft_mc%coefs_b2)
      allocate(cdft_mc%coefs_b1(n_basis, n_occ2))
      allocate(cdft_mc%coefs_b2(n_basis, n_occ2))
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
   call g2g_cdft_set_V(cdft_reg%Vc, cdft_reg%Vs)
end subroutine cdft_mixed_switch

subroutine cdft_mixed_invert_spin()
   ! Changes the sign of the potential Vs.
   use cdft_data, only: cdft_c, cdft_reg

   implicit none
   integer :: ireg

   do ireg = 1, cdft_c%n_regions
      cdft_reg%Vs(ireg) = - cdft_reg%Vs(ireg)
   enddo
   
   call g2g_cdft_set_V(cdft_reg%Vc, cdft_reg%Vs)
end subroutine cdft_mixed_invert_spin

subroutine cdft_mixed_set_coefs(coefs, alpha, state)
   use cdft_data, only: cdft_mc

   implicit none
   LIODBLE, intent(in) :: coefs(:,:)
   logical, intent(in) :: alpha
   integer, intent(in) :: state

   integer :: ntop

   if (alpha) then
      ntop = size(cdft_mc%coefs_a1,2)
      if (state == 1) cdft_mc%coefs_a1 = coefs(:,1:ntop)
      if (state == 2) cdft_mc%coefs_a2 = coefs(:,1:ntop)
   else
      ntop = size(cdft_mc%coefs_b1,2)
      if (state == 1) cdft_mc%coefs_b1 = coefs(:,1:ntop)
      if (state == 2) cdft_mc%coefs_b2 = coefs(:,1:ntop)
   endif
end subroutine cdft_mixed_set_coefs

subroutine cdft_mixed_print(Hab)
   implicit none
   LIODBLE, intent(in) :: Hab(2,2)

   write(*,'(A)') "CDFT - Orthogonalised H for states 1 and 2"
   write(*,*) "  H_11 = ", Hab(1,1), " Eh"
   write(*,*) "  H_22 = ", Hab(2,2), " Eh"
   write(*,*) "  H_21 = ", Hab(2,1), " Eh"
   write(*,*) "  H_12 = ", Hab(1,2), " Eh"

end subroutine cdft_mixed_print

