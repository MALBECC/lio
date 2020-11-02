! Initialises MO coefficient storage for both states in 
! mixed CDFT calculations.
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

! This sub exchanges data between the two states, storing
! the data of the first run in the second state.
subroutine cdft_mixed_switch()
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
   ! Or that is what you do if you want to waste time. Since most 
   ! of the time we will be exchanging configurations, we can just
   ! set the new potentials to opposite signs, getting a better
   ! initial guess.
   do ireg = 1, cdft_c%n_regions
      tmpVc = cdft_reg%Vc(ireg)
      tmpVs = cdft_reg%Vs(ireg)

      cdft_reg%Vc(ireg) = -cdft_reg%Vc(ireg)
      cdft_reg%Vs(ireg) = -cdft_reg%Vs(ireg)

      cdft_reg%Vc2(ireg) = tmpVc
      cdft_reg%Vs2(ireg) = tmpVs
   enddo
   call g2g_cdft_set_V(cdft_reg%Vc, cdft_reg%Vs)
end subroutine cdft_mixed_switch

! Changes the sign of the potential Vs.
subroutine cdft_mixed_invert_spin()
   use cdft_data, only: cdft_c, cdft_reg

   implicit none
   integer :: ireg

   do ireg = 1, cdft_c%n_regions
      cdft_reg%Vs(ireg) = - cdft_reg%Vs(ireg)
   enddo
   
   call g2g_cdft_set_V(cdft_reg%Vc, cdft_reg%Vs)
end subroutine cdft_mixed_invert_spin

! Stores MO coefficients for a given electronic state.
subroutine cdft_mixed_set_coefs(coefs, alpha, state)
   use cdft_data, only: cdft_mc

   implicit none
   LIODBLE, intent(in) :: coefs(:,:)
   logical, intent(in) :: alpha
   integer, intent(in) :: state

   integer :: ntop, nbas

   nbas = size(coefs,1)
   if (alpha) then
      ntop = size(cdft_mc%coefs_a1,2)
      if (state == 1) call copy_coef(coefs, cdft_mc%coefs_a1, nbas, ntop)
      if (state == 2) call copy_coef(coefs, cdft_mc%coefs_a2, nbas, ntop)
   else
      ntop = size(cdft_mc%coefs_b1,2)
      if (state == 1) call copy_coef(coefs, cdft_mc%coefs_b1, nbas, ntop)
      if (state == 2) call copy_coef(coefs, cdft_mc%coefs_b2, nbas, ntop)
   endif
end subroutine cdft_mixed_set_coefs

! Used in the above subroutine, it copies a matrix. 
! This is done so to avoid issues with assumed
! shape arrays, although it should not be a problem
! in modern FORTRAN.
subroutine copy_coef(coefs, mc_coefs, nbas, ntop)
   implicit none
   integer, intent(in)    :: ntop, nbas
   LIODBLE, intent(in)    :: coefs(:,:)
   LIODBLE, intent(inout) :: mc_coefs(:,:)

   integer :: ii, jj

   do ii = 1, nbas
      do jj = 1, ntop
         mc_coefs(ii,jj) = coefs(ii,jj)
      enddo
   enddo
end subroutine copy_coef

! Prints outputs.
subroutine cdft_mixed_print(Hab, Sab)
   implicit none
   LIODBLE, intent(in) :: Hab(2,2), Sab

   write(*,'(A)') "CDFT Summary"
   write(*,'(A)') "------------"
   write(*,'(A)') " > Orthogonalised energies for states 1 and 2"
   write(*,'(A11,F14.7,A3)') "    H_11 = ", Hab(1,1), " Eh"
   write(*,'(A11,F14.7,A3)') "    H_22 = ", Hab(2,2), " Eh"
   write(*,*)
   write(*,'(A)') " > Overlap between states 1 and 2"
   write(*,'(A11,F14.7)') "    S_12 = ", Sab
   write(*,*)
   write(*,'(A)') " > Orthogonalised coupling element Hab in mEh and meV"
   write(*,'(A11,F14.7,A4)') "    H_12 = ", Hab(1,2) * 1000.0D0, " mEh"
   write(*,'(11x,F14.7,A4)') Hab(1,2) * 27211.382543519D0, " meV"

end subroutine cdft_mixed_print