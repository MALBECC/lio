! General setup
subroutine dftd3_setup(n_atoms, atom_z)
   use dftd3_data, only: dftd3, c6_ab, r0_ab, c8_ab, c6_cn, r_cov, c8_coef
   implicit none
   integer, intent(in)  :: n_atoms, atom_z(:)

   if (.not. dftd3) return
   if (.not. allocated(c6_ab  )) allocate(c6_ab(n_atoms,n_atoms))
   if (.not. allocated(r0_ab  )) allocate(r0_ab(n_atoms,n_atoms))
   if (.not. allocated(c8_ab  )) allocate(c8_ab(n_atoms,n_atoms))
   if (.not. allocated(c6_cn  )) allocate(c6_cn(n_atoms,n_atoms,5,5,3))
   if (.not. allocated(r_cov  )) allocate(r_cov(n_atoms))
   if (.not. allocated(c8_coef)) allocate(c8_coef(n_atoms))


   call dftd3_read_c6(c6_cn  , n_atoms, atom_z)
   call dftd3_read_r0(r0_ab  , n_atoms, atom_z)
   call dftd3_read_rc(r_cov  , n_atoms, atom_z)
   call dftd3_read_c8(c8_coef, n_atoms, atom_z)
end subroutine dftd3_setup

subroutine dftd3_finalise()
   use dftd3_data, only: dftd3, c6_ab, r0_ab, c8_ab, c6_cn, r_cov, c8_coef
   implicit none

   if (.not. dftd3) return
   if (allocated(c6_ab)) deallocate(c6_ab)
   if (allocated(r0_ab)) deallocate(r0_ab)
   if (allocated(c8_ab)) deallocate(c8_ab)

   if (allocated(c6_cn))   deallocate(c6_cn)
   if (allocated(r_cov))   deallocate(r_cov)
   if (allocated(c8_coef)) deallocate(c8_coef)
end subroutine dftd3_finalise

! Energy calculations
subroutine dftd3_energy(e_disp, dist_in, n_atoms, dist_sq)
   use dftd3_data, only: dftd3
   implicit none
   integer     , intent(in)    :: n_atoms
   logical     , intent(in)    :: dist_sq
   LIODBLE, intent(in)    :: dist_in(:,:)
   LIODBLE, intent(inout) :: e_disp
   
   LIODBLE :: e_disp2, e_disp3
   LIODBLE, allocatable ::dists(:,:)

   if (.not. dftd3) return
   e_disp2 = 0.0D0
   e_disp3 = 0.0D0
   allocate(dists(n_atoms,n_atoms))
   ! Most of LIO stores interatomic distances as d**2, but for this
   ! part we need d.
   if (dist_sq) then
      dists = sqrt(dist_in)
   else
      dists = dist_in
   endif
   
   call dftd3_set_c6c8(dists, n_atoms)
   call dftd3_2bodies_e(e_disp2, dists, n_atoms)
   call dftd3_3bodies_e(e_disp3, dists, n_atoms)

   ! The 3-body term should be negative (as C9 is negative),
   ! so E3 is added and not substracted.
   e_disp = e_disp - e_disp2 + e_disp3
   deallocate(dists)
end subroutine dftd3_energy

! Gradient calculations. This needs to run the energy routine
! for the same set of positions, in order to keep the C6 and
! C8 coefficients appropriate.
! The gradient calcuation is performed via center forward
! differences since the analytic expression is a hassle to 
! implement. However, this scales as N^4 instead of N^3 when
! taking into account 3-body interactions, so an analytic
! implementation is not entirely out of the menu.
subroutine dftd3_gradients(grad, pos, n_atoms)
   use dftd3_data, only: dftd3
   implicit none
   integer     , intent(in)    :: n_atoms
   LIODBLE, intent(inout) :: grad(:,:), pos(:,:)

   integer      :: iatom, icrd
   LIODBLE :: delta_x, E_gnd, E_new, delta_E
   LIODBLE, allocatable :: new_dist(:,:)
   
   if (.not. dftd3) return
   ! Analytic gradients not available (yet).
   !call dftd3_2bodies_g(grad, dists, pos, n_atoms)
   !call dftd3_3bodies_g(grad, dists, pos, n_atoms)

   allocate(new_dist(n_atoms,n_atoms))

   ! Differences will be calculated using a hundreth of the 
   ! shortest distance between atoms. Needs to be checked
   ! manually since MINVAL (FORTRAN intrinsic) would return
   ! 0.0 due to the diagonal elements.
   call dftd3_calc_dists(pos, n_atoms, new_dist)
   delta_x  = MAXVAL(new_dist)
   do iatom = 1      , n_atoms
   do icrd  = iatom+1, n_atoms
      if (delta_x > new_dist(iatom,icrd)) delta_x = new_dist(iatom,icrd)
   enddo
   enddo
   delta_x = 0.005D0 * delta_x
   
   E_gnd = 0.0D0
   call dftd3_energy(E_gnd, new_dist, n_atoms, .false.)
   
   do iatom = 1, n_atoms
   do icrd  = 1, 3
      ! Forward dE
      pos(iatom,icrd) = pos(iatom,icrd) + delta_x
      call dftd3_calc_dists(pos, n_atoms, new_dist)
      E_new = 0.0D0
      call dftd3_energy(E_new, new_dist, n_atoms, .false.)
      delta_E = E_new - E_gnd

      ! Backward dE
      pos(iatom,icrd) = pos(iatom,icrd) - 2.0D0 * delta_x
      call dftd3_calc_dists(pos, n_atoms, new_dist)
      E_new = 0.0D0
      call dftd3_energy(E_new, new_dist, n_atoms, .false.)
      ! Total dE = dEforward - dEbackward
      ! It is divided by to to keep consistency with dx.
      delta_E = 0.5D0 * (delta_E + (E_gnd - E_new))
      
      pos(iatom,icrd)  = pos(iatom,icrd)  + delta_x
      grad(iatom,icrd) = grad(iatom,icrd) + delta_E / delta_x
   enddo
   enddo

   deallocate(new_dist)
end subroutine dftd3_gradients

subroutine dftd3_calc_dists(pos, n_atoms, dist)
   implicit none
   integer     , intent(in)  :: n_atoms
   LIODBLE, intent(in)  :: pos(:,:)
   LIODBLE, intent(out) :: dist(:,:)
   
   LIODBLE :: tmp
   integer      :: iatom, jatom, icrd

   dist = 0.0D0
   do iatom = 1, n_atoms
   do jatom = iatom+1, n_atoms
       do icrd = 1, 3
           tmp = pos(iatom,icrd) - pos(jatom,icrd)
           dist(iatom,jatom) = dist(iatom,jatom) + tmp * tmp
       enddo
       dist(iatom,jatom) = sqrt(dist(iatom,jatom))
       dist(jatom,iatom) = dist(iatom,jatom)
   enddo
   enddo
end subroutine dftd3_calc_dists
