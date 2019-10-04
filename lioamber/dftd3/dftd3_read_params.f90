subroutine dftd3_read_c6(c6_cn, atom_z, n_atoms)
   implicit none
   integer     , intent(in)  :: n_atoms, atom_z
   real(kind=8), intent(out) :: c6_cn

   integer :: i_line, n_lines, iatom, jatom, cni, cnj, k_ind
   real(kind=8), allocatable :: c6_temp(:,:,:,:,:), pars(:)
   
   allocate(c6_temp(94,94,5,5,3), pars(161925))

   ! n_lines and pars are initialised in pars.f90.
# include "pars.f90"

   k_ind = 1
   do i_line = 1, n_lines
      iat = int(pars(k_ind +1))
      jat = int(pars(k_ind +2))

      cni = int(iat/100) +1
      cnj = int(iat/100) +1
      iat = mod(iat,100)
      jat = mod(jat,100)

      c6_temp(iat, jat, cni, cnj, 1) = pars(k_ind)
      c6_temp(iat, jat, cni, cnj, 2) = pars(k_ind+3)
      c6_temp(iat, jat, cni, cnj, 3) = pars(k_ind+4)

      c6_temp(jat, iat, cnj, cni, 1) = pars(k_ind)
      c6_temp(jat, iat, cnj, cni, 2) = pars(k_ind+4)
      c6_temp(jat, iat, cnj, cni, 3) = pars(k_ind+3)
      k_ind = k_ind + i_line * 5
   enddo

   do iatom = 1, n_atoms
   do jatom = 1, n_atoms
      do cni = 1, 5
      do cnj = 1, 5
         c6_cn(iatom, jatom, cni, cnj, 1) = &
            c6_temp(atom_z(iatom),atom_z(jatom), cni, cnj, 1)
         c6_cn(iatom, jatom, cni, cnj, 2) = &
            c6_temp(atom_z(iatom),atom_z(jatom), cni, cnj, 2)
         c6_cn(iatom, jatom, cni, cnj, 3) = &
            c6_temp(atom_z(iatom),atom_z(jatom), cni, cnj, 3)
      enddo
      enddo
   enddo
   enddo

   deallocate(c6_temp, pars)
end subroutine dftd3_read_c6
