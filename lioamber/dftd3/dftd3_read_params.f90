subroutine dftd3_read_c6(c6_cn, n_atoms, atom_z)
   implicit none
   integer     , intent(in)  :: n_atoms, atom_z(:)
   LIODBLE, intent(out) :: c6_cn(:,:,:,:,:)

   integer :: i_line, n_lines, iatom, jatom, cni, cnj, k_ind
   LIODBLE, allocatable :: c6_temp(:,:,:,:,:), pars(:)
   
   allocate(c6_temp(94,94,5,5,3), pars(161925))

   ! n_lines and pars are initialised in pars.f90.
# include "param_c6.f90"

   k_ind   = 1
   c6_temp = 0.0D0
   do i_line = 1, n_lines
      ! For each k_ind, the first term is the C6 reference, 
      ! second and third are atom identities, and fourth and
      ! fifth are the coordination numbers.
      iatom = int(pars(k_ind +1))
      jatom = int(pars(k_ind +2))

      cni = int(iatom/100) +1
      cnj = int(jatom/100) +1
      iatom = mod(iatom,100)
      jatom = mod(jatom,100)

      c6_temp(iatom, jatom, cni, cnj, 1) = pars(k_ind)
      c6_temp(iatom, jatom, cni, cnj, 2) = pars(k_ind+3)
      c6_temp(iatom, jatom, cni, cnj, 3) = pars(k_ind+4)

      c6_temp(jatom, iatom, cnj, cni, 1) = pars(k_ind)
      c6_temp(jatom, iatom, cnj, cni, 2) = pars(k_ind+4)
      c6_temp(jatom, iatom, cnj, cni, 3) = pars(k_ind+3)
      k_ind = k_ind + 5
   enddo

   c6_cn = 0.0D0
   do iatom = 1, n_atoms
      if (atom_z(iatom) > 0) then
      do jatom = iatom, n_atoms
         if (atom_z(jatom) > 0) then
            do cni = 1, 5
            do cnj = 1, 5
               c6_cn(iatom, jatom, cni, cnj, 1) = &
                  c6_temp(atom_z(iatom),atom_z(jatom), cni, cnj, 1)
               c6_cn(iatom, jatom, cni, cnj, 2) = &
                  c6_temp(atom_z(iatom),atom_z(jatom), cni, cnj, 2)
               c6_cn(iatom, jatom, cni, cnj, 3) = &
                  c6_temp(atom_z(iatom),atom_z(jatom), cni, cnj, 3)

               c6_cn(jatom, iatom, cnj, cni, 1) = &
                                   c6_cn(iatom, jatom, cni, cnj, 1)
               c6_cn(jatom, iatom, cnj, cni, 2) = &
                                   c6_cn(iatom, jatom, cni, cnj, 2)
               c6_cn(jatom, iatom, cnj, cni, 3) = &
                                   c6_cn(iatom, jatom, cni, cnj, 3)
            enddo
            enddo
         endif
      enddo
      endif
   enddo
   deallocate(c6_temp, pars)
end subroutine dftd3_read_c6

subroutine dftd3_read_r0(r0_ab, n_atoms, atom_z)
   implicit none
   integer     , intent(in)  :: n_atoms, atom_z(:)
   LIODBLE, intent(out) :: r0_ab(:,:)

   integer :: iatom, jatom, k_ind
   LIODBLE, allocatable :: r0ab(:), r_temp(:,:)

   allocate(r0ab(4465), r_temp(94,94))
#  include "param_r0.f90"

   k_ind  = 0
   r_temp = 0.0D0
   do iatom = 1, 94
   do jatom = 1, iatom
      k_ind = k_ind +1
      ! Converts Angstrom to AU.
      ! 1 / 0.52917726 = 1.8897259493
      r_temp(iatom,jatom) = r0ab(k_ind) * 1.8897259493D0
      r_temp(jatom,iatom) = r_temp(iatom,jatom)
   enddo
   enddo

   do iatom = 1       , n_atoms
      if (atom_z(iatom) > 0) then
         do jatom = iatom +1, n_atoms
            if (atom_z(jatom) > 0) then
               r0_ab(iatom,jatom) = r_temp(atom_z(iatom),atom_z(jatom))
               r0_ab(jatom,iatom) = r0_ab(iatom,jatom)
            endif
         enddo
      endif
   enddo

   deallocate(r0ab, r_temp)
end subroutine dftd3_read_r0

subroutine dftd3_read_rc(r_cov, n_atoms, atom_z)
   implicit none
   integer     , intent(in)  :: n_atoms, atom_z(:)
   LIODBLE, intent(out) :: r_cov(:)

   LIODBLE, allocatable :: rcov(:)
   integer                   :: iatom

   allocate(rcov(94))
   ! These covalent radii are already in AU and rescaled by 4/3
   ! in order to be used in Eq. 15.
   rcov = (/ &
     0.80628308, 1.15903197, 3.02356173, 2.36845659, 1.94011865, &
     1.88972601, 1.78894056, 1.58736983, 1.61256616, 1.68815527, &
     3.52748848, 3.14954334, 2.84718717, 2.62041997, 2.77159820, &
     2.57002732, 2.49443835, 2.41884923, 4.43455700, 3.88023730, &
     3.35111422, 3.07395437, 3.04875805, 2.77159820, 2.69600923, &
     2.62041997, 2.51963467, 2.49443835, 2.54483100, 2.74640188, &
     2.82199085, 2.74640188, 2.89757982, 2.77159820, 2.87238349, &
     2.94797246, 4.76210950, 4.20778980, 3.70386304, 3.50229216, &
     3.32591790, 3.12434702, 2.89757982, 2.84718717, 2.84718717, &
     2.72120556, 2.89757982, 3.09915070, 3.22513231, 3.17473967, &
     3.17473967, 3.09915070, 3.32591790, 3.30072128, 5.26603625, &
     4.43455700, 4.08180818, 3.70386304, 3.98102289, 3.95582657, &
     3.93062995, 3.90543362, 3.80464833, 3.82984466, 3.80464833, &
     3.77945201, 3.75425569, 3.75425569, 3.72905937, 3.85504098, &
     3.67866672, 3.45189952, 3.30072128, 3.09915070, 2.97316878, &
     2.92277614, 2.79679452, 2.82199085, 2.84718717, 3.32591790, &
     3.27552496, 3.27552496, 3.42670319, 3.30072128, 3.47709584, &
     3.57788113, 5.06446567, 4.56053862, 4.20778980, 3.98102289, &
     3.82984466, 3.85504098, 3.88023730, 3.90543362 /)

   r_cov = 0.0D0
   do iatom = 1, n_atoms
      if (atom_z(iatom) > 0) r_cov(iatom) = rcov(atom_z(iatom))
   enddo
   deallocate(rcov)
end subroutine dftd3_read_rc

subroutine dftd3_read_c8(c8_coef, n_atoms, atom_z)
   implicit none
   integer     , intent(in)  :: n_atoms, atom_z(:)
   LIODBLE, intent(out) :: c8_coef(:)

   LIODBLE, allocatable :: r2r4(:)
   integer                   :: iatom

   allocate(r2r4(94))

   ! r2r4 is already precalculated as sqrt(Qa), for usage in Eq. 6.
   ! The term is sqrt( 0.5 * sqrt(Zi) * <Ri^4> / <Ri^2>)
   r2r4 = (/ &
     2.00734898,  1.56637132,  5.01986934,  3.85379032,  3.64446594, &
     3.10492822,  2.71175247,  2.59361680,  2.38825250,  2.21522516, &
     6.58585536,  5.46295967,  5.65216669,  4.88284902,  4.29727576, &
     4.04108902,  3.72932356,  3.44677275,  7.97762753,  7.07623947, &
     6.60844053,  6.28791364,  6.07728703,  5.54643096,  5.80491167, &
     5.58415602,  5.41374528,  5.28497229,  5.22592821,  5.09817141, &
     6.12149689,  5.54083734,  5.06696878,  4.87005108,  4.59089647, &
     4.31176304,  9.55461698,  8.67396077,  7.97210197,  7.43439917, &
     6.58711862,  6.19536215,  6.01517290,  5.81623410,  5.65710424, &
     5.52640661,  5.44263305,  5.58285373,  7.02081898,  6.46815523, &
     5.98089120,  5.81686657,  5.53321815,  5.25477007, 11.02204549, &
    10.15679528,  9.35167836,  9.06926079,  8.97241155,  8.90092807, &
     8.85984840,  8.81736827,  8.79317710,  7.89969626,  8.80588454, &
     8.42439218,  8.54289262,  8.47583370,  8.45090888,  8.47339339, &
     7.83525634,  8.20702843,  7.70559063,  7.32755997,  7.03887381, &
     6.68978720,  6.05450052,  5.88752022,  5.70661499,  5.78450695, &
     7.79780729,  7.26443867,  6.78151984,  6.67883169,  6.39024318, &
     6.09527958, 11.79156076, 11.10997644,  9.51377795,  8.67197068, &
     8.77140725,  8.65402716,  8.53923501,  8.85024712/)

   c8_coef = 0.0D0
   do iatom = 1, n_atoms
      if (atom_z(iatom) > 0) c8_coef(iatom) = r2r4(atom_z(iatom))
   enddo
   deallocate(r2r4)
end subroutine dftd3_read_c8