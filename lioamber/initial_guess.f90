!##############################################################################!
module initial_guess_data
  integer :: initial_guess = 0
  integer :: atomic_ce(54,3)

contains
   subroutine initialise_ce()
      implicit none

      atomic_ce(:,:) = 0
      atomic_ce(1,1)  = 1 ;                                             ! H
      atomic_ce(6,1)  = 4 ; atomic_ce(6,2)  = 2 ;                       ! C
      atomic_ce(7,1)  = 4 ; atomic_ce(7,2)  = 3 ;                       ! N
      atomic_ce(8,1)  = 4 ; atomic_ce(8,2)  = 4 ;                       ! O
      atomic_ce(26,1) = 8 ; atomic_ce(26,2) = 12; atomic_ce(26,3) = 6 ; !Fe

   end subroutine initialise_ce
end module initial_guess_data

!##############################################################################!
module initial_guess_subs

contains

subroutine initial_guess_aufbau(M, MM, RMM, rhoalpha, rhobeta, natom, NCO, &
                                NCOb, Iz, nshell, Nuc, openshell)
   use initial_guess_data, only: atomic_ce, initialise_ce
   implicit none
   integer, intent(in) :: M, MM, natom, NCO, NCOb, nshell(0:2), Iz(natom),Nuc(M)
   logical, intent(in) :: openshell
   double precision, intent(out) :: RMM(MM), rhoalpha(MM), rhobeta(MM)

   double precision :: start_dens(M,M), start_dens_alpha(M,M), &
                       start_dens_beta(M,M), factor
   integer          :: icount, total_iz, index, jcount
   integer          :: n_elecs(natom,3), atom_id, nunp

   call initialise_ce()

   ! TODO: Crear una rutina que ordene por electronegatividad y luego vaya
   ! sumando o restando electrones a n_elecs, de acuerdo a la carga total y a
   ! los orbitales alpha/beta.

   start_dens(:,:) = 0.0D0
   factor = 1.0D0
   !if (openshell) factor = 0.50D0

   total_iz = 0
   do icount = 1, natom
      n_elecs(icount, :) = atomic_ce(Iz(icount), :)
      total_iz = total_iz + Iz(icount)
   enddo

   do icount = 1, nshell(0)+1
      atom_id = Nuc(icount)
      if (n_elecs(atom_id,1) >= 2) then
         start_dens(icount,icount) = 2.0D0 * factor
         n_elecs(atom_id,1) = n_elecs(atom_id,1) - 2
      else if (n_elecs(atom_id,1) > 0) then
         start_dens(icount,icount) = 1.0D0 * factor
         n_elecs(atom_id,1) = 0
      endif
   enddo

   do icount = nshell(0)+1, nshell(1)+nshell(0)+1, 3
      atom_id = Nuc(icount)
      if (n_elecs(atom_id,2) >= 6) then
         start_dens(icount  ,icount  ) = 2.0D0 * factor
         start_dens(icount+1,icount+1) = 2.0D0 * factor
         start_dens(icount+2,icount+2) = 2.0D0 * factor
         n_elecs(atom_id,2) = n_elecs(atom_id,2) - 6
      else if (n_elecs(atom_id,2) > 0) then
         start_dens(icount, icount)     = dble(n_elecs(atom_id,2))*factor /3.0D0
         start_dens(icount+1, icount+1) = dble(n_elecs(atom_id,2))*factor /3.0D0
         start_dens(icount+2, icount+2) = dble(n_elecs(atom_id,2))*factor /3.0D0
         n_elecs(atom_id,2) = 0
      endif
   enddo

   do icount = nshell(1)+nshell(0)+1, nshell(2)+nshell(0)+nshell(1)+1, 6
      atom_id = Nuc(icount)
      if (n_elecs(atom_id,3) >= 10) then
         start_dens(icount, icount)     = 5.0D0 * factor /3.0D0
         start_dens(icount+1, icount+1) = 5.0D0 * factor /3.0D0
         start_dens(icount+2, icount+2) = 5.0D0 * factor /3.0D0
         start_dens(icount+3, icount+3) = 5.0D0 * factor /3.0D0
         start_dens(icount+4, icount+4) = 5.0D0 * factor /3.0D0
         start_dens(icount+5, icount+5) = 5.0D0 * factor /3.0D0
         n_elecs(atom_id,3) = n_elecs(atom_id,3) - 10
      else if (n_elecs(atom_id,3) > 0) then
         start_dens(icount, icount)     = dble(n_elecs(atom_id,3))*factor /6.0D0
         start_dens(icount+1, icount+1) = dble(n_elecs(atom_id,3))*factor /6.0D0
         start_dens(icount+2, icount+2) = dble(n_elecs(atom_id,3))*factor /6.0D0
         start_dens(icount+3, icount+3) = dble(n_elecs(atom_id,3))*factor /6.0D0
         start_dens(icount+4, icount+4) = dble(n_elecs(atom_id,3))*factor /6.0D0
         start_dens(icount+5, icount+5) = dble(n_elecs(atom_id,3))*factor /6.0D0
         n_elecs(atom_id,3) = 0
      endif
   enddo

   if (openshell) then
      start_dens_alpha(:,:) = start_dens(:,:) * dble(NCO ) / dble(total_iz)
      start_dens_beta(:,:)  = start_dens(:,:) * dble(NCOb) / dble(total_iz)
      start_dens(:,:)       = start_dens_alpha(:,:) + start_dens_beta(:,:)
      call sprepack('L', M, rhoalpha, start_dens_alpha)
      call sprepack('L', M, rhobeta , start_dens_beta)
      print*, "Initial guess Alpha"
      do icount = 1, M
         print*, start_dens_alpha(icount, icount)
      enddo
      print*, "Initial guess Beta"
      do icount = 1, M
         print*, start_dens_beta(icount, icount)
      enddo
   else
      start_dens(:,:) = start_dens(:,:) * dble(NCO*2) / dble(total_iz)
   endif

   call sprepack('L', M, RMM, start_dens)
end subroutine initial_guess_aufbau



end module initial_guess_subs
