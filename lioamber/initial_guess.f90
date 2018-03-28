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
   integer         , intent(in)  :: M, MM, natom, NCO, NCOb, nshell(0:2), &
                                    Iz(natom),Nuc(M)
   logical         , intent(in)  :: openshell
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
      total_iz = total_iz + Iz(icount)
      n_elecs(icount, :) = atomic_ce(Iz(icount), :)
   enddo

   do icount = 1, nshell(0)
      atom_id = Nuc(icount)
      if (n_elecs(atom_id,1) >= 2) then
         start_dens(icount,icount) = 2.0D0 * factor
         n_elecs(atom_id,1) = n_elecs(atom_id,1) - 2
      else if (n_elecs(atom_id,1) > 0) then
         start_dens(icount,icount) = 1.0D0 * factor
         n_elecs(atom_id,1) = 0
      endif
   enddo

   do icount = nshell(0)+1, nshell(1)+nshell(0), 3
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

   do icount = nshell(1)+nshell(0)+1, nshell(2)+nshell(0)+nshell(1), 6
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

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!     CASE OF NO STARTING GUESS PROVIDED, 1 E FOCK MATRIX USED
!     FCe = SCe; (X^T)SX = 1
!     F' = (X^T)FX
!     => (X^-1*C)^-1 * F' * (X^-1*C) = e
!
! This subroutine takes the Hmat as a vector, transforms it into a matrix,
! diagonalizes it, and builds the density from the resulting orbitals.
!
! TODO:
! (*) The starting guess should take care of calculating everything, it
!     shouldn't need to eat Hmat. I do understand the convenience since
!     it is the same for the whole SCF loop. Perhaps make it optional??
! (*) Diagonalizing fock and building the density matrix are things that
!     are or should be performed by dedicated subroutines, and these
!     should be called both here and in the SCF cycle.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine initial_guess_1e(Nmat, Nvec, NCO, ocupF, hmat_vec, Xmat, densat_vec)
   use liosubs_math, only: transform
   use liosubs_dens, only: messup_densmat

   implicit none
   integer         , intent(in)    :: Nmat, Nvec, NCO
   double precision, intent(in)    :: ocupF, Xmat(Nmat,Nmat)
   double precision, intent(inout) :: hmat_vec(Nvec), densat_vec(Nvec)

   double precision, allocatable   :: morb_energy(:), morb_coefon(:,:),   &
                                      morb_coefat(:,:), morb_coefoc(:,:), &
                                      hmat(:,:), dens_mao(:,:)
   double precision, allocatable   :: WORK(:)
   integer                         :: LWORK, info

   call g2g_timer_start('initial guess')
   call g2g_timer_sum_start('initial guess')

   allocate( morb_coefon(Nmat, Nmat), morb_energy(Nvec), dens_mao(Nmat, Nmat) )
   allocate( morb_coefat(Nmat, Nmat), morb_coefoc(Nmat, NCO), hmat(Nmat,Nmat) )

   call spunpack('L', Nmat, hmat_vec, hmat )
   morb_coefon(:,:) = transform( hmat, Xmat )
   morb_energy(:)   = 0.0d0

   LWORK = -1
   if ( allocated(WORK) ) deallocate(WORK)
   allocate( WORK( 1 ) )
   call dsyev('V', 'L', Nmat, morb_coefon, Nmat, morb_energy, WORK, LWORK, info)

   LWORK = NINT( WORK(1) )
   if ( allocated(WORK) ) deallocate(WORK)
   allocate( WORK( LWORK ) )
   call dsyev('V', 'L', Nmat, morb_coefon, Nmat, morb_energy, WORK, LWORK, info)

   morb_coefat = matmul( Xmat, morb_coefon )
   morb_coefoc(1:Nmat,1:NCO) = morb_coefat(1:Nmat,1:NCO)
   dens_mao = ocupF * matmul( morb_coefoc, transpose(morb_coefoc) )
   call messup_densmat( dens_mao )
   call sprepack( 'L', Nmat, densat_vec, dens_mao)

   call g2g_timer_stop('initial guess')
   call g2g_timer_sum_stop('initial guess')
end subroutine initial_guess_1e

end module initial_guess_subs
