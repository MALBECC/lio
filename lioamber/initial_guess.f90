!## INITIAL_GUESS.F90 #########################################################!
! This file contains data and subroutines used to perform the initial guess in !
! SCF. Input variable "initial_guess" is cointained here.                      !
! Modules contained:                                                           !
!       (*) initial_guess_data                                                 !
!       (*) initial_guess_subs                                                 !
!                                                                              !
! Externally accessed subroutines:                                             !
!       (*) get_initial_guess                                                  !
! This subroutine receives the following variables as input:                   !
!   M:         (integer) Number of basis functions.                            !
!   MM:        (integer) The size of vectorised matrices,  M*(M+1)/2.          !
!   NCO:       (integer) Number of occupied orbitals (alpha orbitals in OS).   !
!   NCOb:      (integer) Number of occupied beta orbitals (only used in OS).   !
!   Xmat:      (dp real, size(M,M)) The basechange matrix from atomic to       !
!              canonical basis.                                                !
!   Hvec:      (dp real, size(MM)) Vectorised 1-e integrals (only for 1e       !
!              integral-type guess).                                           !
!   openshell: (logical) Indicates if the calculation is OS.                   !
!   natom:     (integer) Total number of atoms (used in aufbau-type guess).    !
!   Iz:        (integer, size(natom)) Contains the system's atomic numbers.    !
!   nshell:    (integer, size(0:2)) Contains the number of s, p, and d basis.  !
!   Nuc:       (integer, size(M)) For each i basis function, Nuc(i) indicates  !
!              the atom in which that function is centered.                    !
! The following are the outputs:                                               !
!   Rhovec:   (dp real, size(MM)) Vectorised density matrix.                   !
!   Rhoalpha: (dp real, size(MM)) Vectorised alpha density matrix. (only used  !
!             in OS).                                                          !
!   Rhobeta: (dp real, size(MM)) Vectorised beta density matrix (only used in  !
!            OS).                                                              !
!                                                                              !
! Internal subroutines:                                                        !
!        (*) initialise_eec       (initialises EEC array for aufbau-type guess)!
!        (*) initial_guess_aufbau (performs aufbau-type guess)                 !
!        (*) initial_guess_1e     (performs 1e-integral-type guess)            !
!                                                                              !
! Last modified: 28/03/18 Federico N. Pedron                                   !
!##############################################################################!
module initial_guess_data
  integer :: initial_guess = 0
  integer :: atomic_eec(0:54,3)

contains
   subroutine initialise_eec()
      implicit none

      atomic_eec(0:54,:) = 0
      atomic_eec(1,1)  = 1 ;                                               ! H
      atomic_eec(2,1)  = 2 ;                                               ! He
      atomic_eec(3,1)  = 3 ;                                               ! Li
      atomic_eec(4,1)  = 4 ;                                               ! Be
      atomic_eec(5,1)  = 4 ; atomic_eec(5,2)  = 1 ;                        ! B
      atomic_eec(6,1)  = 4 ; atomic_eec(6,2)  = 2 ;                        ! C
      atomic_eec(7,1)  = 4 ; atomic_eec(7,2)  = 3 ;                        ! N
      atomic_eec(8,1)  = 4 ; atomic_eec(8,2)  = 4 ;                        ! O
      atomic_eec(9,1)  = 4 ; atomic_eec(9,2)  = 5 ;                        ! F
      atomic_eec(10,1) = 4 ; atomic_eec(10,2) = 6 ;                        ! Ne
      atomic_eec(11,1) = 5 ; atomic_eec(11,2) = 6 ;                        ! Na
      atomic_eec(12,1) = 6 ; atomic_eec(12,2) = 6 ;                        ! Mg
      atomic_eec(13,1) = 6 ; atomic_eec(13,2) = 7 ;                        ! Al
      atomic_eec(14,1) = 6 ; atomic_eec(14,2) = 8 ;                        ! Si
      atomic_eec(15,1) = 6 ; atomic_eec(15,2) = 9 ;                        ! P
      atomic_eec(16,1) = 6 ; atomic_eec(16,2) = 10;                        ! S
      atomic_eec(17,1) = 6 ; atomic_eec(17,2) = 11;                        ! Cl
      atomic_eec(18,1) = 6 ; atomic_eec(18,2) = 12;                        ! Ar
      atomic_eec(19,1) = 7 ; atomic_eec(19,2) = 12;                        ! K
      atomic_eec(20,1) = 8 ; atomic_eec(20,2) = 12;                        ! Ca
      atomic_eec(21,1) = 8 ; atomic_eec(21,2) = 12; atomic_eec(21,3) = 1 ; ! Sc 
      atomic_eec(22,1) = 8 ; atomic_eec(22,2) = 12; atomic_eec(22,3) = 2 ; ! Ti
      atomic_eec(23,1) = 8 ; atomic_eec(23,2) = 12; atomic_eec(23,3) = 3 ; ! V
      atomic_eec(24,1) = 7 ; atomic_eec(24,2) = 12; atomic_eec(24,3) = 5 ; ! Cr
      atomic_eec(25,1) = 8 ; atomic_eec(25,2) = 12; atomic_eec(25,3) = 5 ; ! Mn
      atomic_eec(26,1) = 8 ; atomic_eec(26,2) = 12; atomic_eec(26,3) = 6 ; ! Fe
      atomic_eec(27,1) = 8 ; atomic_eec(27,2) = 12; atomic_eec(27,3) = 7 ; ! Co
      atomic_eec(28,1) = 8 ; atomic_eec(28,2) = 12; atomic_eec(28,3) = 8 ; ! Ni
      atomic_eec(29,1) = 7 ; atomic_eec(29,2) = 12; atomic_eec(29,3) = 10; ! Cu
      atomic_eec(30,1) = 8 ; atomic_eec(30,2) = 12; atomic_eec(30,3) = 10; ! Zn
      atomic_eec(31,1) = 8 ; atomic_eec(31,2) = 13; atomic_eec(31,3) = 10; ! Ga
      atomic_eec(32,1) = 8 ; atomic_eec(32,2) = 14; atomic_eec(32,3) = 10; ! Ge
      atomic_eec(33,1) = 8 ; atomic_eec(33,2) = 15; atomic_eec(33,3) = 10; ! As
      atomic_eec(34,1) = 8 ; atomic_eec(34,2) = 16; atomic_eec(34,3) = 10; ! Se
      atomic_eec(35,1) = 8 ; atomic_eec(35,2) = 17; atomic_eec(35,3) = 10; ! Br
      atomic_eec(36,1) = 8 ; atomic_eec(36,2) = 18; atomic_eec(36,3) = 10; ! Kr
      atomic_eec(37,1) = 9 ; atomic_eec(37,2) = 18; atomic_eec(37,3) = 10; ! Rb
      atomic_eec(38,1) = 10; atomic_eec(38,2) = 18; atomic_eec(38,3) = 10; ! Sr
      atomic_eec(39,1) = 10; atomic_eec(39,2) = 18; atomic_eec(39,3) = 11; ! Y
      atomic_eec(40,1) = 10; atomic_eec(40,2) = 18; atomic_eec(40,3) = 12; ! Zr
      atomic_eec(41,1) = 9 ; atomic_eec(41,2) = 18; atomic_eec(41,3) = 14; ! Nb
      atomic_eec(42,1) = 9 ; atomic_eec(42,2) = 18; atomic_eec(42,3) = 15; ! Mo
      atomic_eec(43,1) = 10; atomic_eec(43,2) = 18; atomic_eec(43,3) = 15; ! Tc
      atomic_eec(44,1) = 9 ; atomic_eec(44,2) = 18; atomic_eec(44,3) = 17; ! Ru
      atomic_eec(45,1) = 9 ; atomic_eec(45,2) = 18; atomic_eec(45,3) = 18; ! Rh
      atomic_eec(46,1) = 8 ; atomic_eec(46,2) = 18; atomic_eec(46,3) = 20; ! Pd
      atomic_eec(47,1) = 9 ; atomic_eec(47,2) = 18; atomic_eec(47,3) = 20; ! Ag
      atomic_eec(48,1) = 10; atomic_eec(48,2) = 18; atomic_eec(48,3) = 20; ! Cd
      atomic_eec(49,1) = 10; atomic_eec(49,2) = 19; atomic_eec(49,3) = 20; ! In
      atomic_eec(50,1) = 10; atomic_eec(50,2) = 20; atomic_eec(50,3) = 20; ! Sn
      atomic_eec(51,1) = 10; atomic_eec(51,2) = 21; atomic_eec(51,3) = 20; ! Sb
      atomic_eec(52,1) = 10; atomic_eec(52,2) = 22; atomic_eec(52,3) = 20; ! Te
      atomic_eec(53,1) = 10; atomic_eec(53,2) = 23; atomic_eec(53,3) = 20; ! I
      atomic_eec(54,1) = 10; atomic_eec(54,2) = 24; atomic_eec(54,3) = 20; ! Xe
   end subroutine initialise_eec
end module initial_guess_data

!##############################################################################!
module initial_guess_subs

contains

! This subroutine is the interface between SCF and the initial guess choice.
subroutine get_initial_guess(M, MM, NCO, NCOb, Xmat, Hvec, Rhovec, rhoalpha, &
                             rhobeta, openshell, natom, Iz, nshell, Nuc)
   use initial_guess_data, only: initial_guess

   implicit none
   double precision, intent(in) :: Xmat(:,:), Hvec(:)
   logical         , intent(in) :: openshell
   integer         , intent(in) :: M, MM, NCO, NCOb, natom, Iz(natom), Nuc(M), &
                                   nshell(0:2)

   double precision, intent(inout) :: Rhovec(:), rhoalpha(:), rhobeta(:)
   double precision :: ocupF

   call g2g_timer_start('initial guess')
   call g2g_timer_sum_start('initial guess')

   select case (initial_guess)
   case (0)
      if (.not. openshell) then
         ocupF = 2.0D0
         call initial_guess_1e(M, MM, NCO, ocupF, Hvec, Xmat, Rhovec )
      else
         ocupF = 1.0D0
         call initial_guess_1e(M, MM, NCO , ocupF, Hvec, Xmat, rhoalpha)
         call initial_guess_1e(M, MM, NCOb, ocupF, Hvec, Xmat, rhobeta)
         Rhovec   = rhoalpha + rhobeta
      end if
   case (1)
      call initial_guess_aufbau(M, MM, Rhovec, rhoalpha, rhobeta, natom, NCO,&
                                NCOb, Iz, nshell, Nuc, openshell)
   case default
      write(*,*) "ERROR - Initial guess: Wrong value for input initial_guess."
   end select

   call g2g_timer_stop('initial guess')
   call g2g_timer_sum_stop('initial guess')

   return
end subroutine get_initial_guess


! Perfoms the initial guess using a modified aufbau principle.
subroutine initial_guess_aufbau(M, MM, RMM, rhoalpha, rhobeta, natom, NCO, &
                                NCOb, Iz, nshell, Nuc, openshell)
   use initial_guess_data, only: atomic_eec, initialise_eec
   implicit none
   integer         , intent(in)  :: M, MM, natom, NCO, NCOb, nshell(0:2), &
                                    Iz(natom),Nuc(M)
   logical         , intent(in)  :: openshell
   double precision, intent(out) :: RMM(MM), rhoalpha(MM), rhobeta(MM)

   double precision :: start_dens(M,M), start_dens_alpha(M,M), &
                       start_dens_beta(M,M)
   integer          :: icount, total_iz
   integer          :: n_elecs(natom,3), atom_id

   call initialise_eec()
   start_dens(:,:) = 0.0D0
   n_elecs = 0

   total_iz = 0
   do icount = 1, natom
      total_iz = total_iz + Iz(icount)
      n_elecs(icount, :) = atomic_eec(Iz(icount), :)
   enddo

   do icount = 1, nshell(0)
      atom_id = Nuc(icount)
      if (n_elecs(atom_id,1) >= 2) then
         start_dens(icount,icount) = 2.0D0
         n_elecs(atom_id,1) = n_elecs(atom_id,1) - 2
      else if (n_elecs(atom_id,1) > 0) then
         start_dens(icount,icount) = 1.0D0
         n_elecs(atom_id,1) = 0
      endif
   enddo

   do icount = nshell(0)+1, nshell(1)+nshell(0), 3
      atom_id = Nuc(icount)
      if (n_elecs(atom_id,2) >= 6) then
         start_dens(icount  ,icount  ) = 2.0D0
         start_dens(icount+1,icount+1) = 2.0D0
         start_dens(icount+2,icount+2) = 2.0D0
         n_elecs(atom_id,2) = n_elecs(atom_id,2) - 6
      else if (n_elecs(atom_id,2) > 0) then
         start_dens(icount, icount)     = dble(n_elecs(atom_id,2)) / 3.0D0
         start_dens(icount+1, icount+1) = dble(n_elecs(atom_id,2)) / 3.0D0
         start_dens(icount+2, icount+2) = dble(n_elecs(atom_id,2)) / 3.0D0
         n_elecs(atom_id,2) = 0
      endif
   enddo

   do icount = nshell(1)+nshell(0)+1, nshell(2)+nshell(0)+nshell(1), 6
      atom_id = Nuc(icount)
      if (n_elecs(atom_id,3) >= 10) then
         start_dens(icount, icount)     = 5.0D0 / 3.0D0
         start_dens(icount+1, icount+1) = 5.0D0 / 3.0D0
         start_dens(icount+2, icount+2) = 5.0D0 / 3.0D0
         start_dens(icount+3, icount+3) = 5.0D0 / 3.0D0
         start_dens(icount+4, icount+4) = 5.0D0 / 3.0D0
         start_dens(icount+5, icount+5) = 5.0D0 / 3.0D0
         n_elecs(atom_id,3) = n_elecs(atom_id,3) - 10
      else if (n_elecs(atom_id,3) > 0) then
         start_dens(icount, icount)     = dble(n_elecs(atom_id,3)) / 6.0D0
         start_dens(icount+1, icount+1) = dble(n_elecs(atom_id,3)) / 6.0D0
         start_dens(icount+2, icount+2) = dble(n_elecs(atom_id,3)) / 6.0D0
         start_dens(icount+3, icount+3) = dble(n_elecs(atom_id,3)) / 6.0D0
         start_dens(icount+4, icount+4) = dble(n_elecs(atom_id,3)) / 6.0D0
         start_dens(icount+5, icount+5) = dble(n_elecs(atom_id,3)) / 6.0D0
         n_elecs(atom_id,3) = 0
      endif
   enddo

   if (openshell) then
      start_dens_alpha(:,:) = start_dens(:,:) * dble(NCO ) / dble(total_iz)
      start_dens_beta(:,:)  = start_dens(:,:) * dble(NCOb) / dble(total_iz)
      start_dens(:,:)       = start_dens_alpha(:,:) + start_dens_beta(:,:)
      call sprepack('L', M, rhoalpha, start_dens_alpha)
      call sprepack('L', M, rhobeta , start_dens_beta)
   else
      start_dens(:,:) = start_dens(:,:) * dble(NCO*2) / dble(total_iz)
   endif

   call sprepack('L', M, RMM, start_dens)
   return
end subroutine initial_guess_aufbau

! This subroutine performs the 1e-integral guess. It takes the Hmat as a      !
! vector, transforms it into a matrix, diagonalizes it, and builds the        !
! density from the resulting orbitals.                                        !
subroutine initial_guess_1e(Nmat, Nvec, NCO, ocupF, hmat_vec, Xmat, densat_vec)
   use liosubs_math, only: transform
   use liosubs_dens, only: messup_densmat

   implicit none
   integer         , intent(in)    :: Nmat, Nvec, NCO
   double precision, intent(in)    :: ocupF, Xmat(Nmat,Nmat), hmat_vec(Nvec)
   double precision, intent(inout) :: densat_vec(Nvec)

   double precision, allocatable   :: morb_energy(:), morb_coefon(:,:),   &
                                      morb_coefat(:,:), morb_coefoc(:,:), &
                                      hmat(:,:), dens_mao(:,:)
   double precision, allocatable   :: WORK(:)
   integer                         :: LWORK, info

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

   return
end subroutine initial_guess_1e

end module initial_guess_subs
