module subm_int3lu
contains
subroutine int3lu(E2, rho, Fmat_b, Fmat, Gmat, Ginv, Hmat, open_shell)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Integrals subroutines - 2e integrals, 3 index                                !
! Wavefunction and density fitting functions are calculated using the          !
! Obara-Saika recursive method.                                                !
! Inputs: G, F, standard basis and density basis.                              !
! F should already have the 1e part, and here the Coulomb part is added without!
! storing the integrals separately.                                            !
! Output: F updated with Coulomb part, also Coulomb energy.                    !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   use garcha_mod, only: M, Md, cool, cools, kkind, kkinds, kknumd, kknums, af,&
                         B, MEMO
   implicit none
   logical         , intent(in) :: open_shell
   double precision, intent(in) :: rho(:), Gmat(:), Ginv(:), Hmat(:)
   double precision, intent(inout) :: E2, Fmat_b(:), Fmat(:)

   double precision, allocatable :: Rc(:), aux(:)
   double precision :: Ea, Eb, term
   integer          :: M3, M5, M7, M9, M11, MM, MMd, ll(3), iikk, k_ind, &
                       kk_ind, m_ind

   ! 16 loops for all combinations - 1-2: for wavefunction basis, 3 for the
   ! density fitting.
   ! Rc(k) is constructed adding t(i,j,k)*P(i,j), and cf(k), the variationally
   ! obtained fitting coefficient, is obtained by adding R(i)*G-1(i,k)
   ! if t(i,j,k) is not stored, they should be calculated again in order to
   ! evaluate the corresponding part of the Fock matrix.
   ! V(i,j) is obtained by adding af(k_ind) * t(i,j,k).
   allocate(Rc(Md), aux(md))
   Ea = 0.D0 ; Eb = 0.D0

   MM=M*(M+1)/2
   MMd=Md*(Md+1)/2

   M3  = 1  + MM ! Pew
   M5  = M3 + MM ! now S, also F later
   M7  = M5 + MM ! G matrix
   M9  = M7 + MMd ! G inverted
   M11 = M9 + MMd ! Hmat

   if (MEMO) then
      B = 0.0D0
      call g2g_timer_start('int3lu - start')

      do k_ind = 1, 3
         Ll(k_ind) = k_ind * (k_ind-1) / 2
      enddo

      do k_ind = 1, Md
         Rc(k_ind) = 0.D0
      enddo

      do kk_ind = 1, kknumd
         iikk = (kk_ind - 1) * Md
         do k_ind = 1, Md
            Rc(k_ind) = Rc(k_ind) + rho(kkind(kk_ind)) * cool(iikk + k_ind)
         enddo
      enddo

      do kk_ind = 1, kknums
         iikk = (kk_ind - 1) * Md
         do k_ind = 1, Md
            Rc(k_ind) = Rc(k_ind) + rho(kkinds(kk_ind)) * cools(iikk + k_ind)
         enddo
      enddo

      ! Calculation of variational coefficients and fitting coefficients
      do m_ind = 1, Md
         af(m_ind) = 0.0D0
         do k_ind = 1, m_ind-1
            af(m_ind) = af(m_ind) + &
                        Rc(k_ind) * Ginv(m_ind + (2*Md-k_ind)*(k_ind-1)/2)
         enddo
         do k_ind = m_ind, Md
            af(m_ind) = af(m_ind) + &
                        Rc(k_ind) * Ginv(k_ind + (2*Md-m_ind)*(m_ind-1)/2)
         enddo
      enddo

      ! Initialization of Fock matrix elements
      do k_ind = 1, MM
         Fmat(k_ind) = Hmat(k_ind)
      enddo
      if (open_shell) then
      do k_ind = 1, MM
         Fmat_b(k_ind) = Hmat(k_ind)
      enddo
      endif

      do m_ind = 1, Md
         Ea = Ea + af(m_ind)  * Rc(m_ind)
         do k_ind = 1, m_ind
            Eb = Eb + af(k_ind) * af(m_ind) * &
                      Gmat(m_ind + (2*Md-k_ind)*(k_ind-1)/2)
         enddo
         do k_ind = m_ind+1, Md
            Eb = Eb + af(k_ind) * af(m_ind) * &
                      Gmat(k_ind + (2*Md-m_ind)*(m_ind-1)/2)
         enddo
      enddo

      ! Calculation of all integrals again, in order to build the Fock matrix.
      aux = 0.0D0
      if (open_shell) then
         do k_ind = 1, Md
            aux(k_ind) = af(k_ind)
         enddo
      endif

      call g2g_timer_stop('int3lu - start')
      call g2g_timer_start('int3lu')
      if (open_shell) then
         do kk_ind = 1, kknumd
            iikk = (kk_ind - 1) * Md
            do k_ind = 1, Md
               Fmat(kkind(kk_ind)) = Fmat(kkind(kk_ind)) + &
                                     af(k_ind)  * cool(iikk + k_ind)
               Fmat_b(kkind(kk_ind)) = Fmat_b(kkind(kk_ind)) + &
                                       aux(k_ind) * cool(iikk + k_ind)
            enddo
         enddo
         do kk_ind = 1, kknums
            iikk = (kk_ind - 1) * Md
            do k_ind = 1, Md
               Fmat(kkinds(kk_ind)) = Fmat(kkinds(kk_ind)) + &
                                      af(k_ind)  * cools(iikk + k_ind)
               Fmat_b(kkinds(kk_ind)) = Fmat_b(kkinds(kk_ind)) + &
                                        aux(k_ind) * cools(iikk + k_ind)
            enddo
         enddo
      else
         do kk_ind = 1, kknumd
            iikk = (kk_ind - 1) * Md
            term = 0.0D0
            do k_ind = 1, Md
              term = term + af(k_ind) * cool(iikk + k_ind)
            enddo
            Fmat(kkind(kk_ind)) = Fmat(kkind(kk_ind)) + term
         enddo
         do kk_ind = 1, kknums
            iikk = (kk_ind - 1) * Md
            term = 0.0D0
            do k_ind = 1, Md
              term = term + af(k_ind) * cools(iikk + k_ind)
            enddo
            Fmat(kkinds(kk_ind)) = Fmat(kkinds(kk_ind)) + term
         enddo
      endif
      call g2g_timer_stop('int3lu')
   else
      do k_ind = 1, MM
         Fmat(k_ind) = Hmat(k_ind)
         if (open_shell) Fmat_b(k_ind) = Hmat(k_ind)
      enddo

      call aint_coulomb_fock(Ea)
      do m_ind = 1, Md
         do k_ind = 1, m_ind
            Eb = Eb + af(k_ind) * af(m_ind) * &
                      Gmat(m_ind + (2*Md-k_ind)*(k_ind-1)/2)
         enddo
         do k_ind = m_ind+1, Md
            Eb = Eb + af(k_ind) * af(m_ind) * &
                      Gmat(k_ind + (2*Md-m_ind)*(m_ind-1)/2)
         enddo
      enddo
   endif

   E2 = Ea - Eb / 2.D0
   deallocate(Rc, aux)
   return
end subroutine int3lu
end module subm_int3lu
