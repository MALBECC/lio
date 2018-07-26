module subm_int3lu
contains
subroutine int3lu(E2, RMM, M, Md, cool, cools, kkind, kkinds, kknumd, kknums,&
                  af, B, memo, open_shell)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Integrals subroutines - 2e integrals, 3 index                                !
! Wavefunction and density fitting functions are calculated using the          !
! Obara-Saika recursive method.                                                !
! Inputs: G, F, standard basis and density basis.                              !
! F should already have the 1e part, and here the Coulomb part is added without!
! storing the integrals separately.                                            !
! Output: F updated with Coulomb part, also Coulomb energy.                    !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   implicit none
   integer         , intent(in) :: M, Md, kknumd, kknums, kkind(:), kkinds(:)
   logical         , intent(in) :: open_shell, memo
   real            , intent(in) :: cools(:)
   double precision, intent(in) :: cool(:)
   double precision, intent(inout) :: E2, RMM(:), af(:), B(:,:)

   double precision :: Rc(Md), aux(md)
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
   Ea = 0.D0 ; Eb = 0.D0

   MM=M*(M+1)/2
   MMd=Md*(Md+1)/2

   M3=1+MM ! Pew
   M5=M3+MM ! now S, also F later
   M7=M5+MM ! G matrix
   M9=M7+MMd ! G inverted
   M11=M9+MMd ! Hmat

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
            Rc(k_ind) = Rc(k_ind) + RMM(kkind(kk_ind)) * cool(iikk + k_ind)
         enddo
      enddo

      do kk_ind = 1, kknums
         iikk = (kk_ind - 1) * Md
         do k_ind = 1, Md
            Rc(k_ind) = Rc(k_ind) + RMM(kkinds(kk_ind)) * cools(iikk + k_ind)
         enddo
      enddo

      ! Calculation of variational coefficients and fitting coefficients
      do m_ind = 1, Md
         af(m_ind) = 0.0D0
         do k_ind = 1, m_ind-1
            af(m_ind) = af(m_ind) + &
                        Rc(k_ind) * RMM(M9-1 + m_ind + (2*Md-k_ind)*(k_ind-1)/2)
         enddo
         do k_ind = m_ind, Md
            af(m_ind) = af(m_ind) + &
                        Rc(k_ind) * RMM(M9-1 + k_ind + (2*Md-m_ind)*(m_ind-1)/2)
         enddo
      enddo

      ! Initialization of Fock matrix elements
      do k_ind = 1, MM
         RMM(M5-1 + k_ind) = RMM(M11-1 + k_ind)
      enddo
      if (open_shell) then
      do k_ind = 1, MM
         RMM(M3-1 + k_ind) = RMM(M11-1 + k_ind)
      enddo
      endif

      do m_ind = 1, Md
         Ea = Ea + af(m_ind)  * Rc(m_ind)
         do k_ind = 1, m_ind
            Eb = Eb + af(k_ind) * af(m_ind) * &
                      RMM(M7-1 + m_ind + (2*Md-k_ind)*(k_ind-1)/2)
         enddo
         do k_ind = m_ind+1, Md
            Eb = Eb + af(k_ind) * af(m_ind) * &
                      RMM(M7-1 + k_ind + (2*Md-m_ind)*(m_ind-1)/2)
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
               RMM(M5-1 + kkind(kk_ind)) = RMM(M5-1 + kkind(kk_ind)) + &
                                           af(k_ind)  * cool(iikk + k_ind)
               RMM(M3-1 + kkind(kk_ind)) = RMM(M3-1 + kkind(kk_ind)) + &
                                           aux(k_ind) * cool(iikk + k_ind)
            enddo
         enddo
         do kk_ind = 1, kknums
            iikk = (kk_ind - 1) * Md
            do k_ind = 1, Md
               RMM(M5-1 + kkinds(kk_ind)) = RMM(M5-1 + kkinds(kk_ind)) + &
                                             af(k_ind)  * cools(iikk + k_ind)
               RMM(M3-1 + kkinds(kk_ind)) = RMM(M3-1 + kkinds(kk_ind)) + &
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
            RMM(M5-1 + kkind(kk_ind)) = RMM(M5-1 + kkind(kk_ind)) + term
         enddo
         do kk_ind = 1, kknums
            iikk = (kk_ind - 1) * Md
            term = 0.0D0
            do k_ind = 1, Md
              term = term + af(k_ind) * cools(iikk + k_ind)
            enddo
            RMM(M5-1 + kkinds(kk_ind) -1) = RMM(M5-1 + kkinds(kk_ind)) + term
         enddo
      endif
      call g2g_timer_stop('int3lu')
   else
      do k_ind = 1, MM
         RMM(M5-1 + k_ind) = RMM(M11-1 + k_ind)
         if (open_shell) RMM(M3-1 + k_ind) = RMM(M11-1 + k_ind)
      enddo

      call aint_coulomb_fock(Ea)
      do m_ind = 1, Md
         do k_ind = 1, m_ind
            Eb = Eb + af(k_ind) * af(m_ind) * &
                      RMM(M7-1 + m_ind + (2*Md-k_ind)*(k_ind-1)/2)
         enddo
         do k_ind = m_ind+1, Md
            Eb = Eb + af(k_ind) * af(m_ind) * &
                      RMM(M7-1 + k_ind + (2*Md-m_ind)*(m_ind-1)/2)
         enddo
      enddo
   endif

   E2 = Ea - Eb / 2.D0
   return
end subroutine int3lu
end module subm_int3lu
