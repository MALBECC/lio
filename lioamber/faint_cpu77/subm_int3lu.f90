module subm_int3lu
contains
subroutine int3lu(E2)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Integrals subroutines - 2e integrals, 3 index                                !
! Wavefunction and density fitting functions are calculated using the          !
! Obara-Saika recursive method.                                                !
! Inputs: G, F, standard basis and density basis.                              !
! F should already have the 1e part, and here the Coulomb part is added without!
! storing the integrals separately.                                            !
! Output: F updated with Coulomb part, also Coulomb energy.                    !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   use garcha_mod, only: RMM, af, X, B, ngd, md, M, kknumd, kknums, &
                         MEMO, Nang, natom, NCO, NORM, Nunp, OPEN, pi32, &
                         nshell, nshelld, SVD, cool, cools, kkind, kkinds, &
                         ncontd, ad, cd
   implicit none
   double precision, intent(inout) :: E2

   double precision :: Q(3), W(3), Rc(Md), FF(Md), P(Md), Jx(M), aux(ngd)
   double precision  :: sq3, r0, r1, rcond, bda
   double precision  :: t0, ss9, Ex, Ea, Eb, term

   integer :: M1, M2, M3, M5, M7, M9
   integer :: M10, M11
   integer :: MM, MMd, MMp, Md2, Md3, Md5
   integer :: Nel, iconst, irank, info
   integer :: nd, ndd, nk, np, npd, ns, nsd
   integer :: l1, l2, i, j, k1, ll(3)
   integer :: iikk, k_ind, kk_ind, m_ind

   ! 16 loops for all combinations - 1-2: for wavefunction basis, 3 for the
   ! density fitting.
   ! Rc(k) is constructed adding t(i,j,k)*P(i,j), and cf(k), the variationally
   ! obtained fitting coefficient, is obtained by adding R(i)*G-1(i,k)
   ! if t(i,j,k) is not stored, they should be calculated again in order to
   ! evaluate the corresponding part of the Fock matrix.
   ! V(i,j) is obtained by adding af(k_ind) * t(i,j,k).
   ns  = nshell(0) ; np  = nshell(1) ; nd  = nshell(2)
   nsd = nshelld(0); npd = nshelld(1); ndd = nshelld(2)

   Ex = 0.D0 ; Ea = 0.D0 ; Eb = 0.D0

   M2=2*M
   Md2=2*Md
   MM=M*(M+1)/2
   MMd=Md*(Md+1)/2

   M1=1 !first P
   M3=M1+MM ! Pew
   M5=M3+MM ! now S, also F later
   M7=M5+MM ! G matrix
   M9=M7+MMd ! G inverted
   M11=M9+MMd ! Hmat

   sq3 = 1.D0
   if (NORM) sq3 = sqrt(3.D0)

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
      if (OPEN) then
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
      if (OPEN) then
         do k_ind = 1, Md
            aux(k_ind) = af(k_ind)
         enddo
      endif

      call g2g_timer_stop('int3lu - start')
      call g2g_timer_start('int3lu')
      if (open) then
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
         if (OPEN) RMM(M3-1 + k_ind) = RMM(M11-1 + k_ind)
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
