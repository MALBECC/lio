subroutine tsh_probabilities(C,E,Xexc,Eexc,NCO,M,Mlr,Ndim,Nvirt,Etot,Nstat)
use garcha_mod  , only: natom, Pmat_vec, nucvel, atom_mass
use excited_data, only: TSH, root, gamma_old
   implicit none

   integer, intent(in) :: NCO, M, Mlr, Ndim, Nvirt, Nstat
   LIODBLE, intent(in) :: C(M,Mlr), E(Mlr)
   LIODBLE, intent(in) :: Xexc(Ndim,Nstat), Eexc(Nstat)
   LIODBLE, intent(inout) :: Etot

   integer :: ii, jj
   LIODBLE :: Knr
   LIODBLE, allocatable :: Zvec(:), Xvec(:)
   LIODBLE, allocatable :: Xmat(:,:), Zmat(:,:), Zsym(:,:)
   LIODBLE, allocatable :: gammaWS(:,:), gammaXC(:,:), gammaCou(:,:)
   LIODBLE, allocatable :: gammaH(:,:), gammaT(:,:), gammaTot(:,:)
   LIODBLE, allocatable :: rhoG(:,:)

   if ( .not. TSH ) return
   if ( root == 0 ) return

!TODO: for the moment, this is the best place to put Energy
   Etot = Etot + Eexc(root)

   if ( root > Nstat ) then
      print*, "The root variable is bigger than nstates"
      print*, "Please set root <= nstates"
      stop
   endif
   
   print*, "TSH probabilities"
   allocate(Zvec(Ndim),Xvec(Ndim))
   Xvec = Xexc(:,root) / dsqrt(2.0d0)
   Zvec = Xvec / Eexc(root)

   ! Form Matrix in AO.
   allocate(Zmat(M,M),Xmat(M,M),gammaWS(natom,3))
   call VecToMat(Xvec,Xmat,C,Ndim,NCO,M,Mlr)
   call VecToMat(Zvec,Zmat,C,Ndim,NCO,M,Mlr)

   ! Obtain gamma WS: this calculates -gammaWS
   call gammaWS_calc(Xvec,Zvec,Zmat,E,C,gammaWS,NCO,M,Mlr,Ndim,Nvirt,&
                     natom)

   ! Obtain gamma Core
   allocate(Zsym(M,M),gammaH(natom,3))
   Zsym = Zmat + transpose(Zmat)
   call HVgradcalc(Zsym,gammaH,M,natom,.false.)

   ! Obtain gamma Coulomb
   allocate(rhoG(M,M)); call spunpack_rho('L',M,Pmat_vec,rhoG)
   allocate(gammaCou(natom,3)); gammaCou = 0.0d0
   call g2g_calcgammcou(rhoG,Zsym,gammaCou)
   deallocate(Zsym,rhoG)
   
   ! Obtain gamma XC
   allocate(Zsym(M,M),gammaXC(3,natom)); Zsym = 0.0d0; gammaXC = 0.0d0
   call g2g_calcgradxc(Zmat,Zsym,gammaXC,1)
   deallocate(Zsym)

   ! Obtain gamma T
   allocate(gammaT(natom,3)); gammaT = 0.0d0
   call intSG_Exc(gammaT,Xmat,natom,M)

   ! Obtain Total gamma
   allocate(gammaTot(natom,3))
   gammaTot = gammaWS + gammaH + gammaCou + transpose(gammaXC) - gammaT
!  print*, "gamma Tot"
!  do ii=1,natom
!     print*, ii,gammaTot(ii,1),gammaTot(ii,2),gammaTot(ii,3)
!  enddo
   deallocate(gammaWS,gammaH,gammaCou,gammaXC,gammaT)
   deallocate(Zvec,Xvec,Zmat,Xmat)

   ! Norm of NACMEs: When you perform Dynamic only.
   if (allocated(gamma_old)) then
      Knr = 0.0d0
      do ii=1,natom
      do jj=1,3
         Knr = Knr + gammaTot(ii,jj) / atom_mass(ii)
      enddo
      enddo
      print*, "NORM of NACVs", Knr*Knr*4.0d0

      call coef_propagator(gamma_old,nucvel,natom,Eexc(root))
      gamma_old = gammaTot
   endif
   
   deallocate(gammaTot)
end subroutine tsh_probabilities


subroutine coef_propagator(g,v,natom,dE)
use excited_data, only: dE_accum, lambda, tsh_time_dt, B_old, &
                        tsh_Jstate, tsh_Kstate, tsh_coef
   implicit none

   integer, intent(in) :: natom
   LIODBLE, intent(in) :: dE, g(natom,3)
   LIODBLE, intent(in) :: v(3,natom)

   integer :: ii
   LIODBLE :: Q, Gprob, factor, pop, number_random
   complex(kind=8) :: B, B_tot, B_abs, zero, B1, B2, pot, c_j, c_k
   complex(kind=8), allocatable :: Uprop(:,:)

   !  v = Nuclear Velocity [Bohr/au]
   !  g = Non-Adiabatic Coupling Vector [1/Bohr]
   !  Q = v X h [1/au]
   Q = 0.0d0
   do ii=1,natom
     Q = Q + v(1,ii) * g(ii,1)
     Q = Q + v(2,ii) * g(ii,2)
     Q = Q + v(3,ii) * g(ii,3)
   enddo

!  ========== UNITS CONVERTIONS (ha -> au^-1) ==============
!  1 ha = 6579664992709240.0 sec^-1
!  1 sec = 1e+15 femto
!  1 femto = 41.341347575 atomic time units(au)
!  dE(au^-1) = dE(ha) * (6579664992709240.0 / 1e+15) 
!             / 41.341347575 = 0.1591545844211451
!  dE(au^-1) = dE(ha) * 0.1591545844211451
!  =========================================================
   dE_accum = (dE_accum + dE) * 0.1591545844211451d0
   lambda = lambda + dE_accum * tsh_time_dt * 0.5d0

   B = Q * exp(cmplx(0.0d0,-lambda,8))
   B_tot = tsh_time_dt * 0.5d0 * ( B + B_old )

   ! Save old values
   B_old = B 
   dE_accum = dE
   
   zero = (0.0d0,0.0d0)
   B_abs = abs(B_tot)
   if ( .not. (abs(B_abs) > 0.0D0)) then
      B1 = zero
      B2 = zero
   else
      B1 = B_tot / B_abs
      B2 = conjg(B_tot) / B_abs
   endif
   allocate(Uprop(2,2))
   Uprop(1,1)=cos(B_abs)
   Uprop(1,2)=-B2*sin(B_abs)
   Uprop(2,1)=B1*sin(B_abs)
   Uprop(2,2)=cos(B_abs)

   ! Obtain Coeficcients
   tsh_coef = matmul(tsh_coef, Uprop)
   deallocate(Uprop)

   write(*,"(1X,A,1X,I2,A,I2)") "Transition:",tsh_Jstate," ->",tsh_Kstate
   print*, "poblacion1", real(abs(tsh_coef(1))**2.0d0)
   print*, "poblacion2", real(abs(tsh_coef(2))**2.0d0)

   ! Probability Calculate
   c_j = tsh_coef(tsh_Jstate)
   c_k = conjg(tsh_coef(tsh_Kstate))
   if (tsh_Jstate > tsh_Kstate) then
      pot = exp(-cmplx(0.0d0,0.0d0,8))
      Q = -Q
   else
      pot = exp(cmplx(0.0d0,0.0d0,8))
   endif
   factor = real(c_j*c_k*pot)
   Gprob  = factor*Q*(-2.0d0)

   if ( Gprob < 0.0d0 ) Gprob = 0.0d0
   if ( Gprob > 1.0d0 ) print*, "Probabilitie WRONG?"
   pop = real( abs(c_j)*abs(c_j) )

   if ( Gprob < 0.0d0 ) Gprob = 0.0d0
   if ( Gprob > 1.0d0 ) print*, "Probabilitie WRONG?"
   pop = real( abs(c_j)*abs(c_j) )

   Gprob = Gprob * tsh_time_dt / pop
   call random_number(number_random)
   print*, "probability", Gprob
   print*, "random number", number_random
end subroutine coef_propagator


subroutine gammaWS_calc(Xv,Zv,Zm,E,C,gamm,NCO,M,Mlr,Ndim,Nvirt,natom)
use garcha_mod  , only: ntatom, r, d
use faint_cpu   , only: intSG
use excited_data, only: fittExcited, Cocc, Cocc_trans, Coef_trans
use extern_functional_data, only: need_HF
   implicit none

   integer, intent(in) :: NCO, M, Mlr, Ndim, Nvirt, natom
   LIODBLE, intent(in) :: Xv(Ndim), Zv(Ndim), Zm(M,M)
   LIODBLE, intent(in) :: C(M,Mlr), E(Mlr)
   LIODBLE, intent(out) :: gamm(natom,3)

   integer :: ii, jj, NCOc, pos1, MM, ind
   LIODBLE, allocatable :: F2e(:,:), Fxc(:,:), Ftot(:,:)
   LIODBLE, allocatable :: HZIJ(:,:), scratch(:,:)
   LIODBLE, allocatable :: Wmat(:,:), WmatMO(:,:), Wtot(:)

!  CALCULATE TWO ELECTRON PART
   allocate(F2e(M,M))
   if ( .not. fittExcited ) then
      call g2g_timer_start("Fock 2e LR")
      call g2g_calculate2e(Zm,F2e,1)
      F2e = (F2e+transpose(F2e))
      call g2g_timer_stop("Fock 2e LR")
   elseif ( fittExcited .and. (.not. need_HF) ) then
      call g2g_timer_start("Fock 2e LR")
      call calc2eFITT(Zm,F2e,M)
      call g2g_timer_stop("Fock 2e LR")
   else
      print*, "Error in 2 Electron Repulsion Integrals"
      print*, "Check HF in the functional and fittExcited"
      stop
   endif

!  CALCULATE XC PART
   allocate(Fxc(M,M)); Fxc = 0.0d0
   call g2g_calculateg(Zm,Fxc,2)

!  TOTAL FOCK
   allocate(Ftot(M,M)); Ftot = F2e + Fxc + Fxc
   deallocate(F2e,Fxc)

!  CHANGE BASIS of FOCK. AO -> MO(OCC x OCC)
   allocate(scratch(M,NCO),HZIJ(NCO,NCO))
   call dgemm('N','N',M,NCO,M,1.0d0,Ftot,M,Cocc,M,0.0d0,scratch,M)
   call dgemm('N','N',NCO,NCO,M,1.0d0,Cocc_trans,NCO,scratch,M, &
              0.0d0,HZIJ,NCO)
   deallocate(scratch,Ftot)

   allocate(WmatMO(Mlr,Mlr)); WmatMO = 0.0d0
!  FORM OCC x OCC Block
   NCOc = NCO + 1
   do ii=1,NCO
   do jj=1,ii
      WmatMO(NCOc-ii,NCOc-jj) = HZIJ(NCOc-ii,NCOc-jj)
   enddo
   enddo
   do ii=1,NCO
      WmatMO(ii,ii) = WmatMO(ii,ii) * 0.5d0
   enddo
   deallocate(HZIJ)

!  FORM OCC x VIRT Block
   do ii=1,NCO
   do jj=1,Nvirt
      pos1 = (ii-1) * Nvirt + jj
      WmatMO(NCOc-ii,NCO+jj) = E(NCOc-ii) * Zv(pos1) + Xv(pos1)
   enddo
   enddo

   ! CHANGE BASIS of Wmat. MO -> AO
   allocate(scratch(M,Mlr),Wmat(M,M))
   call dgemm('N','N',M,Mlr,Mlr,1.0d0,C,M,WmatMO,Mlr,0.0d0,scratch,M)
   call dgemm('N','N',M,M,Mlr,1.0d0,scratch,M,Coef_trans,Mlr,0.0d0, &
              Wmat,M)
   deallocate(scratch,WmatMO)

   ! NACVs
   MM = M * (M + 1) / 2
   allocate(Wtot(MM))
   ind = 1
   do ii=1,M
      Wtot(ind) = Wmat(ii,ii)
      ind = ind + 1
      do jj=ii+1,M
         Wtot(ind) = Wmat(ii,jj) + Wmat(jj,ii)
         ind = ind + 1
      enddo
   enddo
   Wtot = (-1.0d0) * Wtot
   gamm = 0.0d0
   call intSG(gamm, Wtot, r, d, natom, ntatom)
   deallocate(Wtot,Wmat)
   gamm = 2.0d0 * gamm
end subroutine gammaWS_calc

