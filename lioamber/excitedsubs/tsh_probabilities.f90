subroutine tsh_probabilities(C,E,Xexc,Eexc,NCO,M,Mlr,Ndim,Nvirt,Nstat)
use garcha_mod  , only: natom, Pmat_vec
use excited_data, only: TSH, root
   implicit none

   integer, intent(in) :: NCO, M, Mlr, Ndim, Nvirt, Nstat
   double precision, intent(in) :: C(M,Mlr), E(Mlr)
   double precision, intent(in) :: Xexc(Ndim,Nstat), Eexc(Nstat)

   integer :: ii, jj
   double precision, allocatable :: Zvec(:), Xvec(:)
   double precision, allocatable :: Xmat(:,:), Zmat(:,:), Zsym(:,:)
   double precision, allocatable :: gammaWS(:,:), gammaXC(:,:), gammaCou(:,:)
   double precision, allocatable :: gammaH(:,:), gammaT(:,:), gammaTot(:,:)
   double precision, allocatable :: rhoG(:,:)

   if ( .not. TSH ) return
   if ( root == 0 ) return

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
   print*, "Wgam"
   do ii=1,natom
      print*, ii,gammaWS(ii,1),gammaWS(ii,2),gammaWS(ii,3)
   enddo

   ! Obtain gamma Core
   allocate(Zsym(M,M),gammaH(natom,3))
   Zsym = Zmat + transpose(Zmat)
   call HVgradcalc(Zsym,gammaH,M,natom,.false.)
   print*, "Hgam"
   do ii=1,natom
      print*, ii,gammaH(ii,1),gammaH(ii,2),gammaH(ii,3)
   enddo

   ! Obtain gamma Coulomb
   allocate(rhoG(M,M)); call spunpack_rho('L',M,Pmat_vec,rhoG)
   allocate(gammaCou(natom,3)); gammaCou = 0.0d0
   call g2g_calcgammcou(rhoG,Zsym,gammaCou)
   deallocate(Zsym,rhoG)
   print*, "COUgam"
   do ii=1,natom
      print*, ii,gammaCou(ii,1),gammaCou(ii,2),gammaCou(ii,3)
   enddo
   
   ! Obtain gamma XC
   allocate(Zsym(M,M),gammaXC(3,natom)); Zsym = 0.0d0; gammaXC = 0.0d0
   call g2g_calcgradxc(Zmat,Zsym,gammaXC,1)
   deallocate(Zsym)
   print*, "XCgam"
   do ii=1,natom
      print*, ii,gammaXC(1,ii),gammaXC(2,ii),gammaXC(3,ii)
   enddo

   ! Obtain gamma T
   allocate(gammaT(natom,3)); gammaT = 0.0d0
   call intSG_Exc(gammaT,Xmat,natom,M)
   print*, "Tgam"
   do ii=1,natom
      print*, ii,gammaT(ii,1),gammaT(ii,2),gammaT(ii,3)
   enddo

   ! Obtain Total gamma
   allocate(gammaTot(natom,3))
   gammaTot = gammaWS + gammaH + gammaCou + transpose(gammaXC) - gammaT
   deallocate(gammaWS,gammaH,gammaCou,gammaXC,gammaT)
   deallocate(Zvec,Xvec,Zmat,Xmat)


   


   



   

end subroutine tsh_probabilities


subroutine gammaWS_calc(Xv,Zv,Zm,E,C,gamm,NCO,M,Mlr,Ndim,Nvirt,natom)
use garcha_mod  , only: ntatom, r, d, PBE0
use faint_cpu   , only: intSG
use excited_data, only: fittExcited, Cocc, Cocc_trans, Coef_trans
   implicit none

   integer, intent(in) :: NCO, M, Mlr, Ndim, Nvirt, natom
   double precision, intent(in) :: Xv(Ndim), Zv(Ndim), Zm(M,M)
   double precision, intent(in) :: C(M,Mlr), E(Mlr)
   double precision, intent(out) :: gamm(natom,3)

   integer :: ii, jj, NCOc, pos1, MM, ind
   double precision, allocatable :: F2e(:,:), Fxc(:,:), Ftot(:,:)
   double precision, allocatable :: HZIJ(:,:), scratch(:,:)
   double precision, allocatable :: Wmat(:,:), WmatMO(:,:), Wtot(:)

!  CALCULATE TWO ELECTRON PART
   allocate(F2e(M,M))
   if ( .not. fittExcited ) then
      call g2g_timer_start("Fock 2e LR")
      call g2g_calculate2e(Zm,F2e,1)
      F2e = (F2e+transpose(F2e))
      call g2g_timer_stop("Fock 2e LR")
   elseif ( fittExcited .and. (.not. PBE0) ) then
      call g2g_timer_start("Fock 2e LR")
      call calc2eFITT(Zm,F2e,M)
      call g2g_timer_stop("Fock 2e LR")
   else
      print*, "Error in 2 Electron Repulsion Integrals"
      print*, "Check PBE0 and fittExcited"
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

