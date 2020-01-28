subroutine Wcalculate(Zvec,Dif,Qvec,GxcAO,Vlr,C, &
                      dE,EneSCF,Wmat,Ndim,M,Mlr,NCO)
use garcha_mod  , only: PBE0
use excited_data, only: Cocc, Cocc_trans, Coef_trans, fittExcited
   implicit none

   integer, intent(in) :: M, Mlr, Ndim, NCO
   double precision, intent(in) :: dE
   double precision, intent(in) :: Zvec(Ndim), Vlr(Ndim), Qvec(Ndim) 
   double precision, intent(in) :: Dif(M,M), GxcAO(M,M)
   double precision, intent(in) :: C(M,Mlr), EneSCF(Mlr)
   double precision, intent(out) :: Wmat(M,M)

   integer :: ii, jj, kk
   integer :: Nvirt, pos1, pos2, NCOc
   double precision :: temp1, temp2
   double precision, allocatable :: F2e(:,:), Fxc(:,:), Ftot(:,:)
   double precision, allocatable :: scratch(:,:), HXIJ(:,:), GXCIJ(:,:)
   double precision, allocatable :: WmatMO(:,:)

   Nvirt = Mlr - NCO

!  CALCULATE TWO ELECTRON PART
   allocate(F2e(M,M))
   if ( .not. fittExcited ) then
      call g2g_timer_start("Fock 2e LR")
      call g2g_calculate2e(Dif,F2e,1)
      F2e = (F2e+transpose(F2e))
      call g2g_timer_stop("Fock 2e LR")
   elseif ( fittExcited .and. (.not. PBE0) ) then
      call g2g_timer_start("Fock 2e LR")
      call calc2eFITT(Dif,F2e,M)
      call g2g_timer_stop("Fock 2e LR")
   else
      print*, "Error in 2 Electron Repulsion Integrals"
      print*, "Check PBE0 and fittExcited"
      stop
   endif

!  CALCULATE XC PART
   allocate(Fxc(M,M)); Fxc = 0.0d0
   call g2g_calculateg(Dif,Fxc,2)

!  TOTAL FOCK
   allocate(Ftot(M,M)); Ftot = F2e + Fxc + Fxc
   deallocate(F2e,Fxc)

!  CHANGE BASIS of FOCK. AO -> MO(OCC x OCC)
   allocate(scratch(M,NCO),HXIJ(NCO,NCO))
   call dgemm('N','N',M,NCO,M,1.0d0,Ftot,M,Cocc,M,0.0d0,scratch,M)
   call dgemm('N','N',NCO,NCO,M,1.0d0,Cocc_trans,NCO,scratch,M,0.0d0,HXIJ,NCO)
   deallocate(scratch,Ftot)

!  CHANGE BASIS of GXC. AO -> MO(OCC x OCC)
   allocate(scratch(M,NCO),GXCIJ(NCO,NCO))
   call dgemm('N','N',M,NCO,M,1.0d0,GxcAO,M,Cocc,M,0.0d0,scratch,M)
   call dgemm('N','N',NCO,NCO,M,1.0d0,Cocc_trans,NCO,scratch,M,0.0d0,GXCIJ,NCO)
   deallocate(scratch)

!  FORM ENERGY WEIGTHED DIFFERENCE DENSITY MATRIX
   allocate(WmatMO(Mlr,Mlr)); WmatMO = 0.0d0
   ! BLOCK OCC x OCC
   NCOc = NCO + 1
   do ii=1,NCO
   do jj=1,ii
      temp1 = 0.0D0
      temp2 = 0.0D0
      do kk=1,Nvirt
         pos1 = (ii-1)*Nvirt+kk !ia
         pos2 = (jj-1)*Nvirt+kk !ja
         temp1 = temp1 + Vlr(pos1)*Vlr(pos2)*2.0D0
         temp2 = temp2 + EneSCF(NCO+kk)*(Vlr(pos1)*Vlr(pos2)*2.0D0)
      enddo
      WmatMO(NCOc-ii,NCOc-jj) = dE * temp1 - temp2 &
                  + HXIJ(NCOc-ii,NCOc-jj) + 2.0D0 * GXCIJ(NCOc-ii,NCOc-jj)
      WmatMO(NCOc-jj,NCOc-ii) = WmatMO(NCOc-ii,NCOc-jj)
   enddo
   enddo
   deallocate(HXIJ,GXCIJ)

   ! BLOCK VIR x VIR
   NCOc = NCO + 1
   do ii=1,Nvirt
   do jj=1,ii
      temp1 = 0.0D0
      temp2 = 0.0D0
      do kk=1,NCO
         pos1 = (kk-1)*Nvirt+ii !ki
         pos2 = (kk-1)*Nvirt+jj !kj
         temp1 = temp1 + Vlr(pos1)*Vlr(pos2)*2.0D0
         temp2 = temp2 + EneSCF(NCOc-kk)*(Vlr(pos1)*Vlr(pos2)*2.0D0)
      enddo
      WmatMO(NCO+ii,NCO+jj) = dE * temp1 + temp2
      WmatMO(NCO+jj,NCO+ii) = WmatMO(NCO+ii,NCO+jj)
   enddo
   enddo

   ! BLOCK OCC x VIR
   NCOc = NCO + 1
   do ii=1,NCO
   do jj=1,Nvirt
      pos1 = (ii-1)*Nvirt+jj !ij
      WmatMO(NCO+jj,NCOc-ii) = Qvec(pos1) + EneSCF(NCOc-ii)*Zvec(pos1)
      WmatMO(NCOc-ii,NCO+jj) = WmatMO(NCO+jj,NCOc-ii)
   enddo
   enddo

   ! CHANGE BASIS of Wmat. MO -> AO
   allocate(scratch(M,Mlr))
   call dgemm('N','N',M,Mlr,Mlr,1.0d0,C,M,WmatMO,Mlr,0.0d0,scratch,M)
   call dgemm('N','N',M,M,Mlr,1.0d0,scratch,M,Coef_trans,Mlr,0.0d0,Wmat,M)
   deallocate(scratch,WmatMO)
end subroutine Wcalculate
