subroutine Zvector(C,Ene,X,TundAO,Xmat,Zvec,Qvec,Gxc,NCO,M,Mlr,Ndim,Nvirt)
use excited_data, only: fittExcited
use extern_functional_data, only: libint_inited
   implicit none

   integer, intent(in) :: NCO, M, Mlr, Ndim, Nvirt
   LIODBLE, intent(in) :: C(M,Mlr), Ene(Mlr), TundAO(M,M)
   LIODBLE, intent(in) :: Xmat(M,M), X(Ndim)
   LIODBLE, intent(inout) :: Zvec(Ndim), Qvec(Ndim), Gxc(M,M)

   logical :: is_calc 
   LIODBLE, allocatable :: FX(:,:), FT(:,:), FXAB(:,:)
   LIODBLE, allocatable :: F2e(:,:,:), PA(:,:,:), Rvec(:)
   LIODBLE, allocatable :: FXIJ(:,:), FTIA(:,:), GXCIA(:,:)

   print*, ""
   print*,"======================================="
   print*,"            Z-VECTOR METHOD"
   print*,"======================================="
   print*, ""
   
!  THIRD DERIVATIVE CALCULATED
   Gxc = 0.0d0
   call g2g_calculateg(Xmat,Gxc,3)

!  SECOND DERIVATIVE CALCULATED
   allocate(FX(M,M),FT(M,M)); FX = 0.0d0; FT = 0.0d0
   call g2g_calculateg(Xmat,FX,2)
   call g2g_calculateg(TundAO,FT,2)

!  2-ELECTRON INTEGRALS CALCULATED
   allocate(PA(M,M,2),F2e(M,M,2))
   PA(:,:,1) = Xmat; PA(:,:,2) = TundAO

   is_calc = .false.
   if ( .not. fittExcited ) then
      call g2g_timer_start("Fock 2e LR")
      call g2g_calculate2e(PA,F2e,2)
      F2e = 2.0d0 * F2e
      call g2g_timer_stop("Fock 2e LR")
      is_calc = .true.
   endif

   if ( .not. is_calc ) then
      if ( fittExcited .and. (.not. libint_inited) ) then
         call g2g_timer_start("Fock 2e LR")
         call calc2eFITT(PA(:,:,1),F2e(:,:,1),M)
         call calc2eFITT(PA(:,:,2),F2e(:,:,2),M)
         call g2g_timer_stop("Fock 2e LR")
      else
         print*, "Error in 2 Electron Repulsion Integrals"
         print*, "Check HF in the functional and fittExcited"
         stop
      endif
   endif
   deallocate(PA)

!  Total Focks
   FX = F2e(:,:,1) + 2.0D0 * FX
   FT = F2e(:,:,2) + 2.0D0 * FT
   deallocate(F2e)

!  Change basis of all Fock Matrices
   allocate(FXAB(Nvirt,Nvirt),FXIJ(NCO,NCO))
   allocate(FTIA(NCO,Nvirt),GXCIA(NCO,Nvirt))
   call ChangeBasisF(FX,FT,Gxc,FXAB,FXIJ,FTIA,GXCIA,M,Nvirt,NCO)
   deallocate(FX,FT)

!  Vector R calculate ( R = A * X )
   allocate(Rvec(Ndim))
   call RCalculate(FXAB,FXIJ,FTIA,GXCIA,X,Rvec,Qvec,NCO,Nvirt,Ndim)

!  Solve equation AX=R with PCG Method
   call PCG_solve(Rvec,C,Ene,Zvec,M,Mlr,NCO,Nvirt,Ndim)
   deallocate(Rvec,FXAB,FXIJ,FTIA,GXCIA)
end subroutine Zvector
