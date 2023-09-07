subroutine OscStr(X,Ene,Coef,OsSt,M,Mlr,NCO,Nvirt,Ndim,Nstat)
use excited_data, only: second_LR, Tdip_save, Coef_trans
   implicit none

   integer, intent(in) :: M, Mlr, Ndim, Nstat, NCO, Nvirt
   LIODBLE, intent(in) :: X(Ndim,Nstat),Ene(Nstat)
   LIODBLE, intent(in) :: Coef(M,Mlr)
   LIODBLE, intent(out) :: OsSt(Nstat)

   integer :: ii
   LIODBLE, dimension(:,:), allocatable :: Tdip
   LIODBLE, dimension(:,:), allocatable :: TdensAO,TdensMO

   allocate(Tdip(Nstat,3),TdensMO(Mlr,Mlr),TdensAO(M,M))

   do ii=1,Nstat
      !  FORM TRANSITION DENSITY IN MO BASIS
      TdensMO = 0.0d0
      call vecMOtomatMO(X,TdensMO,Mlr,NCO,Nvirt,Nstat,0,ii,Ndim)

      !  CHANGE BASIS MO -> AO
      call matMOtomatAO(tdensMO,tdensAO,Coef,Coef_trans,M,Mlr,.true.)

      !  CALCULATE TRANSITION DIPOLE
      call TransDipole(TdensAO,Tdip(ii,:),M)
   enddo
   deallocate(TdensMO,TdensAO)

!  CALCULATE THE OSCILATOR STRENGHT
   call ObtainOsc(Tdip,Ene,OsSt,Nstat,0)

!  Saved Transition Dipole
   if ( second_LR ) then
      if(allocated(Tdip_save)) deallocate(Tdip_save)
      allocate(Tdip_save(Nstat,3))
      Tdip_save = Tdip
   endif

   deallocate(Tdip)
end subroutine OscStr

subroutine open_OscStr(Xa,Xb,Ene,CoefA,CoefB,OsSt,M,Mlr,NCOa,NCOb,Nvirta,Nvirtb,Ndim,Nstat)
use excited_data, only: Coef_trans, Coef_transB
   implicit none

   integer, intent(in) :: M, Mlr, Ndim, Nstat, NCOa, NCOb, Nvirta, Nvirtb
   LIODBLE, intent(in) :: Xa(Ndim,Nstat), Xb(Ndim,Nstat), Ene(Nstat)
   LIODBLE, intent(in) :: CoefA(M,Mlr), CoefB(M,Mlr)
   LIODBLE, intent(out) :: OsSt(Nstat)

   integer :: ii
   LIODBLE, dimension(:,:), allocatable :: Tdip
   LIODBLE, dimension(:,:), allocatable :: TdensAOa,TdensMOa,TdensAOb,TdensMOb

   allocate(Tdip(Nstat,3),TdensMOa(Mlr,Mlr),TdensAOa(M,M))
   allocate(TdensMOb(Mlr,Mlr),TdensAOb(M,M))

   do ii=1,Nstat
      !  alpha: Form transition density in MO and then -> AO
      TdensMOa = 0.0d0
      call vecMOtomatMO(Xa,TdensMOa,Mlr,NCOa,Nvirta,Nstat,0,ii,Ndim)
      call matMOtomatAO(tdensMOa,tdensAOa,CoefA,Coef_trans,M,Mlr,.true.)
      !  beta: Form transition density in MO and then -> AO
      TdensMOb = 0.0d0
      call vecMOtomatMO(Xb,TdensMOb,Mlr,NCOb,Nvirtb,Nstat,0,ii,Ndim)
      call matMOtomatAO(tdensMOb,tdensAOb,CoefB,Coef_transB,M,Mlr,.true.)

      ! Total:
      TdensAOa = TdensAOa + TdensAOb
      call TransDipole(TdensAOa,Tdip(ii,:),M)
   enddo
   deallocate(TdensMOa,TdensAOa,TdensMOb,TdensAOb)

!  CALCULATE THE OSCILATOR STRENGHT
   call ObtainOsc(Tdip,Ene,OsSt,Nstat,0)
   OsST = OsST * 0.5d0
   deallocate(Tdip)
end subroutine open_OscStr
