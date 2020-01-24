subroutine OscStr(X,Ene,Coef,OsSt,M,Mlr,NCO,Nvirt,Ndim,Nstat)
   implicit none

   integer, intent(in) :: M, Mlr, Ndim, Nstat, NCO, Nvirt
   double precision, intent(in) :: X(Ndim,Nstat),Ene(Nstat)
   double precision, intent(in) :: Coef(M,Mlr)
   double precision, intent(out) :: OsSt(Nstat)

   integer :: ii
   double precision, dimension(:,:), allocatable :: Tdip
   double precision, dimension(:,:), allocatable :: TdensAO,TdensMO

   allocate(Tdip(Nstat,3),TdensMO(Mlr,Mlr),TdensAO(M,M))

   do ii=1,Nstat
      !  FORM TRANSITION DENSITY IN MO BASIS
      TdensMO = 0.0d0
      call vecMOtomatMO(X,TdensMO,Mlr,NCO,Nvirt,Nstat,0,ii,Ndim)

      !  CHANGE BASIS MO -> AO
      call matMOtomatAO(tdensMO,tdensAO,Coef,M,Mlr,.true.)

      !  CALCULATE TRANSITION DIPOLE
      call TransDipole(TdensAO,Tdip(ii,:),M)
   enddo
   deallocate(TdensMO,TdensAO)

!  CALCULATE THE OSCILATOR STRENGHT
   call ObtainOsc(Tdip,Ene,OsSt,Nstat)
   deallocate(Tdip)
end subroutine OscStr
