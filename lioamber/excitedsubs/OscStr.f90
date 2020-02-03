subroutine OscStr(X,Ene,Coef,OsSt,M,NCO,Nvirt,Ndim,Nstat)
   implicit none

   integer, intent(in) :: M, Ndim, Nstat, NCO, Nvirt
   double precision, intent(in) :: X(Ndim,Nstat),Ene(Nstat)
   double precision, intent(in) :: Coef(M,M)
   double precision, intent(out) :: OsSt(Nstat)

   integer :: ii
   double precision, dimension(:,:), allocatable :: Tdip
   double precision, dimension(:,:), allocatable :: TdensAO,TdensMO

   allocate(Tdip(Nstat,3),TdensMO(M,M),TdensAO(M,M))

   do ii=1,Nstat
      !  FORM TRANSITION DENSITY IN MO BASIS
      TdensMO = 0.0d0
      call vecMOtomatMO(X,TdensMO,M,NCO,Nvirt,Nstat,0,ii,Ndim)

      !  CHANGE BASIS MO -> AO
      call matMOtomatAO(tdensMO,tdensAO,Coef,M,.true.)

      !  CALCULATE TRANSITION DIPOLE
      call TransDipole(TdensAO,Tdip(ii,:),M)
   enddo
   deallocate(TdensMO,TdensAO)

!  CALCULATE THE OSCILATOR STRENGHT
   call ObtainOsc(Tdip,Ene,OsSt,Nstat)
   deallocate(Tdip)
end subroutine OscStr
