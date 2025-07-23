subroutine OscStr(X,Ene,Coef,OsSt,M,Mlr,NCO,Nvirt,Ndim,Nstat)
use excited_data, only: second_LR, Tdip_save
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
      call matMOtomatAO(tdensMO,tdensAO,Coef,M,Mlr,.true.)

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
