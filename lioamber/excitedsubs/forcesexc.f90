subroutine forcesexc(rhoTot,DiffExc,Zvec,Xmat,Qvec,GxcAO,Xexc,Eexc, &
                     Cscf,Escf,M,Mlr,Ndim,NCO,Nstat)
use garcha_mod, only: natom
use excited_data, only: excited_forces, root, for_exc
use extern_functional_data, only: HF, HF_fac
   implicit none

   integer, intent(in) :: M, Mlr, Ndim, NCO, Nstat
   LIODBLE, intent(in) :: Cscf(M,Mlr), Escf(Mlr)
   LIODBLE, intent(in) :: Xexc(Ndim,Nstat), Eexc(Nstat)
   LIODBLE, intent(in) :: rhoTot(M,M), DiffExc(M,M), Xmat(M,M)
   LIODBLE, intent(in) :: Zvec(Ndim), Qvec(Ndim), GxcAO(M,M)

   LIODBLE :: dE
   LIODBLE, allocatable :: Xlr(:), rhoG(:,:)
   LIODBLE, allocatable :: fXC(:,:), fWS(:,:), fHV(:,:)
   LIODBLE, allocatable :: fCOU(:,:), fEE(:,:)

   if ( .not.excited_forces ) return
   if ( root == 0 ) return
   print*, "*Excited State Gradients"

   allocate(Xlr(Ndim)); Xlr = Xexc(:,root) / dsqrt(2.0d0); dE = Eexc(root)

   ! Overlap Gradients
   allocate(fWS(natom,3))
   call WSgradcalc(Zvec,DiffExc,Qvec,GxcAO,Xlr,Cscf,&
                   dE,Escf,fWS,Ndim,M,Mlr,NCO,natom)

   ! Core and Nuclear Gradients
   allocate(fHV(natom,3))
   call HVgradcalc(rhoTot,fHV,M,natom,.true.)

   ! Coulomb Gradients
   allocate(fCOU(natom,3))
   call COUgradcalc(rhoTot,DiffExc,Xmat,fCOU,M,natom)

   ! XC Gradients
   allocate(fXC(3,natom)); fXC = 0.0d0
   call g2g_calcgradxc(DiffExc,Xmat,fXC,0)

   ! Exact Exchange Gradients
   allocate(fEE(natom,3)); fEE = 0.0d0
   if ( HF /= 0 ) then
      allocate(rhoG(M,M)); rhoG = rhoTot - ( DiffExc + transpose(DiffExc) )
      call g2g_exacgrad_excited(rhoG,DiffExc,Xmat,fEE)

      fEE = HF_fac * fEE
      deallocate(rhoG)
   endif

   ! Total Gradients
   if(allocated(for_exc)) deallocate(for_exc)
   allocate(for_exc(natom,3))
   for_exc = fWS + fHV + fCOU + transpose(fXC) + fEE
   deallocate(fWS,fHV,fCOU,fXC,fEE,Xlr)
!  print*, "total f"
!  do ii=1,natom
!     print*, ii, for_exc(ii,1), for_exc(ii,2), for_exc(ii,3)
!  enddo

end subroutine forcesexc
