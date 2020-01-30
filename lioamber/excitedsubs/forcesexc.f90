subroutine forcesexc(rhoTot,DiffExc,Zvec,Xmat,Qvec,GxcAO,Xexc,Eexc, &
                     Cscf,Escf,M,Mlr,Ndim,NCO,Nstat)
use garcha_mod, only: natom, PBE0
use excited_data, only: excited_forces, root
   implicit none

   integer, intent(in) :: M, Mlr, Ndim, NCO, Nstat
   double precision, intent(in) :: Cscf(M,Mlr), Escf(Mlr)
   double precision, intent(in) :: Xexc(Ndim,Nstat), Eexc(Nstat)
   double precision, intent(in) :: rhoTot(M,M), DiffExc(M,M), Xmat(M,M)
   double precision, intent(in) :: Zvec(Ndim), Qvec(Ndim), GxcAO(M,M)

   integer :: ii
   double precision :: dE
   double precision, allocatable :: Xlr(:), rhoG(:,:)
   double precision, allocatable :: fXC(:,:), fWS(:,:), fHV(:,:)
   double precision, allocatable :: fCOU(:,:), fTOT(:,:), fEE(:,:)

   if ( .not.excited_forces ) return
   if ( root == 0 ) return
   print*, "Excited State Gradients"

   allocate(Xlr(Ndim)); Xlr = Xexc(:,root) / dsqrt(2.0d0); dE = Eexc(root)

   ! Overlap Gradients
   allocate(fWS(natom,3))
   call WSgradcalc(Zvec,DiffExc,Qvec,GxcAO,Xlr,Cscf,&
                   dE,Escf,fWS,Ndim,M,Mlr,NCO,natom)

   ! Core and Nuclear Gradients
   allocate(fHV(natom,3))
   call HVgradcalc(rhoTot,fHV,M,natom) !TODO: posible bug cuando llama a int1G para GPU

   ! Coulomb Gradients
   allocate(fCOU(natom,3))
   call COUgradcalc(rhoTot,DiffExc,Xmat,fCOU,M,natom)

   ! XC Gradients
   allocate(fXC(3,natom)); fXC = 0.0d0
   call g2g_calcgradxc(DiffExc,Xmat,fXC)

   ! Exact Exchange Gradients
   allocate(fEE(natom,3)); fEE = 0.0d0
   if ( PBE0 ) then
      allocate(rhoG(M,M)); rhoG = rhoTot - ( DiffExc + transpose(DiffExc) )
      call g2g_exacgrad_excited(rhoG,DiffExc,Xmat,fEE)
      fEE = 0.25d0 * fEE
      deallocate(rhoG)
   endif

   ! Total Gradients
   allocate(fTOT(natom,3))
   fTOT = fWS + fHV + fCOU + transpose(fXC) + fEE
   deallocate(fWS,fHV,fCOU,fXC,fEE,Xlr)
   print*, "total f"
   do ii=1,natom
      print*, ii, fTOT(ii,1), fTOT(ii,2), fTOT(ii,3)
   enddo

end subroutine forcesexc
