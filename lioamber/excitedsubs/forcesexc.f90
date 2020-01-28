subroutine forcesexc(rhoTot,DiffExc,Zvec,Xmat,Qvec,GxcAO,Xexc,Eexc, &
                     Cscf,Escf,M,Mlr,Ndim,NCO,Nstat)
use garcha_mod, only: natom
use excited_data, only: excited_forces, root
   implicit none

   integer, intent(in) :: M, Mlr, Ndim, NCO, Nstat
   double precision, intent(in) :: Cscf(M,Mlr), Escf(Mlr)
   double precision, intent(in) :: Xexc(Ndim,Nstat), Eexc(Nstat)
   double precision, intent(in) :: rhoTot(M,M), DiffExc(M,M), Xmat(M,M)
   double precision, intent(in) :: Zvec(Ndim), Qvec(Ndim), GxcAO(M,M)

   integer :: ii ! borrar
   double precision :: dE
   double precision, allocatable :: Xlr(:)
   double precision, allocatable :: fXC(:,:), fWS(:,:), fHV(:,:)
   double precision, allocatable :: fCOU(:,:), fTOT(:,:)

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



end subroutine forcesexc
