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

   double precision :: dE
   double precision, allocatable :: Xlr(:), Wexc(:,:)

   if ( .not.excited_forces ) return
   if ( root == 0 ) return
   print*, "Excited State Gradients"

   allocate(Xlr(Ndim)); Xlr = Xexc(:,root) / dsqrt(2.0d0); dE = Eexc(root)
   allocate(Wexc(M,M)); Wexc = 0.0d0
   call Wcalculate(Zvec,DiffExc,Qvec,GxcAO,Xlr,Cscf,&
                   dE,Escf,Wexc,Ndim,M,Mlr,NCO)




end subroutine forcesexc
