module lr_data
!Nvirt = number of virtual molecular orbitals
!dim = dimension of matrix Linear Response (Nvirt x NCO)
!cbas,cbasx = Needed for libint (Temporary)
!nstates = number of excited states to calculate
!root = excited state chosen for optimization
   implicit none

   logical :: lresp = .false.
   integer :: Nvirt, dim
   integer :: nstates = 3
   integer :: root = 0
#ifdef TD_SIMPLE
   real, dimension(:,:), allocatable :: eigenvec, cbas, cbasx
   real, dimension(:), allocatable :: eigenval
#else
   real*8, dimension(:,:), allocatable :: eigenvec, cbas, cbasx
   real*8, dimension(:), allocatable :: eigenval
#endif
end module lr_data

module lrtddft
! This is a main subroutine for Linear Response calculate.
! This module performs a Tamm-Dancoff Aproximation (TDA)

contains

   subroutine linear_response(MatCoef,VecEne)
   use lr_data, only: Nvirt,dim,nstates,eigenval,&
                      eigenvec,cbas,root
   use garcha_mod, only: NCO, M, c, a

   implicit none
#ifdef TD_SIMPLE
   real, intent(in) :: MatCoef(M,M)
   real, intent(in) :: VecEne(M)
   real, dimension(:,:) allocatable :: Kfxc, Kc, A_mat
   real, dimension(:,:,:,:), allocatable :: Kc_Int, Kxc_Int
#else
   real*8, intent(in) :: MatCoef(M,M)
   real*8, intent(in) :: VecEne(M)
   real*8, dimension(:,:) allocatable :: Kfxc, Kc, A_mat
   real*8, dimension(:,:,:,:), allocatable :: Kc_Int, Kxc_Int
#endif
   
   call g2g_timer_start('LINEAR RESPONSE')

   Nvirt = M - NCO
   dim = Nvirt * NCO

   allocate(Kc(dim,dim),Kfxc(dim,dim),A_mat(dim,dim))
   allocate(Kxc_Int(M,M,M,M),Kc_Int(M,M,M,M))

   Kfxc = 0.0D0
   Kc = 0.0D0
   Kxc_Int = 0.0D0
   Kc_Int = 0.0D0

   print*, ""
   print*,"#######################################"
   print*,"         LINEAR RESPONSE - TDA"
   print*,"#######################################"
   print*, ""

   call g2g_timer_start('g2g_LinearResponse')
   call g2g_linear_response(MatCoef,Kfxc,Kc,Kxc_Int,Kc_Int,cbas,dim)
   call g2g_timer_stop('g2g_LinearResponse')
   deallocate(Kxc_Int,Kc_Int)

   call g2g_timer_start('ObtainMatrix')
   call ObtainAmatrix(A_mat,Kc,Kfxc,VecEne,dim,M,NCO)
   call g2g_timer_stop('ObtainMatrix')
   deallocate(Kc,Kfxc)

   call g2g_timer_start('Davidson')
   call davidson(A_mat,dim,eigenval,eigenvec,nstates)
   call g2g_timer_stop('Davidson')
   deallocate(A_mat)

   call PrintResults(eigenvec,eigenval,dim,nstates)

   if(root > 0) then
     call UnDensExc(eigenvec(:,root),MatCoef,dim)
   endif
   deallocate(eigenvec,eigenval)

   call g2g_timer_stop('LINEAR RESPONSE')
   end subroutine linear_response

   subroutine ObtainAmatrix(A,Kcou,Kxc,VecE,N,M,NCO)
   implicit none

   integer, intent(in) :: N, M, NCO
#ifdef TD_SIMPLE
   real, intent(in) :: Kcou(N,N),Kxc(N,N),VecE(M)
   real, intent(out) :: A(N,N)
#else
   real*8, intent(in) :: Kcou(N,N),Kxc(N,N),VecE(M)
   real*8, intent(out) :: A(N,N)
#endif

   integer :: k,l,i,j !counters

   k = NCO
   l = NCO + 1
   do i=1,N
     A(i,i) = Kcou(i,i) + Kxc(i,i) + VecE(l) - VecE(k)
     l = l + 1
     if (l > M) then
       l = NCO + 1
       k = k - 1
     endif
   enddo
   do i=1,N
     do j=i+1,N
         A(i,j) = Kcou(i,j) + Kxc(i,j)
         A(j,i) = A(i,j)
     enddo
   enddo
   end subroutine ObtainAmatrix

   subroutine davidson(A,N,eigval,eigvec,Nstat)
   ! This subroutine not performs Davidson Diagonalization, instead performs
   ! DQDS algorithm
   ! TODO: Direct Davidson Iterative Diagonalization
   implicit none

   integer, intent(in) :: N, Nstat
#ifdef TD_SIMPLE
   real, intent(inout) :: A(N,N)
   real, dimension(:,:), allocatable, intent(out) :: eigvec
   real, dimension(:), allocatable, intent(out) :: eigval
#else
   real*8, intent(inout) :: A(N,N)
   real*8, dimension(:,:), allocatable, intent(out) :: eigvec
   real*8, dimension(:), allocatable, intent(out) :: eigval
#endif

   integer :: IL, IU, num, LWORK, LIWORK, info
   integer, dimension(:), allocatable :: IWORK, ISUPPZ
#ifdef TD_SIMPLE
   real :: VL = 1.0 !NOT REFERENCE
   real :: VU = 20.0 !NOT REFERENCE
   real :: ABSTOL = 1.0E-10
   real, dimension(:), allocatable :: WORK
#else
   real*8 :: VL = 1.0 !NOT REFERENCE
   real*8 :: VU = 20.0 !NOT REFERENCE
   real*8 :: ABSTOL = 1.0E-10
   real*8, dimension(:), allocatable :: WORK
#endif

   IL = 1
   IU = Nstat
   num = IU - IL + 1
   LWORK = -1
   LIWORK = -1

   allocate(ISUPPZ(2*N),WORK(1),IWORK(1),eigvec(N,Nstat),eigval(N))
   
   call DSYEVR('V','I','U',N,A,N,VL,VU,IL,IU,ABSTOL,num,eigval,&
               eigvec,N,ISUPPZ,WORK,LWORK,IWORK,LIWORK,info)
   if(info /= 0) then
     print*, "ERROR IN DSYEVR 1"
     stop
   endif

   LWORK = WORK(1)
   LIWORK = IWORK(1)
   deallocate(WORK,IWORK)
   allocate(WORK(LWORK),IWORK(LIWORK))
   
   call DSYEVR('V','I','U',N,A,N,VL,VU,IL,IU,ABSTOL,num,eigval,&
               eigvec,N,ISUPPZ,WORK,LWORK,IWORK,LIWORK,info)
   if(info /= 0) then
     print*, "ERROR IN DSYEVR 2"
     stop
   endif

   deallocate(WORK,IWORK,ISUPPZ)
   end subroutine davidson
  
   subroutine PrintResults(vec,val,N,nstat)
   use garcha_mod, only: M,NCO
   implicit none

   integer, intent(in) :: N, nstat
#ifdef TD_SIMPLE
   real, intent(in) :: vec(N,nstat),val(nstat)
#else
   real*8, intent(in) :: vec(N,nstat),val(nstat)
#endif

   integer :: i,j,from,to

   from = NCO
   to = NCO + 1

   do j=1, nstat
   write(*,100) j,val(j)
   do i=1, N
      if ( abs(vec(i,j)) > 0.1D0 ) then
         write(*,101) from, to, vec(i,j)
      endif
      to = to + 1
      if ( to == M+1 ) then
          from = from - 1
          to = NCO + 1
      endif
   enddo
      print*, " "
      from = NCO
      to = NCO + 1
   enddo

   100 FORMAT(1X,"STATE",I2,3X,"ENERGY=",F14.7," Hartree")
   101 FORMAT(6X,I2," -> ",I2,2X,F14.7)
   end subroutine PrintResults

   subroutine UnDensExc(vec,Coef,N)
!The excited state density matrix is constructed as:
!  PEat = PFat + dPat
!PEat = excited state density matrix in AO basis
!PFat = ground state density matrix in AO basis
!dPat = diference density matrix in AO basis
   use garcha_mod, only: M, NCO, RMM
   implicit none

   integer, intent(in) :: N
#ifdef TD_SIMPLE
   real, intent(in) :: vec(N),Coef(M,M)
   real, dimension(:,:), allocatable :: Pij, Pab, dPat, X, PEat, PFat
   real :: temp
#else
   real*8, intent(in) :: vec(N),Coef(M,M)
   real*8, dimension(:,:), allocatable :: Pij, Pab, dPat, X, PEat, PFat
   real*8 :: temp
#endif

   integer :: i,j,a,b,row,col,NCOc,Nij,Nab,Nvirt,M2

   Nvirt = M - NCO
   allocate(Pij(NCO,NCO),Pab(Nvirt,Nvirt),X(NCO,Nvirt))
   Pij = 0.0D0
   Pab = 0.0D0
   X = 0.0D0

   do row=0,NCO-1
   do col=1,Nvirt
     X(row+1,col) = vec(row*Nvirt+col)
   enddo
   enddo

!Form block Diference Density Matrix occ-occ in MO basis
   temp = 0.0D0
   do i=1,NCO
   do j=1,NCO
     do a=1,Nvirt
       temp = temp + X(i,a)*X(j,a)
     enddo
     Pij(i,j) = temp * (-1.0D0)
     temp = 0.0D0
   enddo
   enddo
   
!Form block Diference Density Matrix nvirt-nvirt in MO basis
   temp = 0.0D0
   do a=1,Nvirt
   do b=1,Nvirt
     do i=1,NCO
       temp = temp + X(i,a)*X(i,b)
     enddo
     Pab(a,b) = temp
     temp = 0.0D0
   enddo
   enddo

   deallocate(X)
   allocate(dPat(M,M))

!Form All Diference Density Matrix in AO basis
   temp = 0.0D0
   NCOc = NCO + 1
   do row=1,M
   do col=1,M
     do i=1,NCO
     do j=1,NCO
       temp = temp + Pij(i,j)*Coef(row,NCOc-i)*Coef(col,NCOc-j)
     enddo
     enddo
     do a=1,Nvirt
     do b=1,Nvirt
       temp = temp + Pab(a,b)*Coef(row,a+NCO)*Coef(col,b+NCO)
     enddo
     enddo
     dPat(row,col) = temp
     temp = 0.0D0
   enddo
   enddo

   deallocate(Pij,Pab)
   allocate(PFat(M,M)PEat(M,M))

   call spunpack_rho('L',M,RMM,PFat)

!Form Unrelaxed Excited State Density Matrix in AO basis
   do i=1,M
   do j=1,M
     PEat(i,j) = PFat(i,j) + dPat(i,j)
   enddo
   enddo
   end subroutine UnDensExc
end module lrtddft
