module lr_data
!Nvirt = number of virtual molecular orbitals
!dim = dimension of matrix Linear Response (Nvirt x NCO)
!cbas = Needed for libint (Temporary)
!nstates = number of excited states to calculate
!root = excited state chosen for optimization
   implicit none

   logical :: lresp = .false.
   integer :: Nvirt, dim, NCOlr, Mlr
   integer :: nstates = 3
   integer :: root = 0
   real*8, dimension(:,:), allocatable :: eigenvec, cbas
   real*8, dimension(:), allocatable :: eigenval
! Use Frozen Core Approximation
! nfo = number of molecular orbitals occupied with lowest energy deleted
! nfv = number of molecular orbiatls virtual with higher energy deleted
   logical :: FCA = .false.
   integer :: nfo = 3
   integer :: nfv = 3
   
end module lr_data

module lrtddft
! This is a main subroutine for Linear Response calculate.
! This module performs a Tamm-Dancoff Aproximation (TDA)

contains

   subroutine linear_response(MatCoef,VecEne)
   use lr_data, only: Nvirt,dim,nstates,eigenval,&
                      eigenvec,cbas,root,FCA,nfo,nfv,&
                      NCOlr, Mlr
   use garcha_mod, only: NCO
   use basis_data, only: c_raw, max_c_per_atom, a, c, M

   implicit none

   real*8, intent(in) :: MatCoef(M,M)
   real*8, intent(in) :: VecEne(M)

   integer :: i
   real*8, dimension(:), allocatable :: Ene_LR
   real*8, dimension(:,:), allocatable :: KMO, A_mat, Coef_LR
   real*8, dimension(:,:,:,:), allocatable :: KAO

   call g2g_timer_start('LINEAR RESPONSE')
   if (.not. allocated(cbas)) allocate(cbas(M, max_c_per_atom))
   cbas = c_raw

! Initialization of variables
   if (FCA .eqv. .true.) then
   print*,"Using Frozen Core Approximation"
   print*,"nfo, nfv", nfo, nfv
      Nvirt = M - NCO - nfv
      NCOlr = NCO - nfo
      Mlr = M - nfo - nfv
      dim = Nvirt * NCOlr
      allocate(Coef_LR(M,Mlr),Ene_LR(Mlr))
      do i=1, NCOlr
        Coef_LR(:,i) = MatCoef(:,i+nfo)
        Ene_LR(i) = VecEne(i+nfo)
      enddo
      do i=1, Nvirt
         Coef_LR(:,NCOlr+i) = MatCoef(:,i+NCO)
         Ene_LR(NCOlr+i) = VecEne(i+NCO)
      enddo
   else
      nfo = 0
      Mlr = M
      NCOlr = NCO
      Nvirt = M - NCO
      dim = Nvirt * NCO
      allocate(Coef_LR(M,M),Ene_LR(M))
      Coef_LR = MatCoef
      Ene_LR = VecEne
   endif

   allocate(KMO(dim,dim),A_mat(dim,dim))
   allocate(KAO(M,M,M,M))

   KAO = 0.0D0
   KMO = 0.0D0

   print*, ""
   print*,"#######################################"
   print*,"         LINEAR RESPONSE - TDA"
   print*,"#######################################"
   print*, ""

   call g2g_timer_start('g2g_LinearResponse')
   call g2g_linear_response(Coef_LR,KMO,KAO,cbas,dim,NCOlr,Nvirt)
   call g2g_timer_stop('g2g_LinearResponse')
   deallocate(KAO)

   call g2g_timer_start('ObtainMatrix')
   call ObtainAmatrix(A_mat,KMO,Ene_LR,dim,Mlr,NCOlr)
   call g2g_timer_stop('ObtainMatrix')
   deallocate(KMO)

   call g2g_timer_start('Davidson')
   call davidson(A_mat,dim,eigenval,eigenvec,nstates)
   call g2g_timer_stop('Davidson')
   deallocate(A_mat)

   call PrintResults(eigenvec,eigenval,dim,nstates)

   if(root > 0) then
     call UnDensExc(eigenvec(:,root),Coef_LR,dim)
   endif
   deallocate(eigenvec,eigenval)
   deallocate(Coef_LR,Ene_LR)

   call g2g_timer_stop('LINEAR RESPONSE')
   end subroutine linear_response

   subroutine ObtainAmatrix(A,Kcou,VecE,N,M,NCO)
   implicit none

   integer, intent(in) :: N, M, NCO
   real*8, intent(in) :: Kcou(N,N),VecE(M)
   real*8, intent(out) :: A(N,N)

   integer :: k,l,i,j !counters

   k = NCO
   l = NCO + 1
   do i=1,N
     A(i,i) = Kcou(i,i) + VecE(l) - VecE(k)
     l = l + 1
     if (l > M) then
       l = NCO + 1
       k = k - 1
     endif
   enddo
   do i=1,N
     do j=1,i-1
         A(j,i) = Kcou(i,j)
         A(i,j) = A(j,i)
     enddo
   enddo
   end subroutine ObtainAmatrix

   subroutine davidson(A,N,eigval,eigvec,Nstat)
   ! This subroutine not performs Davidson Diagonalization, instead performs
   ! DQDS algorithm
   ! TODO: Direct Davidson Iterative Diagonalization
   implicit none

   integer, intent(in) :: N, Nstat
   real*8, intent(inout) :: A(N,N)
   real*8, dimension(:,:), allocatable, intent(out) :: eigvec
   real*8, dimension(:), allocatable, intent(out) :: eigval

   integer :: IL, IU, num, LWORK, LIWORK, info
   integer, dimension(:), allocatable :: IWORK, ISUPPZ
   real*8 :: VL = 1.0 !NOT REFERENCE
   real*8 :: VU = 20.0 !NOT REFERENCE
   real*8 :: ABSTOL = 1.0E-10
   real*8, dimension(:), allocatable :: WORK

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
   use lr_data, only: Mlr, NCOlr, nfo
!   use garcha_mod, only: M,NCO
   implicit none

   integer, intent(in) :: N, nstat
   real*8, intent(in) :: vec(N,nstat),val(nstat)

   integer :: i,j,from,to

   from = NCOlr
   to = NCOlr + 1

   do j=1, nstat
   write(*,100) j,val(j)
   do i=1, N
      if ( abs(vec(i,j)) > 0.1D0 ) then
         write(*,101) from+nfo, to+nfo, vec(i,j)
      endif
      to = to + 1
      if ( to == Mlr+1 ) then
          from = from - 1
          to = NCOlr + 1
      endif
   enddo
      print*, " "
      from = NCOlr
      to = NCOlr + 1
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
   use lr_data, only: Mlr, NCOlr, Nvirt
   use garcha_mod, only: Pmat_vec
   use basis_data, only: M
   implicit none

   integer, intent(in) :: N
   real*8, intent(in) :: vec(N),Coef(M,Mlr)

   real*8, dimension(:,:), allocatable :: Pij, Pab, dPat, X, PEat, PFat,&
                                          Xtrans, Cocc, Cnvirt, CTocc,&
                                          CTnvirt

   integer :: i,j,a,b,row,col,NCOc,Nij,Nab,M2

   allocate(Pij(NCOlr,NCOlr),Pab(Nvirt,Nvirt),X(NCOlr,Nvirt),Xtrans(Nvirt,NCOlr))
   Pij = 0.0D0
   Pab = 0.0D0
   X = 0.0D0

   do row=0,NCOlr-1
   do col=1,Nvirt
     X(row+1,col) = vec(row*Nvirt+col)
     Xtrans(col,row+1) = vec(row*Nvirt+col)
   enddo
   enddo

!Form block Diference Density Matrix occ-occ in MO basis
   Pij = -1.0D0*matmul(X,Xtrans)

!Form block Diference Density Matrix nvirt-nvirt in MO basis
   Pab = matmul(Xtrans,X)

   deallocate(X,Xtrans)
   allocate(dPat(M,M),Cocc(M,NCOlr),Cnvirt(M,Nvirt))
   allocate(CTocc(NCOlr,M),CTnvirt(Nvirt,M))

!Form All Diference Density Matrix in AO basis
   ! Form block occ - occ
   NCOc = NCOlr + 1
   do i=1,NCOlr
      Cocc(:,i) = Coef(:,NCOc-i)
      CTocc(i,:) = Coef(:,NCOc-i)
   enddo
   dPat=matmul(Cocc,matmul(Pij,CTocc))
   deallocate(Cocc,CTocc)

   ! Form block virt - virt 
   do i=1,Nvirt
      Cnvirt(:,i) = Coef(:,NCOlr+i)
      CTnvirt(i,:) = Coef(:,NCOlr+i)
   enddo
   dPat = dPat + matmul(Cnvirt,matmul(Pab,CTnvirt))
   deallocate(Cnvirt,CTnvirt)
   deallocate(Pij,Pab)

! Extract Rho of ground state
   allocate(PFat(M,M),PEat(M,M))
   call spunpack_rho('L',M,Pmat_vec,PFat)

!Form Unrelaxed Excited State Density Matrix in AO basis
   do i=1,M
   do j=1,M
     PEat(i,j) = PFat(i,j) + dPat(i,j)
   enddo
   enddo
   deallocate(PFat,PEat,dPat)
   end subroutine UnDensExc
end module lrtddft
