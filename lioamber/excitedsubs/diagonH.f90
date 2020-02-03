subroutine diagonH(H,N,Val,Vec)
   implicit none

   integer, intent(in) :: N
   double precision, intent(in) :: H(N,N)
   double precision, intent(out) :: Val(N),Vec(N,N)

   double precision, dimension(:), allocatable :: WORK
   integer, dimension(:), allocatable :: IWORK
   integer :: info, LWORK, LIWORK

   Vec = H
   LWORK = -1
   LIWORK = -1

   allocate(WORK(1),IWORK(1))
   call dsyevd('V','U',N,Vec,N,Val,WORK,LWORK,IWORK,LIWORK,info)
   LWORK=int(WORK(1))
   LIWORK=int(IWORK(1))
   deallocate(WORK,IWORK); allocate(WORK(LWORK),IWORK(LIWORK))
   call dsyevd('V','U',N,Vec,N,Val,WORK,LWORK,IWORK,LIWORK,info)
   deallocate(WORK,IWORK)
end subroutine diagonH
