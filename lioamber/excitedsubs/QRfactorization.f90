subroutine QRfactorization(A,N,M)
   implicit none

   integer, intent(in) :: N, M
   double precision, intent(inout) :: A(N,M)

   double precision,dimension(:),allocatable :: WORK, TAU
   integer :: LWORK,INFO

   allocate(WORK(1),TAU(M))
   call dgeqrf(N,M,A,N,TAU,WORK,-1,INFO)
   LWORK = int(WORK(1))
   deallocate(WORK); allocate(WORK(LWORK))
   call dgeqrf(N,M,A,N,TAU,WORK,LWORK,INFO)
   call dorgqr(N,M,M,A,N,TAU,WORK,LWORK,INFO)
   deallocate(WORK)
end subroutine QRfactorization
