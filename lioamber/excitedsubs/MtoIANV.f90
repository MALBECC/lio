subroutine MtoIANV(F,C,Cocc_trans,A,M,Mlr,NCO,Ndim,Sdim,V1,iv)
use extern_functional_data, only: HF

   implicit none

   integer, intent(in) :: M, Mlr, Ndim, NCO, V1, Sdim, iv
   LIODBLE, intent(in) :: F(M,M), C(M,Mlr), Cocc_trans(NCO,M)
   LIODBLE, intent(inout) :: A(Ndim,Sdim)

   integer :: i, j, k, row, Nvirt, NCOc
   character(len=1) :: code 
   LIODBLE, dimension(:,:), allocatable :: B

   LIODBLE :: temp

   Nvirt = Mlr - NCO
   NCOc = NCO + 1
   allocate(B(NCO,M))

   code = 'N'
   do i=1,3
      if ( HF(i)==1 ) code = 'T'
   enddo
   call dgemm('N',code,NCO,M,M,1.0d0,Cocc_trans,NCO,F,M,0.0d0,B,NCO)

   temp = 0.0D0
   do i=1,NCO
   do j=NCOc,Mlr
     do k=1,M
       temp = temp + B(NCOc-i,k) * C(k,j)
     enddo
     row = (i-1) * Nvirt + (j-NCO)
     A(row,V1+iv) = temp
     temp = 0.0D0
   enddo
   enddo

   deallocate(B)
end subroutine MtoIANV
