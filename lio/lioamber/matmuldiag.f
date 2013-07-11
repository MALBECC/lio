            subroutine matmuldiag(A,B,C,M)
!!!!!!!!  Hace C=(A*B) para matrices cuadradas solo la diagonal
       REAL*8 , intent(inout)  :: A(M,M), B(M,M), C(M,M)
       logical ta, tb



         C=0
         do k= 1,M
          do i=1,M
         C(i,i)= C(i,i) + A(i,k)*B(k,i)
         enddo
        enddo

       return


       end












