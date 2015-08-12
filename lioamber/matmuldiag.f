      subroutine matmuldiag(A,B,C,M)
!!!!!!!!  Hace C=(A*B) para matrices cuadradas solo la diagonal
      REAL*8 , intent(in)  :: A(M,M), B(M,M)
      REAL*8 , intent(out)  :: C(M,M)
      C=0
      do k= 1,M
        do i=1,M
          C(i,i)= C(i,i) + A(i,k)*B(k,i)
        enddo
      enddo

      return
      end

