            subroutine matmulnanoc(A,B,C,M)
!!!!!!!!  Hace C=Bt*(A*B)
       REAL*8 , intent(inout)  :: B(M,M)
       COMPLEX*8, intent(inout) :: A(M,M), C(M,M)
       logical ta, tb

       COMPLEX*8, dimension (:,:), ALLOCATABLE :: scratch,scratch1



       allocate (scratch(M,M),scratch1(M,M))

        do i=1,M
        do j=1,M
         scratch(i,j)=A(j,i)
        enddo
        enddo
         scratch1=0
        do i=1,M
        do j=1,M
         do k= 1,M
         scratch1(i,j)= scratch1(i,j) + scratch(k,i)*B(k,j)
         enddo
        enddo
        enddo

         C=0
        do i=1,M
        do j=1,M
         do k= 1,M
        C(i,j)= C(i,j) + B(k,i)*scratch1(k,j)
         enddo
       enddo
       enddo
       return


       end












