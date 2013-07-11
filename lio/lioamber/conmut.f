            subroutine conmut(A,B,C,M)
!!!!!!!!  Hace C=A*B-B*A
       REAL*8 , intent(inout)  :: A(M,M), B(M,M), C(M,M)
       logical ta, tb

       real*8, dimension (:,:), ALLOCATABLE :: scratch,scratch1,scratch2



       allocate (scratch(M,M),scratch1(M,M),scratch2(M,M))

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


        do i=1,M
        do j=1,M
         scratch(i,j)=B(j,i)
        enddo
        enddo

         scratch2=0
         C=0
        do i=1,M
        do j=1,M
         do k= 1,M
         scratch2(i,j)= scratch2(i,j) + scratch(k,i)*A(k,j)
         enddo
        C(i,j)= scratch1(i,j)-scratch2(i,j)
        enddo
        enddo




       deallocate (scratch,scratch1,scratch2)
       return
 

       end












