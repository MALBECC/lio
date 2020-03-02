subroutine residual(Val,Vec,R,A,W,Ndim,Sdim,Nstat)
   implicit none

   integer, intent(in) :: Ndim, Sdim, Nstat
   LIODBLE, intent(in) :: Val(Nstat), Vec(Sdim,Nstat)
   LIODBLE, intent(in) :: R(Ndim,Nstat), A(Ndim,Sdim)
   LIODBLE, intent(out) :: W(Ndim,Nstat)

   integer :: i
   LIODBLE, dimension(:), allocatable :: temp

   allocate(temp(Ndim))
   do i=1,Nstat
     call dgemm('N','N',Ndim,1,Sdim,1.0d0,A,Ndim,Vec(:,i),Sdim,&
                 0.0d0,temp,Ndim)
     W(1:Ndim,i) = temp - Val(i) * R(1:Ndim,i)
   enddo
   deallocate(temp)
end subroutine residual
