subroutine norma(V,N,norm)
  implicit none

  integer, intent(in) :: N
  double precision, intent(in) :: V(N)
  double precision,intent(out) :: norm

  integer :: i

  norm = 0.0d0
  do i=1,N
    norm = norm + V(i)*V(i)
  enddo
end subroutine norma
