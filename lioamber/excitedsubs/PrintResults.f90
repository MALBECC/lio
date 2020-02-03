subroutine PrintResults(vec,val,O,N,nstat,Mlr,NCOlr)
use excited_data, only: nfo
   implicit none

   integer, intent(in) :: N, nstat, Mlr, NCOlr
   double precision, intent(in) :: vec(N,nstat),val(nstat),O(nstat)

   character(len=4) :: j_char, from_char, to_char
   integer :: i,j,from,to
   double precision :: value_X

   from = NCOlr
   to = NCOlr + 1

   do j=1, nstat
   write (j_char, '(i4)') j
   write(*,100) adjustl(j_char), val(j), 45.56335D0/val(j), O(j)
   do i=1, N
      value_X = vec(i,j) / dsqrt(2.0D0)
      if ( abs(value_X) > 0.1D0 ) then
         write (from_char, '(i4)') from+nfo
         write (to_char, '(i4)') to+nfo
         write(*,101) adjustl(from_char), adjustl(to_char), value_X
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

   100 FORMAT(1X,"STATE ",A,3X,"ENERGY=",F8.4," Hartree, ",&
              F12.6," nm"," OSC=",F8.4)
   101 FORMAT(3X,A,"-> ",A,2X,F14.7)
end subroutine PrintResults
