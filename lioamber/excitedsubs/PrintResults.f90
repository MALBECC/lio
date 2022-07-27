subroutine PrintResults(vec,val,O,N,nstat,Mlr,NCOlr)
use garcha_mod  , only: NCO
use excited_data, only: map_occ, map_vir
   implicit none

   integer, intent(in) :: N, nstat, Mlr, NCOlr
   LIODBLE, intent(in) :: vec(N,nstat),val(nstat),O(nstat)

   character(len=4) :: j_char, from_char, to_char
   integer :: i,j,from,to
   LIODBLE :: value_X

   from = NCOlr
   to = NCOlr + 1

   do j=1, nstat
   write (j_char, '(i4)') j
   write(*,100) adjustl(j_char), val(j), 45.56335D0/val(j), O(j)
   do i=1, N
      value_X = vec(i,j) / dsqrt(2.0D0)
      if ( abs(value_X) > 0.1D0 ) then
         write (from_char, '(i4)') map_occ(from)
         write (to_char, '(i4)') map_vir(to-NCOlr)+NCO
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
              F10.4," nm"," OSC=",F12.8)
   101 FORMAT(3X,A,"-> ",A,2X,F14.7)
end subroutine PrintResults

subroutine PrintESA(TD,Ene,OsSt,Ns_slr,Ns_ref)
   implicit none

   integer, intent(in) :: Ns_slr, Ns_ref
   LIODBLE, intent(in) :: TD(Ns_slr,3), Ene(Ns_slr)
   LIODBLE, intent(inout) :: OsSt(Ns_slr)

   integer :: ii

   ! Oscillator Strenght Calculation
   call ObtainOsc(TD,Ene,OsSt,Ns_slr,Ns_ref)
   print*,""
   write(*,"(2X,A,5X,A,5X,A,13X,A)") "STATES","ENERGY[Ha]","LAMBDA[nm]","F. Osc."
   do ii=1,Ns_slr
      write(*,"(1X,I2,A,I2,5X,F8.4,5X,F12.6,2X,F20.10)") Ns_ref,"-> ",Ns_ref+ii,&
                                                       Ene(ii),45.56335D0/Ene(ii),&
                                                       OsSt(ii)
   enddo
   print*,""

end subroutine PrintESA

subroutine open_PrintResults(vecA,vecB,val,O,N,nstat,Mlr,NCOlrA,NCOlrB)
use garcha_mod  , only: NCO, Nunp
use excited_data, only: map_occ, map_vir, map_occb, map_virb
   implicit none

   integer, intent(in) :: N, nstat, Mlr, NCOlrA, NCOlrB
   LIODBLE, intent(in) :: vecA(N,nstat),vecB(N,nstat),val(nstat),O(nstat)

   character(len=4) :: j_char, from_char, to_char
   integer :: i,j,fromA,fromB,toA,toB
   LIODBLE :: value_X

   fromA = NCOlrA
   toA = NCOlrA + 1
   fromB = NCOlrB
   toB = NCOlrB + 1

   do j=1, nstat
      write (j_char, '(i4)') j
      write(*,100) adjustl(j_char), val(j), 45.56335D0/val(j), O(j)
      ! Alpha
      do i=1, N
         value_X = vecA(i,j)
         if ( abs(value_X) > 0.1D0 ) then
            write (from_char, '(i4)') map_occ(fromA)
            write (to_char, '(i4)') map_vir(toA-NCOlrA)+NCO
            write(*,101) adjustl(from_char), adjustl(to_char), "(A)", value_X
         endif
         toA = toA + 1
         if ( toA == Mlr+1 ) then
             fromA = fromA - 1
             toA = NCOlrA + 1
         endif
      enddo
      fromA = NCOlrA
      toA = NCOlrA + 1
      ! Beta
      do i=1, N
         value_X = vecB(i,j)
         if ( abs(value_X) > 0.1D0 ) then
            write (from_char, '(i4)') map_occb(fromB)
            write (to_char, '(i4)') map_virb(toB-NCOlrB)+NCO+Nunp
            write(*,101) adjustl(from_char), adjustl(to_char), "(B)", value_X
         endif
         toB = toB + 1
         if ( toB == Mlr+1 ) then
             fromB = fromB - 1
             toB = NCOlrB + 1
         endif
      enddo
      print*, " "
      fromB = NCOlrB
      toB = NCOlrB + 1
   enddo

   100 FORMAT(1X,"STATE ",A,3X,"ENERGY=",F8.4," Hartree, ",&
              F10.4," nm"," OSC=",F12.8)
   101 FORMAT(3X,A,"-> ",A,1X,A,2X,F14.7)
end subroutine open_PrintResults
