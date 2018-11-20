module SCF_aux
contains

   SUBROUTINE neighbor_list_2e()
! Para hacer lineal la integral de 2 electrone con lista de vecinos. Nano
      USE garcha_mod, ONLY : natom, natomc, r, d, jatc, rmax, nshell, atmin,   &
                           & nnps, nnpp, nnpd, M, nuc
      IMPLICIT NONE
      INTEGER :: i,j, iij, iik, iikk
      REAL*8  :: zij, ti, tj, alf, rexp

      do i=1,natom
         natomc(i)=0
         do j=1,natom
            d(i,j)=(r(i,1)-r(j,1))**2+(r(i,2)-r(j,2))**2+(r(i,3)-r(j,3))**2
            zij=atmin(i)+atmin(j)
            ti=atmin(i)/zij
            tj=atmin(j)/zij
            alf=atmin(i)*tj
            rexp=alf*d(i,j)
            if (rexp.lt.rmax) then
               natomc(i)=natomc(i)+1
               jatc(natomc(i),i)=j
            endif
         enddo
      enddo

      do iij=nshell(0),1,-1
        nnps(nuc(iij))=iij
      enddo

      do iik=nshell(0)+nshell(1),nshell(0)+1,-1
        nnpp(nuc(iik))=iik
      enddo

      do iikk=M,nshell(0)+nshell(1)+1,-1
        nnpd(nuc(iikk))=iikk
      enddo
   END SUBROUTINE neighbor_list_2e


   subroutine seek_nan(vecToTest, vecStart, vecEnd, phrase)
       implicit none
       real*8           , intent(in) :: vecToTest(*)     ! Vector to analize.
       integer          , intent(in) :: vecStart, vecEnd ! Vector range to analize.
       character (len=*), intent(in) :: phrase           ! Output phrase for NaN.
       integer :: iNick

       if (vecStart .gt. vecEnd) then
           write(*,*) "Error: vector start index greater than end index."
           write(*,*) phrase
           stop
       endif

       do iNick = vecStart, vecEnd
           if (vecToTest(iNick) .ne. vecToTest(iNick)) then
               write(*,*) "NaN found in: ", phrase, iNick
               stop
           end if
       enddo
   endsubroutine seek_nan


end module SCF_aux
