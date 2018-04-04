!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module tmpaux_SCF
   implicit none
   contains
!
! TODO:
!
!  neighbor_list_2e: Find out what this is about.
!
!  starting_guess: Generates density starting guess. Find out how it build
!                  and add a better description.
!
!  COPY_VEC: this should go inside of maskrmm, or replaced by a subroutine
!            there,
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
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

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   subroutine xPut( xMat, nCol, rMat )
   real*8 , intent(inout) :: xMat(:,:)
   integer, intent(in)    :: nCol
   real*8 , intent(in)    :: rMat(:,:)
   integer                :: ii, jj

   if ( size(xMat,1) /= size(rMat,1) ) then
      print*, "ERROR IN xPut - wrong number of rows"; stop
   end if

   if ( size(xMat,2) > (nCol + size(rMat,2)) ) then
      print*, "ERROR IN xPut - wrong number of cols"; stop
   end if

   do jj = nCol, size(rMat,2)
   do ii = 1, size(rMat,1)
      xMat(ii,jj) = rMat(ii,jj)
   end do
   end do

   end subroutine xPut

!------------------------------------------------------------------------------!
   subroutine xGet( xMat, nCol, rMat )
   real*8 , intent(in)    :: xMat(:,:)
   integer, intent(in)    :: nCol
   real*8 , intent(inout) :: rMat(:,:)
   integer                :: ii, jj

   if ( size(xMat,1) /= size(rMat,1) ) then
      print*, "ERROR IN xGet - wrong number of rows"; stop
   end if

   if ( size(xMat,2) < (nCol + size(rMat,2)) ) then
      print*, "ERROR IN xGet - wrong number of cols"; stop
   end if

   do jj = nCol, size(rMat,2)
   do ii = 1, size(rMat,1)
      rMat(ii,jj) = xMat(ii,jj)
   end do
   end do

   end subroutine xGet



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   SUBROUTINE COPY_VEC(VEC,DIM_VEC,POINTER_RMM)
!     subrutina temporal para empezar a romper RMM
!     copia el vector VEC a RMM posicion POINTER_RMM

      use garcha_mod, ONLY: RMM
      IMPLICIT NONE
      integer, intent(in) :: DIM_VEC,POINTER_RMM
      real*8, dimension(DIM_VEC), intent(in) :: VEC
      integer :: i

      do i=1, DIM_VEC
         RMM(POINTER_RMM+i-1)=VEC(i)
      end do

   END SUBROUTINE COPY_VEC



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
end module tmpaux_SCF
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
