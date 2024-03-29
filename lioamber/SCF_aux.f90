!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% SCF_AUX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! This module contains auxiliary subroutines used in the SCF cycle and related !
! procedures.                                                                  !
!                                                                              !
! Subroutines included are:                                                    !
!  * fix_densmat       (divides non-diagonal terms of a matrix by 2)           !
!  * messup_densmat    (multiplies non-diagonal terms of a matrix by 2)        !
!  * neighbour_list_2e (makes lists for 2e integrals)                          !
!  * seek_NaN          (searches for NaNs in a vector)                         !
!  * standard_coefs    (makes sure the first element of a matrix is positive)  !
!                                                                              !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module SCF_aux
contains

subroutine fix_densmat(dens_mat)
       ! Fixes density matrix non-diagonal terms.
       implicit none
       double precision, intent(inout) :: dens_mat(:,:)
       integer :: icount, jcount

       do jcount = 1, size(dens_mat, 2)
       do icount = 1, size(dens_mat, 1)
          if (icount /= jcount) &
             dens_mat(icount,jcount) = dens_mat(icount,jcount) / 2.0D0
       enddo
       enddo
end subroutine fix_densmat

subroutine messup_densmat(dens_mat)
   ! Prepares density matrix for placement in RMM and usage in integrals.
   implicit none
   double precision, intent(inout) :: dens_mat(:,:)
   integer :: icount, jcount

   do jcount = 1, size(dens_mat, 2)
   do icount = 1, size(dens_mat, 1)
      if (icount /= jcount) &
         dens_mat(icount,jcount) = 2.0D0 * dens_mat(icount,jcount)
   enddo
   enddo
end subroutine messup_densmat

subroutine neighbour_list_2e(natom, ntatom, r, d)
   ! Makes neighbour list for 2e integrals in order to give it linear
   ! scaling. Also calculates distances (squared) between atoms.
   use basis_data, only: natomc, jatc, rmax, nshell, atmin, nnps, nnpp, nnpd, &
                         M, nuc
   implicit none
   integer         , intent(in)    :: natom, ntatom
   double precision, intent(in)    :: r(ntatom,3)
   double precision, intent(inout) :: d(natom,natom)
   integer          :: icount, jcount
   double precision :: rexp

   do icount = 1, natom
      natomc(icount) = 0
      do jcount = 1, natom
         d(icount,jcount) = (r(icount,1)-r(jcount,1))*(r(icount,1)-r(jcount,1))&
                          + (r(icount,2)-r(jcount,2))*(r(icount,2)-r(jcount,2))&
                          + (r(icount,3)-r(jcount,3))*(r(icount,3)-r(jcount,3))
         rexp = d(icount,jcount) * atmin(icount) * atmin(jcount) &
                / (atmin(icount) + atmin(jcount))
         if (rexp .lt. rmax) then
            natomc(icount) = natomc(icount) +1
            jatc(natomc(icount),icount) = jcount
         endif
      enddo
   enddo

   do icount = nshell(0), 1, -1
     nnps(nuc(icount)) = icount
   enddo
   do icount = nshell(0) + nshell(1), nshell(0) +1, -1
     nnpp(nuc(icount)) = icount
   enddo
   do icount = M, nshell(0) + nshell(1) +1, -1
     nnpd(nuc(icount)) = icount
   enddo
end subroutine neighbour_list_2e

subroutine seek_nan(vecToTest, vecStart, vecEnd, phrase)
    implicit none
    double precision , intent(in) :: vecToTest(*)     ! Vector to analize.
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
end subroutine seek_nan

subroutine standard_coefs(coef_mat)
   ! Standardises coefficient matrix. Essentaly, it makes sure the first value
   ! is always positive.
   implicit none
   double precision, intent(inout) :: coef_mat(:,:)
   integer :: icount, jcount

   do jcount = 1, size(coef_mat, 2)
      if ( coef_mat(1,jcount) < 0.0D0 ) then
         do icount = 1, size(coef_mat, 1)
            coef_mat(icount,jcount) = (-1.0D0) * coef_mat(icount,jcount)
         end do
      end if
   end do
end subroutine standard_coefs

end module SCF_aux
