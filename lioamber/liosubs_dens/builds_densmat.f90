!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine builds_densmat( Msize, Nocup, Focup, coef_mat, dens_mat )
!
! Msize:      number of atomic basis functions.
! Nocup:      number of occupied orbitals.
! Focup:      orbital ocupation factor (2.0d0 for CS, 1.0d0 for OS)
! coef_mat:   matrix containing Fock coefficients.
! dens_mat:   matrix containing output density matrix.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   implicit none
   integer, intent(in)  :: Msize
   integer, intent(in)  :: Nocup
   real*8 , intent(in)  :: Focup
   real*8 , intent(in)  :: coef_mat(Msize, Msize)
   real*8 , intent(out) :: dens_mat(Msize, Msize)

   real*8 , allocatable :: coef_occ(:, :)
   integer              :: ii, jj

!  Copies the occupied orbitals into a truncated matrix
   if ( .not.allocated(coef_occ) ) allocate( coef_occ(Msize,Nocup) )
   do jj = 1, Nocup
   do ii = 1, Msize
      coef_occ(ii, jj) = coef_mat(ii, jj)
   enddo
   enddo

!  Obtains dens_mat as (coef_occ)*(coef_occ^t).
   dens_mat(:,:) = 0.0D0
   call DGEMM( 'N', 'T', Msize, Msize, Nocup, Focup, coef_occ, Msize, coef_occ,&
             & Msize, 0.0D0, dens_mat, Msize)
    
   deallocate(coef_occ)
end subroutine builds_densmat


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
