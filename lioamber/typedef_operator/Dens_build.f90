!carlos: this subroutine builds the density matrix.
subroutine Dens_build(this, Msize, Nocup, Focup, coef_mat)
   implicit none
   class(operator), intent(inout) :: this
   integer, intent(in)            :: Msize
   integer, intent(in)            :: Nocup
   LIODBLE , intent(in)            :: Focup
   LIODBLE , intent(in)            :: coef_mat(Msize, Msize)

   LIODBLE , allocatable :: coef_occ(:, :), dens_mat(:,:)
   integer :: ii, jj

   !  Copies the occupied orbitals into a truncated matrix
   if ( .not.allocated(coef_occ) ) allocate( coef_occ(Msize,Nocup) )
   do jj = 1, Nocup
   do ii = 1, Msize
      coef_occ(ii, jj) = coef_mat(ii, jj)
   enddo
   enddo

   !  Obtains dens_mat as (coef_occ)*(coef_occ^t).
   if ( .not.allocated(dens_mat) ) allocate( dens_mat(Msize,Msize) )
   dens_mat = 0.0D0
   call DGEMM('N', 'T', Msize, Msize, Nocup, Focup, coef_occ, Msize, coef_occ, &
              Msize, 0.0D0, dens_mat, Msize)
   this%data_AO = dens_mat

   deallocate(coef_occ, dens_mat)
end subroutine Dens_build
