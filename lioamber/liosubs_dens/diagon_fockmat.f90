!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine diagon_fockmat( Nsize, fock_mat, Xmat, morb_coefat, morb_energy )

   use linear_algebra, only: matrix_diagon
   implicit none
   integer  , intent(in)  :: Msize
   real*8   , intent(in)  :: fock_mat( Msize, Msize )
   real*8   , intent(in)  :: Xmat( M_in, M_in )
   real*8   , intent(out) :: morb_coefat( M_in, M_in )
   real*8   , intent(out) :: morb_energy( M_in )
   real*8, allocatable    :: eigen_vecs(:,:)


!  Fock(ON) diagonalization, base change of coeficients ( (X^-1)*C )
!------------------------------------------------------------------------------!
      if ( allocated(eigen_vecs) ) deallocate(eigen_vecs)
      allocate( eigen_vecs(M_in,M_in) )


      call g2g_timer_start('SCF - Fock Diagonalization')
      call g2g_timer_sum_start('SCF - Fock Diagonalization (sum)')
      call matrix_diagon( fock_mat, eigen_vecs, morb_energy )
      call g2g_timer_sum_pause('SCF - Fock Diagonalization (sum)')
      call g2g_timer_stop('SCF - Fock Diagonalization')

      call g2g_timer_start('SCF - MOC base change')
      call g2g_timer_sum_start('SCF - MOC base change (sum)')
      call DGEMM( 'N', 'N', Nsize, Nsize, Nsize, 1.0d0, Xmat, Nsize, eigen_vecs, Nsize, 0.0d0, morb_coefat, Nsize )
      call g2g_timer_sum_pause('SCF - MOC base change (sum)')
      call g2g_timer_stop('SCF - MOC base change')

      if ( allocated(eigen_vecs) ) deallocate(eigen_vecs)

end subroutine diagon_fockmat


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
