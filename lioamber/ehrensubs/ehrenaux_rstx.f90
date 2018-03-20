!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
! Save the forces and density matrix of the last position for restart purposes
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine ehrenaux_rsto( rsto_fname, rsto_funit, rsto_nfreq, ndyn_steps,     &
                         & step_number, Natom, forces, nucvel, &
                         & Nbasis, densA, densB )
   implicit none
   character(len=*), intent(in) :: rsto_fname
   integer         , intent(in) :: rsto_funit
   integer         , intent(in) :: rsto_nfreq
   integer         , intent(in) :: ndyn_steps
   integer         , intent(in) :: step_number
   integer         , intent(in) :: Natom
   real*8          , intent(in) :: forces( 3, Natom )
   real*8          , intent(in) :: nucvel( 3, Natom )
   integer         , intent(in) :: Nbasis
   complex*16      , intent(in) :: densA( Nbasis, Nbasis )
   complex*16      , intent(in) :: densB( Nbasis, Nbasis )

   logical :: save_this_step
   integer :: ii, jj

   save_this_step = .false.

   if ( rsto_nfreq > 0 ) then
      if ( modulo(step_number,rsto_nfreq) == 1 ) save_this_step = .true.
   endif
   
!  is the ndyn+1 part right?
   if ( step_number == (ndyn_steps+1) ) save_this_step = .true.


   if ( save_this_step ) then
      open( unit=rsto_funit, file=rsto_fname )

      do jj=1,Natom
      do ii=1,3
         write( unit=rsto_funit, fmt=101 ) ii, jj, forces(ii, jj), nucvel(ii, jj)
      enddo
      enddo

      do jj=1,Nbasis
      do ii=1,Nbasis
         write( unit=rsto_funit, fmt=101 ) ii, jj, densA(ii, jj), densB(ii, jj)
      enddo
      enddo

      close( unit=rsto_funit )
   endif

101 format(2x,I3,2x,I3,4(2x,ES20.12))
end subroutine ehrenaux_rsto
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine ehrenaux_rsti( rsti_fname, rsti_funit, Natom, forces, nucvel,      &
                         & Nbasis, densA, densB )
   implicit none
   character(len=*) , intent(in)    :: rsti_fname
   integer          , intent(in)    :: rsti_funit
   integer          , intent(in)    :: Natom
   real*8           , intent(inout) :: forces( 3, Natom )
   real*8           , intent(inout) :: nucvel( 3, Natom )
   integer          , intent(in)    :: Nbasis
   complex*16       , intent(inout) :: densA( Nbasis, Nbasis )
   complex*16       , intent(inout) :: densB( Nbasis, Nbasis )

   integer :: ii, jj, isc, jsc

   print*,'Using restart'
   open( unit=rsti_funit, file=rsti_fname )

   do jj=1,Natom
   do ii=1,3
      read( unit=rsti_funit, fmt=102 ) isc, jsc, forces(ii, jj), nucvel(ii, jj)
   enddo
   enddo

   do jj=1,Nbasis
   do ii=1,Nbasis
      read( unit=rsti_funit, fmt=102 ) isc, jsc, densA(ii, jj), densB(ii,jj)
   enddo
   enddo

   close( unit=rsti_funit )

102 format(2x,I3,2x,I3,4(2x,ES20.12))
end subroutine ehrenaux_rsti
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
