!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
! Save the forces and density matrix of the last position for restart purposes
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine ehrenaux_rsti &
& ( rsti_fname, Natoms, nucfors, nucvels, Nbasis, densA, densB )
   implicit none
   character(len=*), intent(in)    :: rsti_fname
   integer         , intent(in)    :: Natoms
   real*8          , intent(inout) :: nucfors( 3, Natoms )
   real*8          , intent(inout) :: nucvels( 3, Natoms )
   integer         , intent(in)    :: Nbasis
   complex*16      , intent(inout) :: densA( Nbasis, Nbasis )
   complex*16      , intent(inout) :: densB( Nbasis, Nbasis )

   character(len=*), parameter     :: myfmt="(2x,I3,2x,I3,4(2x,ES24.16))"

   integer :: rsti_funit
   integer :: ii, jj, isc, jsc

   rsti_funit = 123456
   open( unit=rsti_funit, file=rsti_fname )

   do jj=1,Natoms
   do ii=1,3
      read( unit=rsti_funit, fmt=myfmt ) &
      & isc, jsc, nucfors(ii, jj), nucvels(ii, jj)
   enddo
   enddo

   do jj=1,Nbasis
   do ii=1,Nbasis
      read( unit=rsti_funit, fmt=myfmt ) &
      & isc, jsc, densA(ii, jj), densB(ii,jj)
   enddo
   enddo

   close( unit=rsti_funit )

end subroutine ehrenaux_rsti
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine ehrenaux_rsto( rsto_fname, rsto_nfreq, ndyn_steps, step_number,     &
                        & Natoms, nucfors, nucvels, Nbasis, densA, densB )
   implicit none
   character(len=*), intent(in) :: rsto_fname
   integer         , intent(in) :: rsto_nfreq
   integer         , intent(in) :: ndyn_steps
   integer         , intent(in) :: step_number
   integer         , intent(in) :: Natoms
   real*8          , intent(in) :: nucfors( 3, Natoms )
   real*8          , intent(in) :: nucvels( 3, Natoms )
   integer         , intent(in) :: Nbasis
   complex*16      , intent(in) :: densA( Nbasis, Nbasis )
   complex*16      , intent(in) :: densB( Nbasis, Nbasis )

   character(len=*), parameter  :: myfmt="(2x,I3,2x,I3,4(2x,ES24.16))"

   logical :: save_this_step
   integer :: rsto_funit
   integer :: ii, jj

   save_this_step = .false.

   if ( rsto_nfreq > 0 ) then
      if ( modulo(step_number,rsto_nfreq) == 0 ) save_this_step = .true.
   endif
   if ( step_number == (ndyn_steps) ) save_this_step = .true.

   if ( .not. save_this_step ) return
!  ELSE

   rsto_funit = 123456
   open( unit=rsto_funit, file=rsto_fname )

   do jj=1,Natoms
   do ii=1,3
      write( unit=rsto_funit, fmt=myfmt ) &
      & ii, jj, nucfors(ii, jj), nucvels(ii, jj)
   enddo
   enddo

   do jj=1,Nbasis
   do ii=1,Nbasis
      write( unit=rsto_funit, fmt=myfmt ) &
      & ii, jj, densA(ii, jj), densB(ii, jj)
   enddo
   enddo

   close( unit=rsto_funit )

end subroutine ehrenaux_rsto
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
