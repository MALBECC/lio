!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine ehrenaux_writedip( step, freq, time, dipmom, fname )
!------------------------------------------------------------------------------!
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   use liosubs, only: find_free_unit, safeio_open

   implicit none
   integer, intent(in)          :: step
   integer, intent(in)          :: freq
   real*8 , intent(in)          :: time
   real*8 , intent(in)          :: dipmom(3)
   character(len=*), intent(in) :: fname

   real*8  :: dipmod
   integer :: funit
   integer :: step_freq
   integer :: int_stat

   character(len=*), parameter :: sfmt0 = '(A)'
   character(len=*), parameter :: sfmt1 = '(F12.3,3(2x,F12.6),4x,F12.6)'
!
!
!
!------------------------------------------------------------------------------!
   if ( freq == 0 ) return
   int_stat = 0

   funit = 10
   call find_free_unit( funit, 1000, int_stat )
   if ( int_stat /= 0 ) then
      print*,'PROBLEM INSIDE EHRENAUX_OPENFU - find_free_unit'
      stop
   endif
!
!
!
!------------------------------------------------------------------------------!
   step_freq = modulo(step, freq)

   if ( step == 1 ) then
      call safeio_open( funit, fname, 3, int_stat )
      if ( int_stat /= 0 ) then
         print*,'PROBLEM INSIDE EHRENAUX_OPENFU - safeio_open:', int_stat
         stop
      endif

!     Writes header
      write( unit=funit, fmt=sfmt0)
      write( unit=funit, fmt=sfmt0) '# DIPOLE MOMENT: time, x, y, z, norm (fs, debyes)'
      write( unit=funit, fmt=sfmt0)
      close( unit=funit, iostat=int_stat )
   endif
!   else if ( step_freq == 0 ) then
   if ( step_freq == 0 ) then
      call safeio_open( funit, fname, 2, int_stat )
      if ( int_stat /= 0 ) then
         print*,'PROBLEM INSIDE EHRENAUX_OPENFU - safeio_open', int_stat
         stop
      endif

!     Writes dipole
      dipmod = sqrt( dipmom(1)**2 + dipmom(2)**2 + dipmom(3)**2 )
      write( unit=funit, fmt=sfmt1) time, dipmom(1), dipmom(2), dipmom(3), dipmod
      close( unit=funit, iostat=int_stat )
   endif

   if ( int_stat /= 0 ) then
      print*,'PROBLEM INSIDE EHRENAUX_OPENFU - close: ', int_stat
      stop
   endif
!
!
!------------------------------------------------------------------------------!
end subroutine ehrenaux_writedip
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
