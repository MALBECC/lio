!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine find_free_unit( free_unit, last_unit, ext_stat )
!------------------------------------------------------------------------------!
!
!  DESCRIPTION TODO
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!   use errlog TODO
   implicit none
   integer, intent(inout)         :: free_unit
   integer, intent(in)            :: last_unit
   integer, intent(out), optional :: ext_stat
!
   integer :: first_try
   integer :: final_try
   logical :: unit_is_open
   logical :: abort_run
   integer :: int_stat
!
!
!  LOOKING FOR A FREE UNIT
!------------------------------------------------------------------------------!
   int_stat = 0
   first_try = max( free_unit, 10 )
   final_try = last_unit

   unit_is_open = .true.
   do free_unit = first_try, final_try
      inquire( unit=free_unit, opened=unit_is_open, iostat=int_stat )
      if (.not.unit_is_open) exit
   end do
!
!
!  DEALING WITH ERRORS
!------------------------------------------------------------------------------!
   abort_run = .false.

   if ( int_stat .ne. 0 ) then
      if ( present(ext_stat) ) then
         ext_stat = int_stat
      else
         abort_run=.true.
      end if
   end if

   if ( unit_is_open ) then
      int_stat = -1
      if ( present(ext_stat) ) then
         ext_stat = int_stat
      else
         abort_run=.true.
      end if
   end if

   if (abort_run) then
      print*
      print*,"CRITICAL ERROR IN FIND_FREE_UNIT"
      print*,"* First unit tried:   ",first_try
      print*,"* Final unit tried:   ",final_try
      print*,"* Last unit checked:  ",free_unit
      print*,"* That unit was open: ",unit_is_open
      print*,"* Internal status:    ",int_stat
      print*
      print*,"If you don't know what provoked this error, please contact us"
      print*,"through our github webpage //github.com/MALBECC/lio"
      print*,"Aborting run..."; print*; stop
   end if

end subroutine find_free_unit
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
