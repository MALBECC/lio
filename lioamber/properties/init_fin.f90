subroutine properties_initialise(open_shell)
   use fileio         , only: safeio_open
   use properties_data, only: fmulliken, fdipole, flowdin, fbecke,     &
                              ffukui, mulliken, dipole, lowdin, becke, &
                              fukui, UIDs
   implicit none
   logical, intent(in) :: open_shell
   character(len=60)   :: tmp_name
   
   if (mulliken) then
      call safeio_open(UIDs%mul, fmulliken, 3)
      tmp_name = trim(fmulliken)//"_spin"
      if (open_shell) call safeio_open(UIDs%muls, tmp_name, 3)
   endif

   if (lowdin) then
      call safeio_open(UIDs%low, flowdin, 3)
      tmp_name = trim(flowdin)//"_spin"
      if (open_shell) call safeio_open(UIDs%lows, tmp_name, 3)
   endif

   if (becke) then
      call safeio_open(UIDs%bec, fbecke, 3)
      tmp_name = trim(fbecke)//"_spin"
      if (open_shell) call safeio_open(UIDs%becs, tmp_name, 3)
   endif

   if (fukui) call safeio_open(UIDs%fuk, ffukui, 3)
   if (dipole) call safeio_open(UIDs%dip, fdipole, 3)
end subroutine properties_initialise

subroutine properties_finalise()
   use fileio         , only: safeio_close
   use properties_data, only: UIDs
   implicit none

   call safeio_close(UIDs%dip)
   call safeio_close(UIDs%fuk)
   call safeio_close(UIDs%bec)
   call safeio_close(UIDs%becs)
   call safeio_close(UIDs%low)
   call safeio_close(UIDs%lows)
   call safeio_close(UIDs%mul)
   call safeio_close(UIDs%muls)
end subroutine properties_finalise