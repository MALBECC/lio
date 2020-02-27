module precision
! Fixed precision
   integer, parameter :: sp = selected_real_kind(6, 37)
   integer, parameter :: dp = selected_real_kind(15, 307)
   integer, parameter :: qp = selected_real_kind(33, 4931)
! Modifiable precision
#  ifdef  CPU_SIMPLE
   integer, parameter :: xp = selected_real_kind(6, 37)
#  else
#  ifdef  CPU_DOUBLE
   integer, parameter :: xp = selected_real_kind(15, 307)
#  endif
#  endif
end module precision


