!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine check_vecsize_n( nsize, vector, vecname, subname )
   implicit none
   integer         , intent(in) :: nsize
   integer         , intent(in) :: vector(:)
   character(len=*), intent(in) :: vecname
   character(len=*), intent(in) :: subname
#  include "check_vecsize.proced.f90"
end subroutine check_vecsize_n
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine check_vecsize_r( nsize, vector, vecname, subname )
   implicit none
   integer         , intent(in) :: nsize
   real(kind=4)          , intent(in) :: vector(:)
   character(len=*), intent(in) :: vecname
   character(len=*), intent(in) :: subname
#  include "check_vecsize.proced.f90"
end subroutine check_vecsize_r
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine check_vecsize_d( nsize, vector, vecname, subname )
   implicit none
   integer         , intent(in) :: nsize
   LIODBLE          , intent(in) :: vector(:)
   character(len=*), intent(in) :: vecname
   character(len=*), intent(in) :: subname
#  include "check_vecsize.proced.f90"
end subroutine check_vecsize_d
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine check_vecsize_c( nsize, vector, vecname, subname )
   implicit none
   integer         , intent(in) :: nsize
   complex(kind=4)       , intent(in) :: vector(:)
   character(len=*), intent(in) :: vecname
   character(len=*), intent(in) :: subname
#  include "check_vecsize.proced.f90"
end subroutine check_vecsize_c
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine check_vecsize_z( nsize, vector, vecname, subname )
   implicit none
   integer         , intent(in) :: nsize
   complex(kind=8)      , intent(in) :: vector(:)
   character(len=*), intent(in) :: vecname
   character(len=*), intent(in) :: subname
#  include "check_vecsize.proced.f90"
end subroutine check_vecsize_z
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
