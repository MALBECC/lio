!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine read_list_n( file_name, list )
   implicit none
   character(len=*), intent(in)  :: file_name
   integer         , intent(out) :: list(:)
#  include "read_list.proced.f90"
end subroutine read_list_n
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine read_list_r( file_name, list )
   implicit none
   character(len=*), intent(in)  :: file_name
   real(kind=4)          , intent(out) :: list(:)
#  include "read_list.proced.f90"
end subroutine read_list_r
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine read_list_d( file_name, list )
   implicit none
   character(len=*), intent(in)  :: file_name
   LIODBLE          , intent(out) :: list(:)
#  include "read_list.proced.f90"
end subroutine read_list_d
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine read_list_c( file_name, list )
   implicit none
   character(len=*), intent(in)  :: file_name
   complex(kind=4)       , intent(out) :: list(:)
#  include "read_list.proced.f90"
end subroutine read_list_c
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine read_list_z( file_name, list )
   implicit none
   character(len=*), intent(in)  :: file_name
   complex(kind=8)      , intent(out) :: list(:)
#  include "read_list.proced.f90"
end subroutine read_list_z
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
