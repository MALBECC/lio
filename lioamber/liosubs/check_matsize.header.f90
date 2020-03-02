!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine check_matsize_n( nsize1, nsize2, matrix, matname, subname )
   implicit none
   integer         , intent(in) :: nsize1
   integer         , intent(in) :: nsize2
   integer         , intent(in) :: matrix(:,:)
   character(len=*), intent(in) :: matname
   character(len=*), intent(in) :: subname
#  include "check_matsize.proced.f90"
end subroutine check_matsize_n
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine check_matsize_r( nsize1, nsize2, matrix, matname, subname )
   implicit none
   integer         , intent(in) :: nsize1
   integer         , intent(in) :: nsize2
   real*4          , intent(in) :: matrix(:,:)
   character(len=*), intent(in) :: matname
   character(len=*), intent(in) :: subname
#  include "check_matsize.proced.f90"
end subroutine check_matsize_r
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine check_matsize_d( nsize1, nsize2, matrix, matname, subname )
   implicit none
   integer         , intent(in) :: nsize1
   integer         , intent(in) :: nsize2
   LIODBLE          , intent(in) :: matrix(:,:)
   character(len=*), intent(in) :: matname
   character(len=*), intent(in) :: subname
#  include "check_matsize.proced.f90"
end subroutine check_matsize_d
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine check_matsize_c( nsize1, nsize2, matrix, matname, subname )
   implicit none
   integer         , intent(in) :: nsize1
   integer         , intent(in) :: nsize2
   complex*8       , intent(in) :: matrix(:,:)
   character(len=*), intent(in) :: matname
   character(len=*), intent(in) :: subname
#  include "check_matsize.proced.f90"
end subroutine check_matsize_c
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine check_matsize_z( nsize1, nsize2, matrix, matname, subname )
   implicit none
   integer         , intent(in) :: nsize1
   integer         , intent(in) :: nsize2
   complex*16      , intent(in) :: matrix(:,:)
   character(len=*), intent(in) :: matname
   character(len=*), intent(in) :: subname
#  include "check_matsize.proced.f90"
end subroutine check_matsize_z
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
