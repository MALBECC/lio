#include "datatypes/datatypes.fh"


function liocmplx(r_part, i_part) result(c_num)
  implicit none
  LIODBLE, intent(in) :: r_part
  LIODBLE, intent(in) :: i_part
  TDCOMPLEX :: c_num

  c_num = CMPLX(real(r_part, COMPLEX_SIZE/2), &
                real(i_part, COMPLEX_SIZE/2), &
                COMPLEX_SIZE/2)
  return 
end function
