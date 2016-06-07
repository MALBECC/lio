!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
function is_matrix_square_r( matrix ) result( is_it )
  implicit none
  real(kind=SPK), intent(in)     :: matrix(:,:)
  logical                        :: is_it
  is_it = ( size(matrix,1) == size(matrix,2) )
end function
!------------------------------------------------------------------------------!
function is_matrix_square_d( matrix ) result( is_it )
  implicit none
  real(kind=DPK), intent(in)     :: matrix(:,:)
  logical                        :: is_it
  is_it = ( size(matrix,1) == size(matrix,2) )
end function
!------------------------------------------------------------------------------!
function is_matrix_square_c( matrix ) result( is_it )
  implicit none
  complex(kind=SPK), intent(in)  :: matrix(:,:)
  logical                        :: is_it
  is_it = ( size(matrix,1) == size(matrix,2) )
end function
!------------------------------------------------------------------------------!
function is_matrix_square_z( matrix ) result( is_it )
  implicit none
  complex(kind=DPK), intent(in)  :: matrix(:,:)
  logical                        :: is_it
  is_it = ( size(matrix,1) == size(matrix,2) )
end function
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
