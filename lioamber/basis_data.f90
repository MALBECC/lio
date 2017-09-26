!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module basis_data
   implicit none
   integer                              :: basis_size
   integer                              :: basis_size_s
   integer                              :: basis_size_p
   integer                              :: basis_size_d
   integer                              :: maximum_contractions
   integer, allocatable, dimension(:)   :: orbital_contractions
   integer, allocatable, dimension(:)   :: parent_atom
   integer, allocatable, dimension(:,:) :: angular_momentum
   real*8 , allocatable, dimension(:,:) :: gauss_expo
   real*8 , allocatable, dimension(:,:) :: gauss_coef
end module
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
