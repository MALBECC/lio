#include "complex_type.fh"
module typedef_cumat
   implicit none

#  ifdef CUBLAS
   logical            :: cublas_initialised = .false.
   integer, parameter :: SIZE_OF_REAL = 8
   integer, external  :: cublas_alloc
   integer, external  :: cublas_set_matrix
   external           :: cublas_init
   external           :: cublas_shutdown
   external           :: cublas_free
#  endif

   ! For real matrices.
   type cumat_r
      integer         :: mat_size
      integer(kind=8) :: cu_pointer
      real(kind=8), allocatable :: matrix(:,:)
   contains
      procedure, pass :: init        => initialise_r
      procedure, pass :: destroy     => destroy_r
      procedure, pass :: multiply    => multiply_r
      procedure, pass :: change_base => change_base_r
   end type cumat_r

   ! For complex matrices.
   type cumat_x
      integer         :: mat_size
      integer(kind=8) :: cu_pointer
      TDCOMPLEX, allocatable :: matrix(:,:)
   contains
      procedure, pass :: init        => initialise_x
      procedure, pass :: destroy     => destroy_x
      procedure, pass :: multiply    => multiply_x
      procedure, pass :: change_base => change_base_x
   end type cumat_x

   interface change_base_r
      module procedure change_base_rr
      module procedure change_base_rx
   end interface change_base_r

   interface change_base_x
      module procedure change_base_xr
      module procedure change_base_xx
   end interface change_base_x
contains
#include "cumat_init_fin.f90"
#include "cumat_multiply.f90"
#include "cumat_bchange.f90"

end module typedef_cumat
