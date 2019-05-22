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

   type cumat_r
      integer         :: mat_size
      integer(kind=8) :: cu_pointer
      real(kind=8), allocatable :: matrix(:,:)


   contains
      procedure, pass :: init        => initialise_r
      procedure, pass :: exterminate => exterminate_r
      procedure, pass :: multiply    => multiply_r

   end type cumat_r

contains
#include "cumat_init_fin.f90"
#include "cumat_multiply.f90"

end module typedef_cumat
