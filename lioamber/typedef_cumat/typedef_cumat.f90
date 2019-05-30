! Main CUblas MATrix datatype file. Handles CPU/GPU duality for certain
! matrices(mainly basechange and TD magnus). It is also posible to allocate
! a GPU-only matrix, but this only works in GPU compilations (cuda > 0).
!
! WARNING: Unless specifically updated, only the GPU allocation of the matrix
!          contains its most recent values. In case this becomes important,
!          use the %update method, which copies the GPU values into the CPU
!          matrix (if stored).
!
#include "../complex_type.fh"
module typedef_cumat
   implicit none

#  ifdef CUBLAS
   logical            :: cublas_initialised = .false.
   integer, parameter :: REAL_SIZE = 8
   integer, external  :: cublas_alloc
   integer, external  :: cublas_set_matrix
   integer, external  :: cublas_get_matrix
   integer, external  :: cublas_xcopy
   integer, external  :: cublas_dcopy
   external           :: cublas_init
   external           :: cublas_shutdown
   external           :: cublas_free
#  endif

   ! Basic properties.
   type :: cumat_basic
      integer         :: mat_size  = 0
      integer(kind=8) :: cu_pointer
      logical         :: allocated = .false.
      logical         :: gpu_only  = .false.
   contains
      ! Property encapsulation.
      procedure, pass :: get_size
      procedure, pass :: gpu_pointer
      procedure, pass :: is_allocated
      procedure, pass :: is_gpu_only
   endtype cumat_basic

   ! For real-type matrices.
   type, extends(cumat_basic) :: cumat_r
      real(kind=8), allocatable :: matrix(:,:)
   contains
      ! Setup methods and the like.
      procedure, pass :: allocate    => allocate_r
      procedure, pass :: init        => initialise_r
      procedure, pass :: destroy     => destroy_r
      procedure, pass :: set_r
      procedure, pass :: set_rp
      generic         :: set         => set_r
      generic         :: set         => set_rp
      procedure, pass :: get_r
      procedure, pass :: get_rp
      generic         :: get         => get_r
      generic         :: get         => get_rp

      ! External operations.
      procedure, pass :: multiply    => multiply_r
      procedure, pass :: change_base => change_base_rr
   endtype cumat_r

   ! For complex-type matrices.
   type, extends(cumat_basic) :: cumat_x
      TDCOMPLEX, allocatable :: matrix(:,:)
   contains
      ! Setup methods and the like.
      procedure, pass :: allocate    => allocate_x
      procedure, pass :: init        => initialise_x
      procedure, pass :: destroy     => destroy_x
      procedure, pass :: set_x
      procedure, pass :: set_xp
      generic         :: set         => set_x
      generic         :: set         => set_xp
      procedure, pass :: get_x
      procedure, pass :: get_xp
      generic         :: get         => get_x
      generic         :: get         => get_xp

      ! External operations.
      procedure, pass :: multiply    => multiply_x
      procedure, pass :: change_base => change_base_xx
   endtype cumat_x

contains
#include "cumat_init_fin.f90"
#include "cumat_multiply.f90"
#include "cumat_bchange.f90"
#include "cumat_set.f90"
#include "cumat_get.f90"
#include "cumat_misc.f90"
end module typedef_cumat
