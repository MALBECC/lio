! Changes the base of the input matrix, using the cumat object as the
! basechange matrix. Mode indicates whether the change is direct (dir)
! or inverse (inv).
subroutine change_base_rr(this, input_matrix, mode)
#ifdef CUBLAS
   use cublasmath, only: basechange_cublas
#else
   use mathsubs  , only: basechange_gemm
#endif
   use liosubs_math
   implicit none
   class(cumat_r)  , intent(in)    :: this
   character(len=3), intent(in)    :: mode
   real(kind=8)    , intent(inout) :: input_matrix(:,:)

#ifdef CUBLAS
   input_matrix = basechange_cublas(size(input_matrix,1), input_matrix, &
                                    this%cu_pointer, mode)
#else
   input_matrix = basechange_gemm(size(input_matrix,1), input_matrix, &
                                  this%matrix, mode)
#endif

end subroutine change_base_rr

! subroutine change_base_rx(this, input_matrix, mode)
! #ifdef CUBLAS
!    use cublasmath, only: basechange_cublas
! #else
!    use mathsubs  , only: basechange_gemm
! #endif
!    implicit none
!    class(cumat_r)  , intent(in)    :: this
!    character(len=3), intent(in)    :: mode
!    TDCOMPLEX       , intent(inout) :: input_matrix(:,:)

! #ifdef CUBLAS
!    input_matrix = basechange_cublas(size(input_matrix,1), input_matrix, &
!                                     this%cu_pointer, mode)
! #else
!    input_matrix = basechange_gemm(size(input_matrix,1), input_matrix, &
!                                   this%matrix)
! #endif

! end subroutine change_base_rx

! subroutine change_base_xr(this, input_matrix, mode)
! #ifdef CUBLAS
!    use cublasmath, only: basechange_cublas
! #else
!    use mathsubs  , only: basechange_gemm
! #endif
!    implicit none
!    class(cumat_x)  , intent(in)    :: this
!    character(len=3), intent(in)    :: mode
!    real(kind=8)    , intent(inout) :: input_matrix(:,:)
!    TDCOMPLEX       , allocatable   :: aux_matrix(:,:)

!    integer :: ii, jj

!    allocate(aux_matrix(size(input_matrix,1), size(input_matrix,2)))
!    do ii = 1, size(input_matrix,1)
!    do jj = 1, size(input_matrix,2)
!       aux_matrix(ii,jj) = cmplx(input_matrix(ii,jj), 0.0D0)
!    enddo
!    enddo

! #ifdef CUBLAS
!    aux_matrix = basechange_cublas(size(input_matrix,1), aux_matrix, &
!                                   this%cu_pointer, mode)
! #else
!    aux_matrix = basechange_gemm(size(input_matrix,1), aux_matrix, &
!                                 this%matrix)
! #endif

!    do ii = 1, size(input_matrix,1)
!    do jj = 1, size(input_matrix,2)
!       input_matrix(ii,jj) = dble(aux_matrix(ii,jj))
!    enddo
!    enddo

!    deallocate(aux_matrix)
! end subroutine change_base_xr

subroutine change_base_xx(this, input_matrix, mode)
#ifdef CUBLAS
   use cublasmath, only: basechange_cublas
#else
   use mathsubs  , only: basechange_gemm
#endif
   implicit none
   class(cumat_x)  , intent(in)    :: this
   character(len=3), intent(in)    :: mode
   TDCOMPLEX       , intent(inout) :: input_matrix(:,:)

#ifdef CUBLAS
   input_matrix = basechange_cublas(size(input_matrix,1), input_matrix, &
                                    this%cu_pointer, mode)
#else
   input_matrix = basechange_gemm(size(input_matrix,1), input_matrix, &
                                  this%matrix, mode)
#endif

end subroutine change_base_xx