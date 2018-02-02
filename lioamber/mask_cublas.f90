!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module mask_cublas
   implicit none
   logical            :: cublas_not_init = .true.
   integer, parameter :: SIZE_OF_REAL = 8
#  ifdef CUBLAS
      integer, external :: CUBLAS_ALLOC
      integer, external :: CUBLAS_SET_MATRIX
      external          :: CUBLAS_INIT
      external          :: CUBLAS_SHUTDOWN
      external          :: CUBLAS_FREE
#  endif
   contains
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine cublas_setmat( Msize, matrix, dev_ptr )
   implicit none
   integer  , intent(in)    :: Msize
   real*8   , intent(in)    :: matrix(Msize, Msize)
   integer*8, intent(inout) :: dev_ptr
   integer                  :: stat

#  ifdef CUBLAS

   if (cublas_not_init) call CUBLAS_INIT()

   stat = CUBLAS_ALLOC( Msize*Msize, SIZE_OF_REAL, dev_ptr)
   if ( stat /= 0 ) then
      write(*,*) "CUBLAS MEMORY ALLOCATION FAILED"
      call CUBLAS_SHUTDOWN
      stop
   endif

   stat = CUBLAS_SET_MATRIX( Msize, Msize, SIZE_OF_REAL, matrix, Msize, &
                           & dev_ptr, Msize )
   if ( stat /= 0 ) then
      write(*,*) "CUBLAS MATRIX SETTING FAILED"
      call CUBLAS_SHUTDOWN
      stop
   endif

#  endif
end subroutine cublas_setmat


!------------------------------------------------------------------------------!
subroutine cublas_release( dev_ptr )
   implicit none
   integer*8, optional, intent(inout) :: dev_ptr

#  ifdef CUBLAS
   if (present(dev_ptr)) then
      call CUBLAS_FREE(dev_ptr)
   else
      cublas_not_init = .true.
      call CUBLAS_SHUTDOWN
   end if
#  endif

end subroutine cublas_release


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
end module mask_cublas
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
