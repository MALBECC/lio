! Performs a McWeeny purification of rho ON, in order to
! guarantee idempotency.
! This procedure should not be done in the AO basis.

subroutine purify_ON(this)
   implicit none
   class(operator) , intent(inout) :: this

   LIODBLE, allocatable :: temp_r2(:,:), temp_r3(:,:)
   TDCOMPLEX, allocatable :: temp_c2(:,:), temp_c3(:,:)
   

   if (allocated(this%data_ON)) then
      allocate(temp_r2(size(this%data_ON,1), size(this%data_ON,1)) , &
               temp_r3(size(this%data_ON,1), size(this%data_ON,1)))

      temp_r2 = matmul(this%data_ON, this%data_ON)
      temp_r3 = matmul(this%data_ON, temp_r2)

      this%data_ON = 3.0 * temp_r2 - 2.0 * temp_r3

      deallocate(temp_r2, temp_r3)
   endif

   if (allocated(this%dataC_ON)) then
      allocate(temp_c2(size(this%dataC_ON,1), size(this%dataC_ON,1)) , &
               temp_c3(size(this%dataC_ON,1), size(this%dataC_ON,1)))

      temp_c2 = matmul(this%dataC_ON, this%dataC_ON)
      temp_c3 = matmul(this%dataC_ON, temp_c2)

      this%dataC_ON = 3.0 * temp_c2 - 2.0 * temp_c3

      deallocate(temp_c2, temp_c3)
   endif

   

end subroutine purify_ON

