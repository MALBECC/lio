! Performs a McWeeny purification of rho ON, in order to
! guarantee idempotency.
! This procedure should not be done in the AO basis.

subroutine purify_ON(this, is_open_shell)
   implicit none
   logical        , intent(in)    :: is_open_shell
   class(operator), intent(inout) :: this

   LIODBLE  , allocatable :: temp_r2(:,:), temp_r3(:,:)
   TDCOMPLEX, allocatable :: temp_c2(:,:), temp_c3(:,:)
   integer :: iteration
   LIODBLE :: factor_a, factor_b

   
   if (is_open_shell) then
      factor_a = 3.0D0
      factor_b = 2.0D0
   else
      ! Taking into account that in closed shell only rho / 2 is idempotent,
      ! factor_a is 3 * 1 / 4 * 2 and factor_b is 2 * 1 / 8 * 2.
      factor_a = 1.5D0
      factor_b = 0.5D0
   endif

   if (allocated(this%data_ON)) then
      allocate(temp_r2(size(this%data_ON,1), size(this%data_ON,1)) , &
               temp_r3(size(this%data_ON,1), size(this%data_ON,1)))

      do iteration = 1, 20
         temp_r2 = matmul(this%data_ON, this%data_ON)
         temp_r3 = matmul(this%data_ON, temp_r2)

         this%data_ON = factor_a * temp_r2 - factor_b * temp_r3
      enddo
      deallocate(temp_r2, temp_r3)
   endif

   if (allocated(this%dataC_ON)) then
      allocate(temp_c2(size(this%dataC_ON,1), size(this%dataC_ON,1)) , &
               temp_c3(size(this%dataC_ON,1), size(this%dataC_ON,1)))

      do iteration = 1, 20
         temp_c2 = matmul(this%dataC_ON, this%dataC_ON)
         temp_c3 = matmul(this%dataC_ON, temp_c2)

         this%dataC_ON = real(factor_a,COMPLEX_SIZE/2) * temp_c2 - &
                         real(factor_b,COMPLEX_SIZE/2) * temp_c3
      enddo
      deallocate(temp_c2, temp_c3)
   endif

end subroutine purify_ON


! Verifies ON-basis operator idempotency property.
subroutine check_idempotency_ON(this, is_open_shell)
   implicit none
   logical        , intent(in)    :: is_open_shell
   class(operator), intent(inout) :: this

   LIODBLE  , allocatable :: temp_r2(:,:)
   TDCOMPLEX, allocatable :: temp_c2(:,:)
   
   integer   :: i_elem, j_elem
   LIODBLE   :: check_val, factor
   check_val = 0.0D0

   factor = 0.5D0
   if (is_open_shell) factor = 1.0D0

   if (allocated(this%data_ON)) then
      allocate(temp_r2(size(this%data_ON,1), size(this%data_ON,1)))

      temp_r2 = factor * matmul(this%data_ON, this%data_ON)
      
      do i_elem = 1, size(this%data_ON,1)
      do j_elem = i_elem, size(this%data_ON,1)
         check_val = abs(temp_r2(i_elem, j_elem) - this%data_ON(i_elem, j_elem))
         
         if (check_val > 1.0E-6) then
            write(*,*) "Warning: Potential idempotency issue in ON matrix."
            write(*,*) "  Index: ", i_elem, j_elem
            write(*,*) "  Values (normal and squared): ", &
                       this%data_ON(i_elem, j_elem),      &
                       temp_r2(i_elem, j_elem)
         endif
      enddo
      enddo
      deallocate(temp_r2)
   endif

   if (allocated(this%dataC_ON)) then
      allocate(temp_c2(size(this%dataC_ON,1), size(this%dataC_ON,1)))

      temp_c2 = real(factor,COMPLEX_SIZE/2) * &
                matmul(this%dataC_ON, this%dataC_ON)

      do i_elem = 1, size(this%dataC_ON,1)
      do j_elem = i_elem, size(this%dataC_ON,1)
         check_val = abs( real(temp_c2(i_elem, j_elem)) - &
                          real(this%dataC_ON(i_elem, j_elem)) )         
      
         if (check_val > 1.0E-6) then
            write(*,*) "Warning: Potential idempotency issue in ON matrix."
            write(*,*) "  Index: ", i_elem, j_elem
            write(*,*) "  Values (normal and squared): ", &
                       this%dataC_ON(i_elem, j_elem),     &
                       temp_c2(i_elem, j_elem)
         endif

      enddo
      enddo

      deallocate(temp_c2)
   endif
end subroutine check_idempotency_ON