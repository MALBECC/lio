!##############################################################################!
module auxmod_subs
   implicit none
   contains


!##############################################################################!

subroutine safe_open( file_unit, file_name )
   implicit none
   integer         , intent(in) :: file_unit
   character(len=*), intent(in) :: file_name
   integer                      :: mystat

   open( unit=file_unit, file=file_name, iostat=mystat )
   if ( mystat /= 0 ) then
      print*, "ERROR: safe_open"
      print*, "  could not open file ", file_name
      print*, "  iostat = ", mystat
      stop
   end if

end subroutine safe_open

!------------------------------------------------------------------------------!

subroutine safe_close( file_unit )
   implicit none
   integer         , intent(in) :: file_unit
   integer                      :: mystat

   close( unit=file_unit, iostat=mystat )
   if ( mystat /= 0 ) then
      print*, "ERROR: safe_close"
      print*, "  could not close unit ", file_unit
      print*, "  iostat = ", mystat
      stop
   end if

end subroutine safe_close

!------------------------------------------------------------------------------!

subroutine safe_rewind( file_unit )
   implicit none
   integer, intent(in) :: file_unit
   integer             :: mystat

   rewind( unit=file_unit, iostat=mystat )
   if ( mystat /= 0 ) then
      print*, "ERROR: safe_rewind"
      print*, "  could not rewind unit ", file_unit
      print*, "  iostat = ", mystat
      stop
   end if

end subroutine safe_rewind

!------------------------------------------------------------------------------!

subroutine goto_line( funit, nlines )
   implicit none
   integer         , intent(in)  :: funit
   integer         , intent(in)  :: nlines
   integer                       :: nn
   character(len=80)             :: dataline

   do nn = 1, nlines - 1
      read( unit=funit, fmt='(A)' ) dataline
   end do

end subroutine goto_line

!------------------------------------------------------------------------------!

subroutine goto_word( funit, dataword, dataline, wordline )
   implicit none
   integer         , intent(in)  :: funit
   character(len=*), intent(in)  :: dataword
   character(len=*), intent(out) :: dataline
   integer         , intent(out) :: wordline

   logical :: keep_looking

   wordline = 0
   keep_looking = .true.
   do while (keep_looking)
      wordline = wordline + 1
      read( unit=funit, fmt='(A)' ) dataline
      if ( index( dataline, trim(dataword) ) > 0 ) then
         keep_looking = .false.
      end if
   end do

end subroutine goto_word

!------------------------------------------------------------------------------!

subroutine get_ordered_vec( vecsize, vecvals, vecidns )
   implicit none
   integer, intent(in)  :: vecsize
   integer, intent(in)  :: vecvals( vecsize )
   integer, intent(out) :: vecidns( vecsize )

   integer :: vecmax, vecmin, nowval
   integer :: nowidn, newidn

   vecmax = vecvals(1)
   vecmin = vecvals(1)
   do nowidn = 2, vecsize
      if ( vecvals(nowidn) > vecmax ) vecmax = vecvals(nowidn)
      if ( vecvals(nowidn) < vecmin ) vecmin = vecvals(nowidn)
   end do

   newidn = 0
   do nowval = vecmin, vecmax
      do nowidn = 1, vecsize
         if ( vecvals(nowidn) == nowval ) then
            newidn = newidn + 1
            vecidns(newidn) = nowidn
         end if
      end do
   end do

end subroutine get_ordered_vec

!------------------------------------------------------------------------------!

subroutine update_idref( nvec1, nvec2, old_id1, id1_of_id2 )
   implicit none
   integer, intent(in)    :: nvec1
   integer, intent(in)    :: nvec2
   integer, intent(in)    :: old_id1( nvec1 )
   integer, intent(inout) :: id1_of_id2( nvec2 )
   integer                :: id1, id2, this_id1

   do id2 = 1, nvec2

      this_id1 = id1_of_id2(id2)

      do id1 =1, nvec1
         if ( this_id1 == old_id1(id1) ) id1_of_id2(id2) = id1
      end do

   end do

end subroutine update_idref

!------------------------------------------------------------------------------!

subroutine reorder_int1vec( vecsize, vecidns, vecvals )
   implicit none
   integer, intent(in)    :: vecsize
   integer, intent(in)    :: vecidns( vecsize )
   integer, intent(inout) :: vecvals( vecsize )

   integer, allocatable :: tempvec(:)
   integer              :: kk

   allocate( tempvec(vecsize) )

   do kk = 1, vecsize
      tempvec( kk ) = vecvals( vecidns(kk) )
   end do

   vecvals(:) = tempvec(:)
   deallocate( tempvec )

end subroutine reorder_int1vec

!------------------------------------------------------------------------------!

subroutine reorder_dbl3vec( vecsize, vecidns, vecvals )
   implicit none
   integer, intent(in)    :: vecsize
   integer, intent(in)    :: vecidns( vecsize )
   real*8 , intent(inout) :: vecvals( 3, vecsize )

   real*8 , allocatable :: tempvec(:,:)
   integer              :: kk

   allocate( tempvec(3,vecsize) )

   do kk = 1, vecsize
      tempvec( :, kk ) = vecvals( :, vecidns(kk) )
   end do

   vecvals(:,:) = tempvec(:,:)
   deallocate( tempvec )

end subroutine reorder_dbl3vec

!------------------------------------------------------------------------------!

subroutine reorder_dblmat( matsize, vecidns, matvals )
   implicit none
   integer, intent(in)    :: matsize
   integer, intent(in)    :: vecidns( matsize )
   real*8 , intent(inout) :: matvals( matsize, matsize )

   real*8 , allocatable :: tempmat(:,:)
   integer              :: ii, jj

   allocate( tempmat(matsize,matsize) )

   do jj = 1, matsize
   do ii = 1, matsize
      tempmat( ii, jj ) = matvals( vecidns(ii), vecidns(jj) )
   end do
   end do

   matvals(:,:) = tempmat(:,:)
   deallocate( tempmat )

end subroutine reorder_dblmat


!##############################################################################!
end module auxmod_subs
!##############################################################################!
