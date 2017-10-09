!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module find_free_unit_auxmod
implicit none
contains
#include "../find_free_unit.f90"
end module find_free_unit_auxmod

program find_free_unit_ut
   use find_free_unit_auxmod
   implicit none
   logical :: test_wins(4)
   integer :: free_unit
   integer :: last_unit
   integer :: call_stat
!
!
!
! Test 1: Chooses first option
!------------------------------------------------------------------------------!
    free_unit = 20
    last_unit = 22
    call_stat = 0
    call find_free_unit( free_unit, last_unit, call_stat )
    if ( (free_unit == 20) .and. (call_stat == 0) ) then
       print*,"UTEST 1 - OK"
    else
       print*,"UTEST 1 - FAILED!"
       print*,"free_unit (=20): ",free_unit
       print*,"call_stat (=0):  ",call_stat
    end if
!
!
! Test 2: Skips opened file
!------------------------------------------------------------------------------!
    free_unit = 20
    last_unit = 22
    call_stat = 0
    open(file='Output-empty20.o', unit=20)
    call find_free_unit( free_unit, last_unit, call_stat )
    close(unit=20)
    if ( (free_unit == 21) .and. (call_stat == 0) ) then
       print*,"UTEST 2 - OK"
    else
       print*,"UTEST 2 - FAILED!"
       print*,"free_unit (=21): ",free_unit
       print*,"call_stat (=0):  ",call_stat
    end if 

!
!
! Test 3: Error message if all files are open
!------------------------------------------------------------------------------!
    free_unit = 20
    last_unit = 22
    call_stat = 0
    open(file='Output-empty20.o',unit=20)
    open(file='Output-empty21.o',unit=21)
    open(file='Output-empty22.o',unit=22)
    call find_free_unit( free_unit, last_unit, call_stat )
    close(unit=20)
    close(unit=21)
    close(unit=22)
    if ( (free_unit == 23) .and. (call_stat == -1) ) then
       print*,"UTEST 3 - OK"
    else
       print*,"UTEST 3 - FAILED!"
       print*,"free_unit (=22): ",free_unit
       print*,"call_stat (/=0): ",call_stat
    end if 

!
!
! Test 3: Error message if range doesn't make sense
!------------------------------------------------------------------------------!
    free_unit = 22
    last_unit = 20
    call_stat = 0
    call find_free_unit( free_unit, last_unit, call_stat )
    if ( (free_unit == 22) .and. (call_stat == -1) ) then
       print*,"UTEST 4 - OK"
    else
       print*,"UTEST 4 - FAILED!"
       print*,"free_unit (=22): ",free_unit
       print*,"call_stat (/=0): ",call_stat
    end if 
!
!
!------------------------------------------------------------------------------!
end program find_free_unit_ut
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
