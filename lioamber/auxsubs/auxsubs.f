!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       module auxsubs
!--------------------------------------------------------------------!
!
! README
! 
!   This module contains different kinds of miscelaneous
!   auxilliary subroutines and functions used throughout
!   the main program.
!
!   In order to add a new procedure, one needs only to
!   create a file in the containing folder and then include
!   it in the section for procedures using the following
!   syntax:
!
!     > include 'file.f'
!
!   If one has many variations of one procedure, one can
!   group all of them together inside an interface by 
!   creating another file where the interface is specified
!   and then including this file in the section for headers
!   using the same syntax as before.
!
!   One has to remember to be specially carefull when 
!   introducing the parameters of a subroutine that is
!   overloaded (part of an interface). The most advisable
!   thing to do is to always use variables to pass the
!   arguments instead of doing so directly.
!
!--------------------------------------------------------------------!
! INCLUDE FILES WITH HEADERS
       include 'commutate_h.f'
       include 'gaussbell_h.f'
       contains
!--------------------------------------------------------------------!
! INCLUDE FILES WITH PROCEDURES
       include 'commutate.f'
       include 'gaussbell.f'
       end module
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
