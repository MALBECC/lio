!##############################################################################!
program inputgen
!##############################################################################!
!
! This program automates the creation of input files for lio.
!
! Currently its only feature is translation from gaussian and orca
! restart outputs into lio restart inputs. Features to be added in
! the future.
!
! The structure works as follow:
!
!  - keywords module contains the keywords used throughout the code.
!    These should be set by the input and never changed.
!
!  - data_lio stores the data that will be used to make the lio input.
!
!  - data_src stores the data 
!
!
! Code originally written by gonzalodm.
! Adapted by FFR as a stand-alone.
!
!##############################################################################!
   use keywords    , only: reads_keywords, write_keywords
   use data_control, only: reads_data    , write_data    , shape_data

   implicit none
   character(len=80) :: input_fname

   write( unit=6, fmt='(A)' ) &
   & "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
   write( unit=6, fmt='(A)' ) ">  LIO INPUT GENERATOR  <"
   write( unit=6, fmt='(A)' ) &
   & "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
   write( unit=6, fmt='(A)' ) ""

   print*, " > Starting procedure reads_keywords..."; print*

   if ( command_argument_count() > 0 ) then
      call get_command_argument( 1, input_fname )
      call reads_keywords( input_fname )
   else 
      print*, "WARNING: No input provided; using default keywords..."
      print*
   end if
   print*, " > Normal termination of reads_keywords."; print*
   print*, " > Writing keywords..."; print*
   
   call write_keywords(6)

   print*, " > Starting procedure reads_data..."; print*
   call reads_data()
   print*, " > Normal termination of reads_data."; print*

   print*, " > Starting procedure shape_data..."; print*
   call shape_data()
   print*, " > Normal termination of shape_data."; print*

   print*, " > Starting procedure write_data..."; print*
   call write_data()
   print*, " > Normal termination of write_data."; print*

end program inputgen

!##############################################################################!
