!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module build_info
!
! This module stores the data of the compilation options and version of the
! code. This allows to record in the output files the information of which
! version of lio was used and how was it compiled, in order to improve and
! facilitate reproducibility of results.
!
! Information is passed into local constants ( parameters named BDATA_XXX )
! through variables defined during compilation ( BPASS_XXX ).
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
#  ifdef BPASS_VERSION
      character(len=40), parameter :: BDATA_VERSION   = BPASS_VERSION
#  else
      character(len=40), parameter :: BDATA_VERSION   = 'NO DATA'
#     define INCOMPLETE_BDATA
#  endif


#  ifdef BPASS_GITHASH
      character(len=40), parameter :: BDATA_GITHASH   = BPASS_GITHASH
#  else
      character(len=40), parameter :: BDATA_GITHASH   = 'NO DATA'
#     define INCOMPLETE_BDATA
#  endif


#  ifdef BPASS_CUDA
      integer,           parameter :: BDATA_CUDA      = BPASS_CUDA
#  else
      integer,           parameter :: BDATA_CUDA      = -666
#     define INCOMPLETE_BDATA
#  endif


#  ifdef BPASS_INTEL
      integer,           parameter :: BDATA_INTEL     = BPASS_INTEL
#  else
      integer,           parameter :: BDATA_INTEL     = -666
#     define INCOMPLETE_BDATA
#  endif


#  ifdef BPASS_PARALLEL
      integer,           parameter :: BDATA_PARALLEL  = BPASS_PARALLEL
#  else
      integer,           parameter :: BDATA_PARALLEL  = -666
#     define INCOMPLETE_BDATA
#  endif


#  ifdef BPASS_PRECISION
      integer,           parameter :: BDATA_PRECISION = BPASS_PRECISION
#  else
      integer,           parameter :: BDATA_PRECISION = -666
#     define INCOMPLETE_BDATA
#  endif


#  ifdef BPASS_ANALYTICS
      integer,           parameter :: BDATA_ANALYTICS = BPASS_ANALYTICS
#  else
      integer,           parameter :: BDATA_ANALYTICS = -666
#     define INCOMPLETE_BDATA
#  endif


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
contains

   subroutine print_build_info( outunit )
      implicit none
      integer, intent(in) :: outunit
      write( unit=outunit, fmt= '(A)') ''
      write( unit=outunit, fmt= '(2A)')   '  LIO VERSION  : ', BDATA_VERSION
      write( unit=outunit, fmt= '(A)') ''
      write( unit=outunit, fmt= '(2A)')   '  GIT HASHTAG  : ', BDATA_GITHASH
      write( unit=outunit, fmt= '(A)')    ''
      write( unit=outunit, fmt= '(2A)')   '  COMPILATION OPTIONS '
      write( unit=outunit, fmt= '(A,I2)') '  *  cuda      = ', BDATA_CUDA
      write( unit=outunit, fmt= '(A,I2)') '  *  intel     = ', BDATA_INTEL
      write( unit=outunit, fmt= '(A,I2)') '  *  parallel  = ', BDATA_PARALLEL
      write( unit=outunit, fmt= '(A,I2)') '  *  precision = ', BDATA_PRECISION
      write( unit=outunit, fmt= '(A,I2)') '  *  analytics = ', BDATA_ANALYTICS
      write( unit=outunit, fmt= '(A)') ''
#     ifdef INCONSISTENT_BUILD
         write( unit=outunit, fmt= '(A)') &
         '  WARNING! INCONSISTENT COMPILATION:'
         write( unit=outunit, fmt= '(A)') &
         '  Part of the code has been compiled using a certain set'
         write( unit=outunit, fmt= '(A)') &
         '  of options, while another part has been compiled with'
         write( unit=outunit, fmt= '(A)') &
         '  a different set. Please make clean between compilations'
         write( unit=outunit, fmt= '(A)') &
         '  with different options to have a more reliable run.'
         write( unit=outunit, fmt= '(A)') ''
#     endif
#     ifdef INCOMPLETE_BDATA
         write( unit=outunit, fmt= '(A)') &
         '  WARNING! INCOMPLETE DATA:'
         write( unit=outunit, fmt= '(A)') &
         '  There has been a problem with the data passed during'
         write( unit=outunit, fmt= '(A)') &
         '  compilation. This is a bug and needs to be reported'
         write( unit=outunit, fmt= '(A)') &
         '  to us in github.com/MALBECC/lio.'
         write( unit=outunit, fmt= '(A)') ''
#     endif
      write( unit=outunit, fmt= '(A)')  &
      '============================================================'

   end subroutine

end module
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
