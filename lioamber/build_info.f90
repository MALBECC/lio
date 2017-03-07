!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module build_info
!
! This module stores the data of the compilation options and version of the
! code. This allows to record in the output files the information of which
! version of lio was used and how was it compiled, in order to improve and
! facilitate reproducibility of results.
!
! Information is passed into local constants ( parameters named BUILD_XXX )
! through variables defined at compilation ( BPASS_XXX ).
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
#  ifdef BPASS_VERSION
      character(len=10), parameter :: BUILD_VERSION   = BPASS_VERSION
#  else
      character(len=10), parameter :: BUILD_VERSION   = 'NO DATA'
#     define NON_REFERABLE_BUILD
#  endif

!------------------------------------------------------------------------------!
#  ifdef BPASS_GITHASH
      character(len=40), parameter :: BUILD_GITHASH   = BPASS_GITHASH
#  else
      character(len=40), parameter :: BUILD_GITHASH   = 'NO DATA'
#     define NON_REFERABLE_BUILD
#  endif

!------------------------------------------------------------------------------!
#  ifdef BPASS_CUDA
      integer,           parameter :: BUILD_CUDA      = BPASS_CUDA
#  else
      integer,           parameter :: BUILD_CUDA      = -1
#     define NON_REFERABLE_BUILD
#  endif

!------------------------------------------------------------------------------!
#  ifdef BPASS_INTEL
      integer,           parameter :: BUILD_INTEL     = BPASS_INTEL
#  else
      integer,           parameter :: BUILD_INTEL     = -1
#     define NON_REFERABLE_BUILD
#  endif

!------------------------------------------------------------------------------!
#  ifdef BPASS_PROFDEB
      integer,           parameter :: BUILD_PROFDEB   = BPASS_PROFDEB
#  else
      integer,           parameter :: BUILD_PROFDEB   = -1
#     define NON_REFERABLE_BUILD
#  endif

!------------------------------------------------------------------------------!
#  ifdef BPASS_PARALLEL
      integer,           parameter :: BUILD_PARALLEL  = BPASS_PARALLEL
#  else
      integer,           parameter :: BUILD_PARALLEL  = -1
#     define NON_REFERABLE_BUILD
#  endif

!------------------------------------------------------------------------------!
#  ifdef BPASS_PRECISION
      integer,           parameter :: BUILD_PRECISION = BPASS_PRECISION
#  else
      integer,           parameter :: BUILD_PRECISION = -1
#     define NON_REFERABLE_BUILD
#  endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
contains

   subroutine print_build_info( outunit )
      implicit none
      integer, intent(in) :: outunit
      write( unit=outunit, fmt= '(A)') &
      '============================================================'
      write( unit=outunit, fmt= '(2A)')   '  LIO VERSION  : ', BUILD_VERSION
      write( unit=outunit, fmt= '(2A)')   '  GIT HASHTAG  : ', BUILD_GITHASH
      write( unit=outunit, fmt= '(A)')    ''
      write( unit=outunit, fmt= '(2A)')   '  COMPILATION OPTIONS '
      write( unit=outunit, fmt= '(A,I2)') '  *  cuda      = ', BUILD_CUDA
      write( unit=outunit, fmt= '(A,I2)') '  *  intel     = ', BUILD_INTEL
      write( unit=outunit, fmt= '(A,I2)') '  *  profdeb   = ', BUILD_PROFDEB
      write( unit=outunit, fmt= '(A,I2)') '  *  parallel  = ', BUILD_PARALLEL
      write( unit=outunit, fmt= '(A,I2)') '  *  precision = ', BUILD_PRECISION
#     ifdef NON_REFERABLE_BUILD
         write( unit=outunit, fmt= '(A)') ''
         write( unit=outunit, fmt= '(A)') &
         '  WARNING! NON REFERABLE COMPILATION'
         write( unit=outunit, fmt= '(A)') ''
#     endif
      write( unit=outunit, fmt= '(A)')  &
      '============================================================'
   end subroutine

end module
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
