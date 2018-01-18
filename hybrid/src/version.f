module version_info

implicit none

! This file MUST be updated after every self-consistent commit,
! and the PL ("patch level") number increased by one, unless the
! modification involves raising a minor or major version number,
! in which case the PL should be reset to zero.

! A self-consistent commit is a group of changes that fix a bug
! or implement a new feature, in such a way that the program can
! be compiled (no loose ends left). An update to the CHANGES file
! should be an integral part of a commit (the PL number should be
! included for reference.)

! After it is done, this file should be commited.

integer, dimension(3), save  :: num_version = (/2,0,0/)
character(len=80), parameter :: version_str =  &
"HYBRID 2.0 -- [ Lio-hybrid ] (22 Dec 2017)" 

end module version_info
!================================================================

subroutine prversion

! Simple routine to print the version string. Could be extended to
! provide more information, if needed.

! Use free format in file to make more room for long option strings...

use version_info, only: version_str
implicit none

write(6,'(a)') trim(version_str)
!write(6,'(2a)') 'Architecture  : ', &
!"SIESTA_ARCH"
write(6,'(2a)') 'Compiler flags: ', &
"FFLAGS"
!#ifdef MPI
!write(6,'(a)') 'PARALLEL version'
!#else
write(6,'(a)') 'SERIAL version'
!#endif

!#ifdef CDF
!write(6,'(a)') 'NetCDF-capable'
!#endif

end subroutine prversion
!----------------------------------------------------------

subroutine get_version(v)
  use version_info, only: num_version
  implicit none
  integer, intent(out)  :: v(3)
  v = num_version
end subroutine get_version

