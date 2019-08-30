!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
#include "../complex_type.fh"
subroutine movieprint( natoms, mbasis, nstep, nucids, nucpos, elecdens )

   use fileio_data, only: movie_count, movie_nfreq, movie_name0
   implicit none
   integer        , intent(in) :: natoms
   integer        , intent(in) :: mbasis
   integer        , intent(in) :: nstep
   integer        , intent(in) :: nucids(natoms)
   real(kind=8)   , intent(in) :: nucpos(natoms,3)
   complex(kind=8), intent(in) :: elecdens(mbasis,mbasis)

   character(len=3), allocatable :: nucids_name(:)
   character(len=4)  :: charnum
   character(len=51) :: filename_n, filename_e
   integer           :: ndig1, ndig2, ndig3, ndig4
   character(len=1)  :: cdig1, cdig2, cdig3, cdig4

   if ( movie_nfreq == 0 ) return
   if ( mod(nstep,movie_nfreq) /= 0 ) return

   if ( movie_count > 9999 ) then
       print*, "Movie out of range 9999 - not printing anymore..."
       return
   end if

   ndig1 = mod( movie_count, 10 )
   ndig2 = mod( movie_count, 100 )
   ndig3 = mod( movie_count, 1000 )
   ndig4 = movie_count

   ndig2 = ndig2 / 10
   ndig3 = ndig3 / 100
   ndig4 = ndig4 / 1000
   
   write( unit=cdig1, fmt='(I1)' ) ndig1
   write( unit=cdig2, fmt='(I1)' ) ndig2
   write( unit=cdig3, fmt='(I1)' ) ndig3
   write( unit=cdig4, fmt='(I1)' ) ndig4

   charnum = cdig4 // cdig3 // cdig2 // cdig1
   filename_n = trim(adjustl(movie_name0)) // "_nu" // charnum // ".out"
   filename_e = trim(adjustl(movie_name0)) // "_el" // charnum // ".out"

   open(unit = 123456, file=filename_n)
   allocate( nucids_name(natoms) )
   call translate_atomlist( natoms, nucids, nucids_name )
   call write_nucpos( natoms, nucids_name, nucpos, 123456 )
   deallocate( nucids_name )
   close(unit = 123456)

   open(unit = 123456, file=filename_e)
   call write_rho_restart( elecdens, mbasis, 123456)
   close(unit = 123456)

!  Increase movie counter. First one will be the ground state, 0000.
   movie_count = movie_count + 1

end subroutine movieprint
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
