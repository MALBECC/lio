!##############################################################################!
module parser_orca
!##############################################################################!
!
!
!
!
implicit none
contains
!##############################################################################!

subroutine orca_atoms( fname, natoms )
   use auxmod_subs, only: goto_word, safe_open, safe_close, safe_rewind, &
                        & goto_line
   implicit none
   character(len=*), intent(in)  :: fname
   integer         , intent(out) :: natoms
   
   integer            :: datapos0
   character(len=100) :: dataline

   call safe_open( 200, fname )
   call goto_word( 200, 'INTERNAL COORDINATES (ANGSTROEM)', dataline, datapos0 )
   call safe_rewind( 200 )
   goto_line( 200, datapos0-3 )
   read( unit=200, fmt=*) natoms
   natoms = natoms + 1
   call safe_close( 200 )

end subroutine orca_atoms


!##############################################################################!
subroutine orca_basis( fname, mbasis )
   use auxmod_subs, only: goto_word, safe_open, safe_close, safe_rewind, &
                        & goto_line
   implicit none
   character(len=*), intent(in)  :: fname
   integer         , intent(out) :: mbasis

   integer            :: datapos0
   character(len=100) :: dataline
   character(len=20)  :: scrap(4)

   call safe_open( 200, fname )
   call goto_word( 200, 'Basis Dimension', dataline, datapos0 )
   read( unit=dataline, fmt=*) scrap(4), mbasis
   call safe_close( 200 )

end subroutine orca_basis


!##############################################################################!
subroutine orca_atoms_info( fname, natoms, atomnum, atompos )
   use auxmod_subs, only: goto_word, safe_open, safe_close
   implicit none
   character(len=*), intent(in)  :: fname
   integer         , intent(in)  :: natoms
   integer         , intent(out) :: atomnum(natoms)
   real*8          , intent(out) :: atompos(3,natoms)

   character(len=100)  :: dataline
   character(len=2)    :: atomname
   integer             :: scrap0, kk

   call safe_open( 200, fname )

   call goto_word( 200, 'CARTESIAN COORDINATES (ANGSTROEM)', dataline, scrap0 )
   read( unit=200, fmt='(A)' ) dataline
   do kk = 1, natoms
      read( unit=200, fmt=* ) atomname, atompos(:,kk)
   end do

   call goto_word( 200, 'CARTESIAN COORDINATES (A.U.)', dataline, scrap0 )
   read( unit=200, fmt='(A)' ) dataline
   read( unit=200, fmt='(A)' ) dataline
   do kk = 1, natoms
      read( unit=200, fmt=* ) scrap0, atomname, atomnum(kk)
   end do

   call safe_close( 200 )

end subroutine orca_atoms_info


!##############################################################################!
subroutine orca_basis_info( fname, mbasis, atom_ofba, angm_ofba )

   use auxmod_subs, only: goto_word, safe_open, safe_close
   implicit none
   character(len=*), intent(in)  :: fname
   integer         , intent(in)  :: mbasis
   integer         , intent(out) :: atom_ofba(mbasis)
   integer         , intent(out) :: angm_ofba(mbasis)

   character(len=100)   :: dataline, scrap(4)
   character(len=10)    :: basetype
   integer              :: scrap0
   integer, allocatable :: shell_mvec(:), shell_nvec(:)
   integer              :: kline, ki, kf, kk, indx, ki_extra

   call safe_open( 200, fname )
   call goto_word( 200, 'INITIAL GUESS ORBITALS', dataline, scrap0 )

   do kk = 1, 5
      read( unit=200, fmt='(A)' ) dataline
   end do

   do kk = 1, mbasis

      read( unit=200, fmt='(A)' ) dataline
      read( unit=dataline(1:3), fmt=*) atom_ofba(kk)
      read( unit=dataline     , fmt=*) scrap(1), basetype

      select case( trim(basetype) )
         case("1s","2s","3s","4s","5s","6s")
            angm_ofba(kk) = 0

         case default
            print*, "ERROR - Unidentified kind of basis..."
            print*, " > Basis type = "trim(basetype)
            stop

      end select

   end do

   call safe_close( 200 )

end subroutine orca_basis_info


!##############################################################################!
subroutine orca_densmat( fname, mbasis, elecstate, densmat )

   use auxmod_subs, only: goto_word, safe_open, safe_close
   implicit none
   character(len=*), intent(in)  :: fname
   integer         , intent(in)  :: elecstate
   integer         , intent(in)  :: mbasis
   real*8          , intent(out) :: densmat(mbasis,mbasis)

   character(len=100)  :: dataline, scrap0
   integer             :: total_blocks, nblock, nline, ii, jj

   call safe_open( 200, fname )
   select case( elecstate )
      case (0)
         call goto_word( 200, 'DENSITY', dataline )
         call goto_word( 200, 'DENSITY', dataline )

      case (1)
         print*, "ERROR - Not available for orca yet"
         print*, "elecstate = ", elecstate
         stop

      case default
         print*, "ERROR - Unidentified kind of electronic state..."
         print*, "elecstate = ", elecstate
         stop

   end select
  
   total_block = (mbasis/6)
   if ( mod(mbasis,6) /= 0 ) total_block = total_block + 1

   do nblock = 1, total_blocks
   do nline  = 1, mbasis
      jj_i = (nblock-1) * 6
      jj_f = min( nblock*6, mbasis )
      read( unit=200, fmt=* ) scrap0, densmat( nline, jj_i : jj_f )
   end do
   end do

   call safe_close( 200 )

end subroutine orca_densmat


!##############################################################################!
end module parser_orca
!##############################################################################!
