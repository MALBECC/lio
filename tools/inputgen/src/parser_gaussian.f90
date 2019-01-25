!##############################################################################!
module parser_gaussian
!##############################################################################!
!
!
!
!
implicit none
contains
!##############################################################################!

subroutine gaufchk_atoms( fname, natoms )
   use auxmod_subs, only: goto_word, safe_open, safe_close
   implicit none
   character(len=*), intent(in)  :: fname
   integer         , intent(out) :: natoms
   
   character(len=100) :: dataline
   character(len=100) :: crap(4)

   call safe_open( 200, fname )
   call goto_word( 200, 'Atomic numbers', dataline )
   read( unit=dataline, fmt=*) crap(1:4), natoms
   call safe_close( 200 )

end subroutine gaufchk_atoms


!##############################################################################!
subroutine gaufchk_basis( fname, mbasis )
   use auxmod_subs, only: goto_word, safe_open, safe_close
   implicit none
   character(len=*), intent(in)  :: fname
   integer         , intent(out) :: mbasis

   character(len=100) :: dataline
   character(len=100) :: crap(5)

   call safe_open( 200, fname )
   call goto_word( 200, 'Number of basis functions', dataline )
   read( unit=dataline, fmt=*) crap(1:5), mbasis
   call safe_close( 200 )

end subroutine gaufchk_basis


!##############################################################################!
subroutine gaufchk_atoms_info( fname, natoms, atomnum, atompos )
   use auxmod_subs, only: goto_word, safe_open, safe_close
   implicit none
   character(len=*), intent(in)  :: fname
   integer         , intent(in)  :: natoms
   integer         , intent(out) :: atomnum(natoms)
   real*8          , intent(out) :: atompos(3,natoms)

   character(len=100)  :: dataline
   integer             :: kline, ki, kf, kk, ii, jj
   real*8, allocatable :: posvec(:)

   real*8, parameter :: ANG_IN_BOHR = 0.529177d0

   call safe_open( 200, fname )

   call goto_word( 200, 'Atomic numbers', dataline )
   do kline = 1, 1 + (natoms / 6 )
      ki = (kline-1) * 6+1
      kf = min( natoms, kline*6 )
      read( unit=200, fmt=* ) atomnum(ki:kf)
   end do

   call goto_word( 200, 'Current cartesian coordinates', dataline )
   allocate( posvec( 3*natoms ) )
   do kline = 1, 1 + ( 3*natoms / 5 )
      ki = (kline-1) * 5 + 1
      kf = min( 3*natoms, kline*5 )
      read( unit=200, fmt=* ) posvec(ki:kf)
   end do

   kk = 0
   do jj = 1, natoms
   do ii = 1, 3
      kk = kk+1
      atompos(ii,jj) = posvec(kk) * ANG_IN_BOHR
   end do
   end do

   call safe_close( 200 )

end subroutine gaufchk_atoms_info


!##############################################################################!
subroutine gaufchk_basis_info( fname, mbasis, atom_ofba, angm_ofba )

   use auxmod_subs, only: goto_word, safe_open, safe_close
   implicit none
   character(len=*), intent(in)  :: fname
   integer         , intent(in)  :: mbasis
   integer         , intent(out) :: atom_ofba(mbasis)
   integer         , intent(out) :: angm_ofba(mbasis)

   character(len=100)   :: dataline, crap(4)

   integer              :: nshells
   integer, allocatable :: shell_mvec(:), shell_nvec(:)
   integer              :: kline, ki, kf, kk, indx, ki_extra

   call safe_open( 200, fname )
   call goto_word( 200, 'Shell types', dataline )
   read( unit=dataline, fmt=*) crap(1:4), nshells
   allocate( shell_mvec(nshells) )
   allocate( shell_nvec(nshells) )

   do kline = 1, 1 + (nshells / 6 )
      ki = (kline-1) * 6+1
      kf = min( nshells, kline*6 )
      read( unit=200, fmt=* ) shell_mvec(ki:kf)
   end do

   call goto_word( 200, 'Shell to atom map', dataline )
   do kline = 1, 1 + (nshells / 6 )
      ki = (kline-1) * 6+1
      kf = min( nshells, kline*6 )
      read( unit=200, fmt=* ) shell_nvec(ki:kf)
   end do


   indx = 0
   do kk = 1, nshells

      select case( shell_mvec(kk) )
         case (0)
            ki_extra = 1
            angm_ofba(indx+1) = 0

         case (1,-1)
            ki_extra = 4
            angm_ofba(indx+1) = 0
            angm_ofba(indx+2) = 1
            angm_ofba(indx+3) = 1
            angm_ofba(indx+4) = 1

         case (2,-2)
            ki_extra = 10
            angm_ofba(indx+ 1) = 0
            angm_ofba(indx+ 2) = 1
            angm_ofba(indx+ 3) = 1
            angm_ofba(indx+ 4) = 1
            angm_ofba(indx+ 5) = 2
            angm_ofba(indx+ 6) = 2
            angm_ofba(indx+ 7) = 2
            angm_ofba(indx+ 8) = 2
            angm_ofba(indx+ 9) = 2
            angm_ofba(indx+10) = 2

         case default
            print*, "ERROR - Unidentified kind of basis..."
            stop

      end select

      do ki = 1, ki_extra
         atom_ofba(indx+ki) = shell_nvec(kk)
      end do
      indx = indx + ki_extra

      if ( indx > mbasis ) then
         print*, "ERROR: discrepancy with number of basis..."
         stop
      end if

   end do

   call safe_close( 200 )

end subroutine gaufchk_basis_info


!##############################################################################!
subroutine gaufchk_densmat( fname, mbasis, elecstate, densmat )

   use auxmod_subs, only: goto_word, safe_open, safe_close
   implicit none
   character(len=*), intent(in)  :: fname
   integer         , intent(in)  :: elecstate
   integer         , intent(in)  :: mbasis
   real*8          , intent(out) :: densmat(mbasis,mbasis)

   character(len=100)  :: dataline, crap(4)
   integer             :: vecsize
   real*8, allocatable :: vecdens(:)
   integer             :: kline, ki, kf, indx, ii, jj

   call safe_open( 200, fname )
   select case( elecstate )
      case (0)
         call goto_word( 200, 'Total SCF Density', dataline )

      case (1)
         call goto_word( 200, 'Total CI Density', dataline )

      case default
         print*, "ERROR - Unidentified kind of electronic state..."
         print*, "elecstate = ", elecstate
         stop

   end select
  
   vecsize = mbasis * (mbasis+1) / 2
   allocate( vecdens(vecsize) )

   do kline = 1, 1 + (vecsize/5)
      ki = (kline-1) * 5+1
      kf = min( vecsize, kline*5 )
      read( unit=200, fmt=* ) vecdens(ki:kf)
   end do

   do jj = 1, mbasis
   do ii = 1, jj
      indx = ii + ( jj * (jj-1) / 2 )
      densmat(ii,jj) = vecdens(indx)
      densmat(jj,ii) = vecdens(indx)
   end do
   end do

   call safe_close( 200 )

end subroutine gaufchk_densmat


!##############################################################################!
end module parser_gaussian
!##############################################################################!
