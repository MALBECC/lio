!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% WRITE_OUTPUT.F90 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! This file contains several output-file printing routines. Currently includes:!
! * atom_name        (gets atomic number and returns atomic symbol)            !
! * write_dipole_td  (handles dipole moment printing in TD)                    !
! * write_forces     (handles grandient printing to output)                    !
! * write_force_log  (prints forces components to a Forces.log file)           !
! * write_orbitals   (prints orbitals and energies to output)                  !
! * write_orbitals_op(prints orbitals and energies to output, open shell)      !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%% ATOM_NAME %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Takes atomic number Z and translates it to its symbol.                       !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine atom_name(atom_Z, symb)
 implicit none
 integer         , intent(in)  :: atom_Z
 character(LEN=3), intent(out) :: symb

 character(LEN=3) :: name(118)
 name = (/'H  ','HE ','LI ','BE ','B  ','C  ','N  ','O  ','F  ','NE ','NA ', &
          'MG ','AL ','SI ','P  ','S  ','CL ','AR ','K  ','CA ','SC ','TI ', &
          'V  ','CR ','MN ','FE ','CO ','NI ','CU ','ZN ','GA ','GE ','AS ', &
          'SE ','BR ','KR ','RB ','SR ','Y  ','ZR ','NB ','MO ','TC ','RU ', &
          'RH ','PD ','AG ','CD ','IN ','SN ','SB ','TE ','I  ','XE ','CS ', &
          'BA ','LA ','CE ','PR ','ND ','PM ','SM ','EU ','GD ','TB ','DY ', &
          'HO ','ER ','TM ','YB ','LU ','HF ','TA ','W  ','RE ','OS ','IR ', &
          'PT ','AU ','HG ','TL ','PB ','BI ','PO ','AT ','RN ','FR ','RA ', &
          'AC ','TH ','PA ','U  ','NP ','PU ','AM ','CM ','BK ','CF ','ES ', &
          'FM ','MD ','NO ','LR ','RF ','DB ','SG ','BH ','HS ','MT ','DS ', &
          'UUU','UUB','UUT','UUQ','UUP','UUH','UUS','UUO'/)
 symb = name(atom_Z)

end subroutine atom_name
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% WRITE_DIPOLE_TD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Prints the dipole momment vector in TD calculations, where dipxyz is the     !
! dipole moment vector, time is the current time (in fs) in TD, and uid is the !
! output file UID.                                                             !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine write_dipole_td(dipxyz, time, uid)
   implicit none
   LIODBLE, intent(in) :: dipxyz(3), time
   integer         , intent(in) :: uid

   write(UID,100) time, dipxyz(1), dipxyz(2), dipxyz(3)

   return
100 format (e15.8,' ', e15.8,' ',e15.8,' ',e15.8)
end subroutine write_dipole_td
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%% WRITE_DIPOLE_TD_HEADER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Prints the dipole moment file's header in TD calculations, where time_step is!
! the one used in TD, fx-fy-fz are the perturbation field coordinates, and uid !
! is the output UID.                                                           !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine write_dipole_td_header(time_step, fx, fy, fz, uid)
   implicit none
   LIODBLE, intent(in) :: time_step, fx, fy, fz
   integer         , intent(in) :: uid

   write(UID, 100) time_step, fx, fy, fz

   return
100 format ('# ',e12.5,' ', e12.5,' ',e12.5,' ',e12.5)
end subroutine write_dipole_td_header
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%% WRITE_FORCES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Prints the final forces for chosen atoms to output. dxyz is the gradient     !
! vector, natom is the total amount of atoms to print, offset is the starting  !
! atom number, and uid is the output file UID.                                 !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine write_forces(dxyz, natom, offset, uid)
   implicit none
   integer         , intent(in) :: uid, natom, offset
   LIODBLE, intent(in) :: dxyz(3, natom+offset)
   integer :: icount

   do icount = offset+1, offset+natom
      write(uid,100) icount, dxyz(1, icount), dxyz(2, icount), dxyz(3, icount)
   enddo

   return
100 format (I5,2x,f10.6,2x,f10.6,2x,f10.6)
end subroutine write_forces
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%% WRITE_FORCE_LOG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Prints the forces componets to output. ff1G are the 1e gradients, ffSG are   !
! the overlap gradients and ff3G are the 2e gradients (coulomb/XC). natom is   !
! the number of atoms and fileunit is the output file unit. ffT is the total   !
! transposed (because amber).                                                  !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine write_force_log(ffT, ff1G, ffSG, ff3G, ffECPG, natom, fileunit, first_step)
   implicit none
   integer         , intent(in) :: natom, fileunit
   logical         , intent(in) :: first_step
   LIODBLE, intent(in) :: ff1G(natom,3), ffSG(natom,3), ff3G(natom,3),&
                                   ffT(3,natom), ffECPG(natom,3)
   LIODBLE :: ffTall(3)
   character(len=40) :: outfmt
   integer           :: kcount

   outfmt = '(1X, A4, 1X, I4, 3(2X,E14.7))'
   ffTall=0.d0

   if (first_step) then
      open(unit = fileunit, file='Forces.log')
   else
      open(unit = fileunit, file='Forces.log', position='APPEND' )
   endif

   write(fileunit,'(A)') &
      '------------------------------------------------------------'

   do kcount = 1, natom
      write(fileunit, outfmt) 'TOTS', kcount, ffT(1,kcount), ffT(2,kcount), &
                              ffT(3,kcount)
      ffTall(1:3)=ffTall(1:3)+ffT(1:3,kcount)
   enddo
   write(fileunit,'(A)') &
      '------------------------------------------------------------'
      write(fileunit, outfmt) 'SUMT', natom, ffTall(1:3)
   write(fileunit,'(A)') &
      '------------------------------------------------------------'
   write(fileunit,*)

   do kcount = 1, natom
      write(fileunit, outfmt) 'FF1G', kcount, ff1G(kcount,1), ff1G(kcount,2), &
                              ff1G(kcount,3)
      write(fileunit, outfmt) 'FFSG', kcount, ffSG(kcount,1), ffSG(kcount,2), &
                              ffSG(kcount,3)
      write(fileunit, outfmt) 'FF3G', kcount, ff3G(kcount,1), ff3G(kcount,2), &
                              ff3G(kcount,3)
      write(fileunit, outfmt) 'FFEP', kcount, ffECPG(kcount,1), ffECPG(kcount,2), &
                              ffECPG(kcount,3)
      write(fileunit,*)
   enddo
   write(fileunit,*)
   close(fileunit)

end subroutine write_force_log
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%% WRITE_ORBITALS / WRITE_ORBITALS_OPEN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Prints orbital energies and coefficients. M is the total number of basis,    !
! NCO is the number of occupied orbitals, Eorbs is the energy of each MO,      !
! MO_coeff is the MO coefficient matrix, and uid is the output file UID.       !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine write_orbitals(M, NCO, E_orbs, MO_coeff, uid)
   implicit none
   integer, intent(in)          :: M, NCO, uid
   LIODBLE, intent(in) :: E_orbs(M), MO_coeff(M*M)
   integer                      :: icount, jcount

   write(uid,*) 'ORBITAL COEFFICIENTS AND ENERGIES, CLOSED SHELL'
   do icount = 1, NCO
      write(UID, 850) icount, E_orbs(icount)
      do jcount = 1, M
         write(UID, 400) MO_coeff(jcount + icount)
      enddo
   enddo
   do icount = NCO+1, M
      write(uid, 851) icount, E_orbs(icount)
      do jcount = 1, M
         write(uid, 400) MO_coeff(jcount + icount)
      enddo
   enddo

   return
400 format(4(E14.7E2, 2x))
850 format('MOLECULAR ORBITAL #',2x,I3,3x,'ORBITAL ENERGY ',F14.7)
851 format('MOLECULAR ORBITAL #',2x,I3,3x,'ORBITAL ENERGY ',F14.7, '(NON OCC.)')
end subroutine write_orbitals

subroutine write_orbitals_op(M, NCO, NUnp, E_orbs, E_orbs_b, MO_coeff, &
                               MO_coeff_b, uid)
   implicit none
   integer, intent(in)          :: M, NCO, NUnp, uid
   LIODBLE, intent(in) :: E_orbs(M), E_orbs_b(M), MO_coeff(M*M), &
                                   MO_coeff_b(M*M)
   integer                      :: icount, jcount

   write(uid,*) 'ORBITAL COEFFICIENTS AND ENERGIES, OPEN SHELL ALPHA'
   do icount = 1, NCO
      write(UID, 850) icount, E_orbs(icount)
      do jcount = 1, M
         write(UID, 400) MO_coeff(jcount + icount)
      enddo
   enddo
   do icount = NCO+1, M
      write(uid, 851) icount, E_orbs(icount)
      do jcount = 1, M
         write(uid, 400) MO_coeff(jcount + icount)
      enddo
   enddo

   write(uid,*) 'ORBITAL COEFFICIENTS AND ENERGIES, OPEN SHELL BETA'
   do icount = 1, NCO+NUnp
      write(UID, 850) icount, E_orbs_b(icount)
      do jcount = 1, M
         write(UID, 400) MO_coeff_b(jcount + icount)
      enddo
   enddo
   do icount = NCO+NUnp+1, M
      write(uid, 851) icount, E_orbs(icount)
      do jcount = 1, M
         write(uid, 400) MO_coeff_b(jcount + icount)
      enddo
   enddo

   return
400 format(4(E14.7E2, 2x))
850 format('MOLECULAR ORBITAL #',2x,I3,3x,'ORBITAL ENERGY ',F14.7)
851 format('MOLECULAR ORBITAL #',2x,I3,3x,'ORBITAL ENERGY ',F14.7, '(NON OCC.)')
end subroutine write_orbitals_op
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!% WRITE_ORBITAL_POPULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Writes the orbital population to an output file.                             !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine write_orbital_population(ocup_a, ocup_b)
   implicit none
   LIODBLE, intent(in)           :: ocup_a(:)
   LIODBLE, intent(in), optional :: ocup_b(:)

   integer :: ocup_size, icount
   logical :: is_open
   
   ocup_size = size(ocup_a,1)

   inquire(unit = 1444, opened = is_open)
   if (.not. is_open) then
      open(file = 'orb_pops', unit = 1444)
      write(1444,'(A)') "Orbital occupations (density eigenvalues)"
   endif
   
   write(1444,'(A)') ""
   if (present(ocup_b)) then
      do icount = 1, ocup_size
         write(1444,'(I4,F14.7,F14.7,F14.7)')          &
               icount, ocup_a(icount), ocup_b(icount), &
               ocup_a(icount) + ocup_b(icount) 
      enddo
   else
      do icount = 1, ocup_size
         write(1444,'(I4,F14.7)') icount, ocup_a(icount)
      enddo
   endif
end subroutine write_orbital_population
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!