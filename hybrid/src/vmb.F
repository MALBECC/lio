      subroutine temp(iunit,natoms,ma,va,ntcon,tempe)
C *************************************************************************
C Subroutine to calculate instantaneous temperature
C
C Written by P.Ordejon, November'97
C ************************* UNITS ******************************************
C Temperature in Kelvin
C Atomic masses in atomic mass units
C
C Space units depend on input option:
C
C   Option iunit = 1:
C     Distances are in Angstrom
C   Option iunit = 2:
C     Distances are in Bohr
C ************************* INPUT *********************************************
C integer iunit         : Units option: 1 or 2 (see UNITS above)
C integer natoms        : Number of atoms in the simulation cell
C real*8 ma(natoms)     : Atomic masses 
C real*8 va(3,natoms)   : Atomic velocities 
C integer               : Total number of position constraints imposed
C ************************* OUTOPUT *******************************************
C real*8 tempe          : Instantaneous system temperature 
C *****************************************************************************
      implicit none

      integer 
     .   natoms,iunit,ntcon

      double precision
     .  ma(natoms),tempe,va(3,natoms)

C Internal variables ..........................................................
 
      integer
     .  i,ia

      double precision
     .  Ang,eV,fovermp,kin

C ........................

      if (iunit .ne. 1 .and. iunit .ne. 2) then
        write(6,*) 'temp: Wrong iunit option;  must be 1 or 2'
        stop
      endif

C Define constants and conversion factors .....................................

      Ang = 1.d0 / 0.529177d0
c      eV  = 1.d0 / 13.60580d0
      eV     = 1._dp / 27.211396132_dp
      if (iunit .eq. 1) then
C  convert F/m in (eV/Amstrong)/amu  to  Amstrong/fs**2
        fovermp = 0.009579038d0
      else
C  convert F/m in (Ry/Bohr)/amu  to  Bohr/fs**2
        fovermp = 0.009579038d0 *Ang**2 / eV
      endif

C ........................

C Calculate kinetic energy and temperature ...................
C Kinetic energy of atoms
      kin = 0.0d0
      do ia = 1,natoms
        do i = 1,3
          kin = kin + 0.5d0 * ma(ia) * va(i,ia)**2 / fovermp
        enddo
      enddo

C Instantaneous temperature (Kelvin)
      if (iunit .eq. 1) then
        tempe = 2.0d0*kin/(3.0d0*natoms-3.0d0-ntcon)/8.617d-5
      else
        tempe = 2.0d0*kin/(3.0d0*natoms-3.0d0-ntcon)/8.617d-5/eV
      endif
C .....................

      return
      end
    
