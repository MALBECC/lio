c subroutine that asignates masses and specie
      subroutine assign(na_u,nac,atname,iza,izs,masst)

      use precision, only : dp
      use fdf
      use sys
      implicit          none
      integer           i,j,k, na_u, nac, iu
      integer           iza(na_u), izs(na_u+nac)
      real(dp)          masst(na_u+nac), ATMASS, mass
      character         atname(nac)*4,atn4*4,atn1*1,atn2*2
      character*2       sym(na_u+nac), SYMBOL, name

C Assigantes solute atomic masses (masst) and atomic symbol (sym)
      masst=0.0
      do i=1,na_u
        masst(i) = ATMASS(iza(i))
        sym(i) = SYMBOL(iza(i))
        if(masst(i).eq.0.0) then
          call die("assign: There are solute atoms without mass")
        endif
      enddo

C Assigantes solvent atomic masses (masst) and atomic symbol (sym)
        izs = 0
        izs(1:na_u)=iza(1:na_u)
        k=na_u+1
        do i=1,nac
          atn4 = atname(i)
          atn1 = atn4(1:1)
          atn2 = atn4(1:2)
          if (atn1.eq.'H')  izs(k)=1
          if (atn1.eq.'C')  izs(k)=6
          if (atn1.eq.'N')  izs(k)=7
          if (atn1.eq.'O')  izs(k)=8
          if (atn1.eq.'F')  izs(k)=9
          if (atn2.eq.'Na') izs(k)=11
          if (atn2.eq.'Mg') izs(k)=12
          if (atn1.eq.'P')  izs(k)=15
          if (atn1.eq.'S')  izs(k)=16
          if (atn2.eq.'Cl') izs(k)=17
          if (atn1.eq.'K')  izs(k)=19
          if (atn2.eq.'Ca') izs(k)=20
          if (atn2.eq.'Mn') izs(k)=25
          if (atn2.eq.'FE') izs(k)=26
          if (atn2.eq.'Cu') izs(k)=29
          if (atn2.eq.'Zn') izs(k)=30
          if (atn2.eq.'Br') izs(k)=35
          if (atn1.eq.'I')  izs(k)=53
          if (atn2.eq.'Pt') izs(k)=78
 
      if(izs(k).eq.0) then
       call die('assign: There are solvent atoms without atomic number')
      endif
       masst(k) = ATMASS (izs(k))
       sym(k)= SYMBOL(izs(k))
      if(masst(k).eq.0.0) then
       call die('assign: There are solvent atoms without mass')
      endif
       k=k+1
        enddo

c Read atomic masses of different atoms from fdf block 

      if ( fdf_block('NewMasses', iu) ) then
 5     continue
       read(iu,'(A2)',advance='no',err=10,end=10) name 
       if(name.eq.'%e'.or.name.eq.'%E') return 
       read(iu,*,err=10,end=10) mass 
       write(6,"(/,a, A2, a, f6.2)")
     . 'assign: Read atomic mass of:  ',name,'as',mass

c assignates new masses
        do i=1,na_u+nac
        if(name.eq.sym(i)) masst(i) = mass
        enddo

       goto 5
 6     continue
      endif

      return
 10   stop 'assign: Problem reading from "NewMasses" block'
      end subroutine assign
c*************************************************************************
      FUNCTION SYMBOL( Z )
C Given the atomic number, returns the atomic symbol (e.g. 'Na')
C Written by J. Soler
      use sys
      implicit none
      character(len=2)    :: SYMBOL  ! Atomic symbol
      integer, intent(in) :: Z       ! Atomic number

      integer, parameter  :: NZ=103
      character(len=2) :: NAME(NZ)
      DATA NAME /'H' ,'He','Li','Be','B' ,'C' ,'N' ,'O' ,'F' ,'Ne',
     .           'Na','Mg','Al','Si','P' ,'S' ,'Cl','Ar','K' ,'Ca',
     .           'Sc','Ti','V' ,'Cr','Mn','Fe','Co','Ni','Cu','Zn',
     .           'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y' ,'Zr',
     .           'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',
     .           'Sb','Te','I' ,'Xe','Cs','Ba','La','Ce','Pr','Nd',
     .           'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
     .           'Lu','Hf','Ta','W' ,'Re','Os','Ir','Pt','Au','Hg',
     .           'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th',
     .           'Pa','U' ,'Np','Pu','Am','Cm','Bk','Cf','Es','Fm',
     .           'Md','No','Lr'/

      IF (Z.EQ.0 .OR. Z.EQ.-100) THEN
         SYMBOL = 'BS'
      ELSE IF (ABS(Z).LE.NZ) THEN
         SYMBOL = NAME(ABS(Z))
      ELSE
         WRITE(6,*) 'SYMBOL: ERROR: No data for Z =', Z
         SYMBOL = ' '
      ENDIF

      END function symbol

      FUNCTION ATMASS( Z )
C Returns the average atomic mass from the atomic number Z.
C Dta taken from VCH periodic table.
C Written by J.M.Soler. April'97.
      use precision, only : dp
      use sys
      implicit none

      real(dp)             :: ATMASS ! Average atomic mass, in amu
      integer, intent(in)  :: Z      ! Atomic number

      integer, PARAMETER  :: NZ=94
      character(len=50) message

      real AMASS(0:NZ)
      DATA AMASS / 0.00,
     .     1.01,  4.00,  6.94,  9.01, 10.81, 12.01, 14.01, 16.00,
     .    19.00, 20.18, 22.99, 24.31, 26.98, 28.09, 30.97, 32.07,
     .    35.45, 39.95, 39.10, 40.08, 44.96, 47.88, 50.94, 52.00,
     .    54.94, 55.85, 58.93, 58.69, 63.55, 65.39, 69.72, 72.61,
     .    74.92, 78.96, 79.90, 83.80, 85.47, 87.62, 88.91, 91.22,
     .    92.91, 95.94, 98.91,101.07,102.91,106.42,107.87,112.41,
     .   114.82,118.71,121.75,127.60,126.90,131.29,132.91,137.33,
     .   138.91,140.12,140.91,144.24,146.92,150.36,151.97,157.25,
     .   158.93,162.50,164.93,167.26,168.93,173.04,174.97,178.49,
     .   180.95,183.85,186.21,190.2 ,192.22,195.08,196.97,200.59,
     .   204.38,207.2 ,208.98,208.98,209.99,222.02,223.02,226.03,
     .   227.03,232.04,231.04,238.03,237.05,244.06/

      IF (Z.LT.0 .OR. Z.GT.NZ) THEN
         write(message,'(a,i4)') 'ATMASS: ERROR: No data for Z =',Z
         call die(message)
      ELSE
         ATMASS=real(AMASS(Z),kind=dp)
      ENDIF

      END function atmass

c*************************************************************************



