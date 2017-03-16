!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% READECP.F90 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! This file contains subroutines aimed to read Efective Core Potential data.   !
! Subroutines included are:                                                    !
! * lecturaECP                                                                 !
! * dataECPelement                                                             !
! The routines read ECP parameters located in $LIOHOME/dat/ECP/"tipeECP".      !
! V 1.0 september 2015 - Nicolas Foglia.                                       !
!------------------------------------------------------------------------------!
! ECP parameters format:                                                       !
!ATOM_SYMBOL(in capitals)                                                      !
!Junk    LMax_of_ECP   Z_CORE                                                  !
!Junk                                                                          !
!Number_of_functions_for_LMax                                                  !
!n1_ECP b1_ECP a1_ECP                                                          !
!n2_ECP b2_ECP a2_ECP                                                          !
!...                                                                           !
!...                                                                           !
!Number_of_functions_for_L=0                                                   !
!n1_ECP b1_ECP a1_ECP                                                          !
!n2_ECP b2_ECP a2_ECP                                                          !
!...                                                                           !
!...                                                                           !
!Number_of_functions_for_L=1                                                   !
!...                                                                           !
!...                                                                           !
!end                                                                           !
!------------------------------------------------------------------------------!
!             LMax-1   l                                                       !
! V = V_LMax + Σ       Σ |lm> V_l <lm|                                         !
!             l=0     m=-l                                                     !
!                                                                              !
! where Vl = Σ ai_ECP * r^ni_ECP * exp(-bi_ECP*r^2)                            !
!            i                                                                 !
!------------------------------------------------------------------------------!
! Example:                                                                     !
!                                                                              !
!V                                                                             !
!V-ECP     2     10                                                            !
!d-ul potential                                                                !
!  1                                                                           !
!1     15.3000000             -2.5483000                                       !
!s-ul potential                                                                !
!  3                                                                           !
!0      1.8256900              3.5213800                                       !
!2      4.9531400            225.0512300                                       !
!2      2.9512300           -228.5468500                                       !
!p-ul potential                                                                !
!  2                                                                           !
!0     40.5262000              4.5328100                                       !
!2      7.3458990             48.9660120                                       !
!CR     0                                                                      !
!CR-ECP     2     10                                                           !
!d-ul potential                                                                !
!  1                                                                           !
!1     23.8530452             -2.6648250                                       !
!s-ul potential                                                                !
!  3                                                                           !
!0      1.6345600              4.5679900                                       !
!2      6.9382500            271.0586100                                       !
!2      5.4215800           -127.6456800                                       !
!p-ul potential                                                                !
!  2                                                                           !
!0     41.6151600              4.2216200                                       !
!2      8.2232100             55.6485300                                       !
!end                                                                           !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% LecturaECP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Reads atomic symbols and ECP type values.                                    !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine LecturaECP
  use ECP_mod   , only : ecptypes, tipeECP, ZlistECP, asignacion, Zcore, Lmax, &
                         expnumbersECP, nECP, bECP, aECP, verbose_ECP
  use garcha_mod, only : Iz, natom
  
  implicit none
  integer   :: i, jj, ppotat, ppotati, ppotat_acum
  character :: simb*3

  write(*,*)    ; write(*,4160) ; write(*,4161)
  write(*,4162) ; write(*,4163) ; write(*,4164)

  ppotat = 0    ; ppotat_acum = 0    ; Lmax = 0
  Zcore  = 0    ; nECP        = 0 
  bECP   = 0.d0 ; aECP        = 0.d0

  do i = 1, ecptypes
      ppotati = 0 ! Amount of atoms sharing a type of pseudopotential.
      
      ! Takes atomic number (Z) and obtains atom symbol in capital letters.
      call asignacion(ZlistECP(i), simb) 
      ! Reads pseudopotential values for the desired element.
      call dataECPelement(ZlistECP(i), simb)
      
      do jj = 1, natom
         if (Iz(jj).eq.ZlistECP(i)) then
            ppotati=ppotati+1
         endif
      enddo
      ppotat = ppotat + ppotati ! Total amount of atoms with pseudopotentials.      

      if (verbose_ECP .gt. 0) then
         write(*,4165) simb, ppotati,tipeECP
      else
         if (ppotati .gt. 0) write(*,4165) simb, ppotati, tipeECP
      endif
          ppotat_acum = ppotat_acum + ppotati
          ppotati     = 0
  enddo

  write(*,4166)
  write(*,4167) , ppotat_acum
  write(*,4168)
  write(*,*)

4160 FORMAT(4x,"╔═════════════════════════════════════&
      ══════════════╗")
4161 FORMAT(4x,"║    READING EFFECTIVE CORE POTENTIAL PARAMETERS    ║")
4162 FORMAT(4x,"╠══════════════╦═════════════════════╦&
      ══════════════╣")
4163 FORMAT(4x,"║   ATOM TYPE  ║   NUMBER IN SYSTEM  ║   ECP TYPE   ║")
4164 FORMAT(4x,"╠══════════════╬═════════════════════╬&
      ══════════════╣")
4165 FORMAT(4x,"║",a9,5x,"║",i11,10x,"║",5x,a9 ,"║")
4166 FORMAT(4x,"╠══════════════╬═════════════════════╬&
      ══════════════╝")
4167 FORMAT(4x,"║     TOTAL    ║",i11,10x,"║")
4168 FORMAT(4x,"╚══════════════╩═════════════════════╝&
      ")
end subroutine lecturaECP
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% LecturaECP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Reads ECP values for a given element.                                        !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine dataECPelement(Z, elemento)

   use ECP_mod, only : Zcore, Lmax, expnumbersECP , nECP, bECP, aECP, tipeECP
   ! Zcore(Z) is the core's charge for the chosen ECP, given the atom's Z. 
   ! Lmax(Z) is the chosen ECP's maximum L, given the atom's Z.
   ! expnumbersECP(Z, l) is the amount of ECP terms for the given l angular 
   ! momentum and atom's Z.
   ! nECP, bECP y aECP contain the n, b and A coefficients for each kind of ECP
   ! contraction. 
   !  Vl = Σ ai_ECP * r^ni_ECP * exp(-bi_ECP*r^2)
   !       i 
   ! tipeECP contains the type of pseudopotential used.

   implicit none
   integer  , intent(in) :: Z          ! Nuclear atomic charge.
   character, intent(in) :: elemento*3 ! The name of the atomic element to search.
   character :: simbolo*3              ! Symbol used to search in the ECP file.

   ! Auxiliar variables
   character :: boba*6, lio*255
   integer   :: u, l
   logical   :: existeECP, found

   if (tipeECP.eq."USSERDEF") then
      existeECP = .true.
      inquire(file = "ECPPARAM", exist=existeECP)
   else
      call get_environment_variable("LIOHOME", lio)!
      existeECP = .true.
      inquire(file = trim(lio)//"/dat/ECP/"//tipeECP, exist = existeECP)
   endif

   ! Verifies if ECP parameters file exists.
   if (.not. existeECP) then
      write(*,*) "Effective Core Potential parameters file not found"

      if (tipeECP.eq."USSERDEF") then
         write(*,*) "Effective Core potential parameter file have", &
         " to be name ECPPARAM"
      else       
         write(*,*) "check ", trim(lio), "/dat/ECP/"
      endif

      stop
   else
      if (tipeECP.eq."USSERDEF") then
         open(unit = 7, file = "ECPPARAM")
      else
         open(unit = 7, file = trim(lio)//"/dat/ECP/"//tipeECP)
      end if

      found = .true.
      do while(found)
         read(7,*) simbolo ! Reads file line by line until element is found.
         if (simbolo.eq.elemento) then
            found = .false.       
            read(7,*) boba, Lmax(Z), Zcore(Z) 
            read(7,*)
            read(7,*) expnumbersECP(Z, Lmax(Z))
 
            do u = 1, expnumbersECP(Z, Lmax(Z))
               read(7,*) nECP(Z,Lmax(Z),u), bECP(Z,Lmax(Z),u), &
                         aECP(Z,Lmax(Z),u)
            enddo

            do l = 0, Lmax(Z)-1 ! Repeats from l=0 to Lmax-1
               read(7,*)
               read(7,*) expnumbersECP(Z,l)
               do u=1, expnumbersECP(Z,l)
                  read(7,*) nECP(Z,l,u), bECP(Z,l,u), aECP(Z,l,u)
               end do
            end do
         elseif (simbolo .eq. "end") then
            ! Ends program if pseudopotential data is not found.
            write (*,*) "Element ",elemento ," not found in ECP file."
            if (tipeECP.eq."USSERDEF") then
               write(*,*) "Please check the ECP params file."
            else
               write(*,*) "Please check ", trim(lio), "/dat/ECP/", tipeECP
            endif
            stop
         endif
      enddo
      close(7)
   endif
end subroutine dataECPelement
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!


