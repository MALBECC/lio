!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% EXTRAS.F90 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! This file mostly contains subroutines regarding output writing and formats.  !
! Subroutines included are:                                                    !
! * lio_logo        (prints "WELCOME TO LIO")                                  !
! * lio_logo2       (prints logo in fancier characters)                        !
! * new_write_nml   (prints namelist options)                                  !
! * write_zlist_ecp (prints the Z used for ECP)                                !
! * writeEnergies   (prints energy summary)                                    !
! * write_e_step    (prints the energy for each SCF iteration)                 !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

subroutine  lio_logo()
   implicit none
   write(*,*)
   write(*,9001 ); write(*,9002 ); write(*,9003 ); write(*,9004 )
   write(*,9005 ); write(*,9006 ); write(*,9007 ); write(*,9008 )
   write(*,9009 ); write(*,9010 ); write(*,9011 ); write(*,9012 )
   write(*,9013 ); write(*,9014 ); write(*,9015 ); write(*,9016 )
   write(*,9017 ); write(*,9018 ); write(*,9019 ); write(*,9020 )
   write(*,9021 ); write(*,9022 ); write(*,9023 ); write(*,9024 )
   write(*,9025 ); write(*,9026 )
   write(*,*)

9001 format(4x,"╔═════════════════════════════════════&
      ════════════════════════════════════╗")
9002 format(4x,"║                                                                         ║")
9003 format(4x,"║  E     W     E ELCOMEW  C        MEW       ELCO    L       C   WELCOME  ║")
9004 format(4x,"║  O    LCO    L O        E       W   O     C    E   ME     ME   C        ║")
9005 format(4x,"║   E   M W   O  W        L      L     E   M         E      EL   E        ║")
9006 format(4x,"║   O   E C   W  C        M      M         E      E  O E   C M   L        ║")
9007 format(4x,"║   W   O E   C  EWELCOM  E      E         O      O  W L   E E   MEWELCO  ║")
9008 format(4x,"║    O E   C M   L        O      O         W      W  C  E E  O   E        ║")
9009 format(4x,"║    W L   E E   M        W      W         C      C  E  L O  W   O        ║")
9010 format(4x,"║    COM   LCO   E        C      C     L   E      E  L  MEW  C   W        ║")
9011 format(4x,"║     W     E    O        E       W   O     C    E   M   L   E   C        ║")
9012 format(4x,"║     C     L    WELCOME  LCOMEW   OME       WELC    E   M   L   EWELCOM  ║")
9013 format(4x,"║                                                                         ║")
9014 format(4x,"║                                                                         ║")
9015 format(4x,"║             TOTOTOT    OTOT             O      L     LIOL               ║")
9016 format(4x,"║                O      T    O            O      L    O    I              ║")
9017 format(4x,"║                O     O      T           O      L   I      O             ║")
9018 format(4x,"║                O     O      T           O      L   I      O             ║")
9019 format(4x,"║                O     O      T           O      L   I      O             ║")
9020 format(4x,"║                O     O      T           O      L   I      O             ║")
9021 format(4x,"║                O     O      T           O      L   I      O             ║")
9022 format(4x,"║                O     O      T           O      L   I      O             ║")
9023 format(4x,"║                O      T    O            O      L    O    I              ║")
9024 format(4x,"║                O       OTOT             OLIOLI L     LIOL               ║")
9025 format(4x,"║                                                                         ║")
9026 format(4x,"╚═════════════════════════════════════&
      ════════════════════════════════════╝")

end subroutine lio_logo

subroutine lio_logo2()
   write(*,*)
   write(*,1700); write(*,1722); write(*,1701); write(*,1702)
   write(*,1703); write(*,1704); write(*,1705);	write(*,1706)
   write(*,1722); write(*,1708); write(*,1709); write(*,1710)
   write(*,1711); write(*,1712); write(*,1713); write(*,1714)
   write(*,1715); write(*,1716); write(*,1717); write(*,1718)
   write(*,1719); write(*,1720); write(*,1722); write(*,1721)
   write(*,*)

1700 format(4x,"╔════════════════════════════════════"&
     ,"════════════════════════════════════════"&
     ,"══════════════════╗")
1701 format(4x,"║",4x,"██╗    ██╗    ███████╗    ██╗          ██████╗"&
     ,"     ██████╗     ███╗   ███╗    ███████╗",4x,"║")
1702 format(4x,"║",4x,"██║    ██║    ██╔════╝    ██║         ██╔════╝"&
     ,"    ██╔═══██╗    ████╗ ████║    ██╔════╝",4x,"║")
1703 format(4x,"║",4x,"██║ █╗ ██║    █████╗      ██║         ██║     "&
     ,"    ██║   ██║    ██╔████╔██║    █████╗  ",4x,"║")
1704 format(4x,"║",4x,"██║███╗██║    ██╔══╝      ██║         ██║     "&
     ,"    ██║   ██║    ██║╚██╔╝██║    ██╔══╝  ",4x,"║")
1705 format(4x,"║",4x,"╚███╔███╔╝    ███████╗    ███████╗    ╚███"&
     ,"███╗    ╚██████╔╝    ██║ ╚═╝ ██║    ███████╗",4x,"║")
1706 format(4x,"║",4x," ╚══╝╚══╝     ╚══════╝    ╚══════╝     ╚════"&
     ,"═╝     ╚═════╝     ╚═╝     ╚═╝    ╚══════╝",4x,"║")
1708 format(4x,"║",4x,"                                ████████╗     █"&
     ,"█████╗                                 ",4x,"║")
1709 format(4x,"║",4x,"                                ╚══██╔══╝    █"&
     ,"█╔═══██╗                                ",4x,"║")
1710 format(4x,"║",4x,"                                   ██║       █"&
     ,"█║   ██║                                ",4x,"║")
1711 format(4x,"║",4x,"                                   ██║       █"&
     ,"█║   ██║                                ",4x,"║")
1712 format(4x,"║",4x,"                                   ██║       ╚"&
     ,"██████╔╝                                ",4x,"║")
1713 format(4x,"║",4x,"                                   ╚═╝        ╚"&
     ,"═════╝                                 ",4x,"║")
1714 format(4x,"║",4x,"                                               "&
     ,"                                       ",4x,"║")
1715 format(4x,"║",4x,"                             ██╗         ██╗   "&
     ,"  ██████╗                              ",4x,"║")
1716 format(4x,"║",4x,"                             ██║         ██║   "&
     ," ██╔═══██╗                             ",4x,"║")
1717 format(4x,"║",4x,"                             ██║         ██║   "&
     ," ██║   ██║                             ",4x,"║")
1718 format(4x,"║",4x,"                             ██║         ██║   "&
     ," ██║   ██║                             ",4x,"║")
1719 format(4x,"║",4x,"                             ███████╗    ██║   "&
     ," ╚██████╔╝                             ",4x,"║")
1720 format(4x,"║",4x,"                             ╚══════╝    ╚═╝   "&
     ,"  ╚═════╝                              ",4x,"║")
1721 format(4x,"╚════════════════════════════════════"&
     ,"════════════════════════════════════════"&
     ,"══════════════════╝")
1722 format(4x,"║",94x,"║")

end subroutine lio_logo2


subroutine new_write_nml(charge)

   use garcha_mod, only :
   use ECP_mod   , only : ecpmode, ecptypes, tipeECP, ZlistECP, cutECP,      &
                          local_nonlocal, ecp_debug, ecp_full_range_int,     &
                          verbose_ECP, Cnorm, FOCK_ECP_read, FOCK_ECP_write, &
                          Fulltimer_ECP, cut2_0, cut3_0
#ifdef CUBLAS
   use cublasmath
#endif

   integer, intent(in) :: charge

   ! LIO
   write(*,8000)       ; write(*,8100)      ; write(*,8001)

   ! System  
   write(*,8000)       ; write(*,8101)      ; write(*,8002)
   write(*,8210) natom ; write(*,8211) nsol ; write(*,8212) charge
   write(*,8300) Nunp  ; write(*,8003)

   ! Theory level
   write(*,8000)          ; write(*,8102)       ; write(*,8002)
   write(*,8220) open     ; write(*,8221) nmax  ; write(*,8222) basis_set
   write(*,8223) fitting_set                    ; write(*,8224) int_basis
   write(*,8225) diis     ; write(*,8226) ndiis ; write(*,8227) gold
   write(*,8228) told     ; write(*,8229) Etold ; write(*,8230) hybrid_converg
   write(*,8231) good_cut ; write(*,8301) Rmax  ; write(*,8302) RmaxS
   write(*,8003)

   ! Effective Core Potential
   write(*,8000)         ; write(*,8103)          ; write(*,8002)
   write(*,8240) Ecpmode ; write(*,8241) Ecptypes ; write(*,8242) TipeECP

   call write_zlist_ecp(ZlistECP, Ecptypes)
   write(*,8244) Fock_ECP_read  ; write(*,8245) Fock_ECP_write
   write(*,8246) cutECP         ; write(*,8247) cut2_0    
   write(*,8248) cut3_0         ; write(*,8249) Verbose_ECP 
   write(*,8250) ECP_debug      ; write(*,8251) fulltimer_ECP
   write(*,8252) local_nonlocal ; write(*,8253) ECP_full_range_int
   write(*,8003)                
 
   ! TDDFT
   write(*,8000)           ; write(*,8104)            ; write(*,8002)
   write(*,8260) Timedep   ; write(*,8261) Propagator ; write(*,8262) Ntdstep
   write(*,8263) Tdstep    ; write(*,8264) Nbsh       ; write(*,8265) Field
   write(*,8266) Fx        ; write(*,8267) Fy         ; write(*,8268) Fz
   write(*,8269) Tdrestart ; write(*,8314) Exter      ; write(*,8003)

   ! Output Options
   write(*,8000)             ; write(*,8105)          ; write(*,8002)
   write(*,8280) Verbose     ; write(*,8281) writeXYZ ; write(*,8282) writeDens
   write(*,8283) writeForces ; write(*,8003)


   ! Restart Options
   write(*,8000)       ; write(*,8106)            ; write(*,8002)
   write(*,8290) VCinp ; write(*,8291) Frestartin ; write(*,8292) Rho_restart_in
   write(*,8293) Rho_restart_out                  ; write(*,8003)

   ! Miscellaneous
   write(*,8000)           ; write(*,8107)         ; write(*,8002)
   write(*,8303) PredCoef  ; write(*,8304) Idip    ; write(*,8305) IntSolDouble
   write(*,8306) DgTRIG    ; write(*,8307) Iexch   ; write(*,8308) Integ
   write(*,8309) Dens      ; write(*,8310) Igrid   ; write(*,8311) Igrid2
   write(*,8312) A0        ; write(*,8313) Epsilon ; write(*,8315) Int_basis
   write(*,8316) Cubegen_only                      ; write(*,8317) Cube_Res
   write(*,8318) Cube_Dens ; write(*,8319) Cube_Dens_file
   write(*,8320) Cube_Orb  ; write(*,8321) Cube_Sel; write(*,8322) Cube_Orb_File
   write(*,8323) Cube_Elec ; write(*,8324) Cube_Elec_File
   write(*,8003)
 
   return
8000 format(4x,"╔══════════════════════════════════&
 ════════════════╗")
8001 format(4x,"╚══════════════════════════════════&
 ════════════════╝")
8002 format(4x,"╠══════════════════════╦═══════════&
 ════════════════╣")
8003 format(4x,"╚══════════════════════╩═══════════&
 ════════════════╝")

8100 format(4x,"║                     LIO Input                    ║")
8101 format(4x,"║                      System                      ║")
8102 format(4x,"║                   Theory Level                   ║")
8103 format(4x,"║             Effective Core Potential             ║")
8104 format(4x,"║                       TDDFT                      ║")
8105 format(4x,"║                   write Option                   ║")
8106 format(4x,"║                      Restart                     ║")
8107 format(4x,"║                      Others                      ║")

!System
8210 format(4x,"║  Natom               ║  ",7x,i6,12x,"║")
8211 format(4x,"║  Nsol                ║  ",7x,i8,10x,"║")
8212 format(4x,"║  Charge              ║  ",9x,i5,11x,"║")
8300 format(4x,"║  Nunp                ║  ",i5,20x,"║")

!Theory level
8220 format(4x,"║  Open                ║  ",l,23x,"║")
8221 format(4x,"║  Nmax                ║  ",i5,20x,"║")
8222 format(4x,"║  Basis_Set           ║  ",a24," ║")
8223 format(4x,"║  Fitting_Set         ║  ",a24," ║")
8224 format(4x,"║  Int_Basis           ║  ",l,23x,"║")
8225 format(4x,"║  Diis                ║  ",l,23x,"║")
8226 format(4x,"║  Ndiis               ║  ",i3,22x,"║")
8227 format(4x,"║  Gold                ║  ",f14.8,11x,"║")
8228 format(4x,"║  Told                ║  ",f14.8,11x,"║")
8229 format(4x,"║  Etold               ║  ",f14.8,11x,"║")
8230 format(4x,"║  Hybrid_converg      ║  ",l,23x,"║")
8231 format(4x,"║  Good_cut            ║  ",f14.8,11x,"║")
8301 format(4x,"║  Rmax                ║  ",f14.8,11x,"║")
8302 format(4x,"║  RmaxS               ║  ",f14.8,11x,"║")

!Effective Core Potential
8240 format(4x,"║  Ecpmode             ║  ",l,23x,"║")
8241 format(4x,"║  Ecptypes            ║  ", i3,22x,"║")
8242 format(4x,"║  TipeECP             ║  ",a25,"║")
8243 format(4x,"║  Zlistecp            ║  ",i2,i2,i2,i2,i2,i2,"║")
8244 format(4x,"║  Fock_ECP_read       ║  ",l,23x,"║")
8245 format(4x,"║  Fock_ECP_write      ║  ",l,23x,"║")
8246 format(4x,"║  cutECP              ║  ",l,23x,"║")
8247 format(4x,"║  cut2_0              ║  ",f14.8,11x,"║")
8248 format(4x,"║  cut3_0              ║  ",f14.8,11x,"║")
8249 format(4x,"║  Verbose_ECP         ║  ",i2,23x,"║")
8250 format(4x,"║  ECP_debug           ║  ",l,23x,"║")
8251 format(4x,"║  fulltimer_ECP       ║  ",l,23x,"║")
8252 format(4x,"║  local_nonlocal      ║  ",i2,23x,"║")
8253 format(4x,"║  ECP_full_range_int  ║  ",l,23x,"║")

!TDDFT
8260 format(4x,"║  Timedep             ║  ",i2,23x,"║")
8261 format(4x,"║  Propagator          ║  ",i2,23x,"║")
8262 format(4x,"║  Ntdstep             ║  ",i10,15x,"║")
8263 format(4x,"║  Tdstep              ║  ",f14.8,11x,"║")
8264 format(4x,"║  Nbsh                ║  ",i4,21x,"║")
8265 format(4x,"║  Field               ║  ",l,23x,"║")
8266 format(4x,"║  Fx                  ║  ",f14.8,11x,"║")
8267 format(4x,"║  Fy                  ║  ",f14.8,11x,"║")
8268 format(4x,"║  Fz                  ║  ",f14.8,11x,"║")
8269 format(4x,"║  Tdrestart           ║  ",l,23x,"║")
8314 format(4x,"║  Exter               ║  ",l,23x,"║")

!write Options
8280 format(4x,"║  Verbose             ║  ",l,23x,"║")
8281 format(4x,"║  writeXYZ            ║  ",l,23x,"║")
8282 format(4x,"║  writeDens           ║  ",l,23x,"║")
8283 format(4x,"║  writeForces         ║  ",l,23x,"║")
!Restart
8290 format(4x,"║  VCinp               ║  ",l,23x,"║")
8291 format(4x,"║  Frestartin          ║  ",a25,"║")
8292 format(4x,"║  Rho_restart_in      ║  ",l,23x,"║")
8293 format(4x,"║  Rho_restart_out     ║  ",l,23x,"║")

!Others
8303 format(4x,"║  PredCoef            ║  ",l,23x,"║")
8304 format(4x,"║  Idip                ║  ",i4,21x,"║")
8305 format(4x,"║  IntSolDouble        ║  ",l,23x,"║")
8306 format(4x,"║  DgTRIG              ║  ",f14.8,11x,"║")
8307 format(4x,"║  Iexch               ║  ",i5,20x,"║")
8308 format(4x,"║  Integ               ║  ",l,23x,"║")
8309 format(4x,"║  Dens                ║  ",l,23x,"║")
8310 format(4x,"║  Igrid               ║  ",i3,22x,"║")
8311 format(4x,"║  Igrid2              ║  ",i3,22x,"║")
8312 format(4x,"║  A0                  ║  ",f14.7,11x,"║")
8313 format(4x,"║  Epsilon             ║  ",f14.7,11x,"║")
8315 format(4x,"║  Int_basis           ║  ",l,23x,"║")
8316 format(4x,"║  Cubegen_only        ║  ",l,23x,"║")
8317 format(4x,"║  Cube_Res            ║  ",i5,20x,"║")
8318 format(4x,"║  Cube_Dens           ║  ",l,23x,"║")
8319 format(4x,"║  Cube_Dens_file      ║  ",a25,"║")
8320 format(4x,"║  Cube_Orb            ║  ",l,23x,"║")
8321 format(4x,"║  Cube_Sel            ║  ",i5,20x,"║")
8322 format(4x,"║  Cube_Orb_File       ║  ",a25,"║")
8323 format(4x,"║  Cube_Elec           ║  ",l,23x,"║")
8324 format(4x,"║  Cube_Elec_File      ║  ",a25,"║")

end subroutine new_write_nml

subroutine write_zlist_ecp(ZlistECP, D)
   integer, intent(in) :: ZlistECP(128), D
   integer :: i, k, lines, rest

   if (D .lt. 6) then
      if (D .eq. 1) write(*,8538) ZlistECP(1)
      if (D .eq. 2) write(*,8539) ZlistECP(1:2)
      if (D .eq. 3) write(*,8540) ZlistECP(1:3)
      if (D .eq. 4) write(*,8541) ZlistECP(1:4)
      if (D .eq. 5) write(*,8542) ZlistECP(1:5)
   else
      lines = D/6
      rest  = mod(D,6)
      write(*,8543) ZlistECP(1:6)

      do i=1,lines-1
         k = 6*i + 1
         write(*,8544) ZlistECP(k:k+5)
      enddo

      if (rest .eq. 1) write(*,8545) ZlistECP(6*lines+1:D)
      if (rest .eq. 2) write(*,8546) ZlistECP(6*lines+1:D)
      if (rest .eq. 3) write(*,8547) ZlistECP(6*lines+1:D)
      if (rest .eq. 4) write(*,8548) ZlistECP(6*lines+1:D)
      if (rest .eq. 5) write(*,8549) ZlistECP(6*lines+1:D)
   endif

8538 format(4x,"║  Zlistecp            ║ ",i3,"                       ║")
8539 format(4x,"║  Zlistecp            ║ ",i3,i3,"                    ║")
8540 format(4x,"║  Zlistecp            ║ ",i3,i3,i3,"                 ║")
8541 format(4x,"║  Zlistecp            ║ ",i3,i3,i3,i3,"              ║")
8542 format(4x,"║  Zlistecp            ║ ",i3,i3,i3,i3,i3,"           ║")
8543 format(4x,"║  Zlistecp            ║ ",i3,i3,i3,i3,i3,i3,"        ║")
8544 format(4x,"║                      ║ ",i3,i3,i3,i3,i3,i3,"        ║")
8545 format(4x,"║                      ║ ",i3"                        ║")
8546 format(4x,"║                      ║ ",i3,i3,"                    ║")
8547 format(4x,"║                      ║ ",i3,i3,i3,"                 ║")
8548 format(4x,"║                      ║ ",i3,i3,i3,i3"               ║")
8549 format(4x,"║                      ║ ",i3,i3,i3,i3,i3"            ║")

end subroutine write_zlist_ecp


subroutine writeEnergies(E1, E2, En, Eecp, Exc, Es, ecpmode, E_restrain)
   use garcha_mod, only : number_restr, nsol
   implicit none
   real*8 , intent(in) :: E1, E2, En, Eecp, Exc, Es, E_restrain
   logical, intent(in) :: ecpmode

   write(*,*)       ; write(*,7000)         ; write(*,7004)
   write(*,7001)    ; write(*,7005) E1-Eecp ; write(*,7002)
   write(*,7006) E2 ; write(*,7002)         ; write(*,7007) En
   
   ! Adds ECP energy if ECP is active.
   if (ecpmode) then   
      write(*,7002) ; write(*,7008) Eecp
   endif
 
   ! Prints QM/MM energy if solvent is present.  
   if (nsol.gt.0) then
      write(*,7002) ; write(*,7012) Es
   endif
 
   write(*,7002) ; write(*,7009) Exc ;  write(*,7002)

   if (number_restr.gt.0) then
      write(*,7011) E_restrain ; write(*,7002)
   end if

   write(*,7010) E1 + E2 + En + Exc ; write(*,7003) ; write(*,*)

   return
7000 format(4x,"╔══════════════════════════════════&
 ══════════╗")
7001 format(4x,"╠══════════════════╦═══════════════&
 ══════════╣")
7002 format(4x,"╠══════════════════╬═══════════════&
 ══════════╣")
7003 format(4x,"╚══════════════════╩═══════════════&
 ══════════╝")

7004 format(4x,"║        ENERGY CONTRIBUTIONS in A.U.        ║")
7005 format(4x,"║   ONE ELECTRON   ║",4x,F14.7,7x"║")
7006 format(4x,"║   COULOMB        ║",4x,F14.7,7x,"║")
7007 format(4x,"║   NUCLEAR        ║",4x,F14.7,7x,"║")
7008 format(4x,"║   E. CORE POT    ║",4x,F14.7,7x"║")
7009 format(4x,"║   EXC. - CORR.   ║",4x,F14.7,7x"║")
7012 format(4x,"║   E QM-MM        ║",4x,F14.7,7x"║")
7010 format(4x,"║   TOTAL          ║",4x,F14.7,7x"║")
7011 format(4x,"║   E. RESTR.      ║",4x,F14.7,7x"║")

end subroutine writeenergies


subroutine write_e_step(step, energy)
   implicit none
   integer, intent(in) :: step
   real*8 , intent(in) :: energy

   write(6,8500)
   write(6,8501) step, energy
   write(6,8502)

   return

8500 format(4x,"╔════════╦═════════════╦═══════════&
 ╦══════════════════════╗")
8501 format(4x,"║ Iter # ║",2x,I10,1x,"║ QM Energy ║",4x,F14.7,4x,"║")
8502 format(4x,"╚════════╩═════════════╩═══════════&
 ╩══════════════════════╝")
end subroutine write_e_step

subroutine write_conv_status(GOOD, TOLD, EGOOD, ETOLD)
   implicit none
   real*8, intent(in) :: GOOD, TOLD, EGOOD, ETOLD

   write(6,8601)
   write(6,8602)
   write(6,8603)
   write(6,8604) GOOD , TOLD
   write(6,8605) EGOOD, ETOLD
   write(6,8606)
 
   return

8601 format(4x,"           ╔════════════╦═════════════╗")
8602 format(4x,"           ║    Value   ║ Conv. Crit. ║")
8603 format(4x,"╔══════════╬════════════╬══════════&
 ═══╣")
8604 format(4x,"║ Good     ║",1x,E10.3,1x,"║",1x,E10.3,2x,"║")
8605 format(4x,"║ En. Good ║",1x,E10.3,1x,"║",1x,E10.3,2x,"║")
8606 format(4x,"╚══════════╩════════════╩══════════&
 ═══╝")
end subroutine write_conv_status
