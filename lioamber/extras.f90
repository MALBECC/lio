
	SUBROUTINE NEW_WRITE_NML()

      use garcha_mod
      use ECP_mod, only : ecpmode, ecptypes, tipeECP, ZlistECP &
      ,cutECP,local_nonlocal, ecp_debug,ecp_full_range_int &
      ,verbose_ECP,Cnorm,FOCK_ECP_read, FOCK_ECP_write,Fulltimer_ECP &
      ,cut2_0,cut3_0
      use fileio_data, only: verbose

      return;
      if (verbose .lt. 4) return
!LIO
	write(*,8000)
	write(*,8100)
        write(*,8001)

!System
        write(*,8000)
        write(*,8101)
        write(*,8002)
	write(*,8210) natom
        write(*,8211) nsol
        write(*,8212) charge
        write(*,8300) Nunp
	write(*,8003)

!Theory level
	write(*,8000)
        write(*,8102)
        write(*,8002)
	write(*,8220) open
        write(*,8221) nmax
        write(*,8222) basis_set
        write(*,8223) fitting_set
        write(*,8224) int_basis
        write(*,8225) diis
        write(*,8226) ndiis
        write(*,8227) gold
        write(*,8228) told
        write(*,8229) Etold
        write(*,8230) hybrid_converg
        write(*,8231) good_cut
        write(*,8301) Rmax
        write(*,8302) RmaxS
        write(*,8003)

!Effective Core Potential
        write(*,8000)
        write(*,8103)
        write(*,8002)
	write(*,8240) Ecpmode
        write(*,8241) Ecptypes
        write(*,8242) TipeECP
!        write(*,8243) Zlistecp(1:Ecptypes)
        call WriteZlistECP(ZlistECP,Ecptypes)
        write(*,8244) Fock_ECP_read
        write(*,8245) Fock_ECP_write
        write(*,8246) cutECP
        write(*,8247) cut2_0
        write(*,8248) cut3_0
        write(*,8249) Verbose_ECP
        write(*,8250) ECP_debug
        write(*,8251) fulltimer_ECP
        write(*,8252) local_nonlocal
        write(*,8253) ECP_full_range_int
        write(*,8003)

!TDDFT
        write(*,8000)
        write(*,8104)
        write(*,8002)
        write(*,8260) Timedep
        write(*,8261) Propagator
        write(*,8262) Ntdstep
        write(*,8263) Tdstep
        write(*,8264) Nbsh
        write(*,8265) Field
        write(*,8266) Fx
        write(*,8267) Fy
        write(*,8268) Fz
        write(*,8269) Tdrestart
        write(*,8314) Exter
        write(*,8003)

!Write Option
        write(*,8000)
        write(*,8105)
        write(*,8002)
        write(*,8280) Verbose
        write(*,8281) WriteXYZ
        write(*,8282) WriteDens
        write(*,8283) WriteForces
        write(*,8003)


!Restart
        write(*,8000)
        write(*,8106)
        write(*,8002)
	write(*,8290) VCinp
        write(*,8291) Frestartin
        write(*,8292) Rho_restart_in
        write(*,8293) Rho_restart_out
        write(*,8003)

!Others
        write(*,8000)
        write(*,8107)
        write(*,8002)
        write(*,8303) PredCoef
        write(*,8304) Idip
        write(*,8305) IntSolDouble
        write(*,8306) DgTRIG
        write(*,8307) Iexch
        write(*,8308) Integ
        write(*,8309) Dens
        write(*,8310) Igrid
        write(*,8311) Igrid2
        write(*,8312) A0
        write(*,8313) Epsilon
        write(*,8315) Int_basis
        write(*,8316) Cubegen_only
        write(*,8317) Cube_Res
        write(*,8318) Cube_Dens
        write(*,8319) Cube_Dens_file
        write(*,8320) Cube_Orb
        write(*,8321) Cube_Sel
        write(*,8322) Cube_Orb_File
        write(*,8323) Cube_Elec
        write(*,8324) Cube_Elec_File
        write(*,8003)





 8000 FORMAT(4x,"╔═════════════════════════════════", &
 "═════════════════╗")
 8001 FORMAT(4x,"╚═════════════════════════════════", &
 "═════════════════╝")
 8002 FORMAT(4x,"╠══════════════════════╦══════════", &
 "═════════════════╣")
 8003 FORMAT(4x,"╚══════════════════════╩══════════", &
 "═════════════════╝")

 8100 FORMAT(4x,"║                     LIO Input                    ║")
 8101 FORMAT(4x,"║                      System                      ║")
 8102 FORMAT(4x,"║                   Theory Level                   ║")
 8103 FORMAT(4x,"║             Effective Core Potential             ║")
 8104 FORMAT(4x,"║                       TDDFT                      ║")
 8105 FORMAT(4x,"║                   Write Option                   ║")
 8106 FORMAT(4x,"║                      Restart                     ║")
 8107 FORMAT(4x,"║                      Others                      ║")

!System
 8210 FORMAT(4x,"║  Natom               ║  ",7x,i6,12x,"║")
 8211 FORMAT(4x,"║  Nsol                ║  ",7x,i8,10x,"║")
 8212 FORMAT(4x,"║  Charge              ║  ",9x,i5,11x,"║")
 8300 FORMAT(4x,"║  Nunp                ║  ",i5,20x,"║")

!Theory level
 8220 FORMAT(4x,"║  Open                ║  ",l,23x,"║")
 8221 FORMAT(4x,"║  Nmax                ║  ",i5,20x,"║")
 8222 FORMAT(4x,"║  Basis_Set           ║  ",a24," ║")
 8223 FORMAT(4x,"║  Fitting_Set         ║  ",a24," ║")
 8224 FORMAT(4x,"║  Int_Basis           ║  ",l,23x,"║")
 8225 FORMAT(4x,"║  Diis                ║  ",l,23x,"║")
 8226 FORMAT(4x,"║  Ndiis               ║  ",i3,22x,"║")
 8227 FORMAT(4x,"║  Gold                ║  ",f14.8,11x,"║")
 8228 FORMAT(4x,"║  Told                ║  ",f14.8,11x,"║")
 8229 FORMAT(4x,"║  Etold               ║  ",f14.8,11x,"║")
 8230 FORMAT(4x,"║  Hybrid_converg      ║  ",l,23x,"║")
 8231 FORMAT(4x,"║  Good_cut            ║  ",f14.8,11x,"║")
 8301 FORMAT(4x,"║  Rmax                ║  ",f14.8,11x,"║")
 8302 FORMAT(4x,"║  RmaxS               ║  ",f14.8,11x,"║")


!Effective Core Potential
 8240 FORMAT(4x,"║  Ecpmode             ║  ",l,23x,"║")
 8241 FORMAT(4x,"║  Ecptypes            ║  ", i3,22x,"║")
 8242 FORMAT(4x,"║  TipeECP             ║  ",a25,"║")
 8243 FORMAT(4x,"║  Zlistecp            ║  ",i2,i2,i2,i2,i2,i2,"║")
 8244 FORMAT(4x,"║  Fock_ECP_read       ║  ",l,23x,"║")
 8245 FORMAT(4x,"║  Fock_ECP_write      ║  ",l,23x,"║")
 8246 FORMAT(4x,"║  cutECP              ║  ",l,23x,"║")
 8247 FORMAT(4x,"║  cut2_0              ║  ",f14.8,11x,"║")
 8248 FORMAT(4x,"║  cut3_0              ║  ",f14.8,11x,"║")
 8249 FORMAT(4x,"║  Verbose_ECP         ║  ",i2,23x,"║")
 8250 FORMAT(4x,"║  ECP_debug           ║  ",l,23x,"║")
 8251 FORMAT(4x,"║  fulltimer_ECP       ║  ",l,23x,"║")
 8252 FORMAT(4x,"║  local_nonlocal      ║  ",i2,23x,"║")
 8253 FORMAT(4x,"║  ECP_full_range_int  ║  ",l,23x,"║")

!TDDFT
 8260 FORMAT(4x,"║  Timedep             ║  ",i2,23x,"║")
 8261 FORMAT(4x,"║  Propagator          ║  ",i2,23x,"║")
 8262 FORMAT(4x,"║  Ntdstep             ║  ",i10,15x,"║")
 8263 FORMAT(4x,"║  Tdstep              ║  ",f14.8,11x,"║")
 8264 FORMAT(4x,"║  Nbsh                ║  ",i4,21x,"║")
 8265 FORMAT(4x,"║  Field               ║  ",l,23x,"║")
 8266 FORMAT(4x,"║  Fx                  ║  ",f14.8,11x,"║")
 8267 FORMAT(4x,"║  Fy                  ║  ",f14.8,11x,"║")
 8268 FORMAT(4x,"║  Fz                  ║  ",f14.8,11x,"║")
 8269 FORMAT(4x,"║  Tdrestart           ║  ",l,23x,"║")
 8314 FORMAT(4x,"║  Exter               ║  ",l,23x,"║")


!Write Options
 8280 FORMAT(4x,"║  Verbose             ║  ",l,23x,"║")
 8281 FORMAT(4x,"║  WriteXYZ            ║  ",l,23x,"║")
 8282 FORMAT(4x,"║  WriteDens           ║  ",l,23x,"║")
 8283 FORMAT(4x,"║  WriteForces         ║  ",l,23x,"║")

!Restart
 8290 FORMAT(4x,"║  VCinp               ║  ",l,23x,"║")
 8291 FORMAT(4x,"║  Frestartin          ║  ",a25,"║")
 8292 FORMAT(4x,"║  Rho_restart_in      ║  ",l,23x,"║")
 8293 FORMAT(4x,"║  Rho_restart_out     ║  ",l,23x,"║")

!Others
 8303 FORMAT(4x,"║  PredCoef            ║  ",l,23x,"║")
 8304 FORMAT(4x,"║  Idip                ║  ",i4,21x,"║")
 8305 FORMAT(4x,"║  IntSolDouble        ║  ",l,23x,"║")
 8306 FORMAT(4x,"║  DgTRIG              ║  ",f14.8,11x,"║")
 8307 FORMAT(4x,"║  Iexch               ║  ",i5,20x,"║")
 8308 FORMAT(4x,"║  Integ               ║  ",l,23x,"║")
 8309 FORMAT(4x,"║  Dens                ║  ",l,23x,"║")
 8310 FORMAT(4x,"║  Igrid               ║  ",i3,22x,"║")
 8311 FORMAT(4x,"║  Igrid2              ║  ",i3,22x,"║")
 8312 FORMAT(4x,"║  A0                  ║  ",f14.7,11x,"║")
 8313 FORMAT(4x,"║  Epsilon             ║  ",f14.7,11x,"║")
 8315 FORMAT(4x,"║  Int_basis           ║  ",l,23x,"║")
 8316 FORMAT(4x,"║  Cubegen_only        ║  ",l,23x,"║")
 8317 FORMAT(4x,"║  Cube_Res            ║  ",i5,20x,"║")
 8318 FORMAT(4x,"║  Cube_Dens           ║  ",l,23x,"║")
 8319 FORMAT(4x,"║  Cube_Dens_file      ║  ",a25,"║")
 8320 FORMAT(4x,"║  Cube_Orb            ║  ",l,23x,"║")
 8321 FORMAT(4x,"║  Cube_Sel            ║  ",i5,20x,"║")
 8322 FORMAT(4x,"║  Cube_Orb_File       ║  ",a25,"║")
 8323 FORMAT(4x,"║  Cube_Elec           ║  ",l,23x,"║")
 8324 FORMAT(4x,"║  Cube_Elec_File      ║  ",a25,"║")

	END SUBROUTINE

        SUBROUTINE WriteZlistECP(ZlistECP,D)
        INTEGER, INTENT(IN), DIMENSION(128) :: ZlistECP !Z de atomos con ECP
        INTEGER, INTENT(IN) :: D
        INTEGER :: i,k,lines,rest

        IF (D .LT. 6) THEN
          IF (D .eq. 1) WRITE(*,8538) ZlistECP(1)
          IF (D .eq. 2) WRITE(*,8539) ZlistECP(1:2)
          IF (D .eq. 3) WRITE(*,8540) ZlistECP(1:3)
          IF (D .eq. 4) WRITE(*,8541) ZlistECP(1:4)
          IF (D .eq. 5) WRITE(*,8542) ZlistECP(1:5)
        ELSE
          lines=D/6
          rest=mod(D,6)
          WRITE(*,8543) ZlistECP(1:6)
          DO i=1,lines-1
            k=6*i+1
            WRITE(*,8544) ZlistECP(k:k+5)
          END DO
          IF (rest .EQ. 1) WRITE(*,8545) ZlistECP(6*lines+1:D)
          IF (rest .EQ. 2) WRITE(*,8546) ZlistECP(6*lines+1:D)
          IF (rest .EQ. 3) WRITE(*,8547) ZlistECP(6*lines+1:D)
          IF (rest .EQ. 4) WRITE(*,8548) ZlistECP(6*lines+1:D)
          IF (rest .EQ. 5) WRITE(*,8549) ZlistECP(6*lines+1:D)
        END IF

 8538 FORMAT(4x,"║  Zlistecp            ║ ",i3,"                       ║")
 8539 FORMAT(4x,"║  Zlistecp            ║ ",i3,i3,"                    ║")
 8540 FORMAT(4x,"║  Zlistecp            ║ ",i3,i3,i3,"                 ║")
 8541 FORMAT(4x,"║  Zlistecp            ║ ",i3,i3,i3,i3,"              ║")
 8542 FORMAT(4x,"║  Zlistecp            ║ ",i3,i3,i3,i3,i3,"           ║")
 8543 FORMAT(4x,"║  Zlistecp            ║ ",i3,i3,i3,i3,i3,i3,"        ║")
 8544 FORMAT(4x,"║                      ║ ",i3,i3,i3,i3,i3,i3,"        ║")
 8545 FORMAT(4x,"║                      ║ ",i3"                        ║")
 8546 FORMAT(4x,"║                      ║ ",i3,i3,"                    ║")
 8547 FORMAT(4x,"║                      ║ ",i3,i3,i3,"                 ║")
 8548 FORMAT(4x,"║                      ║ ",i3,i3,i3,i3"               ║")
 8549 FORMAT(4x,"║                      ║ ",i3,i3,i3,i3,i3"            ║")

        END SUBROUTINE
