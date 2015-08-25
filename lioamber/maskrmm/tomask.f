!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!  FFR_SETUPS
!------------------------------------------------------------------------------!
!
! General use subroutines
!
! ffr_reset_mem()
! ffr_setup_S(Enn_o,Esve_o,Esvn_o)
! ffr_setup_F(Ecoul_o,Exc_o,Eelec_o,Efld_o)
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       SUBROUTINE ffr_setup_S(Enn_o,Esve_o,Esvn_o)
       use garcha_mod, only:nsol
       implicit none
       real*8                      :: Enn,Esve,Esvn
       real*8,intent(out),optional :: Enn_o,Esve_o,Esvn_o
!------------------------------------------------------------------------------!
       call ffr_reset_mem()
       Enn =0.0d0
       Esve=0.0d0
       Esvn=0.0d0

       call int22()
       call int3mem()
       call int3mems()
       call int1(Enn)
       if(nsol.gt.0) then
         call intsol(Esve,Esvn,.true.)
       endif

       if (present(Enn_o))  Enn_o=Enn
       if (present(Esve_o)) Esve_o=Esve
       if (present(Esvn_o)) Esvn_o=Esvn
       RETURN;END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       SUBROUTINE ffr_setup_F(Ecoul_o,Exc_o,Eelec_o,Efld_o)
       use garcha_mod, only:RMM,M,Md,field
       implicit none
       real*8,intent(out),optional :: Ecoul_o,Exc_o,Eelec_o,Efld_o
       real*8                      :: Ecoul,Exc,Eelec,Efld
       integer                     :: MM,M11,kk
!------------------------------------------------------------------------------!
       MM=M*(M+1)/2
       M11=1+3*MM+Md*(Md+1)
       Ecoul = 0.0d0
       Exc   = 0.0d0
       Eelec = 0.0d0
       Efld  = 0.0d0

       call int3lu(Ecoul)
       call g2g_solve_groups(0,Exc,0)
       if (field) then
!         call pert_gauss(Efld)
       endif
       do kk=1,MM
         Eelec=Eelec+RMM(kk)*RMM(M11+kk-1)
       enddo

       if (present(Ecoul_o)) Ecoul_o = Ecoul
       if (present(Exc_o))   Exc_o   = Exc
       if (present(Eelec_o)) Eelec_o = Eelec
       if (present(Efld_o))  Efld_o  = Efld
       RETURN;END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
