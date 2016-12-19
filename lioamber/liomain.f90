!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% LIOMAIN.F90  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! This file contains the liomain subroutine, which performs several common     !
! procedures before and after calling SCF either by LIO alone or when          !
! performing in tantem with AMBER/GROMACS. Currently the routine:              !
! * Allocates matrices needed for SCF calculations.                            !
! * Calls SCF or SCFOP for closed/open shell calculations.                     !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
 
      subroutine liomain(E, dipxyz)
      use garcha_mod, only : M, Smat, RealRho, OPEN
      use ecp_mod,    only : ecpmode
 
          implicit none
          REAL*8  :: dipxyz(3), E

          if (.not.allocated(Smat))    allocate(Smat(M,M))
          if (.not.allocated(RealRho)) allocate(RealRho(M,M))

          if(OPEN) then
              if (ecpmode) stop "ECP is unavailable for Open Shell systems."
              call SCFOP(E, dipxyz)
          else
              call SCF(E, dipxyz)
          endif

          return
      end subroutine liomain
