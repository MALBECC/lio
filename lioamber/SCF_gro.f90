!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% SCF_GRO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! File SCF_gro.f contains SCF_gro subroutine, which performs the SCF cycles    !
! within the molecular dynamics when called from Gromacs software package for a!
! QM/MM calculation. See also init_gromacs.f, which deals with initialization. !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

      subroutine SCF_gro(E, qmcoords, clcoords, clcharge, nsolin)
      use garcha_mod, only : nsol, ntatom, natom, r, v, Em, Rm, pc, nnat,      &
                             writexyz, rqm, Iz, RealRho, Smat, M, ng0, OPEN

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% Parameter Definition and Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
      implicit none
      REAL*8 , intent(in) :: qmcoords(3*natom), clcoords(3*nsolin), &
                             clcharge(nsolin )
      REAL*8  :: dipxyz(3), E
      integer :: ngDyn, i, j, n, nsolin

      dipxyz = 0d0          ; nsol   = nsolin    ;
      ntatom = nsol + natom ; ngDyn  = natom*ng0 ;

      call g2g_timer_sum_start("Total")

      deallocate (r, v, Em, Rm, pc, nnat)
      allocate (r(ntatom,3), v(ntatom,3), Em(ntatom), Rm(ntatom), pc(ntatom),  &
                nnat(100))

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% Coordinate Conversion and XYZ File Output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! This section converts the coordinates array and partial charges array        !
! received from Gromacs into the r (all positions), rqm (QM region positions)  !
! and pc (MM partial charges) arrays.                                          !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
      if(writexyz) write(18,*) ntatom
      if(writexyz) write(18,*)

      do i = 1, natom
          do j = 1, 3
              r(i,j)  = qmcoords((i-1)*3+j)
              rqm(i,j)= qmcoords((i-1)*3+j)
          enddo
          if(writexyz) write(18,345) Iz(i), rqm(i,:)*0.529
      enddo



      do i = 1, nsol
          n     = natom+i
          pc(n) = clcharge(i)
          do j = 1, 3
              r(n, j) = clcoords((i-1)*3+j)
          enddo
          if(writexyz) write(18,346) pc(n), r(n,:)*0.529
      enddo

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% External Routine Calls %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!      
! Currently liomain routine is empty, but it will contain some common variable !
! allocation routine calls. SCFOP calls for an open-shell SCF calculation,     !
! while SFC deals with closed-shell systems.                                   !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       call liomain()
       if (.not.allocated(Smat))    allocate(Smat(M,M))
       if (.not.allocated(RealRho)) allocate(RealRho(M,M))

      if(OPEN) then 
          call SCFOP(E, dipxyz)
      else
          call SCF(E, dipxyz)
      endif

 345  format(2x,I2,2x,3(f10.6,2x))
 346  format(2x,f10.6,2x,3(f10.6,2x))
      return;
      end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
