      subroutine cgvc( na, xa, fa, cell, stress, volume, dxmax, 
     .                 tp, ftol, strtol, varcel, relaxd, usesavecg )   
c ***************************************************************************
c Variable-cell conjugate-gradient geometry optimization
c
c   Energy minimization including atomic coordinates and unit cell vectors.
c   It allows an external target stress:
c              %block MD.TargetStress
c                  3.5  0.0  0.0  0.0  0.0  0.0
c              %endblock MD.TargetStress
c   corresponding to xx, yy, zz, xy, xz, yz.
c   In units of (-MD.TargetPressure)
c   Default: hydrostatic pressure: -1, -1, -1, 0, 0, 0
c
c   Based on E({xa},stress), with {xa} in fractional coor
c   The gradient of the energy given by {cfa} forces (fractional!)
c   The fractional coordinates are multiplied by the initial cell vectors
c   to get them in Bohr for dxmax and preconditioning.
c      
c Written by E. Artacho. November 1998. 
c ******************************** INPUT ************************************
c integer na            : Number of atoms in the simulation cell
c real*8 fa(3,na)       : Atomic forces
c real*8 stress(3,3)    : Stress tensor components
c real*8 volume         : unit cell volume
c real*8 dxmax          : Maximum atomic (or lattice vector) displacement
c real*8 tp             : Target pressure
c real*8 ftol           : Maximum force tolerance
c real*8 strtol         : Maximum stress tolerance
c logical varcel        : true if variable cell optimization
c *************************** INPUT AND OUTPUT ******************************
c real*8 xa(3,na)       : Atomic coordinates
c                         input: current step; output: next step
c real*8 cell(3,3)      : Matrix of the vectors defining the unit cell 
c                         input: current step; output: next step
c                         cell(ix,ja) is the ix-th component of ja-th vector
c ******************************** OUTPUT ***********************************
c logical relaxd        : True when converged
c ***************************************************************************

C
C  Modules
C
      use precision
      use fdf

      implicit          none

      integer           na

      logical           relaxd, varcel, usesavecg

      double precision  cell(3,3), stress(3,3), 
     .                  xa(3,na), fa(3,na),
     .                  tp, ftol, strtol, volume, dxmax

c ---------------------------------------------------------------------------

c Internal variables and arrays

      logical           tarstr, frstme, found

      integer           iu, ia, ndeg, i, j, n, indi, linmin, Node

      double precision  
     .                  cgcntr(0:20), tstres(3,3), celli(3,3), 
     .                  modcel(3), kBar, precon, strain(3,3), 
     .                  sxx, syy, szz, sxy, sxz, syz, cellin(3,3)

      double precision, dimension(:), allocatable, save :: gxa, gfa
      double precision, dimension(:), allocatable, save :: cgaux

      save              frstme, ndeg, tstres, cgcntr, kBar,
     .                  modcel, precon, strain, cellin, linmin
    

      data              frstme      /.true./,
     .                  tarstr      /.false./

c ---------------------------------------------------------------------------

      Node = 0

C Allocate local memory
      allocate(gfa(na*3 + 6))
      allocate(gxa(na*3 + 6))
      if (.not.allocated(cgaux)) then
        allocate(cgaux(2*(na*3 + 6)))
      endif

c If first call to cgvc, check dim and get target stress --------------------

      if ( frstme ) then
  
c look for target stress and read it if found, otherwise generate it --------

        if ( varcel ) then
     
          kBar = 1.d0/1.47108d5

          if (Node.eq.0) then
            tarstr = fdf_block('MD.TargetStress',iu)

            if (tarstr) then
               write(6,'(/a,a)') 'cgvc: Reading %block MD.TargetStress',
     .                        ' (units of MD.TargetPressure).'
               read(iu,*, end=50) sxx, syy, szz, sxy, sxz, syz
               tstres(1,1) = - sxx * tp
               tstres(2,2) = - syy * tp
               tstres(3,3) = - szz * tp
               tstres(1,2) = - sxy * tp
               tstres(2,1) = - sxy * tp
               tstres(1,3) = - sxz * tp
               tstres(3,1) = - sxz * tp
               tstres(2,3) = - syz * tp
               tstres(3,2) = - syz * tp
   50          continue
            else
              write(6,'(/a,a)') 'cgvc: No target stress found, ',
     .                'assuming hydrostatic MD.TargetPressure.'
              do i = 1, 3
                do j = 1, 3
                  tstres(i,j) = 0.d0
                enddo
                tstres(i,i) = - tp
              enddo
            endif

c write target stress down --------------------------------------------------

            write(6,"(/a)") 'cgvc: Target stress (kBar)'
            write(6,"(a,2x,3f12.3)") 
     .       'cgvc:', tstres(1,1)/kBar, tstres(1,2)/kBar, 
     .       tstres(1,3)/kBar
            write(6,"(a,2x,3f12.3)") 
     .       'cgvc:', tstres(2,1)/kBar, tstres(2,2)/kBar, 
     .       tstres(2,3)/kBar
            write(6,"(a,2x,3f12.3)") 
     .       'cgvc:', tstres(3,1)/kBar, tstres(3,2)/kBar, 
     .       tstres(3,3)/kBar
          endif

c moduli of original cell vectors for fractional coor scaling back to au ---

          do n = 1, 3
             modcel(n) = 0.0d0
             do j = 1, 3
                modcel(n) = modcel(n) + cell(j,n)*cell(j,n)
             enddo
             modcel(n) = dsqrt( modcel(n) )
          enddo

          if (Node.eq.0) then
c scale factor for strain variables to share magnitude with coordinates -----
c ---- (a length in Bohrs typical of bond lengths ..) -----------------------

c           precon = modcel(1)
            precon = fdf_physical('MD.PreconditionVariableCell',
     .                           9.4486344d0,'Bohr')
          endif

c dimension of space where E is minimized -----------------------------------

          ndeg = na*3 + 6

c initialize absolute strain and save initial cell vectors -----------------
c initialization to 1. for numerical reasons, later substracted ------------

          do i = 1, 3
            do j = 1, 3
              strain(i,j) = 1.d0
              cellin(i,j) = cell(i,j)
            enddo
          enddo

        else

          ndeg = na*3

        endif

c initialize and read cgaux and cgcntr if present and wanted ---------------

        if (usesavecg) then
          if (Node.eq.0) then
            call iocg( 'read', ndeg*2, cgaux, cgcntr, relaxd, found )
          endif
          if ( found ) then
            linmin = cgcntr(1)
          else
            if (Node.eq.0) then
              write(6,'(/,a)') 'cgvc: WARNING: CG file not found'
            endif
            relaxd = .false.
            cgcntr(0) = 0
            linmin = 1
          endif
        else
          relaxd = .false.
          cgcntr(0) = 0
          linmin = 1
        endif

        frstme = .false.
      endif

c variable cell -------------------------------------------------------------

      if ( varcel ) then

c inverse of matrix of cell vectors  (transverse of) ------------------------

        call reclat( cell, celli, 0 )

c transform coordinates and forces to fractional ---------------------------- 
c but scale them again to Bohr by using the (fix) moduli of the original ----
c lattice vectors (allows using maximum displacement as before) -------------
c convergence is checked here for input forces as compared with ftol --------

        relaxd = .true.
        do ia = 1, na
          do n = 1, 3
            indi = 3*(ia - 1) + n
            gxa(indi) = 0.0d0
            gfa(indi) = 0.0d0
            do i = 1, 3
              gxa(indi) = gxa(indi) + celli(i,n) * xa(i,ia) * modcel(n)
              gfa(indi) = gfa(indi) + cell(i,n) * fa(i,ia) / modcel(n)
              relaxd = relaxd .and. ( dabs(fa(i,ia)) .lt. ftol )
            enddo
          enddo
        enddo

c symmetrizing the stress tensor --------------------------------------------

        do i = 1, 3
           do j = i+1, 3
              stress(i,j) = 0.5d0*( stress(i,j) + stress(j,i) )
              stress(j,i) = stress(i,j)
           enddo
        enddo

c append stress (substracting target stress) and strain to gxa and gfa ------ 
c preconditioning: scale stress and strain so as to behave similarly to x,f -
 
        indi = 3*na
        do i = 1, 3
           do j = i, 3
              indi = indi + 1
              gfa(indi) = -(stress(i,j) - tstres(i,j))*volume/precon
              gxa(indi) = strain(i,j) * precon
           enddo
        enddo

c check stress convergence --------------------------------------------------

        strtol = dabs(strtol)
        do i = 1, 3
           do j = 1, 3
              relaxd = relaxd .and. 
     .          ( dabs(stress(i,j)-tstres(i,j)) .lt. strtol )
           enddo
        enddo

c call conjugate gradient minimization -------------------------------------- 

        if ( .not. relaxd )
     .     call conjgr( ndeg, gxa, gfa, dxmax, 0.d0, cgcntr, cgaux )


c fixed cell ----------------------------------------------------------------

      else

        call conjgr( 3*na, xa, fa, dxmax, ftol, cgcntr, cgaux )
        relaxd = int(cgcntr(0)) .eq. 0

      endif

c checking line minimizations and convergence -------------------------------

      if (nint(cgcntr(1)) .ne. linmin) then
        if (Node.eq.0) then
          write(6,'(/a,i4,a,f10.4)')
     .      'cgvc: Finished line minimization ', linmin,
     .      '.  Mean atomic displacement =', cgcntr(18)/sqrt(dble(na))
        endif
        linmin = nint(cgcntr(1))
      endif

c transform back if variable cell ------------------------------------------- 

      if ( varcel .and. (.not. relaxd) ) then

c new cell ------------------------------------------------------------------

        indi = 3*na
        do i = 1, 3
           do j = i, 3
              indi = indi + 1
              strain(i,j) = gxa(indi) / precon
              strain(j,i) = strain(i,j)
           enddo
        enddo

        do n = 1, 3
           do i = 1, 3
              cell(i,n) = cellin(i,n)
              do j = 1, 3
                 cell(i,n) = cell(i,n) + (strain(i,j)-1.d0)*cellin(j,n)         
              enddo
           enddo
        enddo

c output fractional coordinates to cartesian Bohr, and copy to xa ----------- 

        do ia = 1, na
          do i = 1, 3
            xa(i,ia) = 0.0d0
            do n = 1, 3
              indi = 3*(ia - 1) + n
              xa(i,ia) = xa(i,ia) + cell(i,n) * gxa(indi) / modcel(n)
            enddo
          enddo
        enddo

      endif

c save cgaux ----------------------------------------------------------------

      if (Node.eq.0) 
     .  call iocg( 'write', ndeg*2, cgaux, cgcntr, relaxd, found )

C Deallocate local memory
      deallocate(gfa)
      deallocate(gxa)

      return
      end

