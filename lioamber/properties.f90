!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% PROPERTIES.F90 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! This file contains several molecular properties calculation. Currently       !
! includes:                                                                    !
! * get_degeneration (gets degeneration and degenerated MOs for a chosen MO)   !
! * do_forces        (calculates forces/gradients)                             !
! * do_dipole        (calculates dipole moment)                                !
! Regarding Electronic Population Analysis:                           [ EPA ]  !
! * do_population_analysis (performs the analysis required)                    !
! * mulliken_calc    (calculates atomic Mulliken population charges)           !
! * lowdin_calc      (calculates atomic Löwdin population charges)             !
! Regarding Reactivity Indexes:                                       [ RXI ]  !
! * do_fukui         (performs Fukui function calculation and printing)        !
! * get_softness     (gets the molecule's global softness)                     !
! * fukui_calc       (calculates CS condensed-to-atoms fukui function)         !
! * fukui_calc_os    (calculates CS condensed-to-atoms fukui function)         !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% GET_DEGENERATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Gets degeneration number (nDeg) for the MO of interest (orb) and the MO      !
! indexes (nDegMO) of those MOs which share degeneration.                      !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine get_degeneration(energ, orb, M, nDeg, nDegMO)

    implicit none
    real*8 , intent(in)  :: energ(M)
    integer, intent(in)  :: orb, M
    integer, intent(out) :: nDeg, nDegMO(M)
  
    integer :: i, iCountU, iCountD
    real*8  :: criterium, ratio
         
    if (M.eq.1) then
        nDeg = 1
        return
    endif

    ratio = 0d0 ; nDeg = 0 ; criterium = 0.000005 ; nDegMO = 0 ;

    ! Assings counting limits in order to avoid exploring all MOs.
    ! Special cases: only one MO, or MO of interest is the first one.
    if (orb.gt.1) then 
        iCountU = orb ; iCountD = orb - 1
    else
        iCountU = 2   ; iCountD = 1
    endif

    ! Performs energy comparisons, upwards and downwards.
    do i=iCountU, M, 1
        ratio = abs((energ(i) - energ(orb)) / ((energ(i) + energ(orb))/2))
        if (ratio.le.criterium) then
            nDeg = nDeg + 1 ; nDegMO(nDeg) = i
        else 
            exit
        endif
    enddo
          
    do i=iCountD, 1, -1
        ratio = abs((energ(i) - energ(orb)) / ((energ(i) + energ(orb))/2))
        if (ratio.le.criterium) then 
            nDeg = nDeg + 1 ; nDegMO(nDeg) = i
        else
            exit
        endif
    enddo

    return
end subroutine get_degeneration
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% DO_FORCES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Calculates forces for QM and MM regions and writes them to output.           !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine do_forces(uid)

    use garcha_mod, only : natom, nsol

    implicit none
    integer, intent(in) :: uid
    integer             :: k
    real*8, allocatable :: dxyzqm(:,:), dxyzcl(:,:)

    open(unit=uid, file='forces')

    allocate ( dxyzqm(3, natom) )
    dxyzqm = 0.0

    call dft_get_qm_forces(dxyzqm)
    if (nsol.gt.0) then
        allocate ( dxyzcl(3, natom+nsol) )
        dxyzcl = 0.0
        call dft_get_mm_forces(dxyzcl, dxyzqm)
    endif

    call write_forces(dxyzqm, natom, 0, uid)
    deallocate (dxyzqm)
    
    if(nsol.gt.0) then
        call write_forces(dxyzcl, nsol, natom, uid)       
        deallocate (dxyzcl)
    endif

    return
end subroutine do_forces
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% DO_DIPOLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Sets variables up and calls dipole calculation.                              !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine do_dipole(dipxyz, uid)
    implicit none
    integer, intent(in)    :: uid
    real*8 , intent(inout) :: dipxyz(3)
    real*8                 :: u

    call g2g_timer_sum_start('dipole')
    call dipnew(dipxyz)
    u = sqrt(dipxyz(1)**2 + dipxyz(2)**2 + dipxyz(3)**2)
 
    call write_dipole(dipxyz, u, uid)
    call g2g_timer_sum_stop('dipole')
 
    return
end subroutine do_dipole
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%% ELECTRONIC POPULATION ANALYSIS [ EPA ] %%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% DO_POPULATION_ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Performs the different population analyisis available.                       !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine do_population_analysis()
   use garcha_mod, only : RMM, Smat, RealRho, M, Enucl, Nuc, Iz, natom, &
                          mulliken, lowdin, sqsm
   use ECP_mod   , only : ecpmode, IzECP

   implicit none
   integer :: M1, M5, IzUsed(natom), kk
   real*8  :: q(natom)

   ! Needed until we dispose of RMM.
   M1=1 ; M5=1+M*(M+1)

   ! Iz used to write the population file.
   IzUsed = Iz
   if (ecpmode) IzUsed = IzECP

   ! Decompresses and fixes S and RealRho matrixes, which are needed for
   ! population analysis.
   call int1(Enucl)
   call spunpack('L',M,RMM(M5),Smat)
   call spunpack('L',M,RMM(M1),RealRho)
   call fixrho(M,RealRho)

   ! Performs Mulliken Population Analysis if required.
   if (mulliken) then
       call mulliken_calc(natom,M,RealRho,Smat,Nuc,Iz,q)
       call write_population(85,natom,IzUsed,q,0)
   endif
   ! Performs Löwdin Population Analysis if required.
   if (lowdin) then 
       do kk=1,natom
           q(kk)=real(Iz(kk))
       enddo
       call lowdin_calc(M,natom,RealRho,sqsm,Nuc,q)
       call write_population(85,natom,IzUsed,q,1)
   endif

   return
endsubroutine do_population_analysis

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% MULLIKEN_CALC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Performs a Mulliken Population Analysis and outputs atomic charges.          !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine mulliken_calc(N, M, RealRho, Smat, NofM, q0, q)

    ! RealRho         : Rho written in atomic orbitals.                !
    ! q0, q           : Starting and final Mulliken charges.           !
    ! M, N, NofM, Smat: N° of basis functions, atoms, Nuclei belonging !
    !                   to each function M, overlap matrix.            !
    implicit none
    integer, intent(in)  :: N, M, NofM(M), q0(N)
    real*8 , intent(in)  :: RealRho(M,M), Smat(M,M)
    real*8 , intent(out) :: q(N)

    integer :: i, j, k
    real*8  :: qe

    call g2g_timer_start('mulliken')

    do k=1,N
        q(k)=real(q0(k))
    enddo

    do i=1,M
        do j=1,M
             qe = RealRho(i, j) * Smat(i, j)
             q(NofM(i)) = q(NofM(i)) - qe
        enddo
    enddo

    call g2g_timer_stop('mulliken')

    return
end subroutine mulliken_calc
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%% LOWDIN_CALC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Performs a Löwdin Population Analysis and outputs atomic charges.            !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine lowdin_calc(M, N, rhomat, sqsmat, atomorb, atomicq)
 
    implicit none
    integer,intent(in)   :: M, N, atomorb(M)
    real*8,intent(in)    :: rhomat(M,M), sqsmat(M,M)
    real*8,intent(inout) :: atomicq(N)

    real*8  :: newterm
    integer :: natom
    integer :: i, j, k

    do k=1, M
        natom=atomorb(k)
        do i=1, M
            do j=1, M
                newterm = sqsmat(k, i) * rhomat(i, j) * sqsmat(j, k)
                atomicq(natom) = atomicq(natom) - newterm
            enddo
        enddo
    enddo

    return
end subroutine lowdin_calc
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%% REACTIVITY INDEXES [ RXI ] %%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% GET_SOFTNESS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Gets the molecule's global softness using the energy of HOMO/LUMO ALPHA/BETA !
! MOs. In a closed-shell context, enAH=enBH and enAL=enBL.                     !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine get_softness(enAH, enAL, enBH, enBL, softness)

    implicit none
    real*8, intent(in)  :: enAH, enAL, enBH, enBL
    real*8, intent(out) :: softness

    softness = 4 / (enAH + enBH - enAL - enBL)

    return
end subroutine get_softness
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% DO_FUKUI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Performs Fukui function calls and printing.                                  !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine do_fukui()
    use garcha_mod, only : X, NCO, M, natom, Nuc, Smat, Eorbs, Iz, OPEN

    implicit none
    real*8  :: fukuim(natom), fukuin(natom), fukuip(natom), softness

    if (OPEN) then
    else
        call fukui_calc(X(1,M*2+1), NCO, M, natom, Nuc, Smat, fukuim, fukuip, &
                        fukuin, Eorbs)
        call get_softness(Eorbs(NCO-1), Eorbs(NCO), Eorbs(NCO-1), Eorbs(NCO), &
                          softness)
        call write_fukui(fukuim, fukuip, fukuin, natom, Iz, softness)
    endif
 
    return
end subroutine do_fukui
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% FUKUI_CALC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Performs a closed-shell Fukui function calculation.                          !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine fukui_calc(coef, nOcc, M, N, NofM, Smat, fukuiNeg, fukuiPos,  &
                      fukuiRad, ener)

    ! coef, nOcc      : MO basis coefficients and N° occupied MOs.     !
    ! enAlpha, enBeta : Alpha and Beta MO energies.                    !
    ! M, N, NofM, Smat: N° of basis functions, atoms, Nuclei belonging !
    !                   to each function M, overlap matrix.            !
    ! FukuiXXX        : Fukui function for radical attack (Rad),       !
    !                   electron loss (Pos) and electron gain (Neg).   !
    ! shapeX          : Shape factor for the MO of interest = Atomic   !
    !                   contributions to HOMO/LUMO MOs.                !
    ! nDegX, degMOX   : Degeneration of such MOs.                      !
    implicit none
    integer, intent(in)  :: M, N, NofM(M), nOcc
    real*8 , intent(in)  :: coef(M,M)  , Smat(M,M)  , ener(M) 
    real*8 , intent(out) :: fukuiRad(N), fukuiPos(N), fukuiNeg(N)

    integer :: i, j, k, nDegH, nDegL
    real*8  :: shapeH(N), shapeL(N), dummy
    integer, allocatable :: degMOH(:), degMOL(:)

    ! Calculates MO degenerations and MOs with same energy.
    allocate(degMOH(M), degMOL(M))
    call get_degeneration(ener, nOcc   , M, nDegH, degMOH)
    call get_degeneration(ener, nOcc +1, M, nDegL, degMOL)

    ! Calculates shape factors and Spin-Polarized Fukui functions.
    shapeH = 0d0  ;  shapeL = 0d0

    do i=1, M
        do j=1, M
            do k=1, nDegH
                dummy           = coef(i, degMOH(k)) *         &
                                  coef(j, degMOH(k)) * Smat(i,j)
                shapeH(NofM(i)) = shapeH(NofM(i)) + dummy
            enddo
            do k=1, nDegL
                dummy           = coef(i, degMOL(k)) *         &
                                  coef(j, degMOL(k)) * Smat(i,j)
                shapeL(NofM(i)) = shapeL(NofM(i)) + dummy
            enddo
        enddo
    enddo

    do i=1, N
        fukuiNeg(i) = shapeH(i) / nDegH
        fukuiPos(i) = shapeL(i) / nDegL
        fukuiRad(i) = ( fukuiNeg(i) + fukuiPos(i) ) / 2
    enddo

    deallocate (degMOH, degMOL)
    
    return
end subroutine fukui_calc
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% FUKUI_CALC_OS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Performs an open-shell Fukui function calculation (Spin-Polarized Fukui).    !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine fukui_calc_os(coefAlp, coefBet, nAlpha, nBeta, M, N, NofM, &
                               Smat, fukuiNeg, fukuiPos, fukuiRad, enAlpha, &
                               enBeta)

    ! coefAlp, coefBet: -> Alpha and Beta coefficients.                !
    ! nAlpha , nBeta  : N° occupied Alpha and Beta orbitals.           !
    ! enAlpha, enBeta : Alpha and Beta MO energies.                    !
    ! M, N, NofM, Smat: N° of basis functions, atoms, Nuclei belonging !
    !                   to each function M, overlap matrix.            !
    ! FukuiXXX        : Fukui function for radical attack (Rad),       !
    !                   electron loss (Pos) and electron gain (Neg).   !
    ! shapeXY         : Shape factor for the MO of interest = Atomic   !
    !                   contributions to HOMO/LUMO MOs.                !
    ! nDegXY, degMOXY : Degeneration of such MOs.                      !
    implicit none
    integer, intent(in)  :: M, N, NofM(M), nAlpha, nBeta
    real*8 , intent(in)  :: coefAlp(M,M), coefBet(M,M), Smat(M,M)
    real*8 , intent(in)  :: enAlpha(M)  , enBeta(M)
    real*8 , intent(out) :: fukuiRad(N) , fukuiPos(N) , fukuiNeg(N)
      
    integer :: i, j, k, nDegAH, nDegAL, nDegBH, nDegBL
    real*8  :: shapeAH(N), shapeAL(N), shapeBH(N), shapeBL(N), dummy
    integer, allocatable :: degMOAH(:), degMOAL(:), degMOBH(:), degMOBL(:)

    ! Calculates MO degenerations and MOs with same energy.
    allocate(degMOAH(M), degMOAL(M), degMOBH(M), degMOBL(M))
    call get_degeneration(enAlpha, nAlpha   , M, nDegAH, degMOAH)
    call get_degeneration(enAlpha, nAlpha +1, M, nDegAL, degMOAL)
    call get_degeneration(enBeta , nBeta    , M, nDegBH, degMOBH)
    call get_degeneration(enBeta , nBeta  +1, M, nDegBL, degMOBL)

    ! Calculates shape factors and Spin-Polarized Fukui functions.
    shapeAH = 0d0  ;  shapeAL = 0d0 ;  shapeBH = 0d0  ;  shapeBL = 0d0

    do i=1, M
        do j=1, M
            do k=1, nDegAH
                dummy            = coefAlp(i, degMOAH(k)) *         &
                                   coefAlp(j, degMOAH(k)) * Smat(i,j)
                shapeAH(NofM(i)) = shapeAH(NofM(i)) + dummy
            enddo
            do k=1, nDegAL
                dummy            = coefAlp(i, degMOAL(k)) *         &
                                   coefAlp(j, degMOAL(k)) * Smat(i,j)
                shapeAL(NofM(i)) = shapeAL(NofM(i)) + dummy
            enddo
            do k=1, nDegBH
                dummy            = coefBet(i, degMOBH(k)) *         &
                                   coefBet(j, degMOBH(k)) * Smat(i,j)
                shapeBH(NofM(i)) = shapeBH(NofM(i)) + dummy
            enddo
            do k=1, nDegBL
                dummy            = coefBet(i, degMOBL(k)) *         &
                                   coefBet(j, degMOBL(k)) * Smat(i,j)
                shapeBL(NofM(i)) = shapeBL(NofM(i)) + dummy
            enddo
        enddo
    enddo

    do i=1, N
        fukuiNeg(i) = ( shapeAH(i)/nDegAH + shapeBH(i)/nDegBH ) / 2
        fukuiPos(i) = ( shapeAL(i)/nDegAL + shapeBL(i)/nDegBL ) / 2
        fukuiRad(i) = ( fukuiNeg(i)       + fukuiPos(i)       ) / 2
    enddo 
         
    deallocate (degMOAH, degMOAL, degMOBH, degMOBL)
    
    return
end subroutine fukui_calc_os
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
