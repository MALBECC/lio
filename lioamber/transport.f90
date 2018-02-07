module transport_data
   implicit none
   logical   :: transport_calc   = .false.
   logical   :: generate_rho0    = .false.
   logical   :: gate_field       = .false.
   integer   :: save_charge_freq = 0
   integer   :: pop_drive        = 1
   integer   :: ngroup           = 0
   integer   :: pop_uid          = 678
   integer   :: drive_uid        = 5555
   real*8    :: scratchgamma     = 0.0D0
   real*8    :: driving_rate     = 0.001
   real*8    :: GammaMagnus      = 0.0D0
   real*8    :: GammaVerlet      = 0.0D0
   real*8    :: re_traza         = 0.0D0

   integer   , allocatable :: mapmat(:,:), group(:)
#ifdef TD_SIMPLE
   complex*8 , allocatable :: rhofirst(:,:)
#else
   complex*16, allocatable :: rhofirst(:,:)
#endif
end module transport_data

module transport_subs
   implicit none
contains

subroutine transport_init(M, natom, Nuc, RMM5, overlap, rho)
   use transport_data, only: GammaMagnus, GammaVerlet, ngroup, group, rhofirst,&
                             driving_rate, pop_drive, pop_uid, drive_uid, mapmat
   implicit none
   integer, intent(in)  :: M, natom, Nuc(M)
   real*8 , intent(in)  :: RMM5(M*(M+1)/2)
   real*8 , allocatable, intent(inout) :: overlap(:,:)
#ifdef TD_SIMPLE
   complex*8 , allocatable, intent(inout) :: rho(:,:)
#else
   complex*16, allocatable, intent(inout) :: rho(:,:)
#endif
   integer :: orb_group(M), icount

   ngroup  = 0
   allocate(rhofirst(M,M))
   allocate(group(natom), mapmat(M,M))
   call transport_read_groups(natom)
   call mat_map(Nuc, M, natom)

   do icount = 1, size(orb_group)
      orb_group(icount) = group( Nuc(icount) )
   enddo

   GammaMagnus = driving_rate
   GammaVerlet = driving_rate*0.1D0

   select case (pop_drive)
      case (1)
         open( unit = pop_uid  , file = "MullikenGroup")
         open( unit = drive_uid, file = 'DriveMul')
      case (2)
         open( unit = pop_uid  , file = "LowdinGroup")
         open( unit = drive_uid, file = 'DriveLowd')
      case default
         write(*,*) 'ERROR - Transport: Wrong value for pop_drive '//&
                    '(transport_init)'
   end select

   if (allocated(overlap)) deallocate(overlap)
   allocate(overlap(M,M))
   call spunpack('L', M, RMM5, overlap)
   call transport_generate_rho(M, rho)

   return
end subroutine transport_init

subroutine transport_read_groups(natom)
   use transport_data, only: group, ngroup
   implicit none
   integer, intent(in)  :: natom
   integer :: icount, group_read
   logical :: file_exists

   inquire(file = 'atomgroup', exist = file_exists)
   if (.not.file_exists) then
      write(*,*) ' ERROR - Transport: Cannot find atomgroup file.'
      stop
   endif

   open( unit = 678 , file = 'atomgroup')
   do icount = 1 , natom
      read(678, *) group_read
      group(icount) = group_read
      if (group_read.gt.ngroup) ngroup = group_read
   enddo
   close(unit=678)

   if(ngroup.gt.3) write(*,*) 'WARNING - Transport: If the number of &
   group is greater than 3, then group 1 should be the donor and 2 the &
   acceptor.'

   return
end subroutine transport_read_groups

subroutine transport_generate_rho(M, rho)
   use transport_data, only: generate_rho0, rhofirst
   implicit none
   integer, intent(in)  :: M
#ifdef TD_SIMPLE
   complex*8 , allocatable, intent(inout) :: rho(:,:)
#else
   complex*16, allocatable, intent(inout) :: rho(:,:)
#endif
   integer   :: icount, jcount

   open( unit = 100000, file = 'rhofirst')
   if (generate_rho0) then
      do icount = 1, M
      do jcount = 1, M
         write(100000,*) rho(icount, jcount)
      enddo
      enddo
      rhofirst = rho
      write(*,*) 'RhoFirst has been written.'
   else
      do icount = 1, M
      do jcount = 1, M
         read(100000,*) rhofirst(icount, jcount)
      enddo
      enddo
      write(*,*) 'RhoFirst has been read.'
   endif

   return
end subroutine transport_generate_rho

subroutine transport_rho_trace(M, rho)
   implicit none
   integer, intent(in) :: M
#ifdef TD_SIMPLE
   complex*8 , intent(in) :: rho(M,M)
#else
   complex*16, intent(in) :: rho(M,M)
#endif
   integer   :: icount
   complex*8 :: traza

   traza = dcmplx(0.0D0, 0.0D0)
   do icount = 1, M
      traza = traza + rho(icount, icount)
   enddo
   write(*,*) 'Transport - Rho Trace = ', real(traza)

   return
end subroutine transport_rho_trace

subroutine transport_propagate(M, natom, Nuc, Iz, propagator, istep, &
                               overlap, sqsm, rho1)
   use transport_data, only: save_charge_freq, pop_drive, GammaMagnus, &
                             GammaVerlet, ngroup, group
   implicit none
   integer   , intent(in)    :: M, natom, Nuc(M), Iz(natom), propagator, istep
   real*8    , intent(inout) :: overlap(M,M), sqsm(M,M)
#ifdef TD_SIMPLE
   complex*8 , intent(inout) :: rho1(M,M)
#else
   complex*16, intent(inout) :: rho1(M,M)
#endif
   real*8  :: scratchgamma, gamma
   integer :: save_freq

   save_freq = save_charge_freq
   select case (propagator)
   case (1)
      save_freq = save_freq*10
      gamma     = GammaVerlet
   case (2)
      gamma     = GammaMagnus
   end select

   call g2g_timer_start('Transport-propagation')
   if ((propagator.eq.2) .and. (istep.gt.999)) Then
      scratchgamma = gamma
   else
      scratchgamma = gamma * exp(-0.0001D0 * (dble(istep-1000))**2)
   endif

   call electrostat(rho1, overlap, scratchgamma, M)
   if (mod( istep-1 , save_freq) == 0) then
      call drive_population(M, natom, Nuc, Iz, rho1, overlap, sqsm, 1)
   endif

   call g2g_timer_stop('Transport-propagation')
   return
end subroutine transport_propagate

#ifdef CUBLAS
subroutine transport_propagate_cu(M, natom, Nuc, Iz, propagator, istep, &
                                  overlap, sqsm, rho1, devPtrY)
   use cublasmath    , only: basechange_cublas
   use transport_data, only: ngroup, group, save_charge_freq, GammaMagnus, &
                             GammaVerlet
   implicit none
   integer   , intent(in)    :: M, natom, Nuc(M), Iz(natom), propagator, istep
   integer*8 , intent(in)    :: devPtrY
   real*8    , intent(inout) :: overlap(M,M), sqsm(M,M)
#ifdef TD_SIMPLE
   complex*8 , intent(inout) :: rho1(M,M)
#else
   complex*16, intent(inout) :: rho1(M,M)
#endif
   real*8  :: scratchgamma, gamma
   integer :: save_freq

   save_freq = save_charge_freq
   select case (propagator)
   case (1)
      save_freq = save_freq*10
      gamma     = GammaVerlet
   case (2)
      gamma     = GammaMagnus
   end select

   call g2g_timer_start('Transport-propagation')
   if ((propagator.eq.2) .and. (istep.gt.999)) Then
      scratchgamma = gamma
   else
      scratchgamma = gamma * exp(-0.0001D0 * (dble(istep-1000))**2)
   endif

   call electrostat(rho1, overlap, scratchgamma, M)
   if (mod( istep-1 , save_freq) == 0) then
      call drive_population(M, natom, Nuc, Iz, rho1, overlap, sqsm, 1)
   endif

   call g2g_timer_start('complex_rho_ao_to_on-cu')
   rho1 = basechange_cublas(M, rho1, devPtrY, 'dir')
   call g2g_timer_stop('complex_rho_ao_to_on-cu')

   call g2g_timer_stop('Transport-propagation')

   return
end subroutine transport_propagate_cu
#endif

subroutine transport_population(M, natom, Nuc, Iz, rho1, overlap, smat, &
                                propagator, is_lpfrg, istep)
   use transport_data, only: save_charge_freq
   implicit none
   integer, intent(in) :: M, natom, Iz(natom), Nuc(natom), istep, propagator
   logical, intent(in) :: is_lpfrg
   real*8 , intent(in) :: overlap(M,M), smat(M,M)
#ifdef TD_SIMPLE
      complex*8 , intent(in) :: rho1(M,M)
#else
      complex*16, intent(in) :: rho1(M,M)
#endif

   if ( ((propagator.gt.1) .and. (is_lpfrg) .and.       &
      (mod((istep-1), save_charge_freq*10) == 0)) .or.  &
      (mod((istep-1), save_charge_freq) == 0) ) then
      call drive_population(M, natom, Nuc, Iz, rho1, overlap, smat, 2)
   endif

   return
end subroutine transport_population

subroutine mat_map(Nuc, M, natom)
   ! This subroutine classifies each index for the evolution of the density
   ! matrix during propagation.
   use transport_data, only: group, mapmat
   implicit none
   integer, intent(in)   :: M, natom, Nuc(M)
   integer               :: i, j, group_1, group_2, group_3

   mapmat = 0
   do i=1,M
   do j=1,M
      if ((group(nuc(i)).eq.1).and.(group(nuc(j)).eq.1)) mapmat(i,j)=1
      if ((group(nuc(i)).eq.1).and.(group(nuc(j)).eq.2)) mapmat(i,j)=2
      if ((group(nuc(i)).eq.1).and.(group(nuc(j)).eq.3)) mapmat(i,j)=3
      if ((group(nuc(i)).eq.2).and.(group(nuc(j)).eq.1)) mapmat(i,j)=4
      if ((group(nuc(i)).eq.2).and.(group(nuc(j)).eq.2)) mapmat(i,j)=5
      if ((group(nuc(i)).eq.2).and.(group(nuc(j)).eq.3)) mapmat(i,j)=6
      if ((group(nuc(i)).eq.3).and.(group(nuc(j)).eq.1)) mapmat(i,j)=7
      if ((group(nuc(i)).eq.3).and.(group(nuc(j)).eq.2)) mapmat(i,j)=8
      if ((group(nuc(i)).eq.3).and.(group(nuc(j)).eq.3)) mapmat(i,j)=9
   end do
   end do

   ! Counting the number of basis functions for each section of transport.
   group_1 = 0 ; group_2 = 0 ; group_3 = 0
   do i=1,M
      if (mapmat(i,i).eq.1) group_1 = group_1 + 1
      if (mapmat(i,i).eq.5) group_2 = group_2 + 1
      if (mapmat(i,i).eq.9) group_3 = group_3 + 1
   end do

    write(*,*) 'Transport - Basis functions from group 1 =', group_1
    write(*,*) 'Transport - Basis functions from group 2 =', group_2
    write(*,*) 'Transport - Basis functions from group 3 =', group_3

    return
end subroutine mat_map

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine electrostat(rho1, overlap, Gamma0, M)
   ! This subroutine modifies the density matrix according to the group
   ! containing the basis function indexes.
   use transport_data, only: mapmat, rhofirst
   implicit none
   integer, intent(in) :: M
   real*8,  intent(in) :: overlap(M,M), Gamma0
   integer :: i, j
   real*8  :: GammaIny, GammaAbs

#ifdef TD_SIMPLE
   complex*8 , intent(inout) :: rho1(M,M)
   complex*8 , allocatable   :: rho_scratch(:,:,:)
#else
   complex*16, intent(inout) :: rho1(M,M)
   complex*16, allocatable   :: rho_scratch(:,:,:)
#endif

   call g2g_timer_start('electrostat')
   allocate(rho_scratch(M,M,2))

   rho_scratch = 0.0D0
   do i = 1, M
   do j = 1, M
      select case (mapmat(i,j))
         case (0, 9)
            rho_scratch(i,j,1) = dcmplx(0.0D0,0.0D0)
            rho_scratch(i,j,2) = dcmplx(0.0D0,0.0D0)
         case (1, 2, 4, 5)
            rho_scratch(i,j,1) = rho1(i,j)
            rho_scratch(i,j,2) = rhofirst(i,j)
         case (3, 6, 7, 8)
            rho_scratch(i,j,1) = 0.50D0*rho1(i,j)
            rho_scratch(i,j,2) = 0.50D0*rhofirst(i,j)
      end select
   enddo
   enddo

   GammaIny = Gamma0*0.5D0
   GammaAbs = GammaIny
   write(*,*) 'Transport - GammaAbs, GammaIny =', GammaAbs, GammaIny

   do i = 1, M
   do j = 1, M
      rho1(i,j)= GammaAbs*rho_scratch(i,j,1) - GammaIny*rho_scratch(i,j,2)
      ! Checks NaNs.
      if (rho1(i,j).ne.rho1(i,j)) then
         stop 'Houston, we have a problem - NaN found in Rho1.'
      end if
   enddo
   enddo

   deallocate(rho_scratch)
   call g2g_timer_stop('electrostat')

   return
end subroutine electrostat

subroutine drive_population(M, natom, Nuc, Iz, rho1, overlap, smat, dvopt)
   use transport_data, only: pop_uid, drive_uid, ngroup, group, pop_drive
   implicit none
   integer, intent(in) :: M, natom, Iz(natom), Nuc(natom)
   integer, intent(in) :: dvopt
   real*8 , intent(in) :: overlap(M,M), smat(M,M)
#ifdef TD_SIMPLE
      complex*8 , intent(in) :: rho1(M,M)
#else
      complex*16, intent(in) :: rho1(M,M)
#endif
   real*8  :: qgr(ngroup), traza, q(natom), rho(M,M)
   integer :: i

   qgr(:) = 0.0D0
   traza  = 0.0D0

   if (dvopt==1) then
      q(:)=0.0d0
   else if (dvopt==2) then
      do i = 1, natom
         q(i) = Iz(i)
      enddo
   end if

   rho = real(rho1)
   select case (pop_drive)
      case (1)
         call mulliken_calc(natom, M, rho, overlap, Nuc, q)
      case (2)
         call lowdin_calc(M, natom, rho, smat, Nuc, q)
      case default
         write(*,*) "ERROR - Transport: Wrong value for Pop (drive_population)."
         stop
   end select

   do i = 1, natom
      qgr(group(i)) = qgr(group(i)) + q(i)
   enddo

   if(dvopt==1) then
      do i = 1, ngroup
         write(drive_uid,*) i, i, qgr(i)
         traza = traza + qgr(i)
      enddo
      write(*,*) "Total trace =", traza
   else if (dvopt==2) then
      do i = 1, ngroup
         write(pop_uid,*) i, i, qgr(i)
         traza = traza + qgr(i)
      enddo
      write(pop_uid,*) "Total trace =", traza
   end if

   return
end subroutine drive_population

end module transport_subs
