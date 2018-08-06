module transport_data
   implicit none
   logical   :: transport_calc   = .false.   !Active transport options
   logical   :: generate_rho0    = .false.   !Option to
   logical   :: gate_field       = .false.
   integer   :: save_charge_freq = 1
   integer   :: pop_drive        = 1
   integer   :: nbias            = 0         !Number of electrodes present
   integer   :: pop_uid          = 678
   integer   :: drive_uid        = 5555
   real*8    :: driving_rate     = 0.001
   real*8    :: GammaMagnus      = 0.0D0
   real*8    :: GammaVerlet      = 0.0D0
   real*8    :: re_traza         = 0.0D0
!charly:
   integer, allocatable :: timestep_init(:)


   integer   , allocatable :: mapmat(:,:), group(:)
#ifdef TD_SIMPLE
   complex*8 , allocatable :: rhofirst(:,:,:)
#else
   complex*16, allocatable :: rhofirst(:,:,:)
#endif
end module transport_data

module transport_subs
   implicit none
contains

subroutine transport_init(M, dim3, natom, Nuc, RMM5, overlap, rho, OPEN)
   use transport_data, only: GammaMagnus, GammaVerlet, nbias, group, rhofirst,&
                             driving_rate, pop_drive, pop_uid, drive_uid,     &
                             mapmat, timestep_init
   implicit none
   logical, intent(in)  :: OPEN
   integer, intent(in)  :: M, natom, Nuc(M), dim3
   real*8 , intent(in)  :: RMM5(M*(M+1)/2)
   real*8 , allocatable, intent(inout) :: overlap(:,:)
#ifdef TD_SIMPLE
   complex*8 , allocatable, intent(inout) :: rho(:,:,:)
#else
   complex*16, allocatable, intent(inout) :: rho(:,:,:)
#endif
   integer :: orb_group(M), icount

   allocate(rhofirst(M,M,dim3),timestep_init(nbias))

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
         open( unit = pop_uid  , file = "Mullikenbias")
         open( unit = drive_uid, file = 'DriveMul')
      case (2)
         open( unit = pop_uid  , file = "Lowdinbias")
         open( unit = drive_uid, file = 'DriveLowd')
      case default
         write(*,*) 'ERROR - Transport: Wrong value for pop_drive '//&
                    '(transport_init)'
   end select

   if (allocated(overlap)) deallocate(overlap)
   allocate(overlap(M,M))
   call spunpack('L', M, RMM5, overlap)
   call transport_generate_rho(M, rho, OPEN)

   return
end subroutine transport_init


subroutine transport_read_groups(natom)
!This subroutine read the transport.in file, the file must have in the first
!line the time step where we want to start the dinamic in each electrode:
!   timestep_init(1), timestep_init(2), .... , timestep_init(n)
!At least two electrodes with oposite polarisation must have the same
!timestep_init.
!Then, there must be a list in the same order than the atoms in the xyz file,
!where it's assign the group which each atom belong:
!     0 - Correspond to the device
!   i>0 - Correspond to the i electrode
   use transport_data, only: group, nbias, timestep_init
   implicit none
   integer, intent(in)  :: natom
   integer :: icount, group_read
   logical :: file_exists

   inquire(file = 'transport.in', exist = file_exists)
   if (.not.file_exists) then
      write(*,*) ' ERROR - Transport: Cannot find transport.in file.'
      stop
   endif
   open( unit = 678 , file = 'transport.in')

   read(678,*) timestep_init
   do icount = 1 , natom
      read(678, *) group_read
      group(icount) = group_read
      if (group_read.gt.nbias) then
         write(*,*) "ERROR - A value greater than nbias declared has been found&
                     in transport.in file"
         stop
      end if
   enddo
   close(unit=678)

   return
end subroutine transport_read_groups

subroutine transport_generate_rho(M, rho, OPEN)
   use transport_data, only: generate_rho0, rhofirst
   implicit none
   logical, intent(in)  :: OPEN
   integer, intent(in)  :: M
#ifdef TD_SIMPLE
   complex*8 , allocatable, intent(inout) :: rho(:,:,:)
#else
   complex*16, allocatable, intent(inout) :: rho(:,:,:)
#endif
   integer   :: icount, jcount

   open( unit = 100000, file = 'rhofirst')

   if (OPEN) then
      if (generate_rho0) then
         do icount = 1, M
         do jcount = 1, M
            write(100000,*) rho(icount, jcount,1),rho(icount, jcount,2)
         enddo
         enddo
         rhofirst = rho
         write(*,*) 'RhoFirst has been written.'
      else
         do icount = 1, M
         do jcount = 1, M
            read(100000,*) rhofirst(icount, jcount,1),rhofirst(icount, jcount,2)
         enddo
         enddo
         write(*,*) 'RhoFirst has been read.'
      endif
   else
      if (generate_rho0) then
         do icount = 1, M
         do jcount = 1, M
            write(100000,*) rho(icount, jcount,1)
         enddo
         enddo
         rhofirst = rho
         write(*,*) 'RhoFirst has been written.'
      else
         do icount = 1, M
         do jcount = 1, M
            read(100000,*) rhofirst(icount, jcount,1)
         enddo
         enddo
         write(*,*) 'RhoFirst has been read.'
      endif
   end if

   return
end subroutine transport_generate_rho

subroutine transport_rho_trace(M, dim3, rho)
   implicit none
   integer, intent(in) :: M, dim3
#ifdef TD_SIMPLE
   complex*8 , intent(in) :: rho(M,M,dim3)
#else
   complex*16, intent(in) :: rho(M,M,dim3)
#endif
   integer   :: icount
   complex*8 :: traza

   traza = dcmplx(0.0D0, 0.0D0)
   do icount = 1, M
      traza = traza + rho(icount, icount,1)
   enddo

   if (dim3==2) then
      write(*,*) 'Transport - Alpha Rho Trace = ', real(traza)
      traza = dcmplx(0.0D0, 0.0D0)
      do icount = 1, M
         traza = traza + rho(icount, icount,2)
      enddo
      write(*,*) 'Transport - Beta Rho Trace = ', real(traza)
   else
      write(*,*) 'Transport - Rho Trace = ', real(traza)
   end if

   return
end subroutine transport_rho_trace

subroutine transport_propagate(M, dim3, natom, Nuc, Iz, propagator, istep, &
                               overlap, sqsm, rho1, Ymat ,OPEN)
   use transport_data, only: save_charge_freq, pop_drive, GammaMagnus, &
                             GammaVerlet, nbias, group, timestep_init, nbias
   use mathsubs,       only: basechange_gemm
   implicit none
   logical   , intent(in)    :: OPEN
   integer   , intent(in)    :: dim3
   integer   , intent(in)    :: M, natom, Nuc(M), Iz(natom), propagator, istep
   real*8    , intent(in)    :: Ymat(M,M)
   real*8    , intent(inout) :: overlap(M,M), sqsm(M,M)
#ifdef TD_SIMPLE
   complex*8 , intent(inout) :: rho1(M,M,dim3)
#else
   complex*16, intent(inout) :: rho1(M,M,dim3)
#endif
   real*8  :: gamma
   real*8  :: scratchgamma(nbias)
   integer :: save_freq
   integer :: ii

   save_freq = save_charge_freq
   select case (propagator)
   case (1)
      save_freq = save_freq*10
      gamma     = GammaVerlet
   case (2)
      gamma     = GammaMagnus
   end select

   call g2g_timer_start('Transport-propagation')

!charly: for the momment we can made transport just with magnus
   do ii=1,nbias
      if ((propagator.eq.2) .and. (istep.gt.(timestep_init(ii)+999))) Then
         scratchgamma(ii) = gamma
      else if ((istep>=timestep_init(ii)).and.                                 &
               (istep<=(timestep_init(ii)+999)))then
         scratchgamma(ii) = gamma * exp(-0.0001D0 *                            &
                            (dble(istep-(1000+timestep_init(ii)))**2))
      else if (istep <= timestep_init(ii)) then
         scratchgamma(ii) = 0.0d0
      endif
   end do


   call electrostat(rho1(:,:,1), overlap, scratchgamma, M,Nuc, 1)
   if (OPEN) call electrostat(rho1(:,:,2), overlap, scratchgamma, M, Nuc, 2)
   if (mod( istep-1 , save_freq) == 0) then
      call drive_population(M, dim3, natom, Nuc, Iz, rho1, overlap, sqsm, 1,   &
                            OPEN)
   endif

   rho1(:,:,1)=basechange_gemm(M,rho1(:,:,1),Ymat)
   if (OPEN) rho1(:,:,2) = basechange_gemm(M,rho1(:,:,2),Ymat)

   call g2g_timer_stop('Transport-propagation')
   return
end subroutine transport_propagate

#ifdef CUBLAS
subroutine transport_propagate_cu(M, dim3, natom, Nuc, Iz, propagator, istep, &
                                  overlap, sqsm, rho1, devPtrY, OPEN)
   use cublasmath    , only: basechange_cublas
   use transport_data, only: nbias, group, save_charge_freq, GammaMagnus, &
                             GammaVerlet, timestep_init
   implicit none
   logical, intent(in)       :: OPEN
   integer, intent(in)       :: dim3
   integer   , intent(in)    :: M, natom, Nuc(M), Iz(natom), propagator, istep
   integer*8 , intent(in)    :: devPtrY
   real*8    , intent(inout) :: overlap(M,M), sqsm(M,M)
#ifdef TD_SIMPLE
   complex*8 , intent(inout) :: rho1(M,M,dim3)
#else
   complex*16, intent(inout) :: rho1(M,M,dim3)
#endif
   real*8  :: gamma
   real*8  :: scratchgamma(nbias)
   integer :: save_freq
   integer :: ii

   save_freq = save_charge_freq
   select case (propagator)
   case (1)
      save_freq = save_freq*10
      gamma     = GammaVerlet
   case (2)
      gamma     = GammaMagnus
   end select

   call g2g_timer_start('Transport-propagation')

!charly: for the momment we can made transport just with magnus
   do ii=1,nbias
      if ((propagator.eq.2) .and. (istep.gt.(timestep_init(ii)+999))) Then
         scratchgamma(ii) = gamma
      else if ((istep>=timestep_init(ii)).and.                                 &
               (istep<=(timestep_init(ii)+999)))then
         scratchgamma(ii) = gamma * exp(-0.0001D0 *                            &
                            (dble(istep-(1000+timestep_init(ii)))**2))
      else if (istep <= timestep_init(ii)) then
         scratchgamma(ii) = 0.0d0
      endif
   end do

   call electrostat(rho1(:,:,1), overlap, scratchgamma, M,Nuc, 1)
   if (OPEN) call electrostat(rho1(:,:,2), overlap, scratchgamma, M,Nuc, 2)

   if (mod( istep-1 , save_freq) == 0) then
      call drive_population(M, dim3, natom, Nuc, Iz, rho1, overlap, sqsm, 1,   &
                            OPEN)
   endif

   call g2g_timer_start('complex_rho_ao_to_on-cu')
   rho1(:,:,1) = basechange_cublas(M, rho1(:,:,1), devPtrY, 'dir')
   if (OPEN) rho1(:,:,2) = basechange_cublas(M, rho1(:,:,2), devPtrY, 'dir')
   call g2g_timer_stop('complex_rho_ao_to_on-cu')

   call g2g_timer_stop('Transport-propagation')

   return
end subroutine transport_propagate_cu
#endif

subroutine transport_population(M, dim3, natom, Nuc, Iz, rho1, overlap, smat, &
                                propagator, is_lpfrg, istep, OPEN)
   use transport_data, only: save_charge_freq
   implicit none
   logical, intent(in) :: OPEN
   integer, intent(in) :: dim3
   integer, intent(in) :: M, natom, Iz(natom), Nuc(natom), istep, propagator
   logical, intent(in) :: is_lpfrg
   real*8 , intent(in) :: overlap(M,M), smat(M,M)
#ifdef TD_SIMPLE
      complex*8 , intent(in) :: rho1(M,M,dim3)
#else
      complex*16, intent(in) :: rho1(M,M,dim3)
#endif

   if ( ((propagator.gt.1) .and. (is_lpfrg) .and.       &
      (mod((istep-1), save_charge_freq*10) == 0)) .or.  &
      (mod((istep-1), save_charge_freq) == 0) ) then
      call drive_population(M, dim3, natom, Nuc, Iz, rho1, overlap, smat, 2,   &
                            OPEN)
   endif

   return
end subroutine transport_population


subroutine mat_map(Nuc, M, natom)
! This subroutine classify each element of the matrix density in three cases
! whichs correspond to the three cases present in the driving term of Drive
! Liouville von Neuman, this cases will be use in electrostat. The cases are:
!         1) if i and j  = 0
!         2) if i or j   = 0
!         3) if i and j /= 0
   use transport_data, only: group, nbias ,mapmat
   implicit none
   integer, intent(in)   :: M, natom, Nuc(M)
   integer               :: i, j, group_n(nbias+1)

   mapmat = 0
   do i=1,M
   do j=1,M
      if ((group(nuc(i)).eq.0).and.(group(nuc(j)).eq.0)) mapmat(i,j)=1
      if ((group(nuc(i)).eq.0).and.(group(nuc(j)).ne.0)) mapmat(i,j)=2
      if ((group(nuc(i)).ne.0).and.(group(nuc(j)).eq.0)) mapmat(i,j)=2
      if ((group(nuc(i)).ne.0).and.(group(nuc(j)).ne.0)) mapmat(i,j)=3
   end do
   end do

   ! Counting the number of basis functions for each section of transport.
   group_n = 0
   do i=1,M
      group_n(group(nuc(i))+1) = group_n(group(nuc(i))+1)+1
   end do

   do i=1,nbias+1
      write(*,*) "Transport - Basis functions from group",i-1,"=", group_n(i)
   end do

    return
end subroutine mat_map

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine electrostat(rho1, overlap, Gamma0, M,Nuc, spin)
   ! This subroutine modifies the density matrix according to the group
   ! containing the basis function indexes.
   use transport_data, only: mapmat, rhofirst, nbias, group
   implicit none
   integer, intent(in) :: M, spin
   integer, intent(in) :: Nuc(M)
   real*8,  intent(in) :: overlap(M,M), Gamma0(nbias)
   integer :: i, j
   real*8  :: GammaIny, GammaAbs, tempgamma, tempgamma2

#ifdef TD_SIMPLE
   complex*8 , intent(inout) :: rho1(M,M)
   complex*8 , allocatable   :: rho_scratch(:,:,:)
#else
   complex*16, intent(inout) :: rho1(M,M)
   complex*16, allocatable   :: rho_scratch(:,:,:)
#endif

   call g2g_timer_start('electrostat')
   allocate(rho_scratch(M,M,2))

   tempgamma   = 0.0D0
   tempgamma2  = 0.0D0
   rho_scratch = 0.0D0

   do i=1,M
   do j=1,M

      if (group(nuc(i))/=0.and.group(nuc(j))==0) then
         tempgamma = Gamma0(group(nuc(i)))
      else if (group(nuc(i))==0.and.group(nuc(j))/=0) then
         tempgamma = Gamma0(group(nuc(j)))
      else if (group(nuc(i))/=0.and.group(nuc(j))/=0) then
         tempgamma = min(Gamma0(group(nuc(i))),                                &
                         Gamma0(group(nuc(j))))
         tempgamma2 = max(Gamma0(group(nuc(i))),                               &
                          Gamma0(group(nuc(j))))
      end if

      select case(mapmat(i,j))
         case (1)
            rho_scratch(i,j,1) = dcmplx(0.0D0,0.0D0)
            rho_scratch(i,j,2) = dcmplx(0.0D0,0.0D0)
         case (2)
            rho_scratch(i,j,1) = tempgamma*0.50D0*rho1(i,j)
            rho_scratch(i,j,2) = tempgamma*0.50D0*rhofirst(i,j,spin)
         case (3)
            rho_scratch(i,j,1) = tempgamma*0.50D0*rho1(i,j) +                  &
                                 tempgamma2*0.50D0*rho1(i,j)
            rho_scratch(i,j,2) = tempgamma*0.50D0*rhofirst(i,j,spin) +         &
                                 tempgamma2*0.50D0*rhofirst(i,j,spin)
      end select
   end do
   end do

   GammaIny = 0.5D0
   GammaAbs = GammaIny
   write(*,*) 'Transport - Gammas =', (Gamma0(i),i=1,nbias)

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

subroutine drive_population(M, dim3, natom, Nuc, Iz, rho1, overlap, smat,      &
                            dvopt, OPEN)
   use transport_data, only: pop_uid, drive_uid, nbias, group, pop_drive
   implicit none
   logical, intent(in) :: OPEN
   integer, intent(in) :: dim3
   integer, intent(in) :: M, natom, Iz(natom), Nuc(natom)
   integer, intent(in) :: dvopt
   real*8 , intent(in) :: overlap(M,M), smat(M,M)
#ifdef TD_SIMPLE
      complex*8 , intent(in) :: rho1(M,M, dim3)
#else
      complex*16, intent(in) :: rho1(M,M, dim3)
#endif
   real*8  :: qgr(nbias+1), traza, q(natom), rho(M,M,dim3)
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
         call mulliken_calc(natom, M, rho(:,:,1), overlap, Nuc, q)
         if (OPEN) call mulliken_calc(natom, M, rho(:,:,2), overlap, Nuc,      &
                                      q)
      case (2)
         call lowdin_calc(M, natom, rho(:,:,1), smat, Nuc, q)
         if (OPEN) call lowdin_calc(M, natom, rho(:,:,2), smat, Nuc, q)
      case default
         write(*,*) "ERROR - Transport: Wrong value for Pop (drive_population)."
         stop
   end select

   do i = 1, natom
      qgr(group(i)+1) = qgr(group(i)+1) + q(i)
   enddo

   if(dvopt==1) then
      do i = 1, nbias+1
         write(drive_uid,*) i-1, i-1, qgr(i)
         traza = traza + qgr(i)
      enddo
      write(*,*) "Total trace =", traza
   else if (dvopt==2) then
      do i = 1, nbias+1
         write(pop_uid,*) i-1, i-1, qgr(i)
         traza = traza + qgr(i)
      enddo
      write(pop_uid,*) "Total trace =", traza
   end if

   return
end subroutine drive_population

end module transport_subs
