module DOS_data
   implicit none
   logical               :: dos_calc     = .false.
   logical               :: pdos_calc    = .false.
   logical               :: pdos_allb    = .false.
   integer               :: min_level    = 1
   integer               :: dos_nsteps   = 200000
   real*8                :: dos_sigma    = 0.0004d0
   real*8                :: dos_Eref     = 0.0d0 !-0.180187401161995 !Ag
   integer               :: pdos_nbases  = 0            !Option for just 1 atom
   integer               :: pdos_natoms  = 0
   integer, allocatable  :: pdos_nuc(:)
   integer, allocatable  :: pdos_base(:)
   real*8,  allocatable  :: pdos(:)
   real*8,  allocatable  :: pdos_b(:,:)


end module DOS_data
