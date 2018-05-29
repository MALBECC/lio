!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module ehrendata
!------------------------------------------------------------------------------!
   implicit none
   integer :: nustep_count = 0
   integer :: elstep_count = 0

   integer :: rsti_funit  = 654321
   integer :: rsto_funit  = 123456

   real*8  :: stored_time   = 0.0d0
   real*8  :: stored_energy = 0.0d0
   real*8  :: stored_dipmom(3) = 0.0d0

   complex*16, allocatable :: stored_densM1(:,:)
   complex*16, allocatable :: stored_densM2(:,:)

   !  Run type information
      integer :: ndyn_steps = 0  ! Number of total nuclear dynamic steps
      integer :: edyn_steps = 0  ! Number of total electronic dynamic steps PER
                                 ! nuclear dynamic step.
      logical :: nullify_forces = .false.
   !  ndyn == 0 & edyn == 0   =>   Single point
   !  ndyn /= 0 & edyn == 0   =>   BO atomistic dynamic (aux mm)
   !  ndyn == 0 & edyn /= 0   =>   TD electron dynamic
   !  ndyn /= 0 & edyn /= 0   =>   Ehrenfest dynamic (aux mm)

   !  Output information
      integer           :: wdip_nfreq = 0
      character(len=80) :: wdip_fname = "liorun.dip"

   !  Restarting information
   !  If (rsti_loads), the program will load the restart from (rsti_fname).
   !  If (rsto_saves), the program will save the restart information overwriting
   !  the file (rsto_fname) every (rsto_nfreq) steps and in the last step.
      logical           :: rsti_loads = .false.
      character(len=80) :: rsti_fname = "liorun.rsti"
      logical           :: rsto_saves = .false.
      integer           :: rsto_nfreq = 0
      character(len=80) :: rsto_fname = "liorun.rsto"

   !  External Electrical Field
   !     If (eefld_on), an external field will be applied to the system. The
   !  amplitude in each direction is given by the (eefld_amp) variables. It
   !  can have an oscilating time modulation of a specific (eefld_wavelen) and
   !  also a gaussian envelop centered in (eefld_timepos), with width given by
   !  (eefld_timeamp). Both (eefld_timegih) and (eefld_timegfh) must be true for
   !  a full gaussian, activating the modulation before and after the center
   !  respectively.
      logical :: eefld_on   = .false.
      real*8  :: eefld_ampx = 0.0d0 ! in au
      real*8  :: eefld_ampy = 0.0d0 ! in au
      real*8  :: eefld_ampz = 0.0d0 ! in au
      logical :: eefld_timegih = .false. ! time gaussian initial half
      logical :: eefld_timegfh = .false. ! time gaussian final half
      real*8  :: eefld_timepos =  1.0d0  ! in ps (currently fs!)
      real*8  :: eefld_timeamp =  0.2d0  ! in ps (currently fs!)
      real*8  :: eefld_wavelen =  0.0d0  ! in nm

end module ehrendata
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
