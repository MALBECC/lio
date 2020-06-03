subroutine fstsh_init( dt )
use fstsh_data  , only: tsh_file, all_states, type_coupling, Sovl_old, Sovl_now, &
                        a_old, c_old, r_old, C_scf_old, tsh_nucStep, tsh_Enstep, &
                        current_state, phases_old, vel_old, first_interp,        &
                        Nesup_now, Nesup_old, tsh_time_dt, tsh_minprob
use excited_data, only: root, nstates
use basis_data  , only: M,  max_c_per_atom
use garcha_mod  , only: natom, ntatom
   implicit none

   LIODBLE, intent(in) :: dt
   ! random variables
   integer :: random_size
   integer, dimension(12) :: random_values
   integer, dimension(:), allocatable :: seed

   ! Set nuclear time step
   ! dt = [femto] -> tsh_time_dt = [a.u.]
   tsh_time_dt = dt * 41.341374575751d0 

   ! Open file and print header
   open(unit=tsh_file,file="tsh_lio.log")
   write(tsh_file,*) "=============================="
   write(tsh_file,*) " TRAJECTORY SURFFACE HOPPING"
   write(tsh_file,*) "=============================="

   ! Set and Allocated Variables
   all_states    = nstates + 1
   current_state = root + 1
   allocate(Sovl_old(M,M),Sovl_now(M,M))
   allocate(a_old(M,max_c_per_atom),c_old(M,max_c_per_atom))
   allocate(r_old(ntatom,3)); r_old = 0.0d0
   allocate(C_scf_old(M,M)); C_scf_old = 0.0d0
   allocate(phases_old(all_states)); phases_old = 1.0d0
   allocate(vel_old(3,natom)); vel_old = 0.0d0
   allocate(Nesup_now(all_states),Nesup_old(all_states))
   Nesup_now = 0.0d0; Nesup_old = 0.0d0
   tsh_nucStep = 0; first_interp = .true.
   Sovl_old = 0.0d0; Sovl_now = 0.0d0; a_old = 0.0d0; c_old = 0.0d0;

   ! Initialization of random values
   call date_and_time(VALUES=random_values)
   call random_seed(size=random_size)
   allocate(seed(random_size))
   seed = random_values
   call random_seed(put=seed)

   ! Print Information
   write(tsh_file,*) " "
   write(tsh_file,*) "Information:"
   write(tsh_file,"(1X,A,I1,A)") "All states= ", all_states, " (1=GS)"
   write(tsh_file,"(1X,A,I1,A)") "Current state= ", current_state
   write(tsh_file,"(1X,A,I1)"  ) "Type of Coupling= ", type_coupling
   write(tsh_file,"(1X,A,F8.4)") "Min. Probability= ", tsh_minprob
   write(tsh_file,"(1X,A,I2,A)") "Electronic Steps= ", tsh_Enstep
   write(tsh_file,"(1X,A,F10.5,A)") "Nuclear Time-Step= ", tsh_time_dt * 0.02418884254d0, " fs."
   write(tsh_file,"(1X,A,F10.5,A)") "Electronic Time-Step= ", (tsh_time_dt/real(tsh_Enstep)) * 0.02418884254d0, " fs."
   write(tsh_file,*) "SEED= ", seed
   deallocate(seed)

end subroutine fstsh_init

subroutine tsh_init(dt)
   use excited_data,only: TSH, tsh_time_dt, tsh_coef, tsh_Jstate, &
                          tsh_Kstate, gamma_old, excited_forces
   use garcha_mod  , only: natom

   LIODBLE, intent(in) :: dt

   ! Random variables
   integer :: random_size
   integer, dimension(12) :: random_values
   integer, dimension(:), allocatable :: seed

   if ( TSH ) then
      ! dt_i = ps
      ! 1 ps = 4.134137d4 au
      ! tsh_time_dt = au
      tsh_time_dt = dt * (41341.3733366d0) ! tsh_time_dt in atomic units
    
      print*, "Init TSH Dynamics"
      ! RANDOM SEED
      call date_and_time(VALUES=random_values)
      call random_seed(size=random_size)
      allocate(seed(random_size))
      seed = random_values
      print*, "SEED:", seed
      call random_seed(put=seed)
      deallocate(seed)
      
      allocate(tsh_coef(2))
      tsh_coef(1) = (0.0d0,0.0d0)
      tsh_coef(2) = (1.0d0,0.0d0)
      tsh_Jstate  = 2
      tsh_Kstate  = 1
      allocate(gamma_old(natom,3))
      gamma_old = 0.0d0
      
      excited_forces = .true.
   endif   

end subroutine tsh_init
