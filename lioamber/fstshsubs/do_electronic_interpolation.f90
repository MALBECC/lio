subroutine do_electronic_interpolation(vel_now,doTSH)
! VARIABLES
! sigma_old = sigma at t-3/2 dt
! sigma_now = sigma at t-1/2 dt
! sigma_0   = sigma at t-dt
! sigma_1   = sigma at t
use garcha_mod, only: natom
use fstsh_data, only: tsh_file, first_interp, vel_old, sigma_old, sigma_now,    &
                      sigma_0, sigma_1, all_states, tsh_nucStep, current_state, &
                      tsh_Enstep, coef_Stat, elec_Coup, elec_Pha, dot_stat,     &
                      Nesup_now, Nesup_old, elec_Ene, after_hopp, tsh_time_dt
   implicit none

   LIODBLE, intent(inout) :: vel_now(3,natom)
   logical, intent(inout) :: doTSH

   integer :: method_sigma, state_before, iter, old_surf, new_surf
   LIODBLE :: dt_elec, tot_time
   LIODBLE, allocatable :: elec_vel(:,:)

   doTSH = .false.

   if ( first_interp ) then
      write(tsh_file,*) "We only save variables. Do not perform electronic interpolation"

      first_interp = .false.
      vel_old = vel_now
      sigma_old = sigma_now; sigma_0 = 0.0d0; sigma_1 = 0.0d0
  
      ! Variables initialization needed at electronic sub-time step
      call interpolation_init( )

      ! Obtain Cdot at actual time
      !   this routine only needs variables at actual electronic time step
      call dot_calculation(coef_Stat, elec_Coup, elec_Pha, dot_Stat, &
                           all_states)

      ! Save electronic variables
      coef_Stat(3,:) = coef_Stat(2,:); coef_Stat(2,:) = coef_Stat(1,:)
      elec_Coup(3,:,:) = elec_Coup(2,:,:); elec_Coup(2,:,:) = elec_Coup(1,:,:)
      elec_Pha(3,:,:) = elec_Pha(2,:,:); elec_Pha(2,:,:) = elec_Pha(1,:,:)
      dot_stat(3,:) = dot_stat(2,:); dot_stat(2,:) = dot_stat(1,:)

      ! Step Nuclear Actualization
      tsh_nucStep = tsh_nucStep + 1
      return
   endif

   method_sigma = 1 ! Normal sigma extrapolation 
   if ( tsh_nucStep == 1 .or. after_hopp ) then
       after_hopp = .false.
       method_sigma = 2 ! Raw sigma extrapolation
   endif

   ! We obtain sigma at actual nuclear time step
   call extrap_sigmaT(sigma_old,sigma_now,sigma_1,all_states,method_sigma)

   write(tsh_file,*) "Starting Electronic Interpolation"
!  call print_sigma(sigma_old,all_states,"t-1.5dt")
!  call print_sigma(sigma_now,all_states,"t-0.5dt")
!  call print_sigma(sigma_0,all_states,  "t-dt   ")
!  call print_sigma(sigma_1,all_states,  "t      ")

   state_before = current_state
   allocate(elec_vel(3,natom))
   dt_elec = tsh_time_dt / real(tsh_Enstep)

   ! Electronic Interpolation
   do iter = 0, tsh_Enstep-1
      tot_time = (dt_elec * iter + tsh_nucStep * tsh_time_dt) * 0.02418884254d0 ! change a.u to femto
      write(tsh_file,"(3X,A,I2,A,F8.4,A)") "Electronic Sub-step= ",iter," Time= ", tot_time, " fs."

      ! Energy, Coupling and Velocities Interpolation
      call interpol_EneCoupVel(Nesup_now,Nesup_old,sigma_1,sigma_0, elec_Ene, &
                               elec_Coup,all_states,vel_now,vel_old,elec_vel, &
                               natom,tsh_time_dt,dt_elec,iter)

      ! Phase Calculation
      call phase_calculation(elec_Ene,elec_Pha,dt_elec,all_states,tsh_nucStep)

      ! Coefficients Evolution
      call coef_evolution(coef_Stat,elec_Coup,elec_Pha,dot_Stat,all_states,dt_elec,tsh_nucStep)

      ! Probabilities of Hopp Calculates
      old_surf = current_state
      call do_probabilities(coef_Stat,elec_Coup,elec_Pha,elec_Ene,all_states,dt_elec,tsh_nucStep)
      new_surf = current_state

      ! Hopping ?
      if ( old_surf /= new_surf ) then
         write(tsh_file,"(4X,A,I2)") "HOPPING in substep: ", iter
      endif
      
      ! Decoherence Term
      call correct_decoherence(coef_Stat,elec_Ene,elec_vel,natom,all_states,dt_elec)

      ! Save Variables
      coef_Stat(3,:) = coef_Stat(2,:); coef_Stat(2,:) = coef_Stat(1,:)
      elec_Ene(3,:) = elec_Ene(2,:); elec_Ene(2,:) = elec_Ene(1,:)
      dot_stat(3,:) = dot_stat(2,:); dot_stat(2,:) = dot_stat(1,:)
      elec_Coup(3,:,:) = elec_Coup(2,:,:); elec_Coup(2,:,:) = elec_Coup(1,:,:)
      elec_Pha(3,:,:) = elec_Pha(2,:,:); elec_Pha(2,:,:) = elec_Pha(1,:,:)

      write(tsh_file,*) " "

   enddo ! END ELEC. INTERP.
   deallocate(elec_vel)

   if ( new_surf /= state_before ) then
      write(tsh_file,"(1X,A,I2,A,I2)") "HOPP Surfaces", state_before, "->", new_surf
      write(tsh_file,"(1X,A)") "Electronic Hopp TRUE"
      call adjust_vel(state_before,new_surf,Nesup_now,vel_now,all_states)
      print*, "POSSIBLE HOPP or FRUSTATED HOPP", state_before, "->", new_surf
      doTSH = .true.
   endif

   ! Sigma Actualization
   sigma_0 = sigma_1; sigma_old = sigma_now

   ! Velocities Actualization
   vel_old = vel_now

   ! Potential Energy Actualization 
   Nesup_old = Nesup_now

   ! Step Nuclear Actualization
   tsh_nucStep = tsh_nucStep + 1

end subroutine do_electronic_interpolation

subroutine extrap_sigmaT(sgmO,sgmN,sgm1,N,method)
use fstsh_data, only: tsh_file
   implicit none

   integer, intent(in)  :: N, method
   LIODBLE, intent(in)  :: sgmO(N,N), sgmN(N,N)
   LIODBLE, intent(out) :: sgm1(N,N)

   integer :: ii, jj

   if ( method == 1 ) then
      write(tsh_file,*) "Sigma Extrapolation"
      sgm1 = 0.5d0 * ( 3.0d0*sgmN - sgmO )
   elseif ( method == 2 ) then
      write(tsh_file,*) "Bad Approximation in sigma extrapolation"
      sgm1 = sgmN
   else 
      write(tsh_file,*) "Bad method in order to sigma extrapolation"
      stop
   endif
end subroutine extrap_sigmaT

subroutine interpolation_init( )
use fstsh_data, only: coef_Stat, dot_Stat, elec_Coup, elec_Ene, elec_Pha, &
                      current_state, all_states, tsh_file
implicit none

   write(tsh_file,*) "Allocatiion of electronic sub-time step variables"

   ! Allocate Variables
   allocate(coef_Stat(3,all_states),dot_Stat(3,all_states))
   allocate(elec_Coup(3,all_states,all_states),elec_Ene(3,all_states))
   allocate(elec_Pha(3,all_states,all_states))

   ! Set initial values
   coef_Stat = cmplx(0.0d0,0.0d0,8)
   coef_Stat(1,current_state) = cmplx(1.0d0,0.0d0,8)
   dot_Stat  = cmplx(0.0d0,0.0d0,8)
   elec_Pha  = 0.0d0; elec_Coup = 0.0d0; elec_Ene = 0.0d0
end subroutine interpolation_init

subroutine dot_calculation(coef, coup, pha, dot, Nsup)
! This routine calculate coeff derivative at actual electronic time
! This only need variables at actual time
   implicit none

   integer, intent(in)     :: Nsup
   TDCOMPLEX, intent(in)   :: coef(3,Nsup)
   LIODBLE, intent(in)     :: coup(3,Nsup,Nsup), pha(3,Nsup,Nsup)
   TDCOMPLEX, intent(inout):: dot(3,Nsup)

   integer    :: kk, ll
   TDCOMPLEX :: res, temp

   temp = cmplx(0.0d0,0.0d0,8)

   do kk=1,Nsup
     res  = cmplx(0.0d0,0.0d0,8)
     do ll=1,Nsup
        temp = exp( cmplx(0.0d0,pha(1,kk,ll),8) )
        res  = res - coef(1,ll) * temp * coup(1,kk,ll)
     enddo
     dot(1,kk) = res
   enddo

end subroutine dot_calculation

subroutine interpol_EneCoupVel(NE,NE_old,NC,NC_old,ene,coup,Nsup,NV,&
                               NV_old,vel,nat,dtN,dtE,itp)
! This routine perform an interpolation of Energy and Coupling between
! two nuclear time steps (t and t-dt) in order to obtain Energy and Coupling at actual
! electronic step
   implicit none

   integer, intent(in) :: Nsup, itp, nat
   LIODBLE, intent(in)  :: dtN, dtE
   LIODBLE, intent(in)  :: NE(Nsup),NE_old(Nsup)
   LIODBLE, intent(in)  :: NC(Nsup,Nsup),NC_old(Nsup,Nsup)
   LIODBLE, intent(in)  :: NV(3,nat), NV_old(3,nat)
   LIODBLE, intent(inout) :: ene(3,Nsup), coup(3,Nsup,Nsup), vel(3,nat)

   ene(1,:) = NE_old(:) + ( NE(:)-NE_old(:) ) * (itp*dtE/dtN)
   coup(1,:,:) = NC_old(:,:) + ( NC(:,:)-NC_old(:,:) ) * (itp*dtE/dtN)
   vel(:,:) = NV_old(:,:) + ( NV(:,:)-NV_old(:,:) ) * (itp*dtE/dtN)
end subroutine interpol_EneCoupVel

subroutine phase_calculation(Ene,Pha,dt,Nsup,inStep)
! This subroutine propagates the phase of coef elec
   implicit none

   integer, intent(in)    :: Nsup, inStep
   LIODBLE,  intent(in)    :: Ene(3,Nsup), dt
   LIODBLE,  intent(inout) :: Pha(3,Nsup,Nsup)

   integer :: ii, jj
   LIODBLE :: dt1, dt2, dt3, temp1, temp2, temp3, tempf

   if ( inStep == 0 ) stop "Something wrong in pahse_calculateion, &
                            step is 0"

   if ( inStep == 1 ) then
      dt1 = dt * 0.5d0
      dt2 = dt * 0.5d0
      dt3 = 0.0d0
   else
      dt1 = dt * 5.0d0 / 12.0d0
      dt2 = dt * 2.0d0 / 3.0d0
      dt3 = dt / 12.0d0
   endif

   do ii=1, Nsup
   do jj=1, Nsup
      temp1 = Ene(1,ii) - Ene(1,jj)
      temp2 = Ene(2,ii) - Ene(2,jj)
      temp3 = Ene(3,ii) - Ene(3,jj)
      tempf = temp1*dt1 + temp2*dt2 - temp3*dt3
      Pha(1,ii,jj) = Pha(2,ii,jj) + tempf
   enddo
   enddo
   Pha(1,:,:) = (-1.0d0) * Pha(1,:,:)

end subroutine phase_calculation

subroutine coef_evolution(coef,coup,pha,dot,Nsup,dt,inStep)
!  This routine performs the electronic propagation
!  if we are in the first nuc time step and first ele time step, 
!  we use euler alghorithm but if the step is greater, we use butcher propagation
   implicit none

   integer,   intent(in)    :: Nsup, inStep
   LIODBLE,   intent(in)    :: pha(3,Nsup,Nsup), dt
   LIODBLE,   intent(inout)    :: coup(3,Nsup,Nsup) ! In reality, this value is not modified
   TDCOMPLEX, intent(inout) :: coef(3,Nsup), dot(3,Nsup)

   LIODBLE, dimension(:,:), allocatable :: mat
   TDCOMPLEX, dimension(:), allocatable :: vec, temp1, temp2, temp3

   if ( inStep == 0 ) stop "Something wrong in coef_evolution, &
                            inStep = 0"

   if ( inStep == 1 ) then
      ! EULER
      coef(1,:) = coef(2,:) + dot(2,:) * dt
      call dot_calculation(coef,coup,pha,dot,Nsup)
      coef(1,:) = coef(2,:) + ( dot(1,:)+dot(2,:) ) * dt * 0.5d0
      call dot_calculation(coef,coup,pha,dot,Nsup)
   else
      ! BUTCHER
      allocate(mat(Nsup,Nsup),temp1(Nsup),temp2(Nsup),temp3(Nsup))
      temp1 = dot(2,:)*9.0d0 + dot(3,:)*3.0d0
      coef(1,:) = coef(3,:) + temp1 * dt / 8.0d0
      mat = coup(1,:,:) ! Save coupling at this time
      coup(1,:,:) = ( mat + coup(2,:,:) ) * 0.5d0
      call dot_calculation(coef,coup,pha,dot,Nsup)
      coup(1,:,:) = mat; deallocate(mat)

      temp1 = ( coef(2,:)*28.0d0 - coef(3,:)*23.0d0 ) / 5.0d0
      temp2 = ( dot(1,:)*32.0d0  - dot(2,:)*60.0d0  ) * dt/15.0d0
      temp3 = ( dot(3,:)*26.0d0 ) * dt/15.0d0
      coef(1,:) = temp1 + temp2 - temp3

      allocate(vec(Nsup)); vec = dot(1,:)
      call dot_calculation(coef,coup,pha,dot,Nsup)
      temp1 = ( coef(2,:)*32.0d0 - coef(3,:) ) / 31.0d0
      temp2 = vec*64.0d0 + dot(1,:)*15.0d0 + dot(2,:)*12.0d0
      temp3 = dot(3,:)
      coef(1,:) = temp1 + ( temp2-temp3 ) * dt/93.0d0
      call dot_calculation(coef,coup,pha,dot,Nsup)
      deallocate(vec,temp1,temp2,temp3)
   endif
end subroutine coef_evolution

subroutine do_probabilities(coef,coup,pha,ene,Nsup,dt,inS)
use fstsh_data, only: current_state, tsh_file, tsh_minprob
! This routine calculate probabilities of HOPP using
! fewest switch algorithm by Tully
! The switch is jj -> kk
   implicit none

   integer,   intent(in) :: Nsup, inS
   TDCOMPLEX, intent(in) :: coef(3,Nsup)
   LIODBLE,   intent(in) :: coup(3,Nsup,Nsup), pha(3,Nsup,Nsup), dt
   LIODBLE,   intent(in) :: ene(3,Nsup)

   integer   :: jj, kk
   TDCOMPLEX :: cj, ck, tmpc
   LIODBLE   :: norm, number_random, tmpr, cj2
   LIODBLE, dimension(:), allocatable :: prob

   if (inS == 0) stop "Something wrong in do_probabilities, &
                       Nuclear Step = 0"

   jj = current_state
   cj  = coef(1,jj)
   allocate(prob(Nsup)); prob=0.0d0
   norm = 0.0d0

   norm = 0.0d0
   do kk=1,Nsup
      tmpc = exp(cmplx(0.0d0,pha(1,kk,jj),8))
      ck = conjg(coef(1,kk))
      tmpr = real(cj*ck*tmpc)
      prob(kk) = tmpr*coup(1,kk,jj)*(-2.0d0)
      norm = norm + abs(coef(1,kk))**2
      if ( prob(kk) < 0.0d0 ) prob(kk) = 0.0d0
   enddo
   cj2 = abs(cj)**2
   prob = prob * dt / cj2
  
   write(tsh_file,"(4X,A)") "Population, Probabilities, Norm of Nacvs"
   do kk=1,Nsup
      write(tsh_file,"(4X,A,I2,F10.5,F10.5,F10.5)") "PPN", kk, abs(coef(1,kk))**2, prob(kk), dabs(coup(1,current_state,kk))**2
   enddo
   write(tsh_file,"(4X,A,F10.5)") "Total Pobl.", norm

   call random_number(number_random)
   write(tsh_file,"(4X,A,F10.5,A,I2,A,F10.5)") "random= ", number_random, " to_state= ", maxloc(prob), " max_prob= ", maxval(prob)

   if ( ene(1,current_state) < ene(1,1) .and. current_state /= 1 ) then
      write(tsh_file,"(4X,A)") "Forcing the system at Ground State"
      prob(1) = 1.0d0
   endif

   if ( maxval(prob) > number_random .and. maxval(prob) > tsh_minprob ) then
      write(tsh_file,"(4X,A,I2,A,I2)") "HOPP= ", current_state, " -> ", maxloc(prob)
      current_state = maxloc(prob,1)
   endif
   deallocate(prob)
end subroutine do_probabilities

subroutine adjust_vel(oldSUP,newSUP,Ene,nucvel,Nsup)
use fstsh_data, only: tsh_file, current_state, after_hopp
use garcha_mod, only: atom_mass, natom
   implicit none

   integer, intent(in)    :: oldSUP, newSUP, Nsup
   LIODBLE, intent(in)    :: Ene(Nsup)
   LIODBLE, intent(inout) :: nucvel(3,natom)

   integer :: ii, jj
   LIODBLE :: kinEold, kinEnew, deltaE, totalDiff, alpha
   LIODBLE :: beta, temp, scalv, factor, mass, sgm1, sgm2, fac
   LIODBLE, dimension(:,:), allocatable :: mom

   ! 1fs = 41.34137457575 a.u.
   fac = 41.34137457575d0

   ! Kinetic Energy before hopp
   kinEold = 0.0d0
   do ii=1,natom
      mass = atom_mass(ii)
      kinEold = kinEold + mass * (nucvel(1,ii)/fac)**2
      kinEold = kinEold + mass * (nucvel(2,ii)/fac)**2
      kinEold = kinEold + mass * (nucvel(3,ii)/fac)**2
   enddo
   kinEold = 0.5d0 * kinEold
   write(tsh_file,"(3X,A,F10.5,A)") "Kinetic Energy before hopp= ",kinEold," ha."

   ! Potential Energy
   deltaE = Ene(newSUP) - Ene(oldSUP)
   write(tsh_file,"(3X,A,F10.5,A)") "Potential Energy HOPP= ",deltaE, " ha."

   totalDiff = deltaE - kinEold
   scalv = 1.0d0 - (deltaE/kinEold)

   write(tsh_file,"(3X,A,F10.5)") "dEp-dEk= ", totalDiff
   if (totalDiff > 0.0d0) then
      write(tsh_file,"(3X,A)") "Forbidden HOPP. Velocities are inverting"
      nucvel = -1.0d0 * nucvel
      current_state = OldSUP
   else
      write(tsh_file,"(3X,A)") "Allowing HOPP. Rescaling Velocities"
      after_hopp = .true. ! The next nuclear step -> raw sigma extrapolation

      ! Momentum calculated
      allocate(mom(3,natom))
      do ii=1,natom
         mass = atom_mass(ii)
         mom(1,ii) = mass * nucvel(1,ii) / fac
         mom(2,ii) = mass * nucvel(2,ii) / fac
         mom(3,ii) = mass * nucvel(3,ii) / fac
      enddo
      ! Coef calculation
      alpha = 0.0d0; beta = 0.0d0
      do ii=1,natom
         mass = atom_mass(ii)
         do jj=1,3
            alpha = alpha + 0.5d0*mom(jj,ii)*mom(jj,ii)/mass
            beta  = beta  + (nucvel(jj,ii)/fac) * mom(jj,ii)
         enddo
      enddo

      temp = beta*beta + 4.0d0*alpha*(-deltaE)
      if ( temp < 0.0d0 ) then
         write(tsh_file,"(3X,A)") "Reserving Velocity"
         scalv = dsqrt(scalv)
         factor = (scalv-beta) / alpha
         do ii=1,natom
            mass = atom_mass(ii)
            nucvel(1,ii) = nucvel(1,ii) + factor * nucvel(1,ii)
            nucvel(2,ii) = nucvel(2,ii) + factor * nucvel(2,ii)
            nucvel(3,ii) = nucvel(3,ii) + factor * nucvel(3,ii)
         enddo
      else
         write(tsh_file,"(3X,A)") "Scaling Moment"
         sgm1 = 1.0d0/(2.0d0*alpha) * (beta+dsqrt(temp))
         sgm2 = 1.0d0/(2.0d0*alpha) * (beta-dsqrt(temp))
         factor = min(dabs(sgm1),dabs(sgm2))
         if (dabs(factor-dabs(sgm1)) < 1.0d-6) factor=sgm1
         if (dabs(factor-dabs(sgm2)) < 1.0d-6) factor=sgm2

         do ii=1,natom
            mass = atom_mass(ii)
            nucvel(1,ii) = nucvel(1,ii) - (factor * mom(1,ii) * fac / mass)
            nucvel(2,ii) = nucvel(2,ii) - (factor * mom(2,ii) * fac / mass)
            nucvel(3,ii) = nucvel(3,ii) - (factor * mom(3,ii) * fac / mass)
         enddo
      endif
      deallocate(mom)
      write(tsh_file,"(3X,A,F10.5)") "Rescaling Factor= ", factor
   endif

   ! Kinetic Energy after hopp
   kinEnew = 0.0d0
   do ii=1,natom
      mass = atom_mass(ii)
      kinEnew = kinEnew + mass * (nucvel(1,ii)/fac)**2
      kinEnew = kinEnew + mass * (nucvel(2,ii)/fac)**2
      kinEnew = kinEnew + mass * (nucvel(3,ii)/fac)**2
   enddo
   kinEnew = 0.5d0 * kinEnew

   write(tsh_file,"(3X,A,F10.5,A)") "Kinetic Energy after hopp= ",kinEnew, " ha."
   write(tsh_file,"(3X,A,F10.5,A)") "Diff. Kinetic Energy= ", kinEnew-kinEold, " ha."
   write(tsh_file,"(3X,A,F10.5,A)") "Diff. Total Energy= ", Ene(newSUP)+kinEnew-Ene(oldSUP)-kinEold, " ha."
end subroutine adjust_vel

subroutine correct_decoherence(coef,ene,vel,natom,Nsup,dt_elec)
use fstsh_data, only: current_state
use garcha_mod, only: atom_mass
   implicit none

   integer, intent(in)      :: Nsup, natom
   TDCOMPLEX, intent(inout) :: coef(3,Nsup)
   LIODBLE, intent(in)      :: ene(3,Nsup), vel(3,natom), dt_elec

   integer :: ii, jj
   LIODBLE :: kin_e, mass, norm, temp1, temp2, fac

   ! fs to au
   fac = 0.02418884254d0 ** 2
   ! Obtain Kinetic Energy
   kin_e = 0.0d0
   do ii=1, natom
      mass = atom_mass(ii)
      kin_e = kin_e + mass * vel(1,ii)**2.0d0 * fac
      kin_e = kin_e + mass * vel(2,ii)**2.0d0 * fac
      kin_e = kin_e + mass * vel(3,ii)**2.0d0 * fac
   enddo

   kin_e = 0.1d0 / kin_e ! 0.1 is the decay damping
   ii = current_state; norm = 0.0d0
   do jj=1, Nsup
      if ( ii /= jj ) then
         temp1 = 1.0d0 / dabs(ene(1,jj)-ene(1,ii))
         temp2 = temp1 * ( 1.0d0 + kin_e )
         temp1 = dsqrt( exp(-dt_elec/temp2) )
         coef(1,jj) = coef(1,jj) * temp1
         norm = norm + abs(coef(1,jj))**2
      endif
   enddo

   temp1 = abs(coef(1,ii))**2.0d0
   temp2 = dsqrt( (1.d0-norm) / temp1 )
   coef(1,ii) = coef(1,ii) * temp2
end subroutine correct_decoherence
