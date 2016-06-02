!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine rmmcalc_fockmao( dens_mao, fock_mao, energy )
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  use garcha_mod, only:M,Md,RMM
  implicit none
  complex*16,intent(in) :: dens_mao(M,M)
  real*8,intent(out)    :: fock_mao(M,M)
  real*8,intent(out)    :: energy

  real*8   :: energy_1e
  real*8   :: energy_coulomb
  real*8   :: energy_xc
  integer  :: MM, MMd, idx0, kk

  call g2g_timer_start('rmmcalc_fockmao')
  energy_1e = 0.0d0
  energy_coulomb = 0.0d0
  energy_xc = 0.0d0

  call rmmput_dens(dens_mao)
  call int3lu(energy_coulomb)
  call g2g_solve_groups(0,energy_xc,0)
  call rmmget_fock(fock_mao)

  MM=M*(M+1)/2
  MMd=Md*(Md+1)/2
  idx0=3*MM+2*MMd
  do kk=1,MM
    energy_1e = energy_1e + RMM(kk)*RMM(idx0+kk)
  enddo

  energy=energy+energy_1e
  energy=energy+energy_Coulomb
  energy=energy+energy_xc
  call g2g_timer_stop('rmmcalc_fockmao')

end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
