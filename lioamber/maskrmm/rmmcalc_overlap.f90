!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine rmmcalc_overlap( overlap_m, energy )
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  use garcha_mod, only: M, Md, nuc, natom, natomc, d, r, atmin, &
                        rmax, jatc, nnps, nnpp, nnpd, nshell, &
                        igrid2, kkind, kkinds, cool, cools
  implicit none
  real*8,intent(out) :: overlap_m(M,M)
  real*8,intent(out) :: energy

  real*8   :: energy_nuc
  real*8   :: energy_solv_t
  real*8   :: energy_solv_f

  real*8   :: alf, rexp
  integer  :: zij, ti, tj
  integer  :: igpu, ii, jj
  logical  :: MEMO

  call g2g_timer_start('rmmcalc_overlap')
  energy_nuc = 0.0d0
  energy_solv_t = 0.0d0
  energy_solv_f = 0.0d0

  do ii=1,natom
    natomc(ii)=0
    do jj=1,natom
      d(ii,jj)=0.0d0
      d(ii,jj)=d(ii,jj)+(r(ii,1)-r(jj,1))**2
      d(ii,jj)=d(ii,jj)+(r(ii,2)-r(jj,2))**2
      d(ii,jj)=d(ii,jj)+(r(ii,3)-r(jj,3))**2
      zij=atmin(ii)+atmin(jj)
      ti=atmin(ii)/zij
      tj=atmin(jj)/zij
      alf=atmin(ii)*tj
      rexp=alf*d(ii,jj)
      if (rexp.lt.rmax) then
        natomc(ii)=natomc(ii)+1
        jatc(natomc(ii),ii)=jj
      endif
    enddo
  enddo

  do ii=nshell(0),1,-1
    nnps(nuc(ii))=ii
  enddo

  do ii=nshell(0)+nshell(1),nshell(0)+1,-1
    nnpp(nuc(ii))=ii
  enddo

  do ii=M,nshell(0)+nshell(1)+1,-1
    nnpd(nuc(ii))=ii
  enddo

  call g2g_timer_start('xc grid setup')
  call g2g_reload_atom_positions(igrid2)
  call g2g_timer_stop('xc grid setup')

  call aint_query_gpu_level(igpu)
  if (igpu.gt.1) call aint_new_step()
  call int1(energy_nuc)
  call rmmget_fock(overlap_m)


  if (allocated(kkind))  deallocate(kkind)
  if (allocated(kkinds)) deallocate(kkinds)
  if (allocated(cool))   deallocate(cool)
  if (allocated(cools))  deallocate(cools)

  call g2g_reload_atom_positions(igrid2)
  call aint_query_gpu_level(igpu)
  if (igpu.gt.1) call aint_new_step()

  if (igpu.le.1) then
    call intsol(energy_solv_f,energy_solv_t,.true.)
  else
    call aint_qmmm_fock(energy_solv_f,energy_solv_t)
  endif

  call int2()

  if (igpu.gt.2) call aint_coulomb_init()
  if (igpu.eq.5) MEMO = .false.
  if (MEMO) then
    call g2g_timer_start('int3mem')
    call int3mem()
    call g2g_timer_stop('int3mem')
  endif

  energy=0.0d0
  energy=energy+energy_nuc
  energy=energy+energy_solv_t
  call g2g_timer_stop('rmmcalc_overlap')

end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
