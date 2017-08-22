!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine RMMcalc0_Init()
!------------------------------------------------------------------------------!
!
! DESCRIPTION
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  use garcha_mod, only: M, igrid2, nuc, natom, natomc, d, r, atmin, &
                        rmax, jatc, nnps, nnpp, nnpd, nshell
  implicit none
  real*8   :: alf,rexp
  integer  :: zij,ti,tj
  integer  :: ii,jj

  call g2g_timer_start('RMMcalc0')
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

  call g2g_reload_atom_positions(igrid2)
  call g2g_timer_stop('RMMcalc0')
end subroutine RMMcalc0_Init
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
