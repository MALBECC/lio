!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine calc_forceDS_dss(natoms,nbasis,pos,vel,Mat0,MatB,fterm)
!------------------------------------------------------------------------------!
!
! DESCRIPTION
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
  use basis_data
  implicit none
  integer,intent(in)     :: natoms          ! Number of atoms
  integer,intent(in)     :: nbasis          ! Number of basis
  real*8,intent(in)      :: pos(3,natoms)
  real*8,intent(in)      :: vel(3,natoms)
  complex*16,intent(in)  :: Mat0(nbasis,nbasis)
  real*8,intent(out)     :: MatB(nbasis,nbasis)
  complex*16,intent(out) :: fterm(3,natoms)

  real*8     :: IMTX(3,4,4)
  real*8     :: ri(3),rj(3)
  real*8     :: ai,aj,cij,ct2,cta
  real*8     :: intx,inty,intz,term1,term2
  integer    :: pi(3),pj(3)
  integer    :: atomi,atomj
  integer    :: kk,ii,jj,ni,nj


  call g2g_timer_start('calc_forceDS_dss')
  fterm=dcmplx(0.0d0,0.0d0)
  MatB(:,:)=0.0d0

  do jj=1,Nbasis
  do ii=1,Nbasis
    atomi = parent_atom(ii)
    atomj = parent_atom(jj)
    ri(:) = pos(:,atomi)
    rj(:) = pos(:,atomj)

    do nj=1,orbital_contractions(jj)
    do ni=1,orbital_contractions(ii)
      ai=gauss_expo(ni,ii)
      aj=gauss_expo(nj,jj)
      call calc_gintmat( 0, 1, ai, aj, ri, rj, IMTX )

      cij=gauss_coef(ni,ii)*gauss_coef(nj,jj)
      ct2=cij*2
      cta=ct2*aj

      do kk=1,3
         pi(:)  = angular_momentum(:,ii)
         pj(:)  = angular_momentum(:,jj)
         pj(kk) = pj(kk)+1

         intx  = IMTX(1,1+pi(1),1+pj(1))
         inty  = IMTX(2,1+pi(2),1+pj(2))
         intz  = IMTX(3,1+pi(3),1+pj(3))
         term1 = cta*intx*inty*intz

         term2=0.0d0
         if (pj(kk).ge.2) then
           pj(kk) = pj(kk)-2
           intx = IMTX(1,1+pi(1),1+pj(1))
           inty = IMTX(2,1+pi(2),1+pj(2))
           intz = IMTX(3,1+pi(3),1+pj(3))
           if (pj(kk).eq.0) term2 = -cij*intx*inty*intz
           if (pj(kk).eq.1) term2 = -ct2*intx*inty*intz
         endif

         MatB(ii,jj)=MatB(ii,jj)+vel(kk,atomj)*(term1+term2)
         fterm(kk,atomj)=fterm(kk,atomj)+Mat0(ii,jj)*(term1+term2)
      enddo

    enddo
    enddo

  enddo
  enddo

  call g2g_timer_stop('calc_forceDS_dss')
end subroutine calc_forceDS_dss
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
