!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine calc_forceDS_dds( natoms, nbasis, pos, vel, Mat0, fterm )
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
  complex*16,intent(out) :: fterm(3,natoms)

  real*8     :: IMTX(3,4,4)
  real*8     :: ri(3),rj(3)
  real*8     :: ai,aj,cij,ci2,cj2
  real*8     :: intx,inty,intz
  real*8     :: term0,term1,term2,term3,term4
  integer    :: pi(3),pj(3)
  integer    :: atomi,atomj
  integer    :: ii,jj,ni,nj,ki,kj


  call g2g_timer_start('calc_forceDS_dds')
  fterm=dcmplx(0.0d0,0.0d0)

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
      cij=gauss_coef(ni,ii)*gauss_coef(nj,jj)
      call calc_gintmat( 0, 1, ai, aj, ri, rj, IMTX )

      do kj=1,3
      do ki=1,3
         pi(:)  = angular_momentum(:,ii)
         pj(:)  = angular_momentum(:,jj)

         pi(ki) = angular_momentum(ki,ii)+1
         pj(kj) = angular_momentum(kj,jj)+1
         intx   = IMTX(1,1+pi(1),1+pj(1))
         inty   = IMTX(2,1+pi(2),1+pj(2))
         intz   = IMTX(3,1+pi(3),1+pj(3))
         term1  = cij*ai*aj*intx*inty*intz

         term2=0.0d0
         if (angular_momentum(kj,jj).gt.0) then
           pi(ki) = angular_momentum(ki,ii)+1
           pj(kj) = angular_momentum(kj,jj)-1
           intx = IMTX(1,1+pi(1),1+pj(1))
           inty = IMTX(2,1+pi(2),1+pj(2))
           intz = IMTX(3,1+pi(3),1+pj(3))
           cj2  = -1.0d0
           if (pj(kj).eq.2) cj2 = -2.0d0
           term2=cij*ai*cj2*intx*inty*intz
         endif

         term3=0.0d0
         if (angular_momentum(ki,ii).gt.0) then
           pi(ki) = angular_momentum(ki,ii)-1
           pj(kj) = angular_momentum(kj,jj)+1
           intx = IMTX(1,1+pi(1),1+pj(1))
           inty = IMTX(2,1+pi(2),1+pj(2))
           intz = IMTX(3,1+pi(3),1+pj(3))
           ci2  = -1.0d0
           if (pi(ki).eq.2) ci2 = -2.0d0
           term3=cij*ci2*aj*intx*inty*intz
         endif

         term4=0.0d0
         if ((angular_momentum(kj,jj).gt.0).and.(angular_momentum(ki,ii).gt.0)) then
           pi(ki) = angular_momentum(ki,ii)-1
           pj(kj) = angular_momentum(kj,jj)-1
           intx = IMTX(1,1+pi(1),1+pj(1))
           inty = IMTX(2,1+pi(2),1+pj(2))
           intz = IMTX(3,1+pi(3),1+pj(3))
           term4=cij*ci2*cj2*intx*inty*intz
         endif

         term0=term1+term2+term3+term4
         fterm(kj,atomj)=fterm(kj,atomj)+Mat0(ii,jj)*vel(ki,atomi)*term0
      enddo
      enddo

    enddo
    enddo

  enddo
  enddo

  call g2g_timer_stop('calc_forceDS_dds')
end subroutine calc_forceDS_dds
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
