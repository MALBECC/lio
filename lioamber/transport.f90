!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
module transport
   implicit none
   logical              :: transport_calc=.false., generate_rho0=.false., gate_field=.false., lpop
   integer              :: save_charge_freq=0
   complex*8            :: traza0, traza
   real*8               :: driving_rate=0.001, scratchgamma, GammaMagnus, GammaVerlet, re_traza

#ifdef TD_SIMPLE
   complex*8,allocatable  :: rhofirst (:,:)
#else
   complex*16,allocatable :: rhofirst (:,:)
#endif
   integer,allocatable,dimension(:,:) :: mapmat

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine mat_map (group, mapmat, Nuc, M, natom)
   implicit none
   integer, intent(in)   :: M, natom
   integer, intent(in)   :: group(natom), Nuc(M)
   integer, intent(out)  :: mapmat (M,M)
   integer               :: i, j, k, nn, n
   write(*,*) 'M =', M
   write(*,*) 'natoms=', natom

   mapmat=0
   i=0
   j=0
   n=0
   nn=0

    DO i=1,M
       DO j=1,M
          IF((group(nuc(i)).eq.1).and.(group(nuc(j)).eq.1)) mapmat(i,j)=1
          IF((group(nuc(i)).eq.1).and.(group(nuc(j)).eq.2)) mapmat(i,j)=2
          IF((group(nuc(i)).eq.1).and.(group(nuc(j)).eq.3)) mapmat(i,j)=3
          IF((group(nuc(i)).eq.2).and.(group(nuc(j)).eq.1)) mapmat(i,j)=4
          IF((group(nuc(i)).eq.2).and.(group(nuc(j)).eq.2)) mapmat(i,j)=5
          IF((group(nuc(i)).eq.2).and.(group(nuc(j)).eq.3)) mapmat(i,j)=6
          IF((group(nuc(i)).eq.3).and.(group(nuc(j)).eq.1)) mapmat(i,j)=7
          IF((group(nuc(i)).eq.3).and.(group(nuc(j)).eq.2)) mapmat(i,j)=8
          IF((group(nuc(i)).eq.3).and.(group(nuc(j)).eq.3)) mapmat(i,j)=9
       ENDDO
    ENDDO

!TEST
       DO i=1,M
             if(mapmat(i,i).eq.1) k=k+1
             if(mapmat(i,i).eq.5) nn=nn+1
             if(mapmat(i,i).eq.9) n=n+1
       ENDDO
       write(*,*) 'Basis from group 1 =', k
       write(*,*) 'Basis from group 2 =', nn
       write(*,*) 'Basis from group 3 =', n


end subroutine mat_map

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine electrostat (rho1,mapmat,overlap,rhofirst,Gamma0, M)
   implicit none
   integer, intent(in) :: M
   integer, intent(in) :: mapmat(M,M)
   integer             :: i, j
   real*8,  intent(in) :: overlap(M,M), Gamma0
   real*8              :: GammaIny, GammaAbs
   
#ifdef TD_SIMPLE
   complex*8, intent(in)     :: rhofirst(M,M)
   complex*8, intent(inout)  :: rho1(M,M)
   complex*8, allocatable    :: rho_scratch (:,:,:)
#else
   complex*16, intent(in)    :: rhofirst(M,M)
   complex*16, intent(inout) :: rho1(M,M)
   complex*16, allocatable   :: rho_scratch (:,:,:)
#endif

   call g2g_timer_start('electrostat')
         
   allocate(rho_scratch(M,M,2))
   rho_scratch=0D0
         
   DO i=1,M
      DO j=1,M
         IF((mapmat(i,j).eq.0).or.(mapmat(i,j).eq.9)) THEN
            rho_scratch(i,j,1)=dcmplx(0.0D0,0.0D0)
            rho_scratch(i,j,2)=dcmplx(0.0D0,0.0D0)
         elseif((mapmat(i,j).eq.1).or.(mapmat(i,j).eq.2).or.(mapmat(i,j).eq.4).or.(mapmat(i,j).eq.5)) THEN
            rho_scratch(i,j,1)=(rho1(i,j))
            rho_scratch(i,j,2)=(rhofirst(i,j))
         elseif((mapmat(i,j).eq.3).or.(mapmat(i,j).eq.7).or.(mapmat(i,j).eq.6).or.(mapmat(i,j).eq.8)) then
            rho_scratch(i,j,1)=(0.50D0*rho1(i,j))
            rho_scratch(i,j,2)=(0.50D0*rhofirst(i,j))
         endif

!         IF(mapmat(i,j).eq.0) THEN
!            rho_scratch(i,j,1)=dcmplx(0.0D0,0.0D0)
!            rho_scratch(i,j,2)=dcmplx(0.0D0,0.0D0)
!         ENDIF
!         IF(mapmat(i,j).eq.9) THEN
!            rho_scratch(i,j,1)=dcmplx(0.0D0,0.0D0)
!            rho_scratch(i,j,2)=dcmplx(0.0D0,0.0D0)
!         ENDIF
!         IF(mapmat(i,j).eq.1) THEN
!            rho_scratch(i,j,1)=(rho1(i,j))
!            rho_scratch(i,j,2)=(rhofirst(i,j))
!         ENDIF
!         IF((mapmat(i,j).eq.3).or.(mapmat(i,j).eq.7)) THEN
!            rho_scratch(i,j,1)=(0.50D0*rho1(i,j))
!            rho_scratch(i,j,2)=(0.50D0*rhofirst(i,j))
!         ENDIF
!         IF((mapmat(i,j).eq.2).or.(mapmat(i,j).eq.4)) THEN
!            rho_scratch(i,j,1)=rho1(i,j)
!            rho_scratch(i,j,2)=(rhofirst(i,j))
!         ENDIF
!         IF(mapmat(i,j).eq.5) THEN
!            rho_scratch(i,j,1)=rho1(i,j)
!            rho_scratch(i,j,2)=rhofirst(i,j)
!         ENDIF
!         IF((mapmat(i,j).eq.6).or.(mapmat(i,j).eq.8)) THEN
!            rho_scratch(i,j,1)=(0.50D0*rho1(i,j))
!            rho_scratch(i,j,2)=(0.50D0*rhofirst(i,j))
!         ENDIF
      ENDDO
   ENDDO

   GammaIny=Gamma0*0.5D0
   GammaAbs=GammaIny
   write(*,*) 'GammaAbs,GammaIny =',GammaAbs,GammaIny
      DO i=1,M
         DO j=1,M
            rho1(i,j)= (GammaAbs*rho_scratch(i,j,1))-(GammaIny*rho_scratch(i,j,2))
         ENDDO
      ENDDO

!-------Stop if NaN-----------!

   DO i=1,M
      DO j=1,M
         IF(rho1(i,j).ne.rho1(i,j)) THEN
            stop 'Huston, we have a problem'
         ENDIF
      ENDDO
   ENDDO

   DEALLOCATE(rho_scratch)

   call g2g_timer_stop('electrostat')

end subroutine electrostat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end module
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
