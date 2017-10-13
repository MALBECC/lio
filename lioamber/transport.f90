!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
module transport
   implicit none
   logical              :: transport_calc=.false., generate_rho0=.false., gate_field=.false., lpop
   integer              :: save_charge_freq=0, Pop_Drive
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

!This subroutine classify each inex for the evolution of density matrix during propagation

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

    do i=1,M
    do j=1,M
       if((group(nuc(i)).eq.1).and.(group(nuc(j)).eq.1)) mapmat(i,j)=1
       if((group(nuc(i)).eq.1).and.(group(nuc(j)).eq.2)) mapmat(i,j)=2
       if((group(nuc(i)).eq.1).and.(group(nuc(j)).eq.3)) mapmat(i,j)=3
       if((group(nuc(i)).eq.2).and.(group(nuc(j)).eq.1)) mapmat(i,j)=4
       if((group(nuc(i)).eq.2).and.(group(nuc(j)).eq.2)) mapmat(i,j)=5
       if((group(nuc(i)).eq.2).and.(group(nuc(j)).eq.3)) mapmat(i,j)=6
       if((group(nuc(i)).eq.3).and.(group(nuc(j)).eq.1)) mapmat(i,j)=7
       if((group(nuc(i)).eq.3).and.(group(nuc(j)).eq.2)) mapmat(i,j)=8
       if((group(nuc(i)).eq.3).and.(group(nuc(j)).eq.3)) mapmat(i,j)=9
    end do
    end do

!Counting the number of basis for each part of transport
    do i=1,M
       if(mapmat(i,i).eq.1) k=k+1
       if(mapmat(i,i).eq.5) nn=nn+1
       if(mapmat(i,i).eq.9) n=n+1
    end do

    write(*,*) 'Basis from group 1 =', k
    write(*,*) 'Basis from group 2 =', nn
    write(*,*) 'Basis from group 3 =', n


end subroutine mat_map

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subroutine electrostat (rho1,mapmat,overlap,rhofirst,Gamma0, M)

!This subroutine modify the density matrix in function of the belong of the indexes

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
         
   do i=1,M
   do j=1,M
      if((mapmat(i,j).eq.0).or.(mapmat(i,j).eq.9)) then
         rho_scratch(i,j,1)=dcmplx(0.0D0,0.0D0)
         rho_scratch(i,j,2)=dcmplx(0.0D0,0.0D0)
      elseif((mapmat(i,j).eq.1).or.(mapmat(i,j).eq.2).or.(mapmat(i,j).eq.4).or.(mapmat(i,j).eq.5)) then
         rho_scratch(i,j,1)=(rho1(i,j))
         rho_scratch(i,j,2)=(rhofirst(i,j))
      elseif((mapmat(i,j).eq.3).or.(mapmat(i,j).eq.7).or.(mapmat(i,j).eq.6).or.(mapmat(i,j).eq.8)) then
         rho_scratch(i,j,1)=(0.50D0*rho1(i,j))
         rho_scratch(i,j,2)=(0.50D0*rhofirst(i,j))
      endif
   end do
   end do

   GammaIny=Gamma0*0.5D0
   GammaAbs=GammaIny
   write(*,*) 'GammaAbs,GammaIny =',GammaAbs,GammaIny
      do i=1,M
         do j=1,M
            rho1(i,j)= (GammaAbs*rho_scratch(i,j,1))-(GammaIny*rho_scratch(i,j,2))
         end do
      end do

!-------Stop if NaN-----------!

   do i=1,M
      do j=1,M
         if(rho1(i,j).ne.rho1(i,j)) then
            stop 'Huston, we have a problem'
         end if
      end do
   end do

   deallocate(rho_scratch)

   call g2g_timer_stop('electrostat')

end subroutine electrostat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine Drive_Population(Pop,ngroup,rho,overlap,&
                            group,smat,q,uid)
        
        use garcha_mod, only : Nuc,natom,M

        implicit none

        integer, intent(in) :: Pop
        integer, intent(in) :: ngroup, uid
        real*8, dimension(M,M), intent(in) :: rho, overlap, smat
        integer, dimension(natom), intent(in) :: group
        real*8, dimension(natom), intent(in) :: q
        real*8, dimension(ngroup) :: qgr
        real*8 :: traza
        integer :: i

        qgr(:) = 0.0D0
        traza= 0.0D0

        if ( Pop == 1) then
             call mulliken_calc(natom,M,rho,overlap,Nuc,q)
        elseif ( Pop == 2 ) then
             call lowdin_calc(M,natom,rho,smat,Nuc,q)
        endif

        do i=1,natom
           qgr(group(i)) = qgr(group(i)) + q(i)
        enddo

        do i=1,ngroup
           write(uid,*) i, i, qgr(i)
           traza = traza + qgr(i)
        enddo
      
        if (uid == 678) then
           write(uid,*) "---------------------------"
        else
           write(uid,*) "tot=", traza
           write(uid,*) "----------------------------"
        endif

end subroutine Drive_Population
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end module
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
