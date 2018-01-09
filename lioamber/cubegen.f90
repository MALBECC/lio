!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module cubegen
   implicit none
   contains
!
! Performs orbital/density plots.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   subroutine cubegen_vecin( Msize, Nocc, coef_vec )
      use garcha_mod, only: cube_dens, cube_orb, cube_elec, cubegen_only, VCINP
      implicit none
      integer, intent(in)  :: Msize
      integer, intent(in)  :: Nocc
      real*8 , intent(in)  :: coef_vec( Msize*(Msize+1)/2 )
      real*8 , allocatable :: coef_mat( :, : )
      integer              :: ii, jj, kk


      if ( .not. (cube_dens.or.cube_orb.or.cube_elec) ) return
 
      if ( cubegen_only .and. (.not.VCINP) ) then
         write(*,*) "cubegen_only CAN ONLY BE USED WITH VCINP"
         stop
      end if

      if (allocated(coef_mat)) deallocate(coef_mat)
      allocate( coef_mat(Msize, Msize) )

!     Can't this copy the whole array?
      kk = 0
      do ii = 1, Nocc
      do jj = 1, Msize
         kk=kk+1
         coef_mat(ii,jj) = coef_vec(kk)
      enddo
      enddo

      call g2g_timer_sum_start('cube gen')
      call cubegen_write( coef_mat )
      call g2g_timer_sum_stop('cube gen')

      deallocate( coef_mat )
   end subroutine cubegen_vecin



!------------------------------------------------------------------------------!
   subroutine cubegen_matin( Msize, ugly_mat )
      use garcha_mod, only: cube_dens, cube_orb, cube_elec, cubegen_only, VCINP
      implicit none
      integer, intent(in)  :: Msize
      real*8 , intent(in)  :: ugly_mat( Msize, 3*Msize )
      real*8 , allocatable :: coef_mat( :, : )
      integer              :: ii, jj


      if ( .not. (cube_dens.or.cube_orb.or.cube_elec) ) return

      if ( cubegen_only .and. (.not.VCINP) ) then
         write(*,*) "cubegen_only CAN ONLY BE USED WITH VCINP"
         stop
      end if

      if (allocated(coef_mat)) deallocate(coef_mat)
      allocate( coef_mat(Msize, Msize) )

      do jj = 1, Msize
      do ii = 1, Msize
         coef_mat(ii,jj) = ugly_mat( ii, 2*Msize + jj )
      enddo
      enddo

      call g2g_timer_sum_start('cube gen')
      call cubegen_write( coef_mat )
      call g2g_timer_sum_stop('cube gen')

      deallocate( coef_mat )
   end subroutine cubegen_matin


!------------------------------------------------------------------------------!
   subroutine cubegen_write( MO_v )
   use garcha_mod

   implicit real*8 (a-h,o-z)

   real*8, dimension(M,M) :: MO_v
   real*8 :: min_exps(120),x0(3),x1(3),origin(3),eval_p(3)
   real*8 :: max_radius, max_dim, vox_dim, p_val, p_func
   real*8 :: p_dist, Morb
   integer :: i,j,k,ii,jj,kk,iii,jjj,kkk
   integer :: ns, np, ni
   integer :: ivox, ivoxx, ivoxy, ivoxz, kk_dens, kk_orb
   integer :: MM, MMd, M1, M2, M3, M5, M7, M9, M11, M13, M15

   parameter(expmax=10)

   MM  = M *  (M+1)  / 2
   MMd = Md * (Md+1) / 2

   M1  = 1         ! first P
   M3  = M1  + MM  ! now Pnew
   M5  = M3  + MM  ! now S, F also uses the same position after S was used
   M7  = M5  + MM  ! now G
   M9  = M7  + MMd ! now Gm
   M11 = M9  + MMd ! now H
   M13 = M11 + MM  ! W ( eigenvalues ), also this space is used in least squares
   M15 = M13 + M   ! aux ( vector for ESSl)

   if (cube_dens) open(unit=4242,file=cube_dens_file)
   if (cube_orb) open(unit=4243,file=cube_orb_file)

   ! first find the rectangular prism that the system lies in
   ivox = cube_res
   min_exps = 999999.D0
   do i=1,M
     j = Iz(Nuc(i))
     do ni=1,ncont(i)
       min_exps(j)=min(min_exps(j),a(i,ni))
     enddo
   enddo

   x0 = 0.D0
   x1 = 0.D0
   do i=1,natom
     ! TODO: not sure about this padding criteria
     max_radius = 2.0D0*sqrt(1.0D0/min_exps(Iz(i)))
     if (i.eq.1) then
       x0 = r(i,:) - max_radius
       x1 = r(i,:) + max_radius
     else
       x0(1) = min(x0(1), r(i,1) - max_radius)
       x0(2) = min(x0(2), r(i,2) - max_radius)
       x0(3) = min(x0(3), r(i,3) - max_radius)

       x1(1) = max(x1(1), r(i,1) + max_radius)
       x1(2) = max(x1(2), r(i,2) + max_radius)
       x1(3) = max(x1(3), r(i,3) + max_radius)
     endif
   enddo
   origin = x0

   ! cube_res gives the # of voxels in the longest dimension
   ! so, next figure out the voxel length based on this criteria
   ! and how many voxels each dimension needs based on voxel length
   max_dim = -1.D0
   do i=1,3
     max_dim = max(max_dim,x1(i)-x0(i))
   enddo
   vox_dim = max_dim / ivox

   ivoxx = ceiling((x1(1)-x0(1)) / vox_dim)
   ivoxy = ceiling((x1(2)-x0(2)) / vox_dim)
   ivoxz = ceiling((x1(3)-x0(3)) / vox_dim)

   if (cube_dens) then
     write(4242,*) "LIO CUBE FILE"
     write(4242,*) "CUBE FORMAT FILE FOR ELECTRON DENSITY"
     write(4242,42) natom,origin(:)
     write(4242,42) ivoxx,vox_dim,0.,0.
     write(4242,42) ivoxy,0.,vox_dim,0.
     write(4242,42) ivoxz,0.,0.,vox_dim
     do i=1,natom
       write(4242,424) Iz(i),0.,r(i,:)
     enddo
   endif
   if (cube_orb) then
     write(4243,*) "LIO CUBE FILE"
     if (cube_sel.eq.0) then
       write(4243,*) "CUBE FORMAT FILE FOR MOLECULAR ORBITALS"
     elseif (cube_sel.gt.0.and.cube_sel.le.NCO) then
       write(4243,*) "CUBE FORMAT FILE FOR SINGLE MOLECULAR ORBITAL"
     else
       write(*,*) "cube_sel VALUE NOT SUPPORTED!"
       stop
     endif
     if (cube_sel.eq.0) then
       write(4243,42) -natom,origin(:)
     else
       write(4243,42) natom,origin(:)
     endif
     write(4243,42) ivoxx,vox_dim,0.,0.
     write(4243,42) ivoxy,0.,vox_dim,0.
     write(4243,42) ivoxz,0.,0.,vox_dim
     do i=1,natom
       write(4243,424) Iz(i),0.,r(i,:)
     enddo
     if (cube_sel.eq.0) then
       write(4243,'(I4)',advance='no') M
       do i=1,M
         write(4243,'(I4)',advance='no') i-NCO
       enddo
       write(4243,*) ""
     endif
   endif

   ns = nshell(0)
   np = nshell(1)
   kk_dens = 1
   kk_orb = 1
   ! for each voxel...
   do i=1,ivoxx
     eval_p(1) = origin(1) + (i-1) * vox_dim
     do j=1,ivoxy
       eval_p(2) = origin(2) + (j-1) * vox_dim
       do k=1,ivoxz
         eval_p(3) = origin(3) + (k-1) * vox_dim
           p_val = 0.D0
         ! calculate function values at this voxel, store in RMM(M15)
         ! s functions
         do ii=1,ns
           p_dist = 0.D0
           do jj = 1,3
             p_dist = p_dist + (eval_p(jj) - r(Nuc(ii),jj))**2
           enddo
           p_func = 0.D0
           do ni=1,ncont(ii)
            if(a(ii,ni)*p_dist.lt.expmax) p_func = p_func + c(ii,ni) * exp(-a(ii,ni)*p_dist)
           enddo

           RMM(M15+ii-1) = p_func
         enddo

         ! p functions
         do ii=ns+1,ns+np,3
           p_dist = 0.D0
           do jj = 1,3
             p_dist = p_dist + (eval_p(jj) - r(Nuc(ii),jj))**2
           enddo
           p_func = 0.D0
           do ni=1,ncont(ii)
            if(a(ii,ni)*p_dist.lt.expmax) p_func = p_func + c(ii,ni) * exp(-a(ii,ni)*p_dist)
           enddo

           do jj = 1,3
             RMM(M15+ii+jj-2) = p_func * (eval_p(jj)-r(Nuc(ii),jj))
           enddo
         enddo

         ! d functions
         do ii=ns+np+1,M,6
           p_dist = 0.D0
           do jj = 1,3
             p_dist = p_dist + (eval_p(jj) - r(Nuc(ii),jj))**2
           enddo
           p_func = 0.D0
           do ni=1,ncont(ii)
            if(a(ii,ni)*p_dist.lt.expmax) p_func = p_func + c(ii,ni) * exp(-a(ii,ni)*p_dist)
           enddo

           kkk = 0
           do jj = 1,3
             do jjj = 1,jj
               kkk = kkk + 1
               RMM(M15+ii+kkk-2)=p_func*(eval_p(jj)-r(Nuc(ii),jj))*(eval_p(jjj)-r(Nuc(ii),jjj))
             enddo
           enddo
         enddo


         if (cube_dens) then
           ! calculate density for this voxel
           kkk = 0
           do ii=1,M
             do jj=ii,M
               kkk = kkk + 1
               p_val=p_val+RMM(kkk)*RMM(M15+ii-1)*RMM(M15+jj-1)
             enddo
           enddo
           write(4242,'(E13.5)',advance='no') p_val
           if (mod(kk_dens,6) .eq. 0) then
             write(4242,*) ""
           endif
           kk_dens = kk_dens + 1
         endif

         if (cube_orb) then
           if (cube_sel.eq.0) then
             do kk=1,M
               ! calculate orbital or orbital^2 if cube_sqrt_orb=treu for this voxel
               p_val = 0.D0
               do ii=1,M
                 do jj=ii+1,M
                   Morb=2.D0*MO_v(kk,ii)*RMM(M15+ii-1)
                   if (cube_sqrt_orb) Morb=Morb*MO_v(kk,ii)*RMM(M15+ii-1)
                   p_val=p_val+Morb
                 enddo
                   Morb=MO_v(kk,ii)*RMM(M15+ii-1)
                   if (cube_sqrt_orb) Morb=Morb*MO_v(kk,ii)*RMM(M15+ii-1)
                   p_val=p_val+Morb
               enddo
               write(4243,'(E13.5)',advance='no') p_val
               if (mod(kk_orb,6) .eq. 0) then
                 write(4243,*) ""
               endif
               kk_orb = kk_orb + 1
             enddo
           else
             ! calculate orbital or orbital^2 if cube_sqrt_orb=treu for this voxel
             p_val = 0.D0
             do ii=1,M
               do jj=ii+1,M
                 Morb=2.D0*MO_v(cube_sel,ii)*RMM(M15+ii-1)
                 if (cube_sqrt_orb) Morb=Morb*MO_v(cube_sel,jj)*RMM(M15+jj-1)
                 p_val=p_val+Morb
               enddo
                 Morb=MO_v(cube_sel,ii)*RMM(M15+ii-1)
                 if (cube_sqrt_orb) Morb=Morb*MO_v(cube_sel,ii)*RMM(M15+ii-1)
                 p_val=p_val+Morb
             enddo
             write(4243,'(E13.5)',advance='no') p_val
             if (mod(kk_orb,6) .eq. 0) then
               write(4243,*) ""
             endif
             kk_orb = kk_orb + 1
           endif
         endif
       enddo
     enddo
   enddo
   if (cube_dens) close(4242)
   if (cube_orb) close(4243)

   if (cube_elec) then 
   call elec(ivoxx,ivoxy,ivoxz,vox_dim,origin(1),origin(2),origin(3))
   endif

  42  format(I5,3(f12.6))
 424  format(I5,4(f12.6))

   end subroutine cubegen_write


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
end module cubegen
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
