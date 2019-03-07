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
   use garcha_mod, only: RMM, natom, r, nco, Iz,  cube_dens, cube_orb, &
                         cube_elec, cube_sel, cube_orb_file, cube_res, &
                         cube_dens_file, cube_sqrt_orb
   use basis_data, only: M, Md, a, c, ncont, nuc, nshell, MM, MMd

   implicit none

   real*8, dimension(M,M) :: MO_v
   real*8 :: min_exps(120),x0(3),x1(3),origin(3),eval_p(3)
   real*8 :: max_radius, max_dim, vox_dim, p_val, p_func
   real*8 :: p_dist, Morb
   integer :: i,j,k,ii,jj,kk,iii,jjj,kkk
   integer :: ns, np, ni
   integer :: ivox, ivoxx, ivoxy, ivoxz, kk_dens, kk_orb
   integer :: M1, M2, M15

   real, parameter :: expmax = 10

   M1  = 1         ! first P
   M15 = 1 + 4*MM + 2*MMd + M   ! aux ( vector for ESSl)

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
   
!##############################################################################!
!## ELEC ######################################################################!
! Calculation of electrical potential for an arbitrary point of space, making  !
! 1e integrals via the Obara-Saika recursive method.                           !
!                                                                              !
! Inputs:                                                                      !
!   · NX, NY, NZ      : Number of x, y and z points.                           !
!   · xMin, yMin, zMin: Minimum value for x, y and z coordinates.              !
!   · deltax          : Increments for x, y and z coordinates (same for all).  !
!                                                                              !
! Outputs:                                                                     !
!   · A file containing the electrostatic potential surface.                   !
!                                                                              !
!##############################################################################!
subroutine elec(NX, NY, NZ, deltax, xMin, yMin, zMin)
   use garcha_mod   , only: r, d, natom, cube_elec_file, RMM, Iz
   use constants_mod, only: PI, PI32
   use basis_data   , only: M, norm, nShell, nCont, nuc, a, c
   use liosubs_math , only: funct

   implicit none
   double precision, intent(in) :: xMin, yMin, zMin, deltax
   integer         , intent(in) :: NX, NY, NZ

   integer          :: iatom, ifunct, i_ind, jatom, jfunct, j_ind, l1, l2, l3, &
                       l4, lij, lk, M2, MM, nnx, nny, nnz, ntotal, ns, np, nd, &
                       nci, ncj, rho_ind

   double precision :: Q(3), xi(3)
   double precision :: V_nuc, tna, temp, SQ3, Z2, Zij, ccoef, rExp, uf, t1, t2,&
                       f1, f2, s0s, s1s, s2s, s3s, s4s, pj0s, pj0p, pj1s, pj1p,&
                       p0s, p1s, p2s, p3s, PI0p, PI1p, pj2s, d0s, d0p, d1s,    &
                       d1p, d2s

   double precision, allocatable :: pote(:)
   double precision, parameter   :: RMAX = 3.0D0
   integer         , parameter   :: LL(3) = (/0, 1, 3/) ! LL is l * (l -1) / 2

   ! Variable initialization.
   SQ3 = 1.0D0
   if (norm) SQ3 = dsqrt(3.D0)

   ns = nShell(0); np = nShell(1); nd = nShell(2)
   MM = M * (M +1) / 2 ; M2 = 2 * M

   ! Calculates squared distance.
   do iatom = 1, natom
   do jatom = 1, natom
      d(iatom,jatom) = (r(iatom,1) - r(jatom,1)) * (r(iatom,1) - r(jatom,1)) + &
                       (r(iatom,2) - r(jatom,2)) * (r(iatom,2) - r(jatom,2)) + &
                       (r(iatom,3) - r(jatom,3)) * (r(iatom,3) - r(jatom,3))
   enddo
   enddo

   ntotal = NX * NY * NZ
   allocate(pote(ntotal))
   pote = 0.0D0

   ! First loop - (s|s)
   do ifunct = 1, ns
   do jfunct = 1, ifunct
      do nci = 1, nCont(ifunct)
      do ncj = 1, nCont(jfunct)
         Zij  = a(ifunct,nci) + a(jfunct,ncj)
         rExp = d(nuc(ifunct),nuc(jfunct)) * a(ifunct,nci) * a(jfunct,ncj) / Zij

         if (rExp .lt. RMAX) then
            t1   = a(ifunct,nci) / Zij
            t2   = a(jfunct,ncj) / Zij
            Q(1) = t1 * r(nuc(ifunct),1) + t2 * r(nuc(jfunct),1)
            Q(2) = t1 * r(nuc(ifunct),2) + t2 * r(nuc(jfunct),2)
            Q(3) = t1 * r(nuc(ifunct),3) + t2 * r(nuc(jfunct),3)

            ccoef   = c(ifunct,nci) * c(jfunct,ncj)
            temp    = 2.0D0 * PI * exp(-rExp) / Zij
            rho_ind = ifunct + (M2 - jfunct) * (jfunct-1) / 2

            ntotal = 0
            do nnx = 1, NX
            do nny = 1, NY
            do nnz = 1, NZ
               ntotal = ntotal +1

               xi(1) = xMin + (nnx -1) * deltax
               xi(2) = yMin + (nny -1) * deltax
               xi(3) = zMin + (nnz -1) * deltax
               uf    = ((Q(1) - xi(1)) * (Q(1) - xi(1))  + &
                        (Q(2) - xi(2)) * (Q(2) - xi(2))  + &
                        (Q(3) - xi(3)) * (Q(3) - xi(3))) * Zij
               s0s   = temp * FUNCT(0,uf)

               pote(ntotal) = pote(ntotal) + ccoef * RMM(rho_ind) * s0s
            enddo
            enddo
            enddo
         endif
      enddo
      enddo
   enddo
   enddo

   ! Second loop - (p|s)
   do ifunct = ns +1, ns + np, 3
   do jfunct = 1, ns
      do nci = 1, nCont(ifunct)
      do ncj = 1, nCont(jfunct)
         Zij  = a(ifunct,nci) + a(jfunct,ncj)
         rExp = d(nuc(ifunct),nuc(jfunct)) * a(ifunct,nci) * a(jfunct,ncj) / Zij

         if (rExp .lt. RMAX) then
            t1   = a(ifunct,nci) / Zij
            t2   = a(jfunct,ncj) / Zij
            Q(1) = t1 * r(nuc(ifunct),1) + t2 * r(nuc(jfunct),1)
            Q(2) = t1 * r(nuc(ifunct),2) + t2 * r(nuc(jfunct),2)
            Q(3) = t1 * r(nuc(ifunct),3) + t2 * r(nuc(jfunct),3)

            ccoef   = c(ifunct,nci) * c(jfunct,ncj)
            temp    = 2.0D0 * PI * exp(-rExp) / Zij

            ntotal = 0
            do nnx = 1, NX
            do nny = 1, NY
            do nnz = 1, NZ
               ntotal = ntotal +1

               xi(1) = xMin + (nnx -1) * deltax
               xi(2) = yMin + (nny -1) * deltax
               xi(3) = zMin + (nnz -1) * deltax

               uf = ((Q(1) - xi(1)) * (Q(1) - xi(1))  + &
                     (Q(2) - xi(2)) * (Q(2) - xi(2))  + &
                     (Q(3) - xi(3)) * (Q(3) - xi(3))) * Zij
               s0s = temp * FUNCT(0,uf)
               s1s = temp * FUNCT(1,uf)

               do l2 = 1, 3
                  tna     = (Q(l2) - r(nuc(ifunct),l2)) * s0s - &
                            (Q(l2) - xi(l2)           ) * s1s
                  rho_ind = ifunct + l2 -1 + ((M2 - jfunct) * (jfunct -1)) / 2

                  pote(ntotal) = pote(ntotal) + ccoef * tna * RMM(rho_ind)
               enddo
            enddo
            enddo
            enddo
         endif
      enddo
      enddo
   enddo
   enddo

   ! Third loop - (p|p)
   do ifunct = ns+1, ns + np, 3
   do jfunct = ns+1, ifunct , 3
      do nci = 1, nCont(ifunct)
      do ncj = 1, nCont(jfunct)

         Zij  = a(ifunct,nci) + a(jfunct,ncj)
         rExp = d(nuc(ifunct),nuc(jfunct)) * a(ifunct,nci) * a(jfunct,ncj) / Zij

         if (rExp .lt. RMAX) then
            Z2   = 2.0D0 * Zij
            t1   = a(ifunct,nci) / Zij
            t2   = a(jfunct,ncj) / Zij
            Q(1) = t1 * r(nuc(ifunct),1) + t2 * r(nuc(jfunct),1)
            Q(2) = t1 * r(nuc(ifunct),2) + t2 * r(nuc(jfunct),2)
            Q(3) = t1 * r(nuc(ifunct),3) + t2 * r(nuc(jfunct),3)

            ccoef   = c(ifunct,nci) * c(jfunct,ncj)
            temp    = 2.0D0 * PI * exp(-rExp) / Zij

            ntotal = 0
            do nnx = 1, NX
            do nny = 1, NY
            do nnz = 1, NZ
               ntotal = ntotal +1

               xi(1) = xMin + (nnx -1) * deltax
               xi(2) = yMin + (nny -1) * deltax
               xi(3) = zMin + (nnz -1) * deltax

               uf = ((Q(1) - xi(1)) * (Q(1) - xi(1))  + &
                     (Q(2) - xi(2)) * (Q(2) - xi(2))  + &
                     (Q(3) - xi(3)) * (Q(3) - xi(3))) * Zij

               s0s = temp * FUNCT(0,uf)
               s1s = temp * FUNCT(1,uf)
               s2s = temp * FUNCT(2,uf)

               do l1 = 1, 3
                  t1  = Q(l1) - r(nuc(ifunct),l1)
                  t2  = Q(l1) - xi(l1)
                  p0s = t1 * s0s - t2 * s1s
                  p1s = t1 * s1s - t2 * s2s

                  lij = 3
                  if (ifunct .eq. jfunct) lij = l1
                  do l2 = 1, lij
                     t1 = Q(l2) - r(nuc(jfunct),l2)
                     t2 = Q(l2) - xi(l2)

                     tna = (Q(l2) - r(nuc(jfunct),l2)) * p0s - &
                           (Q(l2) - xi(l2))            * p1s
                     if (l1 .eq. l2) tna = tna + (s0s - s1s) / Z2

                     i_ind   = ifunct + l1 -1
                     j_ind   = jfunct + l2 -1
                     rho_ind = i_ind + ((M2 - j_ind) * (j_ind -1)) / 2

                     pote(ntotal) = pote(ntotal) + tna * ccoef * RMM(rho_ind)
                  enddo
               enddo
            enddo
            enddo
            enddo
         endif
      enddo
      enddo
   enddo
   enddo

   ! Fourth loop - (d|s)
   do ifunct = ns + np +1, M, 6
   do jfunct = 1, ns
      do nci = 1, nCont(ifunct)
      do ncj = 1, nCont(jfunct)
         Zij  = a(ifunct,nci) + a(jfunct,ncj)
         rExp = d(nuc(ifunct),nuc(jfunct)) * a(ifunct,nci) * a(jfunct,ncj) / Zij

         if (rExp .lt. RMAX) then
            Z2   = 2.0D0 * Zij
            t1   = a(ifunct,nci) / Zij
            t2   = a(jfunct,ncj) / Zij
            Q(1) = t1 * r(nuc(ifunct),1) + t2 * r(nuc(jfunct),1)
            Q(2) = t1 * r(nuc(ifunct),2) + t2 * r(nuc(jfunct),2)
            Q(3) = t1 * r(nuc(ifunct),3) + t2 * r(nuc(jfunct),3)

            ccoef   = c(ifunct,nci) * c(jfunct,ncj)
            temp    = 2.0D0 * PI * exp(-rExp) / Zij

            ntotal = 0
            do nnx = 1, NX
            do nny = 1, NY
            do nnz = 1, NZ
               ntotal = ntotal + 1

               xi(1) = xMin + (nnx -1) * deltax
               xi(2) = yMin + (nny -1) * deltax
               xi(3) = zMin + (nnz -1) * deltax

               uf = ((Q(1) - xi(1)) * (Q(1) - xi(1))  + &
                     (Q(2) - xi(2)) * (Q(2) - xi(2))  + &
                     (Q(3) - xi(3)) * (Q(3) - xi(3))) * Zij

               s0s = temp * FUNCT(0,uf)
               s1s = temp * FUNCT(1,uf)
               s2s = temp * FUNCT(2,uf)

               do l1 = 1, 3
                  t1  = Q(l1) - r(nuc(ifunct),l1)
                  t2  = Q(l1) - xi(l1)
                  p0s = t1 * s0s - t2 * s1s
                  p1s = t1 * s1s - t2 * s2s

                  do l2 = 1, l1
                     t1  = Q(l2) - r(nuc(ifunct),l2)
                     t2  = Q(l2) - xi(l2)
                     tna = t1 * p0s - t2 * p1s

                     f1 = 1.0D0
                     if (l1 .eq. l2) then
                        tna = tna + (s0s - s1s)   / Z2
                        f1  = SQ3
                     endif

                     i_ind   = ifunct + ll(l1) + l2 -1
                     rho_ind = i_ind + ((M2 - jfunct) * (jfunct -1)) / 2

                     pote(ntotal) = pote(ntotal) + tna * RMM(rho_ind) * ccoef/f1
                  enddo
               enddo
            enddo
            enddo
            enddo
         endif
      enddo
      enddo
   enddo
   enddo

   ! Fifth loop - (d|p)
   do ifunct = ns + np +1, M      , 6
   do jfunct = ns +1     , ns + np, 3
      do nci = 1, nCont(ifunct)
      do ncj = 1, nCont(jfunct)
         Zij  = a(ifunct,nci) + a(jfunct,ncj)
         rExp = d(nuc(ifunct),nuc(jfunct)) * a(ifunct,nci) * a(jfunct,ncj) / Zij

         if (rExp .lt. RMAX) then
            Z2   = 2.0D0 * Zij
            t1   = a(ifunct,nci) / Zij
            t2   = a(jfunct,ncj) / Zij
            Q(1) = t1 * r(nuc(ifunct),1) + t2 * r(nuc(jfunct),1)
            Q(2) = t1 * r(nuc(ifunct),2) + t2 * r(nuc(jfunct),2)
            Q(3) = t1 * r(nuc(ifunct),3) + t2 * r(nuc(jfunct),3)

            ccoef   = c(ifunct,nci) * c(jfunct,ncj)
            temp    = 2.0D0 * PI * exp(-rExp) / Zij

            ntotal = 0
            do nnx = 1, NX
            do nny = 1, NY
            do nnz = 1, NZ
               ntotal = ntotal + 1

               xi(1) = xMin + (nnx -1) * deltax
               xi(2) = yMin + (nny -1) * deltax
               xi(3) = zMin + (nnz -1) * deltax

               uf = ((Q(1) - xi(1)) * (Q(1) - xi(1))  + &
                     (Q(2) - xi(2)) * (Q(2) - xi(2))  + &
                     (Q(3) - xi(3)) * (Q(3) - xi(3))) * Zij

               s0s = temp * FUNCT(0,uf)
               s1s = temp * FUNCT(1,uf)
               s2s = temp * FUNCT(2,uf)
               s3s = temp * FUNCT(3,uf)

               do l1 = 1, 3
                  t1  = Q(l1) - r(nuc(ifunct),l1)
                  t2  = Q(l1) - xi(l1)
                  p0s = t1 * s0s - t2 * s1s
                  p1s = t1 * s1s - t2 * s2s
                  p2s = t1 * s2s - t2 * s3s

                  do l2 = 1, l1
                     t1   = Q(l2) - r(nuc(ifunct),l2)
                     t2   = Q(l2) - xi(l2)
                     pj0s = t1 * s0s - t2 * s1s
                     pj1s = t1 * s1s - t2 * s2s
                     d0s  = t1 * p0s - t2 * p1s
                     d1s  = t1 * p1s - t2 * p2s
                     f1   = 1.0D0

                     if (l1 .eq. l2) then
                        f1  = SQ3
                        d0s = d0s + (s0s - s1s) / Z2
                        d1s = d1s + (s1s - s2s) / Z2
                     endif

                     do l3 = 1, 3
                        t1  = Q(l3) - r(nuc(jfunct),l3)
                        t2  = Q(l3) - xi(l3)
                        tna = t1 * d0s - t2 * d1s

                        if (l1 .eq. l3) tna = tna + (pj0s - pj1s) / Z2
                        if (l2 .eq. l3) tna = tna + (p0s  - p1s ) / Z2

                        i_ind   = ifunct + ll(l1) + l2 -1
                        j_ind   = jfunct + l3 -1
                        rho_ind = i_ind + ((M2 - j_ind) * (j_ind -1)) / 2

                        pote(ntotal) = pote(ntotal) + tna * RMM(rho_ind) * &
                                       ccoef / f1
                     enddo
                  enddo
               enddo
            enddo
            enddo
            enddo
         endif
      enddo
      enddo
   enddo
   enddo

   ! Sixth loop - (d|d)
   do ifunct = ns + np +1, M     , 6
   do jfunct = ns + np +1, ifunct, 6
      do nci = 1, nCont(ifunct)
      do ncj = 1, nCont(jfunct)
         Zij  = a(ifunct,nci) + a(jfunct,ncj)
         rExp = d(nuc(ifunct),nuc(jfunct)) * a(ifunct,nci) * a(jfunct,ncj) / Zij

         if (rExp .lt. RMAX) then
            Z2   = 2.0D0 * Zij
            t1   = a(ifunct,nci) / Zij
            t2   = a(jfunct,ncj) / Zij
            Q(1) = t1 * r(nuc(ifunct),1) + t2 * r(nuc(jfunct),1)
            Q(2) = t1 * r(nuc(ifunct),2) + t2 * r(nuc(jfunct),2)
            Q(3) = t1 * r(nuc(ifunct),3) + t2 * r(nuc(jfunct),3)

            ccoef   = c(ifunct,nci) * c(jfunct,ncj)
            temp    = 2.0D0 * PI * exp(-rExp) / Zij

            ntotal  = 0
            do nnx = 1, NX
            do nny = 1, NY
            do nnz = 1, NZ
               ntotal = ntotal + 1

               xi(1) = xMin + (nnx -1) * deltax
               xi(2) = yMin + (nny -1) * deltax
               xi(3) = zMin + (nnz -1) * deltax

               uf = ((Q(1) - xi(1)) * (Q(1) - xi(1))  + &
                     (Q(2) - xi(2)) * (Q(2) - xi(2))  + &
                     (Q(3) - xi(3)) * (Q(3) - xi(3))) * Zij

               s0s = temp * FUNCT(0,uf)
               s1s = temp * FUNCT(1,uf)
               s2s = temp * FUNCT(2,uf)
               s3s = temp * FUNCT(3,uf)
               s4s = temp * FUNCT(4,uf)

               do l1 = 1, 3
                  t1  = Q(l1) - r(nuc(ifunct),l1)
                  t2  = Q(l1) - xi(l1)
                  p0s = t1 * s0s - t2 * s1s
                  p1s = t1 * s1s - t2 * s2s
                  p2s = t1 * s2s - t2 * s3s
                  p3s = t1 * s3s - t2 * s4s

                  do l2 = 1, l1
                     t1   = Q(l2) - r(nuc(ifunct),l2)
                     t2   = Q(l2) - xi(l2)
                     pj0s = t1 * s0s - t2 * s1s
                     pj1s = t1 * s1s - t2 * s2s
                     pj2s = t1 * s2s - t2 * s3s
                     d0s  = t1 * p0s - t2 * p1s
                     d1s  = t1 * p1s - t2 * p2s
                     d2s  = t1 * p2s - t2 * p3s
                     f1   = 1.D0

                     if (l1 .eq. l2) then
                        f1  = SQ3
                        d0s = d0s + (s0s - s1s) / Z2
                        d1s = d1s + (s1s - s2s) / Z2
                        d2s = d2s + (s2s - s3s) / Z2
                     endif

                     lij = 3
                     if (ifunct .eq. jfunct) lij = l1

                     do l3 = 1, lij
                        t1   = Q(l3) - r(nuc(jfunct),l3)
                        t2   = Q(l3) - xi(l3)
                        d0p  = t1 * d0s - t2 * d1s
                        d1p  = t1 * d1s - t2 * d2s
                        PI0p = t1 * p0s - t2 * p1s
                        PI1p = t1 * p1s - t2 * p2s
                        pj0p = t1 * pj0s - t2 * pj1s
                        pj1p = t1 * pj1s - t2 * pj2s

                        if (l1 .eq. l3) then
                           d0p  = d0p  + (pj0s - pj1s) / Z2
                           d1p  = d1p  + (pj1s - pj2s) / Z2
                           PI0p = PI0p + (s0s  - s1s ) / Z2
                           PI1p = PI1p + (s1s  - s2s ) / Z2
                        endif
                        if (l2.eq.l3) then
                           d0p  = d0p  + (p0s - p1s) / Z2
                           d1p  = d1p  + (p1s - p2s) / Z2
                           pj0p = pj0p + (s0s - s1s) / Z2
                           pj1p = pj1p + (s1s - s2s) / Z2
                        endif

                        lk = l3
                        if (ifunct .eq. jfunct) lk = min(l3, Ll(l1) - Ll(l3)+l2)
                        do l4 = 1, lk
                           f2 = 1.0D0
                           tna = (Q(l4) - R(nuc(jfunct),l4)) * d0p - &
                                 (Q(l4) - xi(l4)           ) * d1p

                           if (l4 .eq. l1) tna = tna + (pj0p - pj1p) / Z2
                           if (l4 .eq. l2) tna = tna + (PI0p - PI1p) / Z2
                           if (l4 .eq. l3) then
                              f2  = SQ3
                              tna = tna + (d0s - d1s) / Z2
                           endif

                           i_ind   = ifunct + ll(l1) + l2 -1
                           j_ind   = jfunct + ll(l3) + l4 -1
                           rho_ind = i_ind + ((M2 - j_ind) * (j_ind -1)) / 2

                           pote(ntotal) = pote(ntotal) + RMM(rho_ind) * tna * &
                                          ccoef / (f1 * f2)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
            enddo
            enddo
         endif
      enddo
      enddo
   enddo
   enddo

   ! Electron potential is finished, now core potential is calculated.
   ntotal = 0
   do nnx = 1, NX
   do nny = 1, NY
   do nnz = 1, NZ
      ntotal = ntotal + 1

      xi(1) = xMin + (nnx -1) * deltax
      xi(2) = yMin + (nny -1) * deltax
      xi(3) = zMin + (nnz -1) * deltax

      V_nuc = 0.0D0
      do iatom = 1, natom
         V_nuc = V_nuc - Iz(iatom) / &
                 sqrt( (xi(1) - r(iatom,1)) * (xi(1) - r(iatom,1)) + &
                       (xi(2) - r(iatom,2)) * (xi(2) - r(iatom,2)) + &
                       (xi(3) - r(iatom,3)) * (xi(3) - r(iatom,3)) )
      enddo
      pote(ntotal) = pote(ntotal) + V_nuc
   enddo
   enddo
   enddo

   ! Writes a cubefile similar to those created by GAUSSIAN or ORCA, with the
   ! following format (as of GAUSSIAN98); all coordinates are in atomic units:
   !     LINE   FORMAT      CONTENTS
   !   ===============================================================
   !   1     A           TITLE
   !   2     A           DESCRIPTION OF PROPERTY STORED IN CUBEFILE
   !   3     I5,3F12.6   #ATOMS, X-,Y-,Z-COORDINATES OF ORIGIN
   !   4-6   I5,3F12.6   #GRIDPOINTS, INCREMENT VECTOR
   !   #ATOMS LINES OF ATOM COORDINATES:
   !   ...   I5,4F12.6   ATOM NUMBER, CHARGE, X-,Y-,Z-COORDINATE
   !   REST: 6E13.5      CUBE DATA (WITH Z INCREMENT MOVING FASTEST, THEN
   !                     Y AND THEN X)
   !
   ! For orbital cube files, #ATOMS will be negative (< 0) with one additional
   ! line after the final atom containing the number of orbitales and their
   ! respective numbers. The orbital number will also be the fastest moving
   ! increment.
   open(unit = 188, file = cube_elec_file)

   write(188,*) 'Elfield'
   write(188,*)
   write(188,678) natom, xMin, yMin, zMin
   write(188,678) NX, deltax, 0.0D0 , 0.0D0
   write(188,678) NY, 0.0D0 , deltax, 0.0D0
   write(188,678) NZ, 0.0D0 , 0.0D0 , deltax
   do iatom = 1, natom
      write(188,677) Iz(iatom), 0.0D0, r(iatom,1), r(iatom,2), r(iatom,3)
   enddo
   write(188,679) pote

   close(188)
   deallocate(pote)

677 format(I5,4(F12.6))
678 format(I5,3(F12.6))
679 format(6(E13.5))
end subroutine elec
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
end module cubegen
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
