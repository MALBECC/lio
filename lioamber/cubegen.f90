!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%%% CUBEGEN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Performs orbital, density and electrostatic potential surface plots. These   !
! Can be opened in VMD, and watched by selecting "isosurface" under the        !
! representation menu.                                                         !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module cubegen
   implicit none
contains

subroutine cubegen_vecin(coef_mat)
   use garcha_mod, only: cube_dens, cube_orb, cube_elec, cubegen_only, VCINP
   implicit none
   real(kind=8), intent(in) :: coef_mat(:,:)


   if ( .not. (cube_dens .or. cube_orb .or. cube_elec) ) return

   if ( cubegen_only .and. (.not.VCINP) ) then
      write(*,*) "cubegen_only CAN ONLY BE USED WITH VCINP"
      stop
   end if

   call g2g_timer_sum_start('cube gen')
   call cubegen_write( coef_mat )
   call g2g_timer_sum_stop('cube gen')
end subroutine cubegen_vecin

subroutine cubegen_matin( Msize, ugly_mat )
   use garcha_mod, only: cube_dens, cube_orb, cube_elec, cubegen_only, VCINP
   implicit none
   integer     , intent(in)  :: Msize
   real(kind=8), intent(in)  :: ugly_mat( Msize, 3*Msize )
   real(kind=8), allocatable :: coef_mat( :, : )
   integer                  :: ii, jj


   if ( .not. (cube_dens .or. cube_orb .or. cube_elec) ) return

   if ( cubegen_only .and. (.not. VCINP) ) then
      write(*,*) "cubegen_only CAN ONLY BE USED WITH VCINP"
      stop
   end if

   if (allocated(coef_mat)) deallocate(coef_mat)
   allocate( coef_mat(Msize, Msize) )

   do jj = 1, Msize
   do ii = 1, Msize
      coef_mat(ii, jj) = ugly_mat( ii, 2*Msize + jj )
   enddo
   enddo

   call g2g_timer_sum_start('cube gen')
   call cubegen_write( coef_mat )
   call g2g_timer_sum_stop('cube gen')

   deallocate( coef_mat )
end subroutine cubegen_matin

subroutine cubegen_write( MO_v )
   use garcha_mod, only: natom, r, nco, Iz,  cube_dens, cube_orb, &
                         cube_elec, cube_sel, cube_orb_file, cube_res, &
                         cube_dens_file, cube_sqrt_orb
   use basis_data, only: M, a, ncont, nuc, nshell

   implicit none

   real(kind=8), intent(in) :: MO_v(:,:)

   real(kind=8) :: x0(3), x1(3), origin(3), eval_p(3)
   real(kind=8) :: max_radius, max_dim, vox_dim, p_val, Morb

   integer      :: i, j, k, ii, jj, kk, ns, np, nd, ni
   integer      :: ivox, ivoxx, ivoxy, ivoxz, kk_dens, kk_orb

   integer      :: dim_array
   real(kind=8), allocatable :: min_exps(:), p_array(:)

   if ( .not. (cube_dens .or. cube_orb .or. cube_elec) ) return
   
   ns = nshell(0)
   np = nshell(1)
   nd = nshell(2)
   dim_array = ns + 3*np + 6*nd

   allocate(min_exps(120), p_array(dim_array))
   if (cube_dens) open(unit = 4242, file = cube_dens_file)
   if (cube_orb)  open(unit = 4243, file = cube_orb_file)

   ! First find the rectangular prism that the system lies in.
   ivox = cube_res
   min_exps = 999999.D0
   do ii=1,M
      jj = Iz(Nuc(ii))
      do ni = 1, ncont(ii)
         min_exps(jj) = min(min_exps(jj), a(ii,ni))
      enddo
   enddo

   x0 = 0.D0
   x1 = 0.D0
   do i = 1,natom
      ! TODO: not sure about this padding criteria
      max_radius = 2.0D0 * sqrt(1.0D0 / min_exps(Iz(i)))
      if (i == 1) then
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
   ! cube_res gives the../test/LIO_test/00_agua/orb.cube # of voxels in the longest dimension
   ! So, next figure ou../test/LIO_test/00_agua/orb.cubet the voxel length based on this criteria
   ! and how many voxel../test/LIO_test/00_agua/orb.cubes each dimension needs based on voxel length
   max_dim = -1.D0
   do i = 1, 3
      max_dim = max(max_dim, x1(i) - x0(i))
   enddo
   vox_dim = max_dim / ivox

   ivoxx = ceiling((x1(1) - x0(1)) / vox_dim)
   ivoxy = ceiling((x1(2) - x0(2)) / vox_dim)
   ivoxz = ceiling((x1(3) - x0(3)) / vox_dim)

   if (cube_dens) then
      write(4242,*) "LIO CUBE FILE"
      write(4242,*) "CUBE FORMAT FILE FOR ELECTRON DENSITY"
      write(4242,42) natom, origin(:)
      write(4242,42) ivoxx, vox_dim, 0.0D0  , 0.0D0
      write(4242,42) ivoxy, 0.0D0  , vox_dim, 0.0D0
      write(4242,42) ivoxz, 0.0D0  , 0.0D0  , vox_dim
      do i = 1, natom
         write(4242,424) Iz(i), 0.0D0, r(i,:)
      enddo
   endif

   if (cube_orb) then
      write(4243,*) "LIO CUBE FILE"
      if (cube_sel == 0) then
         write(4243,*) "CUBE FORMAT FILE FOR MOLECULAR ORBITALS"
      elseif ((cube_sel > 0) .and. (cube_sel <= NCO)) then
         write(4243,*) "CUBE FORMAT FILE FOR SINGLE MOLECULAR ORBITAL"
      else
         write(*,*) "cube_sel VALUE NOT SUPPORTED!"
         stop
      endif

      if (cube_sel == 0) then
         write(4243,42) -natom, origin(:)
      else
         write(4243,42) natom, origin(:)
      endif
      write(4243,42) ivoxx, vox_dim, 0.0D0   , 0.0D0
      write(4243,42) ivoxy, 0.0D0  , vox_dim , 0.0D0
      write(4243,42) ivoxz, 0.0D0  , 0.0D0   , vox_dim
      do i = 1,natom
        write(4243,424) Iz(i), 0.0D0, r(i,:)
      enddo

      if (cube_sel == 0) then
         write(4243,'(I4)',advance='no') M
         do i = 1, M
            write(4243,'(I4)',advance='no') i-NCO
         enddo
         write(4243,*) ""
      endif
   endif

   kk_dens = 1
   kk_orb  = 1
   ! For each voxel...
   do i = 1, ivoxx
      eval_p(1) = origin(1) + (i-1) * vox_dim
      do j = 1, ivoxy
         eval_p(2) = origin(2) + (j-1) * vox_dim
         do k = 1, ivoxz
            eval_p(3) = origin(3) + (k-1) * vox_dim
            p_val = 0.D0
           
            ! Calculate function values at this voxel, store in energy
            call evaluate_basis(eval_p(1),eval_p(2),eval_p(3),dim_array,p_array)

            if (cube_dens) then ! Calculate density for this voxel
               p_val = obtainrho(dim_array,p_array)
               write(4242, '(E13.5)', advance='no') p_val
               if (mod(kk_dens,6) == 0) write(4242,*) ""
               kk_dens = kk_dens + 1 
            endif

            if (cube_orb) then
               if (cube_sel == 0) then
                  do kk = 1, M
                  ! Calculate orbital or orbital^2 if cube_sqrt_orb=t
                     p_val = 0.D0
                     do ii = 1, M
                        do jj = ii+1, M
                           Morb = 2.D0 * MO_v(ii,kk) * p_array(ii)
                           if (cube_sqrt_orb) &
                              Morb = Morb * MO_v(jj,kk) * p_array(jj)
                           p_val = p_val + Morb
                        enddo

                        Morb = MO_v(ii,kk) * p_array(ii)
                        if (cube_sqrt_orb) &
                           Morb = Morb * MO_v(ii,kk) * p_array(ii)
                        p_val = p_val + Morb
                     enddo

                     write(4243, '(E13.5)', advance = 'no') p_val
                     if (mod(kk_orb, 6) == 0) write(4243,*) ""
                     kk_orb = kk_orb + 1
                  enddo
               else
                  ! Calculate orbital or orbital^2 if cube_sqrt_orb=t
                  p_val = 0.D0
                  do ii = 1, M
                     do jj = ii+1,M
                        Morb = 2.D0 * MO_v(ii,cube_sel) * p_array(ii)
                        if (cube_sqrt_orb) &
                           Morb = Morb * MO_v(jj,cube_sel) * p_array(jj)
                        p_val = p_val+Morb
                     enddo
                     Morb = MO_v(ii,cube_sel) * p_array(ii)
                     if (cube_sqrt_orb) &
                        Morb  = Morb * MO_v(ii,cube_sel) * p_array(ii)
                     p_val = p_val + Morb
                  enddo
                  write(4243, '(E13.5)', advance = 'no') p_val
                  if (mod(kk_orb,6) == 0) then
                     write(4243,*) ""
                  endif
                  kk_orb = kk_orb + 1
               endif
            endif
         enddo
      enddo
   enddo

   if (cube_dens) close(4242)
   if (cube_orb)  close(4243)

   if (cube_elec) then
      call elec(ivoxx, ivoxy, ivoxz, vox_dim, origin(1), origin(2), origin(3))
   endif

   deallocate(min_exps, p_array)

  42 format(I5,3(f12.6))
 424 format(I5,4(f12.6))

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
   use garcha_mod   , only: r, d, natom, cube_elec_file, Pmat_vec, Iz
   use constants_mod, only: PI
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

               pote(ntotal) = pote(ntotal) + ccoef * Pmat_vec(rho_ind) * s0s
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

                  pote(ntotal) = pote(ntotal) + ccoef * tna * Pmat_vec(rho_ind)
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

                     pote(ntotal) = pote(ntotal) + tna * ccoef * &
                                    Pmat_vec(rho_ind)
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

                     pote(ntotal) = pote(ntotal) + tna * Pmat_vec(rho_ind) *&
                                    ccoef / f1
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

                        pote(ntotal) = pote(ntotal) + tna * Pmat_vec(rho_ind)*&
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

                           pote(ntotal) = pote(ntotal) + Pmat_vec(rho_ind) * &
                                          tna * ccoef / (f1 * f2)
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


subroutine integrate_rho()
!Computes the integral of electronic density in 2 selected directions xy, xz, yz or thetaphi
   use rhoint
   use constants, only : bohr, pi
   use basis_data, only: nshell
   implicit none
   double precision :: xmin, xmax, ymin, ymax, zmin, zmax, rmin, rmax
   double precision :: dx,dy,dz,dr,dtheta,dphi
   double precision :: integral
   double precision :: x,y,z,r, phi, theta
   integer :: steps_x, steps_y, steps_z, steps_r, steps_theta, steps_phi
   integer :: ix, iy, iz, ir, itheta, iphi
   integer :: ns, np, nd, dim_array
   real(kind=8), allocatable :: p_array(:)

   xmin=w_rho_xmin/bohr
   xmax=w_rho_xmax/bohr
   ymin=w_rho_ymin/bohr
   ymax=w_rho_ymax/bohr
   zmin=w_rho_zmin/bohr
   zmax=w_rho_zmax/bohr
   rmin=w_rho_rmin/bohr
   rmax=w_rho_rmax/bohr

   dx=w_rho_dx/bohr
   dy=w_rho_dy/bohr
   dz=w_rho_dz/bohr
   dr=w_rho_dr/bohr
   dtheta=w_rho_dtheta
   dphi=w_rho_dphi


   steps_x=int((xmax-xmin)/dx)
   steps_y=int((ymax-ymin)/dy)
   steps_z=int((zmax-zmin)/dz)
   steps_r=int((rmax-rmin)/dr)
   steps_theta=int(pi/dtheta)
   steps_phi=int(2.d0*pi/dphi)

   ns = nshell(0)
   np = nshell(1)
   nd = nshell(2)
   dim_array = ns + 3*np + 6*nd
   allocate(p_array(dim_array))


   if (write_int_rho == 'z') then

      open(unit=978, file='rho_z.dat')
      z=zmin-dz
      do iz=0, steps_z
         integral=0.d0
         z=z+dz
         x=xmin-dx*0.5d0
         do ix=0, steps_x
            x=x+dx
            y=ymin-dy*0.5d0
            do iy=0, steps_y
               y=y+dy
               call evaluate_basis(x,y,z,dim_array,p_array)
               integral=integral+obtainrho(dim_array,p_array)*dx*dy
            end do
         end do
         write(978,*) z*bohr, integral/bohr !position in Ang, rho in e/Ang
      end do
      close(978)

   elseif (write_int_rho == 'y') then

      open(unit=978, file='rho_y.dat')
      y=ymin-dy
      do iy=0, steps_y
         integral=0.d0
         y=y+dy
         x=xmin-dx*0.5d0
         do ix=0, steps_x
            x=x+dx
            z=zmin-dz*0.5d0
            do iz=0, steps_z
               z=z+dz
               call evaluate_basis(x,y,z,dim_array,p_array)
               integral=integral+obtainrho(dim_array,p_array)*dx*dz
            end do
         end do
         write(978,*) y*bohr, integral/bohr !position in Ang, rho in e/Ang
      end do
      close(978)

   elseif (write_int_rho == 'x') then

      open(unit=978, file='rho_x.dat')
      x=xmin-dx
      do ix=0, steps_x
         integral=0.d0
         x=x+dx
         z=zmin-dz*0.5d0
         do iz=0, steps_z
            z=z+dz
            y=ymin-dy*0.5d0
            do iy=0, steps_y
               y=y+dy
               call evaluate_basis(x,y,z,dim_array,p_array)
               integral=integral+obtainrho(dim_array,p_array)*dz*dy
            end do
         end do
         write(978,*) x*bohr, integral/bohr !position in Ang, rho in e/Ang
      end do
      close(978)

   elseif (write_int_rho == 'r') then
      open(unit=978, file='rho_r.dat')
      r = rmin-dr
      do ir=0, steps_r
         integral=0.d0
         r=r+dr
         theta=0.d0!-dtheta*0.5
         do itheta=0, steps_theta
            theta=theta+dtheta
            phi=-dphi*0.5d0
            do iphi=0, steps_phi
               phi=phi+dphi
               x=r*dsin(theta)*dcos(phi)
               y=r*dsin(theta)*dsin(phi)
               z=r*dcos(theta)
               call evaluate_basis(x,y,z,dim_array,p_array)
               integral=integral+obtainrho(dim_array,p_array)*dtheta*dphi*dsin(theta)
            end do
         end do
         write(978,*) r*bohr, integral/bohr !position in Ang, rho in e/Ang
      end do
      close(978)
   end if

deallocate(p_array)
end subroutine integrate_rho


subroutine evaluate_basis(x,y,z,dim_array,p_array)
! Computes function values of the whole basis set in (x,y,z) position, store in p_array
   use garcha_mod, only: r
   use basis_data, only: M, ncont, nuc, nshell, a, c
   implicit none
   double precision, intent(in) :: x,y,z
   integer, intent(in) :: dim_array
   double precision, intent(out), dimension(dim_array) :: p_array
   double precision, dimension(3) :: eval_p
   double precision, parameter   :: expmax = 10.0D0
   integer :: ns, np, nd
   integer :: ii, jj, ni, jjj, kkk
   double precision :: p_val, p_func, p_dist

   ns = nshell(0)
   np = nshell(1)
   nd = nshell(2)
   eval_p(1)=x
   eval_p(2)=y
   eval_p(3)=z
   p_val = 0.D0
   p_array = 0.d0

! s functions
   do ii = 1, ns
      p_dist = 0.D0
      do jj = 1, 3
         p_dist = p_dist + (eval_p(jj) - r(Nuc(ii),jj))**2
      enddo
   
      p_func = 0.D0
      do ni = 1, ncont(ii)
         if ((a(ii,ni)*p_dist) < expmax) &
            p_func = p_func + c(ii,ni) * exp(-a(ii,ni) * p_dist)
      enddo
      p_array(ii) = p_func
   enddo

! p functions
   do ii = ns+1, ns+np, 3
      p_dist = 0.D0
      do jj = 1, 3
         p_dist = p_dist + (eval_p(jj) - r(Nuc(ii),jj))**2
      enddo

      p_func = 0.D0
      do ni = 1, ncont(ii)
         if ((a(ii,ni)*p_dist) < expmax) &
            p_func = p_func + c(ii,ni) * exp(-a(ii,ni) * p_dist)
      enddo

      do jj = 1,3
         p_array(ii+jj-1) = p_func * (eval_p(jj) - r(Nuc(ii),jj))
      enddo
   enddo

! d functions
   do ii = ns+np+1, M, 6
      p_dist = 0.D0
      do jj = 1,3
         p_dist = p_dist + (eval_p(jj) - r(Nuc(ii),jj))**2
      enddo

      p_func = 0.D0
      do ni = 1, ncont(ii)
         if ((a(ii,ni) * p_dist) < expmax) &
            p_func = p_func + c(ii,ni) * exp(-a(ii,ni) * p_dist)
      enddo

      kkk = 0
      do jj = 1, 3
         do jjj = 1, jj
            kkk = kkk + 1
            p_array(ii+kkk-1) = p_func * (eval_p(jj)  - r(Nuc(ii),jj)) *&
                                   (eval_p(jjj) - r(Nuc(ii),jjj))
         enddo
      enddo
   enddo

end subroutine evaluate_basis


double precision function obtainrho(dim_array,p_array)
! Calculate density value
   use garcha_mod, only: Pmat_vec
   use basis_data, only: M
   implicit none
   double precision :: p_val
   integer, intent(in) :: dim_array
   double precision, intent(in) :: p_array(dim_array)
   integer :: ii, jj, kkk

   p_val=0.d0
   kkk = 0
   do ii = 1 , M
      do jj = ii, M
         kkk   = kkk + 1
         p_val = p_val + Pmat_vec(kkk) * p_array(ii) * p_array(jj)
      enddo
   enddo
   obtainrho = p_val
   return
end function obtainrho



end module cubegen
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
