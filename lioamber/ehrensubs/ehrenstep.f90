!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine ehrenstep( propagator_id, dt, nbasis, natoms,                       &
                    & nucpos, nucvel, qm_forces_ds, Sinv, Uinv, Linv,          &
                    & dens_old, dens_mid, dens_new, fock_mid, energy, dipmom )
!------------------------------------------------------------------------------!
!
! All matrices must be in OM
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   implicit none
   integer, intent(in)       :: propagator_id
   real*8 , intent(in)       :: dt
   integer, intent(in)       :: nbasis
   integer, intent(in)       :: natoms

   real*8, intent(in)        :: nucpos(3, natoms)
   real*8, intent(in)        :: nucvel(3, natoms)
   real*8, intent(inout)     :: qm_forces_ds(3, natoms)

   real*8, intent(in)        :: Sinv(nbasis, nbasis)
   real*8, intent(in)        :: Uinv(nbasis, nbasis)
   real*8, intent(in)        :: Linv(nbasis, nbasis)

   complex*16, intent(in)    :: dens_old(nbasis, nbasis)
   complex*16, intent(in)    :: dens_mid(nbasis, nbasis)
   complex*16, intent(inout) :: dens_new(nbasis, nbasis)
   real*8    , intent(inout) :: fock_mid(nbasis, nbasis)
   real*8    , intent(inout) :: energy
   real*8    , intent(inout) :: dipmom(3)

   real*8    , allocatable   :: Bmat(:,:)
   real*8    , allocatable   :: Dmat(:,:)
   complex*16, allocatable   :: Tmat(:,:)
   complex*16, allocatable   :: dens_mao(:,:)

   integer, parameter :: propagator_id_verlet = 1
   integer, parameter :: propagator_id_magnus = 2

   allocate( Bmat(nbasis,nbasis), Dmat(nbasis,nbasis), Tmat(nbasis,nbasis) )
   allocate( dens_mao(nbasis,nbasis) )

!  Fock and Force calculation (needs density and fock in AO)
!  (this should leave the right Rho in RMM for get_forces)
   dens_mao = dens_mid
   dens_mao = matmul(dens_mao, Linv)
   dens_mao = matmul(Uinv, dens_mao)
   call RMMcalc3_FockMao( dens_mao, fock_mid, dipmom, energy)
   call calc_forceDS( natoms, nbasis, nucpos, nucvel, dens_mao, fock_mid, Sinv,&
                    & Bmat, qm_forces_ds )


!  Set ups propagation cuasi-fock matrix (needs fock in ON)
   fock_mid = matmul(fock_mid, Uinv)
   fock_mid = matmul(Linv, fock_mid)
   Dmat = calc_Dmat( nbasis, Linv, Uinv, Bmat )
   Tmat = DCMPLX(fock_mid) + DCMPLX(0.0d0,1.0d0) * DCMPLX(Dmat)


!  Density Propagation (works in ON)
   if (propagator_id==propagator_id_verlet) then
      call ehren_verlet( nbasis, dt, Tmat, dens_old, dens_mid, dens_new )

   else if (propagator_id==propagator_id_magnus) then
      call ehren_magnus( nbasis, 20, dt, Tmat, dens_old, dens_new )

   else
      print*,'Unidentified substep!'; stop

   endif

   deallocate( Bmat, Dmat, Tmat, dens_mao)
end subroutine ehrenstep
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
