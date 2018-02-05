!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!     CASE OF NO STARTING GUESS PROVIDED, 1 E FOCK MATRIX USED
!     FCe = SCe; (X^T)SX = 1
!     F' = (X^T)FX
!     => (X^-1*C)^-1 * F' * (X^-1*C) = e
!
! This subroutine takes the Hmat as a vector, transforms it into a matrix, 
! diagonalizes it, and builds the density from the resulting orbitals.
!
! TODO:
!
! (*) The starting guess should take care of calculating everything, it
!     shouldn't need to eat Hmat. I do understand the convenience since
!     it is the same for the whole SCF loop. Perhaps make it optional??
!
! (*) Diagonalizing fock and building the density matrix are things that
!     are or should be performed by dedicated subroutines, and these
!     should be called both here and in the SCF cycle.
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine starting_guess( Nmat, Nvec, NCO, hmat_vec, Xmat, densat_vec )
   use liosubs_math  , only: transform

   implicit none
   integer, intent(in)    :: Nmat
   integer, intent(in)    :: Nvec
   integer, intent(in)    :: NCO
   real*8 , intent(in)    :: Xmat(Nmat,Nmat)
   real*8 , intent(inout) :: hmat_vec(Nvec)
   real*8 , intent(inout) :: densat_vec(Nvec)
   
   real*8, allocatable   :: morb_energy(:)
   real*8, allocatable   :: morb_coefon(:,:)
   real*8, allocatable   :: morb_coefat(:,:)
   real*8, allocatable   :: morb_coefoc(:,:)
   real*8, allocatable   :: dens_mao(:,:)
   real*8, allocatable   :: hmat(:,:)
   real*8, allocatable   :: WORK(:)
   integer               :: LWORK
   integer               :: info

   call g2g_timer_start('initial guess')
   call g2g_timer_sum_start('initial guess')

   allocate( morb_energy(Nvec) )
   allocate( morb_coefon(Nmat, Nmat) )
   allocate( morb_coefat(Nmat, Nmat) )
   allocate( morb_coefoc(Nmat, NCO) )
   allocate( hmat(Nmat,Nmat) )
   allocate( dens_mao(Nmat, Nmat) )

   call spunpack('L', Nmat, hmat_vec, hmat )
   morb_coefon(:,:) = transform( hmat, Xmat )
   morb_energy(:)   = 0.0d0

   LWORK = -1
   if ( allocated(WORK) ) deallocate(WORK)
   allocate( WORK( 1 ) )
   call dsyev( 'V', 'L', Nmat, morb_coefon, Nmat, morb_energy, WORK, LWORK, info )

   LWORK = NINT( WORK(1) )
   if ( allocated(WORK) ) deallocate(WORK)
   allocate( WORK( LWORK ) )
   call dsyev( 'V', 'L', Nmat, morb_coefon, Nmat, morb_energy, WORK, LWORK, info )

   morb_coefat = matmul( Xmat, morb_coefon )
   morb_coefoc(1:Nmat,1:NCO) = morb_coefat(1:Nmat,1:NCO)
   dens_mao = (2.0d0) * matmul( morb_coefoc, transpose(morb_coefoc) )
   call messup_densmat( dens_mao )
   call sprepack( 'L', Nmat, densat_vec, dens_mao)
   
   call g2g_timer_stop('initial guess')
   call g2g_timer_sum_stop('initial guess')
end subroutine starting_guess
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
