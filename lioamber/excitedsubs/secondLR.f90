subroutine get_perturbed(C1,C0,X,M,Nv,No,Nd,Ns)
! Inputs:
! C0 : Unperturbed MOs
! X  : Transition Vectors from First LR
! M  : Number of basis
! Nv : Number of virtual MOs
! No : Number of occupied MOs
! Nd : Dimension of excitations
! Ns : Number of states

! Outputs:
! C1 : Perturbed MOs
use garcha_mod  , only: Pmat_vec
use excited_data, only: lambda_LR, state_LR, read_other
use SCF_aux,      only: messup_densmat
   implicit none
   integer, intent(in) :: M, Nv, No, Nd, Ns
   LIODBLE, intent(in) :: C0(:,:), X(Nd,Ns)
   LIODBLE, intent(out):: C1(M,M)

   integer :: NCOc, ii, jj, ele
   LIODBLE :: EE
   LIODBLE, allocatable :: temp(:), Xp(:)
   LIODBLE, allocatable :: D1(:,:)

   allocate(temp(M),Xp(Nd))

   ! Select which state is used in the perturbation
   if ( read_other ) then
      ! aqui pongo el estado que lee
      print*, "reading other state"
      open(unit=456,file="vectors_root")
      read(456,*) EE
      do ii=1,Nd
          read(456,*) Xp(ii)
      enddo
      close(456)
   else
      Xp = X(:,state_LR)
   endif

   ! Form Perturbed MOs: C1
   NCOc = No + 1
   do ii=1,No ! For occupied
      temp = 0.0d0
      do jj=1,Nv
         ele = (ii - 1) * Nv + jj
         temp = temp + Xp(ele) * C0(:,No+jj)
      enddo
      C1(:,NCOc-ii) = C0(:,NCOc-ii) + lambda_LR * temp
   enddo
   do ii=1,Nv ! For virtuals
      temp = 0.0d0
      do jj=1,No
         ele = (jj-1) * Nv + ii
         temp = temp + Xp(ele) * C0(:,NCOc-jj)
      enddo
      C1(:,No+ii) = C0(:,No+ii) - lambda_LR * temp
    enddo
   deallocate(temp,Xp)

   ! Form Perturbed Density of the interested state: D1
   allocate(D1(M,M)); D1 = 0.0d0
   call builds_densmat(M,No,2.0d0,C1,D1)
   call messup_densmat(D1)
   call sprepack('L',M,Pmat_vec,D1)
   deallocate(D1)

end subroutine get_perturbed

subroutine excited_absorption(Xflr,Eflr,Xslr,C1,Dip0,M,No,Nv,Nd,Ns)
! Inputs
! Xflr : Transition vectors from First LR
! Eflr : Excitation energies from First LR
! Xslr : Transition vectors from Second LR
! C1   : Perturbed MOs
! Dip0 : Unperturbed Transition Dipole
! No   : Number of occupied MOs
! Nv   : Number of virtual MOs
! Nd   : Dimension of excitations
! Ns   : Number of states
use excited_data, only: Ctol, Tdip_save, lambda_LR, state_LR, read_other
   implicit none
   
   integer, intent(in)    :: M, No, Nv, Nd, Ns
   LIODBLE, intent(in)    :: Xflr(Nd,Ns), Eflr(Ns), Dip0(Ns,3), C1(M,M)
   LIODBLE, intent(inout) :: Xslr(Nd,Ns)

   logical :: asigned
   integer :: ii, jj, kk, Ns_slr, ele, id
   LIODBLE :: temp, scal, EE
   LIODBLE, allocatable :: Ene(:), OsSt(:), ovlap(:), ss(:)
   LIODBLE, allocatable :: Xtemp(:,:), TD(:,:)

   write(*,"(1X,A)") "EXCITED STATE ABSORPTION"
   write(*,"(1X,A,2X,F8.4)") "Using lambda in SLR:", lambda_LR
   write(*,"(1X,A,2X,I3)") "Using state  in SLR:", state_LR
   write(*,"(1X,A,10X,F8.4,1X,A)") "Using Ctol: ",Ctol,"to projection States"

   ! Get the correct order in Perturbed Transition Vectors
   allocate(Xtemp(Nd,Ns),ovlap(Ns),ss(Ns))
   do ii=1,Ns !Perturbed
      asigned = .false.
      ovlap   = 0.0d0
      do jj=1,Ns !Unperturbed
         temp = 0.0d0
         do kk=1,Nd
            temp = temp + Xslr(kk,ii) * Xflr(kk,jj)
         enddo
         ovlap(jj) = abs(temp)
         ss(jj)    = temp
         if ( (abs(temp) > (1.0d0-Ctol)) .and. (abs(temp) < (1.0d0+Ctol)) ) then
            scal = 1.0d0
            asigned = .true.
            if ( temp < 0.0d0 ) scal = -1.d0
            if ( ii /= jj ) print*, "Switch: ", ii, "->", jj
            Xtemp(:,jj) = Xslr(:,ii) * scal
            exit ! exit of loop jj
         endif
      enddo
      if (.not. asigned ) then
         print*, "The state", ii, "was not asigned"
         write(*,"(1X,A,I2,A,I2,1X,A,F8.4)") "MAX: |X1(",ii,")*X0(",maxloc(ovlap),")^T|=",maxval(ovlap)
         scal = 1.0d0
         id = maxloc(ovlap,1)
         if ( ss(id) < 0.0d0 ) scal = -1.d0
         Xtemp(:,id) = Xslr(:,ii) * scal
      endif
   enddo
   Xslr = Xtemp
   deallocate(Xtemp,ovlap,ss)

   ! Get the perturbed transition dipole
   allocate(Ene(Ns),OsSt(Ns)); Ene = 0.0d0 ! Ene and OsSt are not considered here
   call OscStr(Xslr,Ene,C1,OsSt,M,M,No,Nv,Nd,Ns)
   deallocate(Ene,OsSt)
   
   ! Calculate Transition Dipole of the state
   Ns_slr = Ns - state_LR
   EE = Eflr(state_LR)
   if ( read_other ) then
      print*, "Energy from other state"
      Ns_slr = Ns
      state_LR = 0
      open(unit=456,file="vectors_root")
      read(456,*) EE
      close(456)
   endif
   allocate(TD(Ns_slr,3),Ene(Ns_slr),OsSt(Ns_slr))
   do ii=1, Ns_slr
       ele = state_LR + ii
       TD(ii,:) = ( Tdip_save(ele,:) - Dip0(ele,:) ) / lambda_LR
       Ene(ii)  = Eflr(ele) - EE
   enddo

   ! Oscillator Strenght Calculation and Print Results
   call PrintESA(TD,Ene,OsSt,Ns_slr,state_LR)
   deallocate(TD,Ene,OsSt)

end subroutine excited_absorption

subroutine builds_densmat( Msize, Nocup, Focup, coef_mat, dens_mat )
!
! Msize:      number of atomic basis functions.
! Nocup:      number of occupied orbitals.
! Focup:      orbital ocupation factor (2.0d0 for CS, 1.0d0 for OS)
! coef_mat:   matrix containing Fock coefficients.
! dens_mat:   matrix containing output density matrix.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   implicit none
   integer, intent(in)  :: Msize
   integer, intent(in)  :: Nocup
   LIODBLE , intent(in)  :: Focup
   LIODBLE , intent(in)  :: coef_mat(Msize, Msize)
   LIODBLE , intent(out) :: dens_mat(Msize, Msize)

   LIODBLE , allocatable :: coef_occ(:, :)
   integer               :: ii, jj

!  Copies the occupied orbitals into a truncated matrix
   if ( .not.allocated(coef_occ) ) allocate( coef_occ(Msize,Nocup) )
   do jj = 1, Nocup
   do ii = 1, Msize
      coef_occ(ii, jj) = coef_mat(ii, jj)
   enddo
   enddo

!  Obtains dens_mat as (coef_occ)*(coef_occ^t).
   dens_mat(:,:) = 0.0D0
   call DGEMM( 'N', 'T', Msize, Msize, Nocup, Focup, coef_occ, Msize, coef_occ,&
             & Msize, 0.0D0, dens_mat, Msize)

   deallocate(coef_occ)
end subroutine builds_densmat
