!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Lineal search subroutines
!
! This procedure improve convergency in systems finding best lineal
! combinations of density matrix in steps n and n-1. To do so, it
! evaluates different combinations of rho(n) and rho(n-1), using the
! variable lambda as a weight. The possibilities it evaluates also
! depend on the performance of the algorithm in previous steps,
! which is regulated by the parameters Elast and Pstepsize.
!
! FFR comments: Rho_LS, changed_to_LS and P_oscilation_analysis are
! communication variables that appear here and in some external subs.
! P_oscilation_analysis and changed_to_LS are not even used inside here.
! These should be dealt with differently.
!
!------------------------------------------------------------------------------!
! LOG:
!
! V 1.00 September 2018 Final version - Nicolas Foglia
!
! V 1.01 September 2018 adaptation - FFR
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module converger_ls
   implicit none

!------------------------------!
!  Linear search in rho for convergention:
   integer          :: Rho_LS = 0
!  IF Rho_LS=0 calculate convergence criteria for actual density matrix
!  IF Rho_LS=1 do a lineal search for density matrix only if energy > energy of previus step
!  IF Rho_LS=2 do a lineal search for density matrix in all steeps
   logical          :: changed_to_LS
   logical          :: P_oscilation_analisis = .false.
!------------------------------!

   double precision :: Elast
   double precision :: Pstepsize

   double precision, dimension(:)  , allocatable :: rho_lambda1
   double precision, dimension(:)  , allocatable :: rho_lambda0
   double precision, dimension(:)  , allocatable :: rho_lambda1_alpha
   double precision, dimension(:)  , allocatable :: rho_lambda0_alpha
   double precision, dimension(:)  , allocatable :: rho_lambda1_betha
   double precision, dimension(:)  , allocatable :: rho_lambda0_betha
   double precision, dimension(:,:), allocatable :: P_hist

contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!


subroutine P_conver(niter, En, E1, E2, Ex, good, xnano, rho_a, rho_b)
!  IF Rho_LS=0 calculate convergence criteria for actual density matrix
!  IF Rho_LS=1 do a lineal search for density matrix only if energy > energy of previus step
!  IF Rho_LS=2 do a lineal search for density matrix in all steeps
   use basis_data, only : M 
   implicit none
   integer, intent(in) :: niter !type of lineal search criteria and step number
   double precision, intent(out) :: good ! convergence criteria
   double precision, intent(in) :: En 
   double precision, intent(inout) :: E1, E2, Ex
   real*8, dimension(M,M), intent(inout) :: xnano, rho_a, rho_b ! density matrices of actual steep
   logical :: may_conv ! true if predicted density != density of prvius steep

   may_conv=.true.
   if (Rho_LS.eq.0 .or. Rho_LS.eq.-1) then
   else if (Rho_LS.eq.1 .or. Rho_LS.eq.2) then
     call P_linear_calc(Rho_LS, niter, En, E1, E2, Ex, xnano, may_conv, rho_a, rho_b)
   else
     STOP "wrong Rho_LS value"
   end if
    
   call P_calc_fluctuation(good, xnano)
   if (.not. may_conv) good=-1.d0

end subroutine P_conver



subroutine P_calc_fluctuation(good, xnano)
!  calculates convergence criteria in density matrix, and store new density matrix in RMM(1:MM)
   use garcha_mod, only : RMM
   use basis_data, only : M
   double precision, intent(out) :: good
   real*8, dimension(M,M), intent(in) :: xnano
   integer :: jj,kk, Rposition, M2
   double precision :: del
   M2=2*M
   good = 0.0d0
   do jj=1,M
     do kk=jj,M
       Rposition=kk+(M2-jj)*(jj-1)/2
       del=xnano(jj,kk)-RMM(Rposition)
       del=del*sqrt(2.d0)
       good=good+del**2
       RMM(kk+(M2-jj)*(jj-1)/2)=xnano(jj,kk)
     enddo
   enddo
   good=sqrt(good)/float(M)
end subroutine P_calc_fluctuation



subroutine P_linearsearch_init()
   use garcha_mod, only : OPEN, RMM, rhoalpha, rhobeta
   use basis_data, only : M, Md
   implicit none
   integer :: MM

   MM=M*(M+1)/2
   if (.not. allocated(rho_lambda1)) allocate(rho_lambda1(MM))
   if (.not. allocated(rho_lambda0)) allocate(rho_lambda0(MM))
   
   if(OPEN) then
     if (.not. allocated(rho_lambda1_alpha)) allocate(rho_lambda1_alpha(MM))
     if (.not. allocated(rho_lambda0_alpha)) allocate(rho_lambda0_alpha(MM))
     if (.not. allocated(rho_lambda1_betha)) allocate(rho_lambda1_betha(MM))
     if (.not. allocated(rho_lambda0_betha)) allocate(rho_lambda0_betha(MM))
     rho_lambda0_alpha(1:MM)=rhoalpha(1:MM)
     rho_lambda0_betha(1:MM)=rhobeta(1:MM)
   end if

end subroutine P_linearsearch_init



subroutine P_linearsearch_fin()
   if (allocated(rho_lambda1))       deallocate(rho_lambda1)
   if (allocated(rho_lambda0))       deallocate(rho_lambda0)
   if (allocated(rho_lambda1_alpha)) deallocate(rho_lambda1_alpha)
   if (allocated(rho_lambda0_alpha)) deallocate(rho_lambda0_alpha)
   if (allocated(rho_lambda1_betha)) deallocate(rho_lambda1_betha)
   if (allocated(rho_lambda0_betha)) deallocate(rho_lambda0_betha)
end subroutine P_linearsearch_fin



subroutine P_linear_calc(Rho_LS, niter, En, E1, E2, Ex, xnano,  &
may_conv, rho_a, rho_b)
   use garcha_mod, only : RMM, rhoalpha, rhobeta, OPEN
   use basis_data, only : M
   use liosubs, only: line_search
   implicit none
   integer, intent(in) :: Rho_LS, niter !type of lineal search criteria and step number
   double precision, intent(in) :: En
   double precision, intent(inout) :: E1, E2, Ex 
   real*8,dimension(M,M), intent(inout) :: xnano, rho_a, rho_b ! density matrices of actual steep
   double precision :: dlambda, Blambda !values for combination of density matrices 
   double precision, dimension(0:10) :: E_lambda !energy array in lineal search
   logical, intent(inout) :: may_conv !true if predicted density != previus step density
   double precision, dimension(:), allocatable :: RMM_temp ! auxiliar
   integer :: M2, MM, jj, kk, Rposition, ilambda !auxiliars

   M2=2*M
   MM=M*(M+1)/2
   if (niter .eq. 1)  then
     Pstepsize=1.d0
     rho_lambda0(1:MM)=RMM(1:MM)
     if (OPEN) rho_lambda0_alpha(1:MM)=rhoalpha(1:MM)
     if (OPEN) rho_lambda0_betha(1:MM)=rhobeta(1:MM)
   end if
   
   allocate(RMM_temp(1:MM))
   
   do jj=1,M
     do kk=jj,M
       Rposition=kk+(M2-jj)*(jj-1)/2
       rho_lambda1(Rposition)=xnano(jj,kk)
       if (OPEN) rho_lambda1_alpha(Rposition)=rho_a(jj,kk)
       if (OPEN) rho_lambda1_betha(Rposition)=rho_b(jj,kk)
     enddo
   enddo
   
   RMM(1:MM)=rho_lambda1(1:MM)
   if (OPEN)rhoalpha(1:MM)=rho_lambda1_alpha(1:MM)
   if (OPEN)rhobeta(1:MM)=rho_lambda1_betha(1:MM)
   
   call give_me_energy(E_lambda(10), En, E1, E2, Ex)
   
   if (Elast.lt. E_lambda(10) .or. Rho_LS.eq.2) then
     write(*,*) "This step ", E_lambda(10), "last steep ", Elast
     write(*,*) "doing lineal interpolation in Rho"
     do ilambda=0, 10
       dlambda=Pstepsize*dble(ilambda)/10.d0
       if (dlambda .gt. 1.d0) STOP "dlambda > 1.d0"
       RMM(1:MM)=rho_lambda0(1:MM)*(1.d0-dlambda)+rho_lambda1(1:MM)*dlambda 
       if(OPEN) rhoalpha(1:MM)=rho_lambda0_alpha(1:MM)*(1.d0-dlambda)+rho_lambda1_alpha(1:MM)*dlambda
       if(OPEN) rhobeta(1:MM)=rho_lambda0_betha(1:MM)*(1.d0-dlambda)+rho_lambda1_betha(1:MM)*dlambda
       call give_me_energy(E_lambda(ilambda), En, E1, E2, Ex)
   
       write(*,*) "step ",ilambda, "energy ", E_lambda(ilambda)
     end do
   
     call line_search(11,E_lambda, 1d0, Blambda )
     if (Blambda .ge. 1.d0) Blambda=Blambda-1.0d0
     write(*,*) "Best lambda", Blambda
     Blambda=Blambda*Pstepsize/10.d0
     write(*,*) "Fluctuation: ", Blambda
   else
     Blambda=Pstepsize
   end if
   
   RMM(1:MM)=rho_lambda0(1:MM)*(1.d0-Blambda)+rho_lambda1(1:MM)*Blambda
   if(OPEN) rhoalpha(1:MM)=rho_lambda0_alpha(1:MM)*(1.d0-Blambda)+rho_lambda1_alpha(1:MM)*Blambda
   if(OPEN) rhobeta(1:MM)=rho_lambda0_betha(1:MM)*(1.d0-Blambda)+rho_lambda1_betha(1:MM)*Blambda
   
   do jj=1,M
     do kk=jj,M
       Rposition=kk+(M2-jj)*(jj-1)/2
       xnano(jj,kk)=RMM(Rposition)
       if(OPEN) rho_a(jj,kk)=rhoalpha(Rposition)
       if(OPEN) rho_b(jj,kk)=rhobeta(Rposition)
     enddo
   enddo
   call give_me_energy(Elast, En, E1, E2, Ex)
   
   
   RMM_temp(1:MM)=RMM(1:MM)
   RMM(1:MM)=rho_lambda0(1:MM)
   rho_lambda0(1:MM)=RMM_temp(1:MM)
   
   if(OPEN) then
     RMM_temp(1:MM)=rhoalpha(1:MM)
     rhoalpha(1:MM)=rho_lambda0_alpha(1:MM)
     rho_lambda0_alpha(1:MM)=RMM_temp(1:MM)
   
     RMM_temp(1:MM)=rhobeta(1:MM)
     rhobeta(1:MM)=rho_lambda0_betha(1:MM)
     rho_lambda0_betha(1:MM)=RMM_temp(1:MM)
   end if
   
   deallocate(RMM_temp)
   if (Blambda .le. 4.d-1*Pstepsize) Pstepsize=Pstepsize*0.5d0 
   if (Blambda .ge. 8.d-1*Pstepsize) Pstepsize=Pstepsize*1.2d0
   if (Pstepsize .gt. 1.d0) Pstepsize=1.d0
   if (Blambda .le. 2.d-1*Pstepsize .and. Pstepsize .gt. 1d-4) may_conv=.false.

end subroutine P_linear_calc



subroutine give_me_energy(E, En, E1, E2, Ex)
!  return Energy components for a density matrix stored in RMM(1:MM)
   use garcha_mod, only : RMM, OPEN, MEMO
   use basis_data, only : M, Md
   use faint_cpu, only: int3lu
   implicit none
   double precision, intent(out) :: E
   double precision, intent(in) :: En
   double precision, intent(out) :: E1, E2, Ex
   integer :: kk, MM, MMd, M1, M3, M5, M7, M9, M11
   MM=M*(M+1)/2
   MMd=Md*(Md+1)/2
   
   M1=1 ! first P
   M3=M1+MM ! now Pnew
   M5=M3+MM! now S, F also uses the same position after S was used
   M7=M5+MM! now G
   M9=M7+MMd ! now Gm
   M11=M9+MMd! now H
   
   
   E=0.d0
   E1=0.D0
   E2=0.d0
   Ex=0.d0
   M11=1+3*MM+2*MMd
   
   do kk=1,MM
     E1=E1+RMM(kk)*RMM(M11+kk-1) !Computes 1e energy
   enddo
   ! Computes Coulomb part of Fock, and energy on E2
   call int3lu(E2, RMM(1:MM), RMM(M3:M3+MM), RMM(M5:M5+MM), RMM(M7:M7+MMd), &
        RMM(M9:M9+MMd), RMM(M11:M11+MMd), OPEN, MEMO)
   call g2g_solve_groups(0,Ex,0) ! Computes XC integration / Fock elements
   E=E1+E2+En+Ex

end subroutine give_me_energy


end module converger_ls
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
