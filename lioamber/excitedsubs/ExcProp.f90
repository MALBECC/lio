subroutine ExcProp(CoefA,CoefB,EneA,EneB,Etot)
! This is the main routine to calculate excited states properties
! For now this perform:
! - Linear Response Calculation


! Inputs
! - CoefA: Molecular Orbitals Coefficient of alpha
! - CoefB: Molecular Orbitals COefficient of beta
! - EneA: Molecular Orbitals Energy of alpha
! - EneB: Molecular Orbitals Energy of beta
use garcha_mod, only: OPEN
use excited_data, only: lresp
   implicit none

   double precision, intent(in) :: CoefA(:,:), CoefB(:,:)
   double precision, intent(in) :: EneA(:), EneB(:)
   double precision, intent(inout) :: Etot

   if (lresp .eqv. .false.) return
   if (OPEN  .eqv. .true. ) then 
      print*, "Linear Response doesn't work in Open shell"
      stop
   endif





end subroutine ExcProp
