subroutine basis_initLR(Coef,M,Mlr,NCO,Nvirt,NCOb,Nvirtb,CoefB)
!  This subroutine initializes the matrices needed for 
!  the change basis in LR, Zvector and Forces calculations.
use garcha_mod  , only: OPEN
use excited_data, only: Coef_trans, Cocc, Cocc_trans, Cvir, Cvir_trans, &
                        Coef_transB, CoccB, Cocc_transB, CvirB, Cvir_transB
   implicit none
   integer, intent(in) :: M, Mlr, NCO, Nvirt, NCOb, Nvirtb
   LIODBLE, intent(in) :: Coef(M,Mlr)
   LIODBLE, intent(in), optional :: CoefB(M,Mlr)

   integer :: jj

   ! For Closed Shell
   allocate(Coef_trans(Mlr,M),Cocc(M,NCO),Cocc_trans(NCO,M), &
            Cvir(M,Nvirt),Cvir_trans(Nvirt,M))

   do jj=1,NCO
      Cocc(:,jj) = Coef(:,jj)
      Cocc_trans(jj,:) = Coef(:,jj)
   enddo

   do jj=1,Nvirt
      Cvir(:,jj) = Coef(:,NCO+jj)
      Cvir_trans(jj,:) = Coef(:,NCO+jj)
   enddo
   Coef_trans = transpose(Coef)

   ! For Open Shell
   if ( OPEN ) then
      allocate(Coef_transB(Mlr,M),CoccB(M,NCOb),Cocc_transB(NCOb,M), &
               CvirB(M,Nvirtb),Cvir_transB(Nvirtb,M))

      do jj=1,NCOb
         CoccB(:,jj) = CoefB(:,jj)
         Cocc_transB(jj,:) = CoefB(:,jj)
      enddo

      do jj=1,Nvirtb
         CvirB(:,jj) = CoefB(:,NCOb+jj)
         Cvir_transB(jj,:) = CoefB(:,NCOb+jj)
      enddo
      Coef_transB = transpose(CoefB)
   endif
end subroutine basis_initLR

subroutine basis_deinitLR()
use garcha_mod  , only: OPEN
use excited_data, only: Coef_trans, Cocc, Cocc_trans, Cvir, Cvir_trans, &
                        Coef_transB, CoccB, Cocc_transB, CvirB, Cvir_transB
   implicit none
   deallocate(Coef_trans, Cocc, Cocc_trans, Cvir, Cvir_trans)
   if ( OPEN ) deallocate(Coef_transB, CoccB, Cocc_transB, CvirB, Cvir_transB)
end subroutine basis_deinitLR
