subroutine basis_initLR(Coef,M,Mlr,NCO,Nvirt)
!  This subroutine initializes the matrices needed for 
!  the change basis in LR, Zvector and Forces calculations.
use excited_data, only: Coef_trans, Cocc, Cocc_trans, &
                  Cvir, Cvir_trans
   implicit none
   integer, intent(in) :: M, Mlr, NCO, Nvirt
   double precision, intent(in) :: Coef(M,Mlr)

   integer :: jj

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
end subroutine basis_initLR

subroutine basis_deinitLR()
use excited_data, only: Coef_trans, Cocc, &
                   Cocc_trans, Cvir, Cvir_trans
   implicit none
   deallocate(Coef_trans, Cocc, Cocc_trans, Cvir, Cvir_trans)
end subroutine basis_deinitLR
