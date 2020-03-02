subroutine calc2eFITT(Mat,Fock,M)
use faint_cpu, only: int3lu
use garcha_mod,only: Ginv_vec, Gmat_vec, MEMO, OPEN
   implicit none

   integer, intent(in) :: M
   LIODBLE, intent(in) :: Mat(M,M)
   LIODBLE, intent(out) :: Fock(M,M)

   integer :: MM, ii, jj
   LIODBLE :: E
   LIODBLE, dimension(:), allocatable :: Vvec, Fvec_a
   LIODBLE, dimension(:), allocatable :: Hvec, Fvec_b
   LIODBLE, dimension(:,:), allocatable :: Tmat

   MM = M*(M+1)/2
   allocate(Vvec(MM),Fvec_a(MM),Hvec(MM),Fvec_b(MM))
   Hvec = 0.0d0; Fvec_b = 0.0d0 ! These variables aren't referenced here

   allocate(Tmat(M,M))
   do ii=1,M
      Tmat(ii,ii) = Mat(ii,ii)
      do jj=1,ii-1
         Tmat(ii,jj) = Mat(ii,jj) + Mat(jj,ii)
      enddo
   enddo
   call sprepack('L',M,Vvec,Tmat); deallocate(Tmat)
   call int3lu(E,Vvec,Fvec_b,Fvec_a,Gmat_vec,Ginv_vec,Hvec,OPEN,MEMO)
   Fvec_a = 2.0d0 * Fvec_a
   call spunpack('L',M,Fvec_a,Fock)

   deallocate(Vvec,Fvec_a,Hvec,Fvec_b)

end subroutine calc2eFITT
