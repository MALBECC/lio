            subroutine cumpxt_ZD(A,devPtrX,C,M)
!===============================================================================!
!!!!!!!!  Hace C=A*Xtrans para matrices cuadradas
!===============================================================================!
      implicit none
      integer sizeof_complex
      integer*8 devPtrA
      integer*8 devPtrScratch
      integer,intent(in) :: M
      integer*8,intent(in) :: devPtrX
      integer i,j,stat
      external CUBLAS_INIT, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX
      external CUBLAS_SHUTDOWN, CUBLAS_ALLOC,CUBLAS_CGEMM,CUBLAS_ZGEMM
      integer CUBLAS_ALLOC, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX
      integer CUBLAS_INIT,CUBLAS_CGEMM,CUBLAS_ZGEMM
      COMPLEX*16 , intent(in) :: A(M,M)
      COMPLEX*16, intent(out) :: C(M,M)
      parameter(sizeof_complex=16)
      COMPLEX*16 alpha,beta
!-------------------------------------------------------------------------------!
      alpha=cmplx(1.0D0,0.0D0)
      beta=cmplx(0.0D0,0.0D0)
!--------------------------------------------------------------------!
      stat= CUBLAS_ALLOC(M*M, sizeof_complex, devPtrA)
      if (stat.NE.0) then
        write(*,*) "Allocation failed -cumpxt-1"
        call CUBLAS_SHUTDOWN
        stop
      endif
!--------------------------------------------------------------------!
      stat= CUBLAS_ALLOC(M*M, sizeof_complex, devPtrScratch)
      if (stat.NE.0) then
        write(*,*) "Allocation failed -cumpxt-2"
        call CUBLAS_SHUTDOWN
        stop
      endif
!--------------------------------------------------------------------!
      stat = CUBLAS_SET_MATRIX(M,M,sizeof_complex,A,M,
     > devPtrA,M)
      if (stat.NE.0) then
        write(*,*) "Matrix setting failed -cumpxt-1"
        call CUBLAS_SHUTDOWN
        stop
      endif
!-----------------------GEMM-----------------------------------------!
      stat=CUBLAS_ZGEMM ('N','T',M,M,M,alpha,devPtrA
     > ,M ,devPtrX,M, beta, devPtrScratch,M)
      if (stat.NE.0) then
        write(*,*) "ZGEMM(1) failed -cumpxt"
        call CUBLAS_SHUTDOWN
        stop
      endif
!--------------------------------------------------------------------!
      stat = CUBLAS_GET_MATRIX(M,M,sizeof_complex,devPtrScratch,M,
     > C,M)
      if (stat.NE.0) then
        write(*,*) "matrix copy failed -cumpxt"
        call CUBLAS_SHUTDOWN
        stop
      endif
!--------------------------------------------------------------------!
      call CUBLAS_FREE ( devPtrA )
      call CUBLAS_FREE ( devPtrScratch )
      return
      end subroutine
!===================================================================!
            subroutine cumpxt_CD(A,devPtrX,C,M)
!-------------------------------------------------------------------!
!!!!!!!!  Hace C=A*Xtrans para matrices cuadradas
!-------------------------------------------------------------------!
      implicit none
      integer sizeof_complex
      integer*8 devPtrA
      integer*8 devPtrScratch
      integer,intent(in) :: M
      integer*8,intent(in) :: devPtrX
      integer i,j,stat
      external CUBLAS_INIT, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX
      external CUBLAS_SHUTDOWN, CUBLAS_ALLOC,CUBLAS_CGEMM,CUBLAS_ZGEMM
      integer CUBLAS_ALLOC, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX
      integer CUBLAS_INIT,CUBLAS_CGEMM,CUBLAS_ZGEMM
      COMPLEX*8 , intent(in) :: A(M,M)
      COMPLEX*8, intent(out) :: C(M,M)
      parameter(sizeof_complex=8)
      COMPLEX*8 alpha,beta
!-------------------------------------------------------------------------------!
      alpha=cmplx(1.0D0,0.0D0)
      beta=cmplx(0.0D0,0.0D0)
!--------------------------------------------------------------------!
      stat= CUBLAS_ALLOC(M*M, sizeof_complex, devPtrA)
      if (stat.NE.0) then
        write(*,*) "Allocation failed -cumpxt-1"
        call CUBLAS_SHUTDOWN
        stop
      endif
!--------------------------------------------------------------------!
      stat= CUBLAS_ALLOC(M*M, sizeof_complex, devPtrScratch)
      if (stat.NE.0) then
        write(*,*) "Allocation failed -cumpxt-2"
        call CUBLAS_SHUTDOWN
        stop
      endif
!--------------------------------------------------------------------!
      stat = CUBLAS_SET_MATRIX(M,M,sizeof_complex,A,M,
     > devPtrA,M)
      if (stat.NE.0) then
        write(*,*) "Matrix setting failed -cumpxt-1"
        call CUBLAS_SHUTDOWN
        stop
      endif
!-----------------------GEMM-----------------------------------------!
      stat=CUBLAS_CGEMM ('N','T',M,M,M,alpha,devPtrA
     > ,M ,devPtrX,M, beta, devPtrScratch,M)
      if (stat.NE.0) then
        write(*,*) "CGEMM(1) failed -cumpxt"
        call CUBLAS_SHUTDOWN
        stop
      endif
!--------------------------------------------------------------------!
      stat = CUBLAS_GET_MATRIX(M,M,sizeof_complex,devPtrScratch,M,
     > C,M)
      if (stat.NE.0) then
        write(*,*) "matrix copy failed -cumpxt"
        call CUBLAS_SHUTDOWN
        stop
      endif
!--------------------------------------------------------------------!
      call CUBLAS_FREE ( devPtrA )
      call CUBLAS_FREE ( devPtrScratch )
      return
      end subroutine
!=====================================================================!
