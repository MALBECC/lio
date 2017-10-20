            subroutine cumxtp_DZ(A,devPtrX,C,M)
!========================================================================!
!!!!!!!!  Hace C=A*X para matrices cuadradas
!========================================================================!
      implicit none
      integer*8 devPtrA
      integer*8 devPtrScratch
      integer,intent(in) :: M
      integer*8,intent(in) :: devPtrX
      REAL*8, dimension (:,:), ALLOCATABLE :: scratch1,scratch2
      integer i,j,stat
      external CUBLAS_INIT, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX
      external CUBLAS_SHUTDOWN, CUBLAS_ALLOC, CUBLAS_CGEMM,CUBLAS_ZGEMM
      integer CUBLAS_ALLOC, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX
      integer CUBLAS_INIT,CUBLAS_CGEMM,CUBLAS_ZGEMM
      integer sizeof_complex
      COMPLEX*16 , intent(in) :: A(M,M)
      COMPLEX*16, intent(out) :: C(M,M)
      COMPLEX*16 :: alpha,beta
      parameter(sizeof_complex=16)
!---------------------------------------------------------------------!
      alpha=cmplx(1.0D0,0.0D0)
      beta=cmplx(0.0D0,0.0D0)
!-------------------------ALLOCATION---------------------------------!
      stat= CUBLAS_ALLOC(M*M, sizeof_complex, devPtrA)
      if (stat.NE.0) then
        write(*,*) "Matrix allocation failed -cumpx-1"
        call CUBLAS_SHUTDOWN
        stop
      endif
      stat= CUBLAS_ALLOC(M*M, sizeof_complex, devPtrScratch)
      if (stat.NE.0) then
        write(*,*) "Matrix allocation failed -cumpx-2"
        call CUBLAS_SHUTDOWN
        stop
      endif
!--------------------------------------------------------------------!
      stat = CUBLAS_SET_MATRIX(M,M,sizeof_complex,A,M,
     > devPtrA,M)
      if (stat.NE.0) then
        write(*,*) "Matrix setting on device -cumpx-1"
        call CUBLAS_SHUTDOWN
        stop
      endif
!-----------------------GEMM-----------------------------------------!
      stat=CUBLAS_ZGEMM ('T','N',M,M,M,alpha,devPtrX
     > ,M ,devPtrA,M, beta, devPtrScratch,M)
      if (stat.NE.0) then
        write(*,*) "ZGEMM failed -cumpx"
        call CUBLAS_SHUTDOWN
        stop
      endif
!-----------------GETTING MATRIX FROM DEVICE-------------------------!
      stat = CUBLAS_GET_MATRIX(M,M,sizeof_complex,devPtrScratch,M,
     > C,M)
      if (stat.NE.0) then
        write(*,*) "Getting matrix from device failed -cumpx"
        call CUBLAS_SHUTDOWN
        stop
      endif
!--------------------------------------------------------------------!
      call CUBLAS_FREE ( devPtrA )
      call CUBLAS_FREE ( devPtrScratch )
      return
      end subroutine
!====================================================================!
            subroutine cumxtp_DC(A,devPtrX,C,M)
!========================================================================!
!!!!!!!!  Hace C=A*X para matrices cuadradas
!========================================================================!
      implicit none
      integer*8 devPtrA
      integer*8 devPtrScratch
      integer,intent(in) :: M
      integer*8,intent(in) :: devPtrX
      REAL*8, dimension (:,:), ALLOCATABLE :: scratch1,scratch2
      integer i,j,stat
      external CUBLAS_INIT, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX
      external CUBLAS_SHUTDOWN, CUBLAS_ALLOC, CUBLAS_CGEMM,CUBLAS_ZGEMM
      integer CUBLAS_ALLOC, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX
      integer CUBLAS_INIT,CUBLAS_CGEMM,CUBLAS_ZGEMM
      integer sizeof_complex
      COMPLEX*8 , intent(in) :: A(M,M)
      COMPLEX*8, intent(out) :: C(M,M)
      COMPLEX*8 :: alpha,beta
      parameter(sizeof_complex=8)
!---------------------------------------------------------------------!
      alpha=cmplx(1.0D0,0.0D0)
      beta=cmplx(0.0D0,0.0D0)
!-------------------------ALLOCATION---------------------------------!
      stat= CUBLAS_ALLOC(M*M, sizeof_complex, devPtrA)
      if (stat.NE.0) then
        write(*,*) "Matrix allocation failed -cumpx-1"
        call CUBLAS_SHUTDOWN
        stop
      endif
      stat= CUBLAS_ALLOC(M*M, sizeof_complex, devPtrScratch)
      if (stat.NE.0) then
        write(*,*) "Matrix allocation failed -cumpx-2"
        call CUBLAS_SHUTDOWN
        stop
      endif
!--------------------------------------------------------------------!
      stat = CUBLAS_SET_MATRIX(M,M,sizeof_complex,A,M,
     > devPtrA,M)
      if (stat.NE.0) then
        write(*,*) "Matrix setting on device -cumpx-1"
        call CUBLAS_SHUTDOWN
        stop
      endif
!-----------------------GEMM-----------------------------------------!
      stat=CUBLAS_CGEMM ('T','N',M,M,M,alpha,devPtrX
     > ,M ,devPtrA,M, beta, devPtrScratch,M)
      if (stat.NE.0) then
        write(*,*) "CGEMM failed -cumpx"
        call CUBLAS_SHUTDOWN
        stop
      endif
!-----------------GETTING MATRIX FROM DEVICE-------------------------!
      stat = CUBLAS_GET_MATRIX(M,M,sizeof_complex,devPtrScratch,M,
     > C,M)
      if (stat.NE.0) then
        write(*,*) "Getting matrix from device failed -cumpx"
        call CUBLAS_SHUTDOWN
        stop
      endif
!--------------------------------------------------------------------!
      call CUBLAS_FREE ( devPtrA )
      call CUBLAS_FREE ( devPtrScratch )
      return
      end subroutine

