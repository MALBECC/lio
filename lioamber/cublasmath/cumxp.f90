subroutine cumxp_DZ(A, devPtrX, C, M)
   ! Hace C=A*X para matrices cuadradas
   implicit none
   integer        , intent(in)  :: M
   integer(kind=8), intent(in)  :: devPtrX
   complex(kind=8), intent(in)  :: A(M, M)
   complex(kind=8), intent(out) :: C(M, M)

   integer(kind=8) :: devPtrA, devPtrScratch
   complex(kind=8) :: alpha, beta
   integer         :: stat, sizeof_complex
   parameter(sizeof_complex=16)

   external CUBLAS_INIT, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX, &
            CUBLAS_SHUTDOWN, CUBLAS_ALLOC, CUBLAS_ZGEMM
   integer  CUBLAS_ALLOC, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX, &
            CUBLAS_INIT, CUBLAS_ZGEMM

   alpha = cmplx(1.0D0,0.0D0)
   beta  = cmplx(0.0D0,0.0D0)

   stat = CUBLAS_ALLOC(M*M, sizeof_complex, devPtrA)
   if (stat /= 0) then
      write(*,*) "Matrix allocation failed -cumpx-1"
      call CUBLAS_SHUTDOWN
      stop
   endif
   stat = CUBLAS_ALLOC(M*M, sizeof_complex, devPtrScratch)
   if (stat /= 0) then
      write(*,*) "Matrix allocation failed -cumpx-2"
      call CUBLAS_SHUTDOWN
      stop
   endif

   stat = CUBLAS_SET_MATRIX(M, M, sizeof_complex, A, M, devPtrA, M)
   if (stat /= 0) then
      write(*,*) "Matrix setting on device -cumpx-1"
      call CUBLAS_SHUTDOWN
      stop
   endif

   stat = CUBLAS_ZGEMM ('N','N', M, M, M, alpha, devPtrX, M , devPtrA,&
                         M, beta, devPtrScratch, M)
   if (stat /= 0) then
      write(*,*) "ZGEMM failed -cumpx"
      call CUBLAS_SHUTDOWN
      stop
   endif

   stat = CUBLAS_GET_MATRIX(M, M, sizeof_complex, devPtrScratch, M, C, M)
   if (stat /= 0) then
      write(*,*) "Getting matrix from device failed -cumpx"
      call CUBLAS_SHUTDOWN
      stop
   endif

   call CUBLAS_FREE ( devPtrA )
   call CUBLAS_FREE ( devPtrScratch )
end subroutine

!====================================================================!
subroutine cumxp_DC(A, devPtrX, C, M)
   ! Hace C=A*X para matrices cuadradas
   implicit none
   integer        , intent(in)  :: M
   integer(kind=8), intent(in)  :: devPtrX
   complex(kind=4), intent(in)  :: A(M, M)
   complex(kind=4), intent(out) :: C(M, M)

   integer(kind=8) :: devPtrA, devPtrScratch
   complex(kind=4) :: alpha, beta
   integer         :: stat, sizeof_complex
   parameter(sizeof_complex=8)

   external CUBLAS_INIT, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX, &
            CUBLAS_SHUTDOWN, CUBLAS_ALLOC, CUBLAS_CGEMM
   integer  CUBLAS_ALLOC, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX, &
            CUBLAS_INIT, CUBLAS_CGEMM

   alpha = cmplx(1.0E0,0.0E0)
   beta  = cmplx(0.0E0,0.0E0)

   stat = CUBLAS_ALLOC(M*M, sizeof_complex, devPtrA)
   if (stat /= 0) then
      write(*,*) "Matrix allocation failed -cumpx-1"
      call CUBLAS_SHUTDOWN
      stop
   endif
   stat = CUBLAS_ALLOC(M*M, sizeof_complex, devPtrScratch)
   if (stat /= 0) then
      write(*,*) "Matrix allocation failed -cumpx-2"
      call CUBLAS_SHUTDOWN
      stop
   endif

   stat = CUBLAS_SET_MATRIX(M, M, sizeof_complex, A, M, devPtrA, M)
   if (stat /= 0) then
      write(*,*) "Matrix setting on device -cumpx-1"
      call CUBLAS_SHUTDOWN
      stop
   endif

   stat = CUBLAS_CGEMM ('N','N', M, M, M, alpha, devPtrX, M , devPtrA, M,&
                        beta, devPtrScratch, M)
   if (stat /= 0) then
      write(*,*) "CGEMM failed -cumpx"
      call CUBLAS_SHUTDOWN
      stop
   endif
   stat = CUBLAS_GET_MATRIX(M, M, sizeof_complex, devPtrScratch, M, C, M)
   if (stat /= 0) then
      write(*,*) "Getting matrix from device failed -cumpx"
      call CUBLAS_SHUTDOWN
      stop
   endif

   call CUBLAS_FREE ( devPtrA )
   call CUBLAS_FREE ( devPtrScratch )
end subroutine
