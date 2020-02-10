!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!  CONMUTATOR - gemm version -
!
!  CONMUTATOR(MA,MB)=MC=[MA*MB-MB*MA]
!====================================================================!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       function commutator_cublas_dd(MA,MB)
     > result(MC)
       implicit none
       real*8,intent(in)      :: MA(:,:)
       real*8,intent(in)      :: MB(:,:)
       real*8,allocatable     :: MC(:,:)
       integer                :: nn
       integer sizeof_real
       integer*8 devPtrA
       integer*8 devPtrB
       integer*8 devPtrC
       parameter(sizeof_real=8)
       REAL*8 alpha,beta
       integer i,j,stat
       external CUBLAS_INIT, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX
       external CUBLAS_SHUTDOWN, CUBLAS_ALLOC, CUBLAS_DGEMM
       integer CUBLAS_ALLOC, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX
       integer CUBLAS_INIT, CUBLAS_DGEMM
       stat=CUBLAS_INIT()
       nn=size(MA,1)
       allocate(MC(nn,nn))
       alpha=1.0000000000
       beta=0.00000000000
       MC=0.0D0
       stat= CUBLAS_ALLOC(nn*nn, sizeof_real, devPtrA)
       if (stat.NE.0) then
        write(*,*) "device memory allocation failed -conmutator_dd"
        call CUBLAS_SHUTDOWN
        stop
       endif
       stat= CUBLAS_ALLOC(nn*nn, sizeof_real, devPtrB)
       if (stat.NE.0) then
        write(*,*) "device memory allocation failed -conmutator_dd"
        call CUBLAS_SHUTDOWN
        stop
       endif
       stat= CUBLAS_ALLOC(nn*nn, sizeof_real, devPtrC)
       if (stat.NE.0) then
        write(*,*) "device memory allocation failed -conmutator_dd"
        call CUBLAS_SHUTDOWN
        stop
       endif
       stat = CUBLAS_SET_MATRIX(nn,nn,sizeof_real,MA,nn,devPtrA,nn)
       if (stat.NE.0) then
        write(*,*) "matrix setting failed -conmutator_dd"
        call CUBLAS_SHUTDOWN
        stop
       endif
       stat = CUBLAS_SET_MATRIX(nn,nn,sizeof_real,MB,nn,devPtrB,nn)
       if (stat.NE.0) then
        write(*,*) "matrix setting failed -conmutator_dd"
        call CUBLAS_SHUTDOWN
        stop
       endif
       stat = CUBLAS_SET_MATRIX(nn,nn,sizeof_real,MC,nn,devPtrC,nn)
       if (stat.NE.0) then
         call CUBLAS_FREE( devPtrC )
         write(*,*) "matrix setting failed -commutator_dd"
         call CUBLAS_SHUTDOWN
         stop
       endif
       stat=CUBLAS_DGEMM ('N','N',nn,nn,nn,alpha,devPtrA
     > ,nn ,devPtrB,nn, beta, devPtrC,nn)
       if (stat.NE.0) then
         call CUBLAS_FREE( devPtrC )
         write(*,*) "DGEMM - 1 - failed -commutator_dd"
         call CUBLAS_SHUTDOWN
         stop
       endif
       beta=(-1.00000000000,0.00000000000)
      stat=CUBLAS_DGEMM ('N','N',nn,nn,nn,alpha,devPtrB
     > ,nn ,devPtrA,nn, beta, devPtrC,nn)
      if (stat.NE.0) then
         call CUBLAS_FREE( devPtrC )
         write(*,*) "DGEMM - 2 - failed -commutator_dd"
         call CUBLAS_SHUTDOWN
         stop
      endif
      stat = CUBLAS_GET_MATRIX(nn, nn, sizeof_real, devPtrC, nn, MC, nn)
      if (stat.NE.0) then
      write(*,*) "data download failed -cuconmutc_r"
      call CUBLAS_FREE ( devPtrA )
      call CUBLAS_FREE ( devPtrB )
      call CUBLAS_FREE ( devPtrC )
      call CUBLAS_SHUTDOWN
      stop
      endif
      call CUBLAS_FREE ( devPtrA )
      call CUBLAS_FREE ( devPtrB )
      call CUBLAS_FREE ( devPtrC )
      return;end function
!
!--------------------------------------------------------------------!
       function commutator_cublas_zd(MA,MB)
     > result(MC)
       implicit none
       complex*16,intent(in)  :: MA(:,:)
       real*8,intent(in)      :: MB(:,:)
       complex*16,allocatable :: MC(:,:)
       complex*16 :: alpha,beta
       integer                :: nn
       integer sizeof_complex
       integer*8 devPtrA
       integer*8 devPtrB
       integer*8 devPtrC
       integer i,j,stat
       external CUBLAS_INIT, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX
       external CUBLAS_SHUTDOWN, CUBLAS_ALLOC, CUBLAS_ZGEMM
       integer CUBLAS_ALLOC, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX
       integer CUBLAS_INIT,CUBLAS_ZGEMM
       COMPLEX*16, dimension (:,:), ALLOCATABLE :: scratch
       parameter(sizeof_complex=16)
!
       stat=CUBLAS_INIT()
       if (stat.NE.0) then
           write(*,*) "initialization failed -commutator_zd"
           call CUBLAS_SHUTDOWN
           stop
       endif
       nn=size(MA,1)
       allocate(MC(nn,nn))
       allocate(scratch(nn,nn))
       alpha=(1.0D0,0.0D0)
       beta=(0.0D0,0.0D0)
       MC=(0.0D0,0.0D0)
       do i=1,nn
          do j=1,nn
             scratch(i,j)=CMPLX(MB(i,j),0.0D0)
          enddo
      enddo
      stat= CUBLAS_ALLOC(nn*nn, sizeof_complex, devPtrA)
      if (stat.NE.0) then
         write(*,*) "allocation failed -commutator_zd"
         call CUBLAS_SHUTDOWN
         stop
      endif
      stat= CUBLAS_ALLOC(nn*nn, sizeof_complex, devPtrB)
      if (stat.NE.0) then
         write(*,*) "allocation failed -commutator_zd"
         call CUBLAS_SHUTDOWN
         stop
      endif
      stat= CUBLAS_ALLOC(nn*nn, sizeof_complex, devPtrC)
      if (stat.NE.0) then
        write(*,*) "allocation failed -conmutator_zd"
        call CUBLAS_SHUTDOWN
        stop
      endif
      stat = CUBLAS_SET_MATRIX(nn,nn,sizeof_complex,MA,nn,devPtrA,nn)
      if (stat.NE.0) then
        write(*,*) "matrix setting failed -conmutator_zd"
        call CUBLAS_SHUTDOWN
        stop
      endif
      stat=CUBLAS_SET_MATRIX(nn,nn,sizeof_complex,scratch,nn,devPtrB,nn)
      if (stat.NE.0) then
        write(*,*) "matrix setting failed -conmutator_zd"
        call CUBLAS_SHUTDOWN
        stop
      endif
      stat = CUBLAS_SET_MATRIX(nn,nn,sizeof_complex,MC,nn,devPtrC,nn)
      if (stat.NE.0) then
        call CUBLAS_FREE( devPtrA )
        write(*,*) "matrix setting failed -conmutator_zd"
        call CUBLAS_SHUTDOWN
        stop
      endif
      stat=CUBLAS_ZGEMM ('N','N',nn,nn,nn,alpha,devPtrB
     > ,nn ,devPtrA,nn, beta, devPtrC,nn)
      if (stat.NE.0) then
        write(*,*) "ZGEMM failed -conmutator_zd"
        call CUBLAS_SHUTDOWN
        stop
      endif
      beta=(-1.0D0,0.0D0)
      stat=CUBLAS_ZGEMM ('N','N',nn,nn,nn,alpha,
     > devPtrA,nn ,devPtrB,nn, beta, devPtrC,nn)
      if (stat.NE.0) then
        write(*,*) "ZGEMM failed -conmutator_zd"
        call CUBLAS_SHUTDOWN
        stop
      endif
      stat=CUBLAS_GET_MATRIX(nn, nn, sizeof_complex, devPtrC, nn, MC,nn)
      if (stat.NE.0) then
         write(*,*) "data upload failed"
         call CUBLAS_FREE ( devPtrA )
         call CUBLAS_FREE ( devPtrB )
         call CUBLAS_FREE ( devPtrC )
         call CUBLAS_SHUTDOWN
         stop
      endif
      call CUBLAS_FREE ( devPtrA )
      call CUBLAS_FREE ( devPtrB )
      call CUBLAS_FREE ( devPtrC )
      call CUBLAS_SHUTDOWN
      DEALLOCATE(scratch)
      return;end function
!
!--------------------------------------------------------------------!
       function commutator_cublas_dz(MA,MB)
     > result(MC)
       implicit none
       real*8,intent(in)  :: MA(:,:)
       complex*16,intent(in)      :: MB(:,:)
       complex*16,allocatable :: MC(:,:)
       complex*16 :: alpha,beta
       integer                :: nn
       integer sizeof_complex
       integer*8 devPtrA
       integer*8 devPtrB
       integer*8 devPtrC
       integer i,j,stat
       external CUBLAS_INIT, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX
       external CUBLAS_SHUTDOWN, CUBLAS_ALLOC,CUBLAS_ZGEMM
       integer CUBLAS_ALLOC, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX
       integer CUBLAS_INIT,CUBLAS_ZGEMM
       COMPLEX*16, dimension (:,:), ALLOCATABLE :: scratch
       parameter(sizeof_complex=16)
!
       stat=CUBLAS_INIT()
       if (stat.NE.0) then
           write(*,*) "initialization failed -commutator_dz"
           call CUBLAS_SHUTDOWN
           stop
       endif
       nn=size(MA,1)
       allocate(MC(nn,nn))
       allocate(scratch(nn,nn))
       alpha=(1.0D0,0.0D0)
       beta=(0.0D0,0.0D0)
       MC=(0.0D0,0.0D0)
       do i=1,nn
          do j=1,nn
             scratch(i,j)=CMPLX(MA(i,j),0.0D0)
          enddo
      enddo
      stat= CUBLAS_ALLOC(nn*nn, sizeof_complex, devPtrA)
      if (stat.NE.0) then
         write(*,*) "allocation failed -commutator_dz"
         call CUBLAS_SHUTDOWN
         stop
      endif
      stat= CUBLAS_ALLOC(nn*nn, sizeof_complex, devPtrB)
      if (stat.NE.0) then
         write(*,*) "allocation failed -commutator_dz"
         call CUBLAS_SHUTDOWN
         stop
      endif
      stat= CUBLAS_ALLOC(nn*nn, sizeof_complex, devPtrC)
      if (stat.NE.0) then
        write(*,*) "allocation failed -conmutator_dz"
        call CUBLAS_SHUTDOWN
        stop
      endif
      stat=CUBLAS_SET_MATRIX(nn,nn,sizeof_complex,scratch,nn,devPtrA,nn)
      if (stat.NE.0) then
        write(*,*) "matrix setting failed -conmutator_dz"
        call CUBLAS_SHUTDOWN
        stop
      endif
      stat=CUBLAS_SET_MATRIX(nn,nn,sizeof_complex,MB,nn,devPtrB,nn)
      if (stat.NE.0) then
        write(*,*) "matrix setting failed -conmutator_dz"
        call CUBLAS_SHUTDOWN
        stop
      endif
      stat = CUBLAS_SET_MATRIX(nn,nn,sizeof_complex,MC,nn,devPtrC,nn)
      if (stat.NE.0) then
        call CUBLAS_FREE( devPtrA )
        write(*,*) "matrix setting failed -conmutator_dz"
        call CUBLAS_SHUTDOWN
        stop
      endif
      stat=CUBLAS_ZGEMM ('N','N',nn,nn,nn,alpha,devPtrB
     > ,nn ,devPtrA,nn, beta, devPtrC,nn)
      if (stat.NE.0) then
        write(*,*) "ZGEMM failed -conmutator_dz"
        call CUBLAS_SHUTDOWN
        stop
      endif
      beta=(-1.0D0,0.0D0)
      stat=CUBLAS_ZGEMM ('N','N',nn,nn,nn,alpha,
     > devPtrA,nn ,devPtrB,nn, beta, devPtrC,nn)
      if (stat.NE.0) then
        write(*,*) "ZGEMM failed -conmutator_dz"
        call CUBLAS_SHUTDOWN
        stop
      endif
      stat=CUBLAS_GET_MATRIX(nn, nn, sizeof_complex, devPtrC, nn, MC,nn)
      if (stat.NE.0) then
         write(*,*) "data download failed -conmutator_dz"
         call CUBLAS_FREE ( devPtrA )
         call CUBLAS_FREE ( devPtrB )
         call CUBLAS_FREE ( devPtrC )
         call CUBLAS_SHUTDOWN
         stop
      endif
      call CUBLAS_FREE ( devPtrA )
      call CUBLAS_FREE ( devPtrB )
      call CUBLAS_FREE ( devPtrC )
      call CUBLAS_SHUTDOWN
      DEALLOCATE(scratch)
      return;end function
!
!--------------------------------------------------------------------!
       function commutator_cublas_zz(MA,MB)
     > result(MC)
       implicit none
       complex*16,intent(in)  :: MA(:,:)
       complex*16,intent(in)      :: MB(:,:)
       complex*16,allocatable :: MC(:,:)
       complex*16 :: alpha,beta
       integer                :: nn
       integer sizeof_complex
       integer*8 devPtrA
       integer*8 devPtrB
       integer*8 devPtrC
       integer i,j,stat
       external CUBLAS_INIT, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX
       external CUBLAS_SHUTDOWN, CUBLAS_ALLOC, CUBLAS_ZGEMM
       integer CUBLAS_ALLOC, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX
       integer CUBLAS_INIT,CUBLAS_ZGEMM
       COMPLEX*16, dimension (:,:), ALLOCATABLE :: scratch
       parameter(sizeof_complex=16)
!
       stat=CUBLAS_INIT()
       if (stat.NE.0) then
           write(*,*) "initialization failed -commutator_zz"
           call CUBLAS_SHUTDOWN
           stop
       endif
       nn=size(MA,1)
       allocate(MC(nn,nn))
       allocate(scratch(nn,nn))
       alpha=(1.0D0,0.0D0)
       beta=(0.0D0,0.0D0)
       MC=(0.0D0,0.0D0)
      stat= CUBLAS_ALLOC(nn*nn, sizeof_complex, devPtrA)
      if (stat.NE.0) then
         write(*,*) "allocation failed -commutator_zz"
         call CUBLAS_SHUTDOWN
         stop
      endif
      stat= CUBLAS_ALLOC(nn*nn, sizeof_complex, devPtrB)
      if (stat.NE.0) then
         write(*,*) "allocation failed -commutator_zz"
         call CUBLAS_SHUTDOWN
         stop
      endif
      stat= CUBLAS_ALLOC(nn*nn, sizeof_complex, devPtrC)
      if (stat.NE.0) then
        write(*,*) "allocation failed -conmutator_zz"
        call CUBLAS_SHUTDOWN
        stop
      endif
      stat=CUBLAS_SET_MATRIX(nn,nn,sizeof_complex,MA,nn,devPtrA,nn)
      if (stat.NE.0) then
        write(*,*) "matrix setting failed -conmutator_zz"
        call CUBLAS_SHUTDOWN
        stop
      endif
      stat=CUBLAS_SET_MATRIX(nn,nn,sizeof_complex,MB,nn,devPtrB,nn)
      if (stat.NE.0) then
        write(*,*) "matrix setting failed -conmutator_zz"
        call CUBLAS_SHUTDOWN
        stop
      endif
      stat = CUBLAS_SET_MATRIX(nn,nn,sizeof_complex,MC,nn,devPtrC,nn)
      if (stat.NE.0) then
        call CUBLAS_FREE( devPtrA )
        write(*,*) "matrix setting failed -conmutator_zz"
        call CUBLAS_SHUTDOWN
        stop
      endif
      stat=CUBLAS_ZGEMM ('N','N',nn,nn,nn,alpha,devPtrB
     > ,nn ,devPtrA,nn, beta, devPtrC,nn)
      if (stat.NE.0) then
        write(*,*) "ZGEMM failed -conmutator_zz"
        call CUBLAS_SHUTDOWN
        stop
      endif
      beta=(-1.0D0,0.0D0)
      stat=CUBLAS_ZGEMM ('N','N',nn,nn,nn,alpha,
     > devPtrA,nn ,devPtrB,nn, beta, devPtrC,nn)
      if (stat.NE.0) then
        write(*,*) "ZGEMM failed -conmutator_zz"
        call CUBLAS_SHUTDOWN
        stop
      endif
      stat=CUBLAS_GET_MATRIX(nn, nn, sizeof_complex, devPtrC, nn, MC,nn)
      if (stat.NE.0) then
         write(*,*) "data download failed -conmutator_zz"
         call CUBLAS_FREE ( devPtrA )
         call CUBLAS_FREE ( devPtrB )
         call CUBLAS_FREE ( devPtrC )
         call CUBLAS_SHUTDOWN
         stop
      endif
      call CUBLAS_FREE ( devPtrA )
      call CUBLAS_FREE ( devPtrB )
      call CUBLAS_FREE ( devPtrC )
      call CUBLAS_SHUTDOWN
      return;end function
!
!--------------------------------------------------------------------!
       function commutator_cublas_cd(MA,MB)
     > result(MC)
       implicit none
       complex*8,intent(in)  :: MA(:,:)
       real*8,intent(in)      :: MB(:,:)
       complex*8,allocatable :: MC(:,:)
       complex*8 :: alpha,beta
       integer                :: nn
       integer sizeof_complex
       integer*8 devPtrA
       integer*8 devPtrB
       integer*8 devPtrC
       integer i,j,stat
       external CUBLAS_INIT, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX
       external CUBLAS_SHUTDOWN, CUBLAS_ALLOC,CUBLAS_CGEMM
       integer CUBLAS_ALLOC, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX
       integer CUBLAS_INIT, CUBLAS_CGEMM
       COMPLEX*8, dimension (:,:), ALLOCATABLE :: scratch
       parameter(sizeof_complex=8)
!
       stat=CUBLAS_INIT()
       if (stat.NE.0) then
           write(*,*) "initialization failed -commutator_zd"
           call CUBLAS_SHUTDOWN
           stop
       endif
       nn=size(MA,1)
       allocate(MC(nn,nn))
       allocate(scratch(nn,nn))
       alpha=(1.0E0,0.0E0)
       beta=(0.0E0,0.0E0)
       MC=(0.0E0,0.0E0)
       do i=1,nn
          do j=1,nn
             scratch(i,j)=CMPLX(MB(i,j),0.0E0)
          enddo
      enddo
      stat= CUBLAS_ALLOC(nn*nn, sizeof_complex, devPtrA)
      if (stat.NE.0) then
         write(*,*) "allocation failed -commutator_zd"
         call CUBLAS_SHUTDOWN
         stop
      endif
      stat= CUBLAS_ALLOC(nn*nn, sizeof_complex, devPtrB)
      if (stat.NE.0) then
         write(*,*) "allocation failed -commutator_zd"
         call CUBLAS_SHUTDOWN
         stop
      endif
      stat= CUBLAS_ALLOC(nn*nn, sizeof_complex, devPtrC)
      if (stat.NE.0) then
        write(*,*) "allocation failed -conmutator_zd"
        call CUBLAS_SHUTDOWN
        stop
      endif
      stat = CUBLAS_SET_MATRIX(nn,nn,sizeof_complex,MA,nn,devPtrA,nn)
      if (stat.NE.0) then
        write(*,*) "matrix setting failed -conmutator_zd"
        call CUBLAS_SHUTDOWN
        stop
      endif
      stat=CUBLAS_SET_MATRIX(nn,nn,sizeof_complex,scratch,nn,devPtrB,nn)
      if (stat.NE.0) then
        write(*,*) "matrix setting failed -conmutator_zd"
        call CUBLAS_SHUTDOWN
        stop
      endif
      stat = CUBLAS_SET_MATRIX(nn,nn,sizeof_complex,MC,nn,devPtrC,nn)
      if (stat.NE.0) then
        call CUBLAS_FREE( devPtrA )
        write(*,*) "matrix setting failed -conmutator_zd"
        call CUBLAS_SHUTDOWN
        stop
      endif
      stat=CUBLAS_CGEMM ('N','N',nn,nn,nn,alpha,devPtrB
     > ,nn ,devPtrA,nn, beta, devPtrC,nn)
      if (stat.NE.0) then
        write(*,*) "ZGEMM failed -conmutator_zd"
        call CUBLAS_SHUTDOWN
        stop
      endif
      beta=(-1.0E0,0.0E0)
      stat=CUBLAS_CGEMM ('N','N',nn,nn,nn,alpha,
     > devPtrA,nn ,devPtrB,nn, beta, devPtrC,nn)
      if (stat.NE.0) then
        write(*,*) "ZGEMM failed -conmutator_zd"
        call CUBLAS_SHUTDOWN
        stop
      endif
      stat=CUBLAS_GET_MATRIX(nn, nn, sizeof_complex, devPtrC, nn, MC,nn)
      if (stat.NE.0) then
         write(*,*) "data upload failed"
         call CUBLAS_FREE ( devPtrA )
         call CUBLAS_FREE ( devPtrB )
         call CUBLAS_FREE ( devPtrC )
         call CUBLAS_SHUTDOWN
         stop
      endif
      call CUBLAS_FREE ( devPtrA )
      call CUBLAS_FREE ( devPtrB )
      call CUBLAS_FREE ( devPtrC )
      call CUBLAS_SHUTDOWN
      DEALLOCATE(scratch)
      return;end function
!
!--------------------------------------------------------------------!
       function commutator_cublas_dc(MA,MB)
     > result(MC)
       implicit none
       real*8,intent(in)  :: MA(:,:)
       complex*8,intent(in)      :: MB(:,:)
       complex*8,allocatable :: MC(:,:)
       complex*8 :: alpha,beta
       integer                :: nn
       integer sizeof_complex
       integer*8 devPtrA
       integer*8 devPtrB
       integer*8 devPtrC
       integer i,j,stat
       external CUBLAS_INIT, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX
       external CUBLAS_SHUTDOWN, CUBLAS_ALLOC, CUBLAS_CGEMM
       integer CUBLAS_ALLOC, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX
       integer CUBLAS_INIT, CUBLAS_CGEMM
       COMPLEX*8, dimension (:,:), ALLOCATABLE :: scratch
       parameter(sizeof_complex=8)
!
       stat=CUBLAS_INIT()
       if (stat.NE.0) then
           write(*,*) "initialization failed -commutator_dz"
           call CUBLAS_SHUTDOWN
           stop
       endif
       nn=size(MA,1)
       allocate(MC(nn,nn))
       allocate(scratch(nn,nn))
       alpha=(1.0E0,0.0E0)
       beta=(0.0E0,0.0E0)
       MC=(0.0E0,0.0E0)
       do i=1,nn
          do j=1,nn
             scratch(i,j)=CMPLX(MA(i,j),0.0E0)
          enddo
      enddo
      stat= CUBLAS_ALLOC(nn*nn, sizeof_complex, devPtrA)
      if (stat.NE.0) then
         write(*,*) "allocation failed -commutator_dC"
         call CUBLAS_SHUTDOWN
         stop
      endif
      stat= CUBLAS_ALLOC(nn*nn, sizeof_complex, devPtrB)
      if (stat.NE.0) then
         write(*,*) "allocation failed -commutator_dc"
         call CUBLAS_SHUTDOWN
         stop
      endif
      stat= CUBLAS_ALLOC(nn*nn, sizeof_complex, devPtrC)
      if (stat.NE.0) then
        write(*,*) "allocation failed -conmutator_dc"
        call CUBLAS_SHUTDOWN
        stop
      endif
      stat=CUBLAS_SET_MATRIX(nn,nn,sizeof_complex,scratch,nn,devPtrA,nn)
      if (stat.NE.0) then
        write(*,*) "matrix setting failed -conmutator_dc"
        call CUBLAS_SHUTDOWN
        stop
      endif
      stat=CUBLAS_SET_MATRIX(nn,nn,sizeof_complex,MB,nn,devPtrB,nn)
      if (stat.NE.0) then
        write(*,*) "matrix setting failed -conmutator_dc"
        call CUBLAS_SHUTDOWN
        stop
      endif
      stat = CUBLAS_SET_MATRIX(nn,nn,sizeof_complex,MC,nn,devPtrC,nn)
      if (stat.NE.0) then
        call CUBLAS_FREE( devPtrA )
        write(*,*) "matrix setting failed -conmutator_dc"
        call CUBLAS_SHUTDOWN
        stop
      endif
      stat=CUBLAS_CGEMM ('N','N',nn,nn,nn,alpha,devPtrB
     > ,nn ,devPtrA,nn, beta, devPtrC,nn)
      if (stat.NE.0) then
        write(*,*) "CGEMM failed -conmutator_dc"
        call CUBLAS_SHUTDOWN
        stop
      endif
      beta=(-1.0D0,0.0D0)
      stat=CUBLAS_CGEMM ('N','N',nn,nn,nn,alpha,
     > devPtrA,nn ,devPtrB,nn, beta, devPtrC,nn)
      if (stat.NE.0) then
        write(*,*) "CGEMM failed -conmutator_dc"
        call CUBLAS_SHUTDOWN
        stop
      endif
      stat=CUBLAS_GET_MATRIX(nn, nn, sizeof_complex, devPtrC, nn, MC,nn)
      if (stat.NE.0) then
         write(*,*) "data download failed -conmutator_dc"
         call CUBLAS_FREE ( devPtrA )
         call CUBLAS_FREE ( devPtrB )
         call CUBLAS_FREE ( devPtrC )
         call CUBLAS_SHUTDOWN
         stop
      endif
      call CUBLAS_FREE ( devPtrA )
      call CUBLAS_FREE ( devPtrB )
      call CUBLAS_FREE ( devPtrC )
      call CUBLAS_SHUTDOWN
      DEALLOCATE(scratch)
      return;end function
!
!--------------------------------------------------------------------!
       function commutator_cublas_cc(MA,MB)
     > result(MC)
       implicit none
       complex*8,intent(in)  :: MA(:,:)
       complex*8,intent(in)      :: MB(:,:)
       complex*8,allocatable :: MC(:,:)
       complex*8 :: alpha,beta
       integer                :: nn
       integer sizeof_complex
       integer*8 devPtrA
       integer*8 devPtrB
       integer*8 devPtrC
       integer i,j,stat
       external CUBLAS_INIT, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX
       external CUBLAS_SHUTDOWN, CUBLAS_ALLOC, CUBLAS_CGEMM
       integer CUBLAS_ALLOC, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX
       integer CUBLAS_INIT, CUBLAS_CGEMM
       COMPLEX*8, dimension (:,:), ALLOCATABLE :: scratch
       parameter(sizeof_complex=8)
!
       stat=CUBLAS_INIT()
       if (stat.NE.0) then
           write(*,*) "initialization failed -commutator_cc"
           call CUBLAS_SHUTDOWN
           stop
       endif
       nn=size(MA,1)
       allocate(MC(nn,nn))
       allocate(scratch(nn,nn))
       alpha=(1.0E0,0.0E0)
       beta=(0.0E0,0.0E0)
       MC=(0.0E0,0.0E0)
      stat= CUBLAS_ALLOC(nn*nn, sizeof_complex, devPtrA)
      if (stat.NE.0) then
         write(*,*) "allocation failed -commutator_cc"
         call CUBLAS_SHUTDOWN
         stop
      endif
      stat= CUBLAS_ALLOC(nn*nn, sizeof_complex, devPtrB)
      if (stat.NE.0) then
         write(*,*) "allocation failed -commutator_cc"
         call CUBLAS_SHUTDOWN
         stop
      endif
      stat= CUBLAS_ALLOC(nn*nn, sizeof_complex, devPtrC)
      if (stat.NE.0) then
        write(*,*) "allocation failed -conmutator_cc"
        call CUBLAS_SHUTDOWN
        stop
      endif
      stat=CUBLAS_SET_MATRIX(nn,nn,sizeof_complex,MA,nn,devPtrA,nn)
      if (stat.NE.0) then
        write(*,*) "matrix setting failed -conmutator_cc"
        call CUBLAS_SHUTDOWN
        stop
      endif
      stat=CUBLAS_SET_MATRIX(nn,nn,sizeof_complex,MB,nn,devPtrB,nn)
      if (stat.NE.0) then
        write(*,*) "matrix setting failed -conmutator_cc"
        call CUBLAS_SHUTDOWN
        stop
      endif
      stat = CUBLAS_SET_MATRIX(nn,nn,sizeof_complex,MC,nn,devPtrC,nn)
      if (stat.NE.0) then
        call CUBLAS_FREE( devPtrA )
        write(*,*) "matrix setting failed -conmutator_cc"
        call CUBLAS_SHUTDOWN
        stop
      endif
      stat=CUBLAS_CGEMM ('N','N',nn,nn,nn,alpha,devPtrB
     > ,nn ,devPtrA,nn, beta, devPtrC,nn)
      if (stat.NE.0) then
        write(*,*) "ZGEMM failed -conmutator_cc"
        call CUBLAS_SHUTDOWN
        stop
      endif
      beta=(-1.0D0,0.0D0)
      stat=CUBLAS_CGEMM ('N','N',nn,nn,nn,alpha,
     > devPtrA,nn ,devPtrB,nn, beta, devPtrC,nn)
      if (stat.NE.0) then
        write(*,*) "ZGEMM failed -conmutator_cc"
        call CUBLAS_SHUTDOWN
        stop
      endif
      stat=CUBLAS_GET_MATRIX(nn, nn, sizeof_complex, devPtrC, nn, MC,nn)
      if (stat.NE.0) then
         write(*,*) "data download failed -conmutator_cc"
         call CUBLAS_FREE ( devPtrA )
         call CUBLAS_FREE ( devPtrB )
         call CUBLAS_FREE ( devPtrC )
         call CUBLAS_SHUTDOWN
         stop
      endif
      call CUBLAS_FREE ( devPtrA )
      call CUBLAS_FREE ( devPtrB )
      call CUBLAS_FREE ( devPtrC )
      call CUBLAS_SHUTDOWN
      return;end function
!====================================================================!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
