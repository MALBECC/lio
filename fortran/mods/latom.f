       module latom
       
      integer, dimension(:), ALLOCATABLE :: natomc,nnps,nnpp,nnpd,nns
      integer, dimension(:), ALLOCATABLE :: nnd,nnp
      real*8, dimension (:), ALLOCATABLE :: atmin
      integer, dimension(:,:), ALLOCATABLE :: jatc
      integer kknums,kknump,kknumd
      integer, dimension (:), ALLOCATABLE :: kkind
      parameter (rmax=16.D0)
      real*8, dimension (:), ALLOCATABLE :: cool
      end module latom
