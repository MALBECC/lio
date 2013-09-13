       module latom
       
      integer, dimension(:), ALLOCATABLE :: natomc,nnps,nnpp,nnpd,nns
      integer, dimension(:), ALLOCATABLE :: nnd,nnp
      real*8, dimension (:), ALLOCATABLE :: atmin
      integer, dimension(:,:), ALLOCATABLE :: jatc
      integer kknums,kknumd
      integer, dimension (:), ALLOCATABLE :: kkind,kkinds
      real*8     rmax, rmaxs
      real*8, dimension (:), ALLOCATABLE :: cool
      real*4, dimension (:), ALLOCATABLE :: cools
      end module latom
