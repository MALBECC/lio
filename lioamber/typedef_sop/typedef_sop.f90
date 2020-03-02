!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
#include "../datatypes/datatypes.fh"
module typedef_sop
!------------------------------------------------------------------------------!
!
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   implicit none

   type sop
      integer              :: Nbasis = 0
      LIODBLE , allocatable :: Smat(:,:)
      LIODBLE , allocatable :: Lmat(:,:)
      LIODBLE , allocatable :: Umat(:,:), Utrp(:,:)
      LIODBLE , allocatable :: Vmat(:,:), Vtrp(:,:)
      LIODBLE , allocatable :: Gmat(:,:), Ginv(:,:)
!                             G totally looks like a sigma

      integer              :: method_id = 0
      LIODBLE , allocatable :: Xmat(:,:), Xtrp(:,:)
      LIODBLE , allocatable :: Ymat(:,:), Ytrp(:,:)

   contains
      procedure, pass   :: Sets_smat
      procedure, pass   :: Gets_smat
      procedure, pass   :: Sets_orthog
      procedure, pass   :: Gets_orthog_2m
      procedure, pass   :: Gets_orthog_4m
      procedure, pass   :: Gets_eigens_m
      procedure, pass   :: Gets_eigens_v

      procedure, pass   :: Rebase_fockon
      procedure, pass   :: Rebase_denson
      procedure, pass   :: Rebase_fockao
      procedure, pass   :: Rebase_densao

      procedure, nopass :: Drop_ldvals
      procedure, nopass :: Calc_fulldcmp_4m
      procedure, nopass :: Calc_fulldcmp_7m

   end type sop

contains
#  include "Sets_smat.f90"
#  include "Gets_smat.f90"
#  include "Sets_orthog.f90"
#  include "Gets_orthog_2m.f90"
#  include "Gets_orthog_4m.f90"
#  include "Gets_eigens_m.f90"
#  include "Gets_eigens_v.f90"
#  include "Rebase_bothxx.f90"
#  include "Drop_ldvals.f90"
#  include "Calc_fulldcmp_4m.f90"
#  include "Calc_fulldcmp_7m.f90"

end module typedef_sop
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
