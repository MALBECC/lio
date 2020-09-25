#include "../datatypes/datatypes.fh"
module typedef_operator

   type operator
     LIODBLE, allocatable     :: data_AO(:,:)
     LIODBLE, allocatable     :: data_ON(:,:)
     TDCOMPLEX, allocatable  :: dataC_AO(:,:)
     TDCOMPLEX, allocatable  :: dataC_ON(:,:)

   contains
     procedure, pass :: Sets_data_AO
     procedure, pass :: Gets_data_AO
     procedure, pass :: Sets_dataC_AO
     procedure, pass :: Gets_dataC_AO
     procedure, pass :: Sets_data_ON
     procedure, pass :: Gets_data_ON
     procedure, pass :: Sets_dataC_ON
     procedure, pass :: Gets_dataC_ON
     procedure, pass :: Diagon_datamat
     procedure, pass :: Dens_build
     procedure, pass :: Commut_data_r
     procedure, pass :: Commut_data_c
     procedure, pass :: Shift_diag_ON
     procedure, pass :: Shift_diag_AO
     procedure, pass :: BChange_AOtoON_r
     procedure, pass :: BChange_ONtoAO_r
     procedure, pass :: BChange_AOtoON_x
     procedure, pass :: BChange_ONtoAO_x
     procedure, pass :: purify_ON
     generic         :: BChange_AOtoON => BChange_AOtoON_r
     generic         :: BChange_AOtoON => BChange_AOtoON_x
     generic         :: BChange_ONtoAO => BChange_ONtoAO_r
     generic         :: BChange_ONtoAO => BChange_ONtoAO_x
   end type operator

contains
#include  "Sets_datamat.f90"
#include  "Gets_datamat.f90"
#include  "Diagon_datamat.f90"
#include  "Dens_build.f90"
#include  "Commut_data.f90"
#include  "BChange_data.f90"
#include  "Shift_diag.f90"
#include  "Idempotency_purify.f90"

end module typedef_operator
