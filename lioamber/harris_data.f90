#include "datatypes/datatypes.fh"
module harris_data 
    implicit none 

    LIODBLE, allocatable    :: Tmat_vec(:)              ! Vectorized form of the kinetic energy matrix
    LIODBLE                 :: kinE = 0.D0              ! Kinetic energy
end module harris_data 
