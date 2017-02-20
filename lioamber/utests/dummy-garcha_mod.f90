module garcha_mod

    integer              :: natom, nsol, M, NCO
    integer, allocatable :: nuc(:), Iz(:)
 
    real*8               :: Enucl
    real*8, allocatable  :: rmm(:), Eorbs(:)
    real*8, allocatable  :: x(:,:), smat(:,:), sqsm(:,:), realrho(:,:)

    logical              :: mulliken, lowdin, open

endmodule garcha_mod
