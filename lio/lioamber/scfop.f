c SCF subroutine ----------------------------------
c DIRECT VERSION
c calls all integrals generator subroutines : 1 el integrals,
c 2 el integrals, exchange fitting , so it gets S matrix, F matrix
c and P matrix in lower storage mode ( symmetric matrices)
c
c Dario Estrin, 1992
c---------------------------------------------------
       subroutine SCFOP(coords,E)
       use garcha_mod
c      use qmmm_module, only : qmmm_struct, qmmm_nml

c
      implicit real*8 (a-h,o-z)
       dimension q(natom)
       REAL*8 , intent(inout)  :: coords(ntatom*3)
       real*8, dimension (:,:), ALLOCATABLE ::xnano
       real*8, dimension (:), ALLOCATABLE :: rmm5,rmm15,rmm13
c      dimension d(ntq,ntq)

c
c------------------------------------------------------------------
        write(*,*) 'no anda con capa abierta!!!!'
        return
        end
c
c Pointers
c
