c function calculating density functionals
c for version writing the values of functions in the grid
c on disk ( ONLY for local density functionals)
c
      SUBROUTINE DNS2OP(Densa,Densb,M,M18,NCOa,NCOb,RMM)
c
      implicit real*8 (a-h,o-z)
      INCLUDE 'param'
      dimension F(ng),W(ng),RMM(*)
c
      common /Nc/ Ndens
c now we should evaluate all same loops as the ones used for
c 1 electron matrix elements, but doing only products
c then, the particular density functional wanted is calculated
c
      M18b=M18+M*NCOa
c
      DENSA=0.D0
      DENSB=0.D0
c basis functions evaluated at r are read from disk
c
c
        call CDE(F,M,9)
c
c---------------------------------------------------
c
c now calculation of vector W : density matrix scalar F
c
c alpha spin
      kk=M18-1
      do 51 jj=1,NCOa
      W(jj)=0.0D0
      do 51 ii=1,M
       kk=kk+1
 51    W(jj)=W(jj)+RMM(kk)*F(ii)
c
      do 61 ii=1,NCOa
 61    DENSA=DENSA+W(ii)**2
c
c beta spin
      kk=M18b-1
      do 151 jj=1,NCOb
      W(jj)=0.0D0
      do 151 ii=1,M
       kk=kk+1
 151    W(jj)=W(jj)+RMM(kk)*F(ii)
c
      do 161 ii=1,NCOb
 161    DENSB=DENSB+W(ii)**2
c
      return
      end
c---------------------------------------------------
