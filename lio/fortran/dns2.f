c function calculating density functionals
c for version writing the values of functions in the grid
c on disk ( ONLY for local density functionals)
c
      SUBROUTINE DNS2(Dens,M,M18,NCO,RMM)
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
c
      DENS=0.D0
c basis functions evaluated at r are read from disk
c
c
        call CDE(F,M,9)
c
c---------------------------------------------------
c
c now calculation of vector W : density matrix scalar F
c
      do 1 i=1,M
        W(i)=0.D0
 1    continue
c
      if (Ndens.eq.1) then
      k=0
      do 50 j=1,M
      do 50 i=j,M
       k=k+1
 50   W(j)=W(j)+RMM(k)*F(i)
c
      do 60 i=1,M
       DENS=DENS+F(i)*W(i)
  60  continue
c
      return
c 
      else
c
      kk=0
      do 51 jj=1,NCO
      do 51 ii=1,M
       kk=kk+1
 51    W(jj)=W(jj)+RMM(M18+kk-1)*F(ii)
c
      do 61 ii=1,NCO
 61    DENS=DENS+W(ii)**2
c
      DENS=DENS*2.D0
c
      return
      endif
      end
c---------------------------------------------------
c subroutine that reads vectors from disk
c
      subroutine CDE(FF,MM,i)
       implicit real*8 (a-h,o-z)
       dimension FF(MM)
c
       read(i) FF
       return
       end
c---------------------------------------------------- 
