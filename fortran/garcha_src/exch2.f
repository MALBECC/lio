
c This subroutine deals with the exchange-correlation energy term
c VERSION , USING WITH WRITE.F,
c SAVING BASIS IN DISK
c Output: the vector containing all linear coefficients for the
c potential ( Fock matrix part) and energy (It is used in program
c where the integrals are evaluated)
c-----------------------------------------------------------------
       subroutine exch2(OPEN,Iz,natom,RMM,nshelld,M,Md,M17,
     >                  NCOa,NCOb,B1)
      implicit real*8 (a-h,o-z)
      logical dens1,SVD,integ,OPEN
      integer iconst,igrid,igrid2
      INCLUDE 'param'
      dimension AFUNC(ngd),RMM(*),Nr(0:54),nshelld(0:3)
c
c Number of shells for Least-Squares Fit
      data Nr /20,20,20,25,25,25,25,25,25,25,25,
     >  30,30,30,30,30,30,30,30,
     >  35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,
     >  40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40,40/
c
c output
c fit for exchange correlation density and potential if different
      dimension B1(ngd,3),Iz(natom)
c B1(i,1) exchange correlation density
c B1(i,2) exchange correlation potential
c
      common /fit/ Nang,dens1,integ,Iexch,igrid,igrid2
      common /Sys/ SVD,iconst
      common /Nc/ Ndens
c
      open(unit=8,file='scratch/fit1',form='unformatted')
      open(unit=9,file='scratch/fit2',form='unformatted')
c
c
c initialization
        do 1 k=1,Md
         B1(k,3)=0.D0
         B1(k,1)=0.D0
 1       B1(k,2)=0.D0
c
c
        M18=M17+Md*(Md+1)/2
        NCO=NCOa
c
c loop 12 , over all grid  -----------------------------
c GRID WAS GENERATED BEFORE (write subroutine and is not changed
c during SCF iterations)
c yi density functional for xi ( all necessary data in common)
c
      DO 12 na=1,natom
c
       do 16 n=1,Nr(Iz(na))
c
       do 15 nag=1,Nang
c
         if (OPEN) then
c
         call DNS2OP(DENSA,DENSB,M,M18,NCO,NCOb,RMM)
         call potop(Iexch,DENSA,DENSB,yiex,yiec,y2a,y2b)
         dxi=DENSA+DENSB
         DENS=dxi
c
        else
        call DNS2(DENS,M,M18,NCO,RMM)
        dxi=DENS
c
        call pot(Iexch,dxi,yiex,yiec,y2i)
        endif
c   
        yi = yiex + yiec
*
c reads from disk all fitting functions at xi ------
c
        call CDE(AFUNC,Md,8)
c
c---------------------------------------------------
c -- THIS PART IS THE NORMAL EQUATION METHOD -----
c
        if (OPEN) then
        DO 118 j=1,Md
c   tt = AFUNC but weight was included in write subroutine
         tt=AFUNC(j)
         B1(j,1)=B1(j,1)+yi*tt
         B1(j,2)=B1(j,2)+y2a*tt
         B1(j,3)=B1(j,3)+y2b*tt
c
 118     CONTINUE
c
        else
        DO 11 j=1,Md
c   tt = AFUNC but weight was included in write subroutine
         tt=AFUNC(j)
         B1(j,1)=B1(j,1)+yi*tt
         B1(j,2)=B1(j,2)+y2i*tt
c
 11      CONTINUE
c
         endif
 15    CONTINUE
c
 16   continue
 12   continue
c-------------------------------------------------------
c
c
c ESSL OPTION
#ifdef essl
c
      CALL DPPS(RMM(M17),Md,B1(1,1),1)
      CALL DPPS(RMM(M17),Md,B1(1,2),1)
c
      if (OPEN) then
      CALL DPPS(RMM(M17),Md,B1(1,3),1)
      endif
c
#endif
c LAPACK OPTION
c
#ifdef pack
      call dppsl(RMM(M17),Md,B1(1,1))
      call dppsl(RMM(M17),Md,B1(1,2))
c
      if (OPEN) then
      call dppsl(RMM(M17),Md,B1(1,3))
      endif
c
#endif
c -------------------------------------------
c-------------------------------------------------------
      rewind(8)
      rewind(9)
c
c
      return
      END
c-------------------------------------------------------------
