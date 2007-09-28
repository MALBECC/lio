c MC Subroutine ----------------------------------
c
c
c calls SCF subroutine for calculating energy at each geometry ,
c After that, a move is done, and the process is continued
c
c all specifications given in the namelist SCF are necessary
c
c Dario Estrin, 1992
c---------------------------------------------------
       subroutine MC(MEMO,NORM,natom,Iz,r,v,Nuc,M,ncont,nshell,c,a,
     >     Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,E,
     >    nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write,
          NMCSTEPS,DELT,TEMP)
c
      implicit real*8 (a-h,o-z)
      integer nopt,iconst,igrid,iforce
      logical NORM,ATRHO,VCINP,DIRECT,EXTR,dens,write,SVD
      logical SHFT,OPEN,GRAD,integ,field,sol,free,MEMO
      INCLUDE 'param'
      dimension r(nt,3),nshelld(0:3),nshell(0:3)
      dimension cd(ngd,nl),ad(ngd,nl),Nucd(Md),ncontd(Md)
      dimension c(ng,nl),a(ng,nl),Nuc(M),ncont(M),Iz(nt)
      dimension Em(ntq+nss),Rm(ntq+nss),pc(nss),alpha(nss)
c
c
c
c
c nuclear MD
      dimension Pm(nt)
      dimension RMM(*)
     
c
      dimension xold(nt),yold(nt),zold(nt)
      dimension xnew(nt),ynew(nt),znew(nt)
c auxiliars
c X scratch space in matrices
      dimension X(M,3*M),XX(Md,Md)
c
c
      COMMON /TABLE/ STR(880,0:21)
      common /fit/ Nang,dens,integ,Iexch,igrid,igrid2
      common /dyn/ h,nsteps,Pm,Nscale
      common /Sys/ SVD,iconst
      common /index/ index(ng)
      common /force/ iforce
      common/Ngeom/ ngeo
      common /ENum/ GRAD
      common /cav/ a0,epsilon,field
c------------------------------------------------------------------
c 
c
c Pointers
c
      NCOa=NCO
      NCOb=NCO+Nunp
      iforce=0
c
      Nel=NCOa+NCOb
      GRAD=.false.
      MM=M*(M+1)/2 
      MM2=M**2
      MMd=Md*(Md+1)/2
      M2=2*M
c first P
      M1=1
c now Pnew
      M3=M1+MM
c now S, F also uses the same position after S was used
      M5=M3+MM
c now G
      M7=M5+MM
c now Gm
      M9=M7+MMd
c now H
      M11=M9+MMd
c W ( eigenvalues ), also this space is used in least squares
      M13=M11+MM
c aux ( vector for ESSl)
      M15=M13+M
cLeast squares
      M17=M15+MM
c
c vectors of MO alpha
      M18=M17+MMd
c vectors of MO beta
      M18b=M18+M*NCOa
c

c---------aca empieza MC

c        poner coordenadas de entrada en xold
c           
             do i=1,natom
              xold(i)=r(i,1)
              yold(i)=r(i,2)
              zold(i)=r(i,3)
             enddo
cccccccccccccccccccccccccccccccccccccccccccccccc
c              LEER LA CONFIGURACION INICAL Y DETERMINAR SI LA 
c              CORRIDA ES NUEVA, (LATTICE), VIEJA DESORDENADA O
c              CONTINUACION.
c
c              SI ES CONTINUACION:  GUARDAR PREVIAMENTE: AC
c             # de pasos hechos, IACC,IACC1,IREJ,IREJC,DELT
c              DELT= MAXIMO DESPLAZAMIENTO POSIBLE DE CADA NUCLEO EN CADA
c                    COORDENADA

               ONE = 1.D0
               TWO = 2.D0
               BETA=ONE/(BOLTZ * TEMP)
c
c seteo de ENGOLD
c              ENGOLD=0.0D0
c                 CUANDO LEE LAS CONFIGURACIONES INICIALES GENERA UN 
c                 ENGOLD
c
c aca inicia ciclo MONTE CARLO cccccccccccccccccccccccccccccccc
c
              DO 1000  I = 1, NMCSTEPS
                 DO J = 1, NATOM

                  XNEW(J) = XOLD(J) + DELT*(TWO*RANF()-ONE)
                  YNEW(J) = YOLD(J) + DELT*(TWO*RANF()-ONE)
                  ZNEW(J) = ZOLD(J) + DELT*(TWO*RANF()-ONE)

                 ENDDO
c pasa a coordenadas necesarias para calculo SCF
c
             do i=1,natom
              r(i,1)=xnew(i)
              r(i,2)=ynew(i)
              r(i,3)=znew(i)
             enddo
cccccccccccccccccccccccccccccc

c   llamado a calculo de energia, mediante calculo SCF
      if (.not.OPEN) then
      call SCF(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >         Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,ENGNEW,
     > nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write)
      else
      call SCFop(MEMO,NORM,natom,Iz,r,Nuc,M,ncont,nshell,c,a,
     >         Nucd,Md,ncontd,nshelld,cd,ad,RMM,X,XX,ENGNEW,
     > nopt,OPEN,NMAX,NCO,ATRHO,VCINP,SHFT,Nunp,GOLD,told,write)
      endif
c


                  DELTA = BETA*(ENGNEW-ENGOLD)

          
                  IF (DELTA.LT.ZERO) THEN

22                 IACC = IACC + 1
                   IACC1 = IACC1 + 1

                    ENGOLD = ENGNEW


                  DO J = 1, NATOM

                    XOLD(J) = XNEW(J)
                    YOLD(J) = YNEW(J)
                    ZOLD(J) = ZNEW(J)

                  ENDDO


                   ELSEIF (RANF().LT.DELTA) GOTO 22

                     IREJ = IREJ + 1
                     IREJ1 = IREJ1 + 1

                   ENDIF


                    CALL AVERAGES ( XOLD, YOLD, ZOLD, ENGOLD)  

                    IF (MOD(I,ICHNG).EQ.0) THEN

                       RATIO = PTFIVE-DBLE(IACC1)/DBLE(ICHNG)

                       IF (RATIO.LT.0.4) THEN
                       DELT = DELT * 0.95
                       ELSEIF (RATION.GT.0.6) THEN
                       DELT = DELT * 1.05
                       ENDIF

                    IACC1 = 0

                    ENDIF
                      
C          CADA TANTO GRABAR LOS AC


1000           CONTINUE

C          CADA TANTO GRABAR LOS AC

C          ESCRIBIR CUALQUIER COSA
           return
            END
