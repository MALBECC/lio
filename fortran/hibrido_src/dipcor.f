       SUBROUTINE DIPCOR(HISTO,PMAX,PZMAX,DELP,DELPZ)
       INCLUDE 'param'
        INCLUDE 'COMM'
       DIMENSION HISTO(100,100),PX(130),PY(130),PZ(130)


C=================
       I1=4
       I2=6
C=================

C---POSICION DE OXIGENO DONOR: I1
       RCX = X(I1) 
       RCY = Y(I1) 
       RCZ = Z(I1)

C---POSICION DE OXIGENO ACEPTOR: I2
       RZX = X(I2)
       RZY = Y(I2)
       RZZ = Z(I2)

C---EJE I2-a-I1
       RVX = (RCX - RZX) 
       RVY = (RCY - RZY)
       RVZ = (RCZ - RZZ)

       RVMOD =  DSQRT(RVX**2+RVY**2+RVZ**2)

C---VERSOR EN LA DIRECCION DE TRANSFERENCIA

       RVX = RVX  / RVMOD 
       RVY = RVY  / RVMOD 
       RVZ = RVZ  / RVMOD 

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


C---ORIGEN DE COORDENADAS

       RCX = PTFIVE*(X(I1) + X(I2))
       RCY = PTFIVE*(Y(I1) + Y(I2))
       RCZ = PTFIVE*(Z(I1) + Z(I2))



C---DISTANCIA DEL OXIGENO I DEL SOLVENTE AL OXIGENO DONOR


       DO 1000 I = 1+NATOM,NWAT+NATOM
       IN1 = I + NWAT
       IN2 = IN1 + NWAT
       RDX = X(I) - RCX
       RDY = Y(I) - RCY 
       RDZ = Z(I) - RCZ

       RADII = DSQRT(RDX**2+RDY**2+RDZ**2)

C---COORDENADA Z (CILINDRICA) DEL ATOMO I

       ZPZ = RDX *RVX + RDY * RVY + RDZ * RVZ 

C---COORDENADA RADIAL (CILINDRICA) DEL ATOMO I

       RADIAL = DSQRT( RADII**2 - ZPZ**2)

        IZ = NINT(ZPZ/DELPZ)+100
        IR = NINT(RADIAL/DELP)

	  IZ = MIN(200,IZ)
	  IZ = MAX(1,IZ)

	  IR = MIN(100,IR)
	  IR = MAX(1,IR)

C---DIPOLO
       PX(I) = PC(IN1)*(X(IN1)-X(I)) + PC(IN2)*(X(IN2)-X(I))
       PY(I) = PC(IN1)*(Y(IN1)-Y(I)) + PC(IN2)*(Y(IN2)-Y(I))
       PZ(I) = PC(IN1)*(Z(IN1)-Z(I)) + PC(IN2)*(Z(IN2)-Z(I))
       AMOD = ONE / 
     &    DSQRT(PX(I)**2 + PY(I)**2 + PZ(I)**2)

C---VERSOR EN LA DIRECCION DEL DIPOLO

       PX(I) = PX(I) * AMOD
       PY(I) = PY(I) * AMOD
       PZ(I) = PZ(I) * AMOD


       PRDV = PX(I)*RVX + PY(I)* RXY + PZ(I)*RVZ
c       write(*,*)'dipcor ',prdv


C---COS DEL ANGULO QUE FORMA EL DIPOLO DEL AGUA I CON
C   LA DIRECCION DE TRANSFERENCIA.

       HISTO(IZ,IR) = HISTO(IZ,IR) + PRDV


1000   CONTINUE


       RETURN
       END


	 



