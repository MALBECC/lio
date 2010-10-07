      SUBROUTINE UNGLES(NATSOL,IT,NIN,NDIP)
      INCLUDE 'param'
        INCLUDE 'COMM'
      DIMENSION DTOT (100,100)
      DIMENSION MSGF (100,100)
      REAL*8 BUG
      dimension Tzo(3),Xri(20)
      dimension Xam(50),Yam(50),Zam(50)

c---Calcula CM, Racm=masa total        
       Xacm=0.0D0
       Yacm=0.0D0
       Zacm=0.0D0
       Racm=0.0D0

       nofset=0
       do  i=1,natom+natsol
       if(i.gt.1) nofset=nofset+nnat(i-1)
       do  j=1+nofset,nnat(i)+nofset
       Racm=Racm+wwm(i)
       Xacm=Xacm+wwm(i)*X(j)
       Yacm=Yacm+wwm(i)*Y(j)
       Zacm=Zacm+wwm(i)*Z(j)
       enddo
       enddo
       
       Xacm=Xacm/Racm
       Yacm=Yacm/Racm
       Zacm=Zacm/Racm

       do i=1,10
       Xri(i)=0.0D0
       enddo

c---Pone origen en CM  
c---CM = punto fijo
       nofset=0
       do  i=1,natom+natsol
       if(i.gt.1) nofset=nofset+nnat(i-1)
       do  j=1+nofset,nnat(i)+nofset
       
       xa1=X(j)-Xacm
       xa2=xa1**2
       ya1=Y(j)-Yacm
       ya2=ya1**2
       za1=Z(j)-Zacm
       za2=za1**2
       r2=xa2+ya2+za2

c---Coeficientes del tensor de inercia
       Xri(1)=Xri(1)+wwm(i)*(r2-xa2)
       Xri(4)=Xri(4)+wwm(i)*(r2-ya2)
       Xri(6)=Xri(6)+wwm(i)*(r2-za2)
 
       Xri(2)=Xri(2)-wwm(i)*xa1*ya1
       Xri(3)=Xri(3)-wwm(i)*xa1*za1
       Xri(5)=Xri(5)-wwm(i)*ya1*za1
       enddo
       enddo

       
c       write(*,*)'oxigenos'
c       write(*,*)xri(1),xri(2),xri(3)
c       write(*,*)'     ',xri(4),xri(5)
c       write(*,*)'     ','     ',xri(6)
c       pause 

c---Autovectores
c---Si ITEL=0, esto no importa, los lee despues
c---de diagonalizar ;
c---Si ITEL.NE.0, son los calculados anteriormente      
       X11 = Xmi(1,1) 
       X12 = Xmi(2,1) 
       X13 = Xmi(3,1) 
       X21 = Xmi(1,2) 
       X22 = Xmi(2,2) 
       X23 = Xmi(3,2) 
       X31 = Xmi(1,3) 
       X32 = Xmi(2,3)
       X33 = Xmi(3,3)


c---Diagonalizacion
c---(dspev es = eig.f)
       call dspev('V','L',3,Xri,Tzo,Xmi,3,Xri(7),info)

       if (ITEL.EQ.ZERO) then
       X11 = Xmi(1,1) 
       X12 = Xmi(2,1) 
       X13 = Xmi(3,1) 
       X21 = Xmi(1,2) 
       X22 = Xmi(2,2) 
       X23 = Xmi(3,2) 
       X31 = Xmi(1,3) 
       X32 = Xmi(2,3)
       X33 = Xmi(3,3)
       endif
       
c       write(*,*)'despues de ungles'
c       do i=1,3
c       do j=1,3
cc       write(*,*)i,j,xmi(i,j)
c       enddo
c       enddo
c       pause

c---Producto escalar del autovector viejo con el nuevo
       Vec1= Xmi(1,1)*X11+Xmi(2,1)*X12+Xmi(3,1)*X13
       Vec2= Xmi(1,2)*X21+Xmi(2,2)*X22+Xmi(3,2)*X23
       Vec3= Xmi(1,3)*X31+Xmi(2,3)*X32+Xmi(3,3)*X33

c       write(*,*)'Productos escalares'
c       write(*,*)Vec1,Vec2,Vec3
c       pause


       if (Vec1.LT.ZERO) then
       Xmi(1,1)= -Xmi(1,1)
       Xmi(2,1)= -Xmi(2,1)
       Xmi(3,1)= -Xmi(3,1)
c       write(*,*)'si entra aca es porque Vec1 .lt.zero'
       else
       Xmi(1,1)= Xmi(1,1)
       Xmi(2,1)= Xmi(2,1)
       Xmi(3,1)= Xmi(3,1)
c       write(*,*)' auvc1 ',xmi(1,1),xmi(2,1),xmi(3,1) 
       endif
       if (Vec2.LT.ZERO) then
       Xmi(1,2)= -Xmi(1,2)
       Xmi(2,2)= -Xmi(2,2)
       Xmi(3,2)= -Xmi(3,2)
c       write(*,*)'si entra aca es porque Vec2 .lt.zero'
       else
       Xmi(1,2)= Xmi(1,2)
       Xmi(2,2)= Xmi(2,2)
       Xmi(3,2)= Xmi(3,2)
c       write(*,*)' auvc2 ',xmi(1,2),xmi(2,2),xmi(3,2) 
       endif
       if (Vec3.LT.ZERO) then
       Xmi(1,3)= -Xmi(1,3)
       Xmi(2,3)= -Xmi(2,3)
       Xmi(3,3)= -Xmi(3,3)
c       write(*,*)'si entra aca es porque Vec3 .lt.zero'
       else
       Xmi(1,3)= Xmi(1,3)
       Xmi(2,3)= Xmi(2,3)
       Xmi(3,3)= Xmi(3,3)
c       write(*,*)' auvc1 ',xmi(1,3),xmi(2,3),xmi(3,3) 
       endif


c---Histogramas   
c       INER1= NINT((Tzo(1)-500)/10)+1
c       INER1= MIN(INER1,150)
c       HSINR1(INER1)= HSINR1(INER1)+1

c       INER2= NINT((Tzo(2)-500)/10)+1
c       INER2= MIN(INER2,150)
c       HSINR2(INER2)= HSINR2(INER2)+1

c       INER3= NINT((Tzo(3)-400)/4)+1
c       INER3= MIN(INER3,100)
c       HSINR3(INER3)= HSINR3(INER3)+1

c---Escribe autovectores
       if (mod((it-nin),ndip).EQ.0.and.ndip.gt.0) then
c       WRITE (30,9665) RMK,Tzo(1),Tzo(2),Tzo(3)
       write (30,*)itel
       write (30,*)'V1=',Xmi(1,1),Xmi(2,1),Xmi(3,1)
       write (30,*)'V2=',Xmi(1,2),Xmi(2,2),Xmi(3,2)
       write (30,*)'V3=',Xmi(1,3),Xmi(2,3),Xmi(3,3)
       write (30,*)
       endif

      RETURN

950   format(F10.5,10x,I3)
960   format(F10.5,10x,F10.6,4x,F10.6,4x,F10.6,6x,F10.6)
970   format(1F8.5)
990   format(12x,5F8.5)
9665  format(F18.5,3x,F9.5,2x,F9.5,2x,F9.5)
7743  format(F5.3)

      END



