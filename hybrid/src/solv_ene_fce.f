c this subroutine calculates the solvent energy

      subroutine solv_ene_fce(natot,na_u,nac,ng1,rclas,Em,Rm,pc,
     .    Etot_amber,fcetot_amber,attype,
     .    nbond,nangle,ndihe,nimp,multidihe, multiimp,kbond,bondeq,
     .    kangle,angleeq,kdihe,diheeq,kimp,impeq,perdihe,perimp,
     .    bondtype,angletype,dihetype,imptype,
     .    bondxat,angexat,angmxat,dihexat,dihmxat,impxat,atange,
     .    atangm,atdihe,atdihm,atimp,
     .    ng1type,angetype,angmtype,evaldihe,evaldihm,
     .    dihety,dihmty,impty,nonbonded,scale,nonbondedxat,scalexat,
     .    evaldihelog,evaldihmlog,paso,nparm,
     .    actualiz,listcut,
     .    noat,noaa,sfc,timestep,
     .    water,masst,radblommbond)         

        implicit none
c      vbles grales
       integer   i,j,k,l,na_u,natot,nac,ng1(nac,6),paso 
       double precision  pcA(nac),rclas(3,natot),ramber(3,nac),
     .    EmA(nac),RmA(nac),pc(1:nac),Em(natot),Rm(natot),
c agregue 0 a pc, Nick
     .    Etot_amber,Ebond_amber,Eangle_amber,Edihe_amber,Eimp_amber,
     .    Elj_amber,Eelec_amber,Elj_amber14,Eelec_amber14 
       double precision fcetot_amber(3,nac),
     .   fcebond_amber(3,nac),fceangle_amber(3,nac),
     .   fcedihe_amber(3,nac),fceimp_amber(3,nac),
     .   fcelj_amber(3,nac),fceelec_amber(3,nac),
     .   fcelj_amber14(3,nac),fceelec_amber14(3,nac),
     .   fcetotaxes_amber(3)
       character  attype(nac)*4,noat(nac)*4,noaa(nac)*4
       double precision listcut,sfc,timestep,masst(natot),
     .   ewat,fwat(3,nac)
       logical water
c      vbles de los params de union
       integer   nbond,nangle,ndihe,nimp,nparm,multidihe(nparm),
     .    multiimp(nparm)
       double precision   kbond(nparm),bondeq(nparm),kangle(nparm),
     .   angleeq(nparm),kdihe(nparm),diheeq(nparm),kimp(nparm),
     .   impeq(nparm),perdihe(nparm),perdihe2(nparm),perimp(nparm)
       character   bondtype(nparm)*5,angletype(nparm)*8,
     .   dihetype(nparm)*11,imptype(nparm)*11
c      vbles de bond, angle, dihe e imp
       integer   bondxat(nac),angexat(nac),angmxat(nac),
     .   dihexat(nac),dihmxat(nac),impxat(nac)
       integer   atange(nac,25,2),atangm(nac,25,2),
     .   atdihe(nac,100,3),atdihm(nac,100,3),atimp(nac,25,4)
c	vbles q faltaban
       integer  ng1type(nac,6),angetype(nac,25),angmtype(nac,25),
     .          evaldihe(nac,100,5),evaldihm(nac,100,5),
     .          dihety(nac,100),dihmty(nac,100),impty(nac,25),
     .          nonbonded(nac,100),scale(nac,100),
     .          nonbondedxat(nac),scalexat(nac)
       logical  evaldihelog(nac,100),evaldihmlog(nac,100),actualiz
c parche para omitir interaccion entre extremos terminales
      double precision :: radblommbond


c inicializa las energias y fuerzas
      Etot_amber=0.d0
      Ebond_amber=0.d0
      Eangle_amber=0.d0
      Edihe_amber=0.d0
      Eimp_amber=0.d0
      Elj_amber=0.d0
      Eelec_amber=0.d0
      Elj_amber14=0.d0
      Eelec_amber14=0.d0
      fcetot_amber=0.d0
      fcebond_amber=0.d0
      fceangle_amber=0.d0
      fcedihe_amber=0.d0
      fceimp_amber=0.d0
      fcelj_amber=0.d0
      fceelec_amber=0.d0
      fcetotaxes_amber=0.d0   
      ewat=0.d0
      fwat=0.d0

c cambia variables
      k=1
      do j=1,nac
      pcA(k)=pc(j)
      k=k+1
      enddo
      k=1
      do i=na_u+1,natot
      ramber(1:3,k)=rclas(1:3,i)
      k=k+1
      enddo
      k=1
      do i=na_u+1,natot
      EmA(k)=Em(i)
      RmA(k)=Rm(i)
      k=k+1
      enddo
 
c  pasa a las unidades de Amber
      do i=1,nac
      RmA(i) = RmA(i)*(2.d0**(1.d0/6.d0))*0.529177d0/2.d0
      EmA(i) = EmA(i)*627.5108d0
      ramber(1:3,i)=ramber(1:3,i)*0.529177d0
      enddo

c  llama a subrutina q calcula la energia y fuerzas de bonds
        call amber_bonds(nac,ng1,ramber,Ebond_amber,attype,
     .       nbond,kbond,bondeq,bondtype,bondxat,fcebond_amber,
     .       ng1type,paso,nparm,radblommbond)

c  llama a subrutina q calcula la energia y fuerzas de angles
        call amber_angles(nac,ng1,ramber,Eangle_amber,attype,
     .       nangle,kangle,angleeq,angletype,angexat,angmxat,atange,
     .       atangm,fceangle_amber,angetype,angmtype,paso,nparm)

c  llama a subrutina q calcula la energia y fuerzas de dihedros     
	perdihe2=perdihe  
        call amber_dihes(nac,ng1,ramber,Edihe_amber,
     .            attype,ndihe,kdihe,diheeq,perdihe2,multidihe,
     .            dihetype,dihexat,dihmxat,atdihe,atdihm,
     .            fcedihe_amber,evaldihelog,evaldihe,dihety,
     .            evaldihmlog,evaldihm,dihmty,paso,nparm) 

c  llama a subrutina q calcula la energia y fuerzas de impropers
        call amber_improper(nac,ng1,ramber,Eimp_amber,attype,
     .       nimp,kimp,impeq,perimp,multiimp,imptype,impxat,atimp,
     .       fceimp_amber,impty,paso,nparm)

c  llama a subrutina q calcula la energia y fuerzas de terminos non-bonded
        call amber_nonbonded(nac,ng1,ramber,Elj_amber,
     .       Eelec_amber,Elj_amber14,Eelec_amber14,attype,
     .       EmA,RmA,pcA,bondxat,
     .       angexat,angmxat,atange,atangm,
     .       dihexat,dihmxat,atdihe,atdihm,
     .       fceelec_amber,fcelj_amber,nonbonded,scale,
     .       nonbondedxat,scalexat,paso,actualiz,
     .       listcut,
     .       noat,noaa,sfc,timestep,
     .       na_u,natot,rclas,water,masst,ewat,fwat)       


c  calcula la energia total de amber
       Etot_amber=Ebond_amber+Eangle_amber+Edihe_amber+Eimp_amber+
     .    Elj_amber+Eelec_amber+Elj_amber14+Eelec_amber14+ewat

c  calcula la fuerza total de amber
	do i=1,nac

	fcetot_amber(1:3,i)=fcebond_amber(1:3,i)+fceangle_amber(1:3,i)
     .  +fcedihe_amber(1:3,i)+fceimp_amber(1:3,i)+
     .  fcelj_amber(1:3,i)+fceelec_amber(1:3,i)+fwat(1:3,i)
       fcetot_amber(1:3,i)=(-1.d0)*fcetot_amber(1:3,i)     
       fcetotaxes_amber(1:3)=fcetotaxes_amber(1:3)+fcetot_amber(1:3,i)   
	enddo

      return
      end
c****************************************************************
c subrutina q calcula la energia y fuerzas de bonds

        subroutine amber_bonds(nac,ng1,ramber,Ebond_amber,
     .             attype,nbond,kbond,bondeq,bondtype,bondxat,
     .             fcebond_amber,ng1type,paso,nparm,radblommbond)

       implicit none      
c      vbles grales 
       integer   nac,ng1(nac,6),i,j,k,l,m,n,paso
       double precision   ramber(3,nac),Ebond_amber,
     .                    fcebond_amber(3,nac)
       character   attype(nac)*4
c      vbles de los params de union
       integer   nbond,nparm
       double precision   kbond(nparm),bondeq(nparm)
       character   bondtype(nparm)*5
c      vbles de bond, angle, dihe e imp
       integer   bondxat(nac)
c      vbles de asignacion
       character*4 ty1,ty2
       character*5 tybond
       integer ng1type(nac,6)
       double precision   rij,dist,dx1,dx2,dy1,dy2,dz1,dz2,
     .                    ftotbond(3)
      logical           first                                                                  
      save              first                                                                  
      data              first /.true./    
      logical           ST
       double precision :: radblommbond

	ST=.false.

c asignacion de tipos de union
      if(first) then
       do i=1,nac
c barre at clasicos
        do j=1,bondxat(i)
c barre bonds del atomo i
         do k=1,nbond
c barre todos los bonds leidos en el amber.parm
          tybond=bondtype(k)
          ty1=tybond(1:2)
          ty2=tybond(4:5)

!	write(*,*) "asigno tipo", attype(i), attype(ng1(i,j)), ng1(i,j)

          if(attype(i).eq.ty1.and.attype(ng1(i,j)).eq.ty2) then
            ng1type(i,j)=k
          elseif(attype(i).eq.ty2.and.attype(ng1(i,j)).eq.ty1) then
            ng1type(i,j)=k
          endif

         enddo
        enddo
       enddo
      first=.false.
      endif !asignacion

!	if (.false.) then !poner una variable para prender luego
       do i=1,nac
!check params, Nick
        do j=1,bondxat(i)
	  if (ng1type(i,j) .eq.0) then
	    write(*,*) "bond ", attype(i),"-", attype(ng1(i,j)),
     .  " not defined check amber.parm"
	  write(*,*) i,j, "bonds", bondxat(i), ng1(i,j)
	  ST=.true.
	  end if
	end do
	end do
	if (ST) STOP
!	end if

c calculo de la E y F de bond
       do i=1,nac
        do j=1,bondxat(i)
         rij=dist(ramber(1,i),ramber(2,i),ramber(3,i),
     .            ramber(1,ng1(i,j)),ramber(2,ng1(i,j)),
     .            ramber(3,ng1(i,j)))


          if (rij .lt. radblommbond) then
          Ebond_amber= Ebond_amber+kbond(ng1type(i,j))*
     .    (rij-bondeq(ng1type(i,j)))**2d0       
          dx1=(1.d0/rij)*(ramber(1,i)-ramber(1,ng1(i,j)))
          dx1=2.d0*kbond(ng1type(i,j))*(rij-bondeq(ng1type(i,j)))*dx1
          dy1=(1.d0/rij)*(ramber(2,i)-ramber(2,ng1(i,j))) 
          dy1=2.d0*kbond(ng1type(i,j))*(rij-bondeq(ng1type(i,j)))*dy1  
          dz1=(1.d0/rij)*(ramber(3,i)-ramber(3,ng1(i,j)))
          dz1=2.d0*kbond(ng1type(i,j))*(rij-bondeq(ng1type(i,j)))*dz1
          fcebond_amber(1,i)=fcebond_amber(1,i)+dx1
          fcebond_amber(2,i)=fcebond_amber(2,i)+dy1
          fcebond_amber(3,i)=fcebond_amber(3,i)+dz1
          else
c parche para omitir bonds entre extremos terminales
            Write(*,*) "WARNING, omiting bond ", i, j,
     .  "bond distance ", rij
          end if
        enddo
       enddo

       Ebond_amber= Ebond_amber/2d0
       ftotbond=0.d0

      end
c******************************************************
c  subrutina q calcula la energia y fuerzas de angles
 
        subroutine  amber_angles(nac,ng1,ramber,
     .   Eangle_amber,attype,nangle,kangle,angleeq,angletype,
     .   angexat,angmxat,atange,atangm,fceangle_amber,
     .   angetype,angmtype,paso,nparm)

        implicit none
c      vbles grales
       integer   nac,ng1(nac,6),i,j,k,l,m,n,paso
       double precision   ramber(3,nac),Eangle_amber,
     .                    fceangle_amber(3,nac)
       character   attype(nac)*4
c      vbles de los params de union
       integer   nangle,nparm
       double precision kangle(nparm),angleeq(nparm)
       character angletype(nparm)*8
c      vbles de bond, angle, dihe e imp
       integer   angexat(nac),angmxat(nac)
       integer   atange(nac,25,2),atangm(nac,25,2)
c      vbles de asignacion
       character*4 ty1,ty2,ty3
       character*8 tyangle
       integer angetype(nac,25),angmtype(nac,25)
       double precision  angulo,angle,pi,dx,dy,dz,scal,r12,r32,
     .                   scalar,ftotangle(3),dr12r32,dscalar,dist,
     .                   fesq(3,nac),fmedio(3,nac)
      logical           first                                                                  
      save              first                                                                  
      data              first /.true./    
      
       pi=DACOS(-1.d0)
       fesq=0.d0
       fmedio=0.d0
       ftotangle=0.d0

c asignacion de los tipos de angulos 
      if(first) then
       do i=1,nac
        do j=1,angexat(i)
         do k=1,nangle 
          tyangle=angletype(k)
          ty1=tyangle(1:2)
          ty2=tyangle(4:5)
          ty3=tyangle(7:8)
          if(attype(i).eq.ty1.and.attype(atange(i,j,1)).eq.ty2.and.
     .       attype(atange(i,j,2)).eq.ty3) then
          angetype(i,j)=k
          elseif(attype(i).eq.ty3.and.attype(atange(i,j,1)).eq.ty2.and.
     .       attype(atange(i,j,2)).eq.ty1) then
          angetype(i,j)=k
          endif
         enddo
        enddo
       enddo

       do i=1,nac
        do j=1,angmxat(i)
         do k=1,nangle
          tyangle=angletype(k)
          ty1=tyangle(1:2)
          ty2=tyangle(4:5)
          ty3=tyangle(7:8)
          if(attype(i).eq.ty2.and.attype(atangm(i,j,1)).eq.ty1.and.
     .       attype(atangm(i,j,2)).eq.ty3) then
          angmtype(i,j)=k
          elseif(attype(i).eq.ty2.and.attype(atangm(i,j,1)).eq.ty3.and.
     .       attype(atangm(i,j,2)).eq.ty1) then
          angmtype(i,j)=k
          endif
         enddo
        enddo
       enddo
      first=.false.
      endif ! asignacion

c calcula la E y F de angles
c para el angulo en la esquina
       do i=1,nac
        do j=1,angexat(i)
        angulo=angle(ramber(1,i),ramber(2,i),ramber(3,i),
     .          ramber(1,atange(i,j,1)),ramber(2,atange(i,j,1)),
     .          ramber(3,atange(i,j,1)),ramber(1,atange(i,j,2)),
     .          ramber(2,atange(i,j,2)),ramber(3,atange(i,j,2)))
       Eangle_amber = Eangle_amber + kangle(angetype(i,j))*
     .                ((angulo-angleeq(angetype(i,j)))*
     .                (pi/180d0))**2d0
      scal=scalar(ramber(1,i),ramber(2,i),ramber(3,i),
     .          ramber(1,atange(i,j,1)),ramber(2,atange(i,j,1)),
     .          ramber(3,atange(i,j,1)),ramber(1,atange(i,j,2)),
     .          ramber(2,atange(i,j,2)),ramber(3,atange(i,j,2)))
      r12=dist(ramber(1,i),ramber(2,i),ramber(3,i),
     .          ramber(1,atange(i,j,1)),ramber(2,atange(i,j,1)),
     .          ramber(3,atange(i,j,1)))
      r32=dist(ramber(1,atange(i,j,2)),
     .          ramber(2,atange(i,j,2)),ramber(3,atange(i,j,2)),
     .          ramber(1,atange(i,j,1)),ramber(2,atange(i,j,1)),
     .          ramber(3,atange(i,j,1)))
      dscalar=(ramber(1,atange(i,j,2))-ramber(1,atange(i,j,1)))
      dr12r32=r32*(ramber(1,i)-ramber(1,atange(i,j,1)))/(r12)
      dx=(dscalar*r12*r32-scal*dr12r32)/(r12*r32)**(2.d0)
      dx=-1.d0/(sqrt(1.d0-(scal/(r12*r32))**2.d0))*dx
      dx=2.d0*kangle(angetype(i,j))*(angulo-angleeq(angetype(i,j)))*
     .                (pi/180d0)*dx
      dscalar=(ramber(2,atange(i,j,2))-ramber(2,atange(i,j,1)))
      dr12r32=r32*(ramber(2,i)-ramber(2,atange(i,j,1)))/(r12)
      dy=(dscalar*r12*r32-scal*dr12r32)/(r12*r32)**(2.d0)
      dy=-1.d0/(sqrt(1.d0-(scal/(r12*r32))**2.d0))*dy
      dy=2.d0*kangle(angetype(i,j))*(angulo-angleeq(angetype(i,j)))*
     .                (pi/180d0)*dy
      dscalar=(ramber(3,atange(i,j,2))-ramber(3,atange(i,j,1)))
      dr12r32=r32*(ramber(3,i)-ramber(3,atange(i,j,1)))/(r12)
      dz=(dscalar*r12*r32-scal*dr12r32)/(r12*r32)**(2.d0)
      dz=-1.d0/(sqrt(1.d0-(scal/(r12*r32))**2.d0))*dz
      dz=2.d0*kangle(angetype(i,j))*(angulo-angleeq(angetype(i,j)))*
     .                (pi/180d0)*dz
      fesq(1,i)=fesq(1,i)+dx
      fesq(2,i)=fesq(2,i)+dy   
      fesq(3,i)=fesq(3,i)+dz  
        enddo
       enddo  

c para el angulo en el medio   
       do i=1,nac
        do j=1,angmxat(i)
         angulo=angle(ramber(1,atangm(i,j,1)),ramber(2,atangm(i,j,1)),
     .          ramber(3,atangm(i,j,1)),ramber(1,i),ramber(2,i),
     .          ramber(3,i),ramber(1,atangm(i,j,2)),
     .          ramber(2,atangm(i,j,2)),ramber(3,atangm(i,j,2)))
       Eangle_amber = Eangle_amber + kangle(angmtype(i,j))*
     .                ((angulo-angleeq(angmtype(i,j)))*
     .                (pi/180d0))**2d0
      scal=scalar(ramber(1,atangm(i,j,1)),ramber(2,atangm(i,j,1)),
     .          ramber(3,atangm(i,j,1)),ramber(1,i),ramber(2,i),
     .          ramber(3,i),ramber(1,atangm(i,j,2)),
     .          ramber(2,atangm(i,j,2)),ramber(3,atangm(i,j,2)))
      r12=dist(ramber(1,atangm(i,j,1)),ramber(2,atangm(i,j,1)),
     .          ramber(3,atangm(i,j,1)),ramber(1,i),ramber(2,i),
     .          ramber(3,i))
      r32=dist(ramber(1,atangm(i,j,2)),
     .          ramber(2,atangm(i,j,2)),ramber(3,atangm(i,j,2)),
     .          ramber(1,i),ramber(2,i),ramber(3,i))
      dscalar=2.d0*ramber(1,i)-ramber(1,atangm(i,j,1))-
     .        ramber(1,atangm(i,j,2))
      dr12r32=(r32*(-ramber(1,atangm(i,j,1))+ramber(1,i))/r12)+
     .        (r12*(-ramber(1,atangm(i,j,2))+ramber(1,i))/r32)
      dx=(dscalar*r12*r32-scal*dr12r32)/(r12*r32)**(2.d0)
      dx=-1.d0/(sqrt(1.d0-(scal/(r12*r32))**2.d0))*dx
      dx=2.d0*kangle(angmtype(i,j))*(angulo-angleeq(angmtype(i,j)))*
     .                (pi/180d0)*dx
      dscalar=2.d0*ramber(2,i)-ramber(2,atangm(i,j,1))-
     .        ramber(2,atangm(i,j,2))
      dr12r32=(r32*(-ramber(2,atangm(i,j,1))+ramber(2,i))/r12)+
     .        (r12*(-ramber(2,atangm(i,j,2))+ramber(2,i))/r32)
      dy=(dscalar*r12*r32-scal*dr12r32)/(r12*r32)**(2.d0)
      dy=-1.d0/(sqrt(1.d0-(scal/(r12*r32))**2.d0))*dy
      dy=2.d0*kangle(angmtype(i,j))*(angulo-angleeq(angmtype(i,j)))*
     .                (pi/180d0)*dy
      dscalar=2.d0*ramber(3,i)-ramber(3,atangm(i,j,1))-
     .        ramber(3,atangm(i,j,2))
      dr12r32=(r32*(-ramber(3,atangm(i,j,1))+ramber(3,i))/r12)+
     .        (r12*(-ramber(3,atangm(i,j,2))+ramber(3,i))/r32)
      dz=(dscalar*r12*r32-scal*dr12r32)/(r12*r32)**(2.d0)
      dz=-1.d0/(sqrt(1.d0-(scal/(r12*r32))**2.d0))*dz
      dz=2.d0*kangle(angmtype(i,j))*(angulo-angleeq(angmtype(i,j)))*
     .                (pi/180d0)*dz
      fmedio(1,i)=fmedio(1,i)+dx
      fmedio(2,i)=fmedio(2,i)+dy
      fmedio(3,i)=fmedio(3,i)+dz
        enddo
       enddo
    
       Eangle_amber= Eangle_amber/3d0
 
       do i=1,nac
       fceangle_amber(1,i)=fesq(1,i)+fmedio(1,i)
       fceangle_amber(2,i)=fesq(2,i)+fmedio(2,i)    
       fceangle_amber(3,i)=fesq(3,i)+fmedio(3,i)    
       enddo
        
	end
c****************************************************************
c  subrutina q calcula la energia y fuerzas de dihedros
 
      subroutine  amber_dihes(nac,ng1,ramber,Edihe_amber,
     .            attype,ndihe,kdihe,diheeq,perdihe,multidihe,
     .            dihetype,dihexat,dihmxat,atdihe,atdihm,
     .            fcedihe_amber,evaldihelog,evaldihe,dihety,
     .            evaldihmlog,evaldihm,dihmty,paso,nparm)
	
        implicit none
 
c      vbles grales
       integer   nac,ng1(nac,6),i,j,k,l,m,n,paso
       double precision   ramber(3,nac),Edihe_amber,
     .                    fcedihe_amber(3,nac)
       character   attype(nac)*4
c      vbles de los params de union
       integer   ndihe,nparm,multidihe(nparm)
       double precision kdihe(nparm),diheeq(nparm),perdihe(nparm)
       character dihetype(nparm)*11
c      vbles de bond, angle, dihe e imp
       integer   dihexat(nac),dihmxat(nac)
       integer   atdihe(nac,100,3),atdihm(nac,100,3)
c      vbles de asignacion
       character*4 ty1,ty2,ty3,ty4
       character*11 tydihe
       integer dihety(nac,100),dihmty(nac,100),evaldihe(nac,100,5),
     .         evaldihm(nac,100,5)
       double precision  pi,dihedro,dihedral,E1,dist,
     .                   dx,dy,dz,ftotdihe(3),
     .                   fesq(3,nac),fmedio(3,nac),
     .			 fce(12)
       logical evaldihelog(nac,100),evaldihmlog(nac,100)
      logical           first                                                                  
      save              first                                                                  
      data              first /.true./    

       pi=DACOS(-1.d0)
       ftotdihe=0.d0
       fesq=0.d0
       fmedio=0.d0

c asignacion de los tipos de dihedros
      if(first) then
       evaldihelog=.false.
       evaldihmlog=.false.

       do i=1,nac
        do j=1,dihexat(i)
	dihety(i,j)=1
         m=0
         do k=1,ndihe
         tydihe=dihetype(k)
         ty1=tydihe(1:2)
         ty2=tydihe(4:5)   
         ty3=tydihe(7:8)    
         ty4=tydihe(10:11)   
         if(ty1.eq.'X ') then
	 if(attype(atdihe(i,j,1)).eq.ty2.and.
     .      attype(atdihe(i,j,2)).eq.ty3)  then
         dihety(i,j)=k
         elseif(attype(atdihe(i,j,1)).eq.ty3.and.
     .      attype(atdihe(i,j,2)).eq.ty2)  then
         dihety(i,j)=k
         endif
         elseif(ty1.ne.'X ') then
         if(attype(i).eq.ty1.and.attype(atdihe(i,j,1)).eq.
     .   ty2.and.attype(atdihe(i,j,2)).eq.ty3.and.
     .   attype(atdihe(i,j,3)).eq.ty4) then
         dihety(i,j)=k            
	  if(perdihe(k).lt.0) then
          evaldihelog(i,j)=.true.
          m=m+1
 	  evaldihe(i,j,m)=k
          evaldihe(i,j,m+1)=k+1
          endif
         elseif(attype(i).eq.ty4.and.attype(atdihe(i,j,1)).eq.
     .   ty3.and.attype(atdihe(i,j,2)).eq.ty2.and.
     .   attype(atdihe(i,j,3)).eq.ty1) then
         dihety(i,j)=k
          if(perdihe(k).lt.0) then
          evaldihelog(i,j)=.true.
          m=m+1
          evaldihe(i,j,m)=k
          evaldihe(i,j,m+1)=k+1    
          endif
         endif
         endif
         enddo

	if(dihety(i,j).eq.1) then
C	write(*,*) 'dihety: ',i,j,'sin parametro'
	endif

        enddo
       enddo

       do i=1,nac
        do j=1,dihmxat(i)
	dihmty(i,j)=1
         m=0
         do k=1,ndihe
         tydihe=dihetype(k)
         ty1=tydihe(1:2)
         ty2=tydihe(4:5)
         ty3=tydihe(7:8)
         ty4=tydihe(10:11)
         if(ty1.eq.'X ') then
         if(attype(i).eq.ty2.and.
     .      attype(atdihm(i,j,2)).eq.ty3)  then
         dihmty(i,j)=k
         elseif(attype(i).eq.ty3.and.
     .      attype(atdihm(i,j,2)).eq.ty2)  then
         dihmty(i,j)=k
         endif
         elseif(ty1.ne.'X ') then
         if(attype(atdihm(i,j,1)).eq.ty1.and.attype(i).eq.
     .   ty2.and.attype(atdihm(i,j,2)).eq.ty3.and.
     .   attype(atdihm(i,j,3)).eq.ty4) then
         dihmty(i,j)=k
          if(perdihe(k).lt.0.d0) then
          evaldihmlog(i,j)=.true.
          m=m+1
          evaldihm(i,j,m)=k
          evaldihm(i,j,m+1)=k+1    
          endif
         elseif(attype(atdihm(i,j,1)).eq.ty4.and.attype(i).eq.
     .   ty3.and.attype(atdihm(i,j,2)).eq.ty2.and.
     .   attype(atdihm(i,j,3)).eq.ty1) then
         dihmty(i,j)=k
          if(perdihe(k).lt.0.d0) then
          evaldihmlog(i,j)=.true.
          m=m+1
          evaldihm(i,j,m)=k
          evaldihm(i,j,m+1)=k+1    
          endif
         endif
         endif
         enddo
        if(dihmty(i,j).eq.1) then
C        write(*,*) 'dihmty: ',i,j,'sin parametro'
        endif
        enddo
       enddo
      first=.false.
      endif !asignacion

c calcula la E y F de dihedros 
c para los dihes en la esquina
        do i=1,ndihe
        if(perdihe(i).lt.0.d0) then
        perdihe(i)=-perdihe(i)
        endif
        enddo
        do i=1,nac
         do j=1,dihexat(i)
         dihedral=dihedro(ramber(1,i),ramber(2,i),ramber(3,i),
     .   ramber(1,atdihe(i,j,1)),ramber(2,atdihe(i,j,1)),
     .   ramber(3,atdihe(i,j,1)),
     .   ramber(1,atdihe(i,j,2)),ramber(2,atdihe(i,j,2)),
     .   ramber(3,atdihe(i,j,2)), 
     .   ramber(1,atdihe(i,j,3)),ramber(2,atdihe(i,j,3)),
     .   ramber(3,atdihe(i,j,3)))

c si el dihedro es 500 es error y sale
        if(dihedral.lt.500.1.and.dihedral.gt.499.9) then
        write(*,*) 'ERROR EN EL DIHEDRO(esq)',i,j
        STOP
        endif

       if(evaldihelog(i,j)) then
	do m=1,5
	 if(evaldihe(i,j,m).ne.0) then

       k=evaldihe(i,j,m)
       E1=(kdihe(k)/multidihe(k))*
     .  (1+dCOS((pi/180d0)*(perdihe(k)*dihedral-diheeq(k))))
	Edihe_amber=Edihe_amber+E1
	call diheforce(nac,ramber,
     .                 i,atdihe(i,j,1),atdihe(i,j,2),atdihe(i,j,3),1,
     .		       kdihe(k),diheeq(k),perdihe(k),multidihe(k),fce)    
      dx=fce(1)
      dy=fce(2)
      dz=fce(3)
      fesq(1,i)=fesq(1,i)+dx
      fesq(2,i)=fesq(2,i)+dy
      fesq(3,i)=fesq(3,i)+dz
	endif
	enddo
	else
       k=dihety(i,j)
       E1=(kdihe(k)/multidihe(k))*
     .  (1+dCOS((pi/180d0)*(perdihe(k)*dihedral-diheeq(k))))
	Edihe_amber=Edihe_amber+E1
        call diheforce(nac,ramber,
     .                 i,atdihe(i,j,1),atdihe(i,j,2),atdihe(i,j,3),1,
     .                 kdihe(k),diheeq(k),perdihe(k),multidihe(k),fce)
	dx=fce(1)
        dy=fce(2)
        dz=fce(3)
      fesq(1,i)=fesq(1,i)+dx
      fesq(2,i)=fesq(2,i)+dy
      fesq(3,i)=fesq(3,i)+dz
        endif
        enddo
       enddo       

c para los dihes en el medio
       do i=1,nac
         do j=1,dihmxat(i)
         dihedral=dihedro(
     .   ramber(1,atdihm(i,j,1)),ramber(2,atdihm(i,j,1)),
     .   ramber(3,atdihm(i,j,1)),
     .   ramber(1,i),ramber(2,i),ramber(3,i),
     .   ramber(1,atdihm(i,j,2)),ramber(2,atdihm(i,j,2)),
     .   ramber(3,atdihm(i,j,2)),
     .   ramber(1,atdihm(i,j,3)),ramber(2,atdihm(i,j,3)),
     .   ramber(3,atdihm(i,j,3)))

c si el dihedro es 500 es error y sale
	if(dihedral.lt.500.1.and.dihedral.gt.499.9) then
	write(*,*) 'ERROR EN EL DIHEDRO(dihm)',i,j
	STOP
	endif

      if(evaldihmlog(i,j)) then
        do m=1,5
         if(evaldihm(i,j,m).ne.0) then
       k=evaldihm(i,j,m)
       E1=(kdihe(k)/multidihe(k))*
     .  (1+dCOS((pi/180d0)*(perdihe(k)*dihedral-diheeq(k))))
	Edihe_amber=Edihe_amber+E1
        call diheforce(nac,ramber,
     .                 atdihm(i,j,1),i,atdihm(i,j,2),atdihm(i,j,3),2,
     .                 kdihe(k),diheeq(k),perdihe(k),multidihe(k),fce)
      dx=fce(4)
      dy=fce(5)
      dz=fce(6)	
      fmedio(1,i)=fmedio(1,i)+dx
      fmedio(2,i)=fmedio(2,i)+dy
      fmedio(3,i)=fmedio(3,i)+dz
        endif
        enddo
        else

       k=dihmty(i,j)
       E1=(kdihe(k)/multidihe(k))*
     .  (1+dCOS((pi/180d0)*(perdihe(k)*dihedral-diheeq(k))))
	Edihe_amber=Edihe_amber+E1
        call diheforce(nac,ramber,
     .                 atdihm(i,j,1),i,atdihm(i,j,2),atdihm(i,j,3),2,
     .                 kdihe(k),diheeq(k),perdihe(k),multidihe(k),fce)
      dx=fce(4)
      dy=fce(5)
      dz=fce(6)
      fmedio(1,i)=fmedio(1,i)+dx
      fmedio(2,i)=fmedio(2,i)+dy
      fmedio(3,i)=fmedio(3,i)+dz
        endif
        enddo
       enddo

      Edihe_amber=Edihe_amber/4
      do i=1,nac
      fcedihe_amber(1,i)=fesq(1,i)+fmedio(1,i)
      fcedihe_amber(2,i)=fesq(2,i)+fmedio(2,i)
      fcedihe_amber(3,i)=fesq(3,i)+fmedio(3,i)
      enddo 

       end
c******************************************************************
c  subrutina q calcula la energia y fuerzas de impropers 

       subroutine amber_improper(nac,ng1,ramber,Eimp_amber,attype,
     .            nimp,kimp,impeq,perimp,multiimp,imptype,impxat,
     .            atimp,fimp,impty,paso,nparm)

        implicit none
c      vbles grales
       integer   nac,ng1(nac,6),i,j,k,l,m,n,paso
       double precision   ramber(3,nac),Eimp_amber
       character   attype(nac)*4
c      vbles de los params de union
       integer   nimp,nparm,multiimp(nparm)
       double precision kimp(nparm),impeq(nparm),perimp(nparm)
       character imptype(nparm)*11 
c      vbles de bond, angle, dihe e imp
       integer   impxat(nac)
       integer   atimp(nac,25,4)
c      vbles de asignacion
       character*4 ty1,ty2,ty3,ty4
       character*11 tyimp
       integer impty(nac,25)
       double precision  pi,dihedro,dihedral
c	varianles de fza
	double precision fimp(3,nac),fce(12),fimptot(3),dx,dy,dz
	integer atom
      logical           first                                                                  
      save              first                                                                  
      data              first /.true./
    
       pi=DACOS(-1.d0)
       fimptot=0.d0

c asignacion de los tipos de impropers
      if(first) then
      do i=1,nac
       do j=1,impxat(i)
	impty(i,j)=1
        do k=1,nimp
        tyimp=imptype(k)
        ty1=tyimp(1:2)
        ty2=tyimp(4:5)
        ty3=tyimp(7:8)
        ty4=tyimp(10:11)
	if(ty1.eq.'X '.and.ty2.eq.'X '.and.ty4.eq.'X ') then
	       if(attype(atimp(i,j,3)).eq.ty3) then
        impty(i,j)=k
		elseif(attype(atimp(i,j,2)).eq.ty3) then
        impty(i,j)=k
		endif
        elseif(ty1.eq.'X '.and.ty2.eq.'X ') then
        if(attype(atimp(i,j,3)).eq.ty3.and.
     .     attype(atimp(i,j,4)).eq.ty4) then
        impty(i,j)=k
        elseif(attype(atimp(i,j,1)).eq.ty4.and.
     .     attype(atimp(i,j,2)).eq.ty3) then
        impty(i,j)=k   
        endif
        elseif(ty1.eq.'X ') then
        if(attype(atimp(i,j,2)).eq.ty2.and.
     .     attype(atimp(i,j,3)).eq.ty3.and.
     .     attype(atimp(i,j,4)).eq.ty4) then
        impty(i,j)=k   
        elseif(attype(atimp(i,j,3)).eq.ty2.and.
     .     attype(atimp(i,j,2)).eq.ty3.and.
     .     attype(atimp(i,j,1)).eq.ty4) then
        impty(i,j)=k
        endif
        else
        if(attype(atimp(i,j,1)).eq.ty1.and.
     .     attype(atimp(i,j,2)).eq.ty2.and.
     .     attype(atimp(i,j,3)).eq.ty3.and.
     .     attype(atimp(i,j,4)).eq.ty4) then
        impty(i,j)=k
        elseif(attype(atimp(i,j,4)).eq.ty1.and.
     .     attype(atimp(i,j,3)).eq.ty2.and.
     .     attype(atimp(i,j,2)).eq.ty3.and.
     .     attype(atimp(i,j,1)).eq.ty4) then
        impty(i,j)=k
        endif
        endif
        enddo
       if(impty(i,j).eq.1) then
C        write(*,*) 'impty: ',i,j,'sin parametro'
        endif
       enddo
      enddo
      first=.false.
      endif !asignacion

c calcula la E y F de impropers
      do i=1,nimp
        if(perimp(i).lt.0.d0) then
        perimp(i)=-perimp(i)
        endif
        enddo
        do i=1,nac
         do j=1,impxat(i)
         dihedral=dihedro(ramber(1,atimp(i,j,1)),
     .   ramber(2,atimp(i,j,1)),ramber(3,atimp(i,j,1)),
     .   ramber(1,atimp(i,j,2)),ramber(2,atimp(i,j,2)),
     .   ramber(3,atimp(i,j,2)),
     .   ramber(1,atimp(i,j,3)),ramber(2,atimp(i,j,3)),
     .   ramber(3,atimp(i,j,3)),
     .   ramber(1,atimp(i,j,4)),ramber(2,atimp(i,j,4)),
     .   ramber(3,atimp(i,j,4)))

c si el dihedro es 500 es error y sale
        if(dihedral.lt.500.1.and.dihedral.gt.499.9) then
        write(*,*) 'ERROR EN EL IMPROP',i,j
        STOP
        endif

       k=impty(i,j)
       Eimp_amber=Eimp_amber+(kimp(k)/multiimp(k))*
     .  (1+COS((pi/180)*(perimp(k)*dihedral-impeq(k))))

c busca que atomo del impropio es el i	
	if (i.eq.atimp(i,j,1)) atom=1
	if (i.eq.atimp(i,j,2)) atom=2 
	if (i.eq.atimp(i,j,3)) atom=3 
	if (i.eq.atimp(i,j,4)) atom=4 
	       call diheforce(nac,ramber,
     .         atimp(i,j,1),atimp(i,j,2),atimp(i,j,3),atimp(i,j,4),atom,
     .                 kimp(k),impeq(k),perimp(k),multiimp(k),fce)
      dx=fce((atom-1)*3+1)
      dy=fce((atom-1)*3+2)
      dz=fce((atom-1)*3+3)
      fimp(1,i)=fimp(1,i)+dx
      fimp(2,i)=fimp(2,i)+dy
      fimp(3,i)=fimp(3,i)+dz
        enddo
       enddo
      Eimp_amber=Eimp_amber/4d0

      end
c *******************************************************
c subrutina que calcula la energia y fuerzas de terminos non-bonded

       subroutine amber_nonbonded(nac,ng1,ramber,Elj_amber,
     .       Eelec_amber,Elj_amber14,Eelec_amber14,attype,
     .       Em,Rm,pc,bondxat,
     .       angexat,angmxat,atange,atangm,
     .       dihexat,dihmxat,atdihe,atdihm,
     .       felec,flj,
     .       nonbonded,scale,nonbondedxat,scalexat,paso,
     .       actualiz,listcut,noat,noaa,
     .       sfc,timestep,
     .       na_u,natot,rclas,water,masst,ewat,fwat)       

        implicit none
c      vbles grales
       integer   nac,ng1(nac,6),i,j,k,l,m,n,paso,dimvec,
     . n_pointer,na_u,natot
       double precision   ramber(3,nac),Em(nac),Rm(nac),pc(1:nac),
     .    Elj_amber,Eelec_amber,Elj_amber14,Eelec_amber14
       character*4 attype(nac),noat(nac),noaa(nac)
       double precision ewat,fwat(3,nac),masst(natot),rclas(3,natot)
       logical water
c      vbles de bond, angle, dihe e imp
       integer   bondxat(nac),angexat(nac),angmxat(nac),
     .   dihexat(nac),dihmxat(nac)
       integer   atange(nac,25,2),atangm(nac,25,2),
     .   atdihe(nac,100,3),atdihm(nac,100,3)
c      vbles de asignacion      
       integer nonbonded(nac,100),scale(nac,100),scalexat(nac),
     . nonbondedxat(nac),acs
       double precision A,B,Rij,Eij,dist,distancia,unidades,
     .                  factorlj,factorelec,E1,E2,pi,epsilon,
     .                  dist2,distancia2,rcut,Rij6,distancia2_3,
     .                  fac,dfac,ca,cb,cc,cd,sfc,timestep,e1f,e2f,
     .                  x0,x1,rinn,rout
       logical scale14,nonbond,actualiz 
c variables asoc a las fzas
       double precision dx1,dy1,dz1,dx2,dy2,dz2,felec(3,nac),
     .  felectot(3),fel,flj(3,nac),fljtot(3),listcut
c variables de la lista de vecinos
       integer, dimension(:), allocatable ::         
     .  veclist, veclistemp, veclistxat
       save veclist, veclistxat
      logical           first                                        
      save              first                                    
      data              first /.true./    

	pi=DACOS(-1.d0)
        unidades=((1.602177d-19**2d0)*6.02d23)/(1.0d-10*4d0*
     .     pi*8.8541878d-12*4.184d0*1000d0)
        fel=0.d0
	flj=0.d0
	felec=0.d0
        felectot=0.d0
        fljtot=0.d0
        epsilon=1.d0
        factorlj=2.d0
        factorelec=1.2d0
	acs=200	
	ewat=0.d0
	fwat=0.d0
c	dimvec=nac*3000
c cambio dim vec, Nick  antes alocateaba mas chico de lo q podeia llegar a necesitar
	dimvec=nac*(nac+1)/2
c	if(dimvec.gt.100000000) stop 
c     .  'Solvent Energy and Forces: "dimvec" too large!'


c asigna los atomos nonbonded y su tipo
         if(first) then       
           do i=1,nac 
             m=1
             do j=1,bondxat(i)
	       if(i.lt.ng1(i,j))then
                 nonbonded(i,m)=ng1(i,j)
                 m=m+1
	       endif
             enddo

             do j=1,angexat(i)
	       if(i.lt.atange(i,j,2)) then
	         nonbonded(i,m)=atange(i,j,2)
	         m=m+1
               endif
	     enddo
	     nonbondedxat(i)=m-1

c se fija los que estan 1-4(scaled)
             m=1
             do j=1,dihexat(i)
	       if(i.lt.atdihe(i,j,3)) then
c se fija si ya lo puso (por si hay dihedros repetidos)
	         do n=1,m-1	  
	           if(atdihe(i,j,3).eq.scale(i,n)) goto 10  
	         enddo
c se fija si no lo puso en nonbonded (se no es tambien 1-3)
	         do n=1,nonbondedxat(i)
	           if(atdihe(i,j,3).eq.nonbonded(i,n)) goto 10	
	         enddo	

                 scale(i,m)=atdihe(i,j,3)
                 m=m+1
	       endif
 10          enddo
	     scalexat(i)=m-1
	   enddo
          actualiz=.true.
	  first=.false.
	endif !asignacion

c ahora crea la lista de vecinos si estan dentro de listcut+sfc, 
c se debe actualizar cada 20 fs. 
	acs=int(20/timestep)
	if(mod(paso,acs).eq.0.or.paso.eq.1) actualiz=.true.

c llama a la sub q agrega el water restrain potential
	if(water.and.actualiz) then
        call waters(na_u,nac,natot,rclas,masst,noaa,noat,ewat,fwat)
	endif

c si listcut > 100A, la lista se hace SOLO el 1er paso 
	if(listcut.ge.100.and.paso.ne.1) actualiz=.false.

c crea la lista veclist y el indice veclistxat
	if(actualiz) then
	actualiz=.false.
c	write(6,*) 'Neighbour list actualization in step:',paso

        if(allocated(veclistxat)) deallocate(veclistxat)
        allocate(veclistxat(nac))
        if(allocated(veclist)) deallocate(veclist)
        allocate(veclistemp(dimvec))

	veclistemp=0
	veclistxat=0
 	rcut=(listcut+sfc+2.d0)**2d0
	n=1
	do i=1,nac
	do j=i+1,nac
	   do k=1,nonbondedxat(i)
             if(nonbonded(i,k).eq.j) goto 30
           enddo
           do m=1,scalexat(i)
             if(scale(i,m).eq.j) goto 30
           enddo
	if(noaa(j).eq.'HOH'.and.noat(j).eq.'O') then
		distancia2=dist2(ramber(1,i),ramber(2,i),ramber(3,i),
     .          ramber(1,j),ramber(2,j),ramber(3,j))
		if(distancia2.le.rcut) then
			veclistemp(n)=j
			n=n+1
			veclistemp(n)=j+1
			n=n+1
			veclistemp(n)=j+2
			n=n+1
		endif
	else if(noaa(j).ne.'HOH') then
		distancia2=dist2(ramber(1,i),ramber(2,i),ramber(3,i),
     .          ramber(1,j),ramber(2,j),ramber(3,j))
		if(distancia2.le.rcut) then
!aca pincha con cuts grandes, Nick
!la dimension maxima es nac*(nac+1)/2
                        veclistemp(n)=j
                        n=n+1
		endif
	endif

c fin loop de la lista para cada atomo
 30  	enddo
	veclistxat(i)=n-1
c se fija si la dimension de veclist es suficiente
        if((n-1).gt.dimvec) then
      write(6,*)'Dimension Neighbour list (required, used)=',n-1,dimvec
      write(6,*) 'Solvent Energy and Forces: Stopping Program'
        STOP
        endif
c fin loop de lista para todos atomos
	enddo !nac

c alocatea veclist posta
        allocate(veclist(n-1))
        veclist(1:n-1)=veclistemp(1:n-1)
        deallocate(veclistemp)
	endif !actualiz

c  calculo de la E y F de terminos non-bonded
c loop scale nonbonded

c        write(*,*) "flag conect"
c        do i=1,nac 
c                write(*,*) i, scalexat(i)
c        end do
c	STOP
c test Nick
!	WRITE(*,*) "SCALE, FLAG"
!	do i=1,nac
!	  write(*,*) i, scalexat(i), scale(i,1:scalexat(i))
!	end do

	do i=1,nac
	 do k=1,scalexat(i)
	 j=scale(i,k)	
                dx1=ramber(1,i)-ramber(1,j)
                dy1=ramber(2,i)-ramber(2,j)
                dz1=ramber(3,i)-ramber(3,j)
                distancia2 = dx1*dx1 + dy1*dy1  + dz1*dz1
                distancia = sqrt(distancia2)
                Rij=Rm(i)+Rm(j)
                Eij=sqrt(Em(i)*Em(j))
                Rij6 = Rij**6d0
                distancia2_3 = distancia2**3d0
                E1 = Eij*Rij6/distancia2_3*((Rij6/distancia2_3)-2.d0)
		E1=E1/factorlj
                Elj_amber14=Elj_amber14+E1
	fel = -12.d0*Eij*Rij6/distancia2**4d0*(Rij6/distancia2_3 - 1.d0)
                fel = fel/factorlj
                dx2=-dx1
                dy2=-dy1
                dz2=-dz1
                flj(1,i)=flj(1,i)+dx1*fel
                flj(2,i)=flj(2,i)+dy1*fel
                flj(3,i)=flj(3,i)+dz1*fel
                flj(1,j)=flj(1,j)+dx2*fel
                flj(2,j)=flj(2,j)+dy2*fel
                flj(3,j)=flj(3,j)+dz2*fel
c		write(*,*) "E2 callc", i,j,pc(i),pc(j)
                E2=((pc(i)*pc(j))/distancia)*unidades/epsilon
              	E2=E2/factorelec 
		Eelec_amber14=Eelec_amber14+E2
                fel=-E2/distancia2
		dx2=-dx1
                dy2=-dy1
                dz2=-dz1
                felec(1,i)=felec(1,i)+dx1*fel
                felec(2,i)=felec(2,i)+dy1*fel
                felec(3,i)=felec(3,i)+dz1*fel
                felec(1,j)=felec(1,j)+dx2*fel
                felec(2,j)=felec(2,j)+dy2*fel
                felec(3,j)=felec(3,j)+dz2*fel
        enddo
	enddo


c fin scaled nonbonden
c loop nonscaled nonbonded
        n_pointer=1      
        x0=listcut
        x1=listcut+sfc
	cb=-sfc/x0**2d0
	ca=-cb/2.d0
	cc=ca
	cd=cc-1.d0/x0
	rinn=x0**2d0
	rout=x1**2d0

!	write(*,*) "nonbonded flag"
!	do i=1,nac
!	  write(*,*) i, veclistxat(i), veclist(1:veclistxat(i))
!	end do

        do i=1,nac
	 do k=n_pointer,veclistxat(i)
          j=veclist(k)
                dx1=ramber(1,i)-ramber(1,j)
                dy1=ramber(2,i)-ramber(2,j)
                dz1=ramber(3,i)-ramber(3,j)
		distancia2 = dx1*dx1 + dy1*dy1  + dz1*dz1
                if(distancia2.le.rinn) then
		distancia = sqrt(distancia2)
		Rij=Rm(i)+Rm(j)
                Eij=sqrt(Em(i)*Em(j))
		Rij6 = Rij**6d0
		distancia2_3 = distancia2**3d0
		E1 = Eij*Rij6/distancia2_3*((Rij6/distancia2_3)-2.d0)
                Elj_amber=Elj_amber+E1
		fel = -12.d0*Eij*Rij6/distancia2**4d0*(Rij6/distancia2_3 - 1.d0)
                dx2=-dx1
                dy2=-dy1
                dz2=-dz1
                flj(1,i)=flj(1,i)+dx1*fel
                flj(2,i)=flj(2,i)+dy1*fel
                flj(3,i)=flj(3,i)+dz1*fel
                flj(1,j)=flj(1,j)+dx2*fel
                flj(2,j)=flj(2,j)+dy2*fel
                flj(3,j)=flj(3,j)+dz2*fel
                E2=(pc(i)*pc(j)*unidades/epsilon)*(1.d0/distancia)
		E2F=E2+(pc(i)*pc(j)*unidades/epsilon*cd)
		Eelec_amber=Eelec_amber+E2F
		fel=-E2/distancia2
                dx2=-dx1
                dy2=-dy1
                dz2=-dz1
                felec(1,i)=felec(1,i)+dx1*fel
                felec(2,i)=felec(2,i)+dy1*fel
                felec(3,i)=felec(3,i)+dz1*fel
                felec(1,j)=felec(1,j)+dx2*fel
                felec(2,j)=felec(2,j)+dy2*fel
                felec(3,j)=felec(3,j)+dz2*fel
		elseif(distancia2.gt.rinn.and.distancia2.lt.rout) then
                distancia = sqrt(distancia2)     
		Rij=Rm(i)+Rm(j)
                Eij=sqrt(Em(i)*Em(j))
                Rij6 = Rij**6d0
                distancia2_3 = distancia2**3d0
                E1 = Eij*Rij6/distancia2_3*((Rij6/distancia2_3)-2.d0)
                Elj_amber=Elj_amber+E1
        fel = -12.d0*Eij*Rij6/distancia2**4d0*(Rij6/distancia2_3 - 1.d0)
                dx2=-dx1
                dy2=-dy1
                dz2=-dz1
                flj(1,i)=flj(1,i)+dx1*fel
                flj(2,i)=flj(2,i)+dy1*fel
                flj(3,i)=flj(3,i)+dz1*fel
                flj(1,j)=flj(1,j)+dx2*fel
                flj(2,j)=flj(2,j)+dy2*fel
                flj(3,j)=flj(3,j)+dz2*fel
		fac=(distancia-listcut)/sfc
                E2=pc(i)*pc(j)*unidades/epsilon
                E2F=E2*(ca*fac**2d0+cb*fac+cc)
                Eelec_amber=Eelec_amber+E2F
                fel=E2/distancia
		fel=fel*(2d0*ca*fac+cb)/sfc
                dx2=-dx1
                dy2=-dy1
                dz2=-dz1
                felec(1,i)=felec(1,i)+dx1*fel
                felec(2,i)=felec(2,i)+dy1*fel
                felec(3,i)=felec(3,i)+dz1*fel
                felec(1,j)=felec(1,j)+dx2*fel
                felec(2,j)=felec(2,j)+dy2*fel
                felec(3,j)=felec(3,j)+dz2*fel
		endif

c fin del loop sobre bondlist
          enddo
	   n_pointer = veclistxat(i) + 1
c fin del loop sobre todos los atomos
         enddo
        end
c******************************************************************
c subrutina q calcula el water restarin potential

	subroutine waters(na_u,nac,natot,rclas,masst,noaa,noat,ewat,fwat)
	implicit none
	integer, intent(in) :: na_u,nac,natot

        integer watlist(2000),watlistnum
        double precision rclas(3,natot),ewat,rwat,masscenter(3),
     .  rt(3),rij,E,fwat(3,nac),dx,dy,dz,masst(natot),
     .  kte,ramber(3,natot),dist,dist2,mdist
        character noat(nac)*4,noaa(nac)*4
        integer i,j,k,l,m,n
        double precision pi

        pi=DACOS(-1.d0)
        kte=200.d0
        dist2=0.d0
        mdist=0.d0

c calcula el masscenter del sistema y el rwat
        ramber(1:3,1:natot)=rclas(1:3,1:natot)*0.529177d0
        masscenter=0.d0 
        do i=1,natot
        masscenter(1:3)=masscenter(1:3)+masst(i)*ramber(1:3,i)
        enddo
        masscenter=masscenter/natot
        do i=1,natot
        dist2=(ramber(1,i)-masscenter(1))**2d0+
     .        (ramber(2,i)-masscenter(2))**2d0+
     .        (ramber(3,i)-masscenter(3))**2d0
        if(dist2.gt.mdist) mdist=dist2
        enddo

	rwat=sqrt(mdist) - 2.5d0        
	write(6,'(a,2x,f8.4)') 'Water Cutoff Radius:', rwat

c calculo la matrix con las aguas xa los at MM
        ramber=0.d0
        ramber(1:3,1:nac)=rclas(1:3,na_u+1:natot)*0.529177d0
        k=1
        do i=1,nac
        if(noaa(i).eq.'HOH'.and.noat(i).eq.'O') then
        rij=dist(ramber(1,i),ramber(2,i),ramber(3,i),
     .  masscenter(1),masscenter(2),masscenter(3))
        if(rij.gt.rwat) then
        watlist(k)=i
        k=k+1
c fin de si esta en la zona buffer
        endif
c fin de si es agua
        endif
        enddo
        watlistnum=k-1
        ewat=0.d0
        fwat=0.d0
c calculo de la Ene y la fza para el potencial de las aguas Ewat = 0.0
        do j=1,watlistnum
        i = watlist (j)
        rij = dist(ramber (1,i),ramber(2,i),ramber(3,i),
     .  masscenter(1),masscenter(2), masscenter(3))
        if(rij.gt.rwat) then
        E = kte*((rij-rwat)**2)
        ewat = ewat + E
        dx = (1.d0/rij)*(ramber(1,i)-masscenter(1))
        dx = 2.d0*kte*(rij-rwat)*dx
        dy = (1.d0/rij)*(ramber(2,i)-masscenter(2))
        dy =  2.d0*kte*(rij-rwat)*dy
        dz = (1.d0/rij)*(ramber(3,i)-masscenter(3))
        dz =  2.d0*kte*(rij-rwat)*dz
       fwat(1,i)=fwat(1,i)-dx
       fwat(2,i)=fwat(2,i)-dy
       fwat(3,i)=fwat(3,i)-dz
        endif
        enddo
 
        end
c********************************************************************************

