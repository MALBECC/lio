       module garcha_mod

       implicit real*8 (a-h,o-z)

      INCLUDE 'param'
      logical verbose 
      integer M,Md,natom,ntatom,NMAX,NCO,NUNP,igrid,igrid2
     >  ,Iexch,nsol,npas,npasw,idip,watermod,noconverge,
     > converge,ndiis,NGEO,nang,timedep,ntdstep,propagator,NBCH 
      real*8 GOLD, TOLD, qmmmcut, dgtrig
      parameter (nng=100)
      character*65 title
      character*20 basis,whatis,stdbas
      character*4 date
      character*20 output,fcoord,fmulliken,frestart,frestartin,solv,
     > solv2
      character*4 ctype
      logical exists,MEMO,predcoef
      logical done(ntq),used,NORM,OPEN,ATRHO,DIRECT,VCINP,SHFT,DIIS
      logical TMP1,TMP2,dens,EXTR,write1,SVD,ANG,field1
      logical Coul,Scf1,Prop,GRAD,BSSE,integ,SVD1,sol,tipe
      logical exter,exter1,resp1,popf,primera,writexyz,intsoldouble
      logical OPEN1
      logical dens1,integ1,sol1,free,free1, field, extern
      dimension OCC(40),oc2(400),ATCOEF(100*ng0),ighost(ntq),
     > ighost1(ntq)
      dimension ncf(nng),lt(nng)
      real*8 e_(50,3),wang(50),e_2(116,3),wang2(116),e3(194,3), ! intg1 e intg2
     > wang3(194),Ex                                               !
      integer Nr(0:54),Nr2(0:54)
      real*8 Fx, Fy, Fz, epsilon, a0,tdstep

      real*8, dimension (:,:), ALLOCATABLE :: r,v,rqm,d
      real*8, dimension (:), ALLOCATABLE ::  Em, Rm, pc
       integer, dimension (:), ALLOCATABLE :: Iz, nnat

      dimension isotop(54)!,Pm(nt)
      dimension Rm2(0:54), STR(880,0:21)
      dimension alpha(nss)
c Everything is dimensioned for 2 basis, normal and density
c ncf, lt,at,ct parameters for atomic basis sets
      dimension at(nng),ct(nng),nshell(0:4)
      dimension Num(0:3),nlb(ng),nld(ngd),nshelld(0:4)
       integer iconst1,idip1,ipop1,ispin1,
     > icharge1,Nsol1,natsol1,Ll(3)
      
        real*8, dimension (:), ALLOCATABLE :: af
       real*8, dimension (:,:), ALLOCATABLE :: c,a,cx,ax,cd,ad,B
       integer, dimension (:), ALLOCATABLE :: Nuc,ncont,Nucx,ncontx,Nucd
     >  ,ncontd
        integer, dimension (:), ALLOCATABLE :: indexii, indexiid

c
       
       real*8, dimension (:), ALLOCATABLE :: RMM,RMM1,RMM2,RMM3
       real*8, dimension (:), ALLOCATABLE :: rhoalpha,rhobeta
       real*8, dimension (:,:), ALLOCATABLE :: X, XX
       real*8, dimension (:), ALLOCATABLE :: old1,old2,old3 

c
      parameter(pi32=5.56832799683170698D0,pi=3.14159265358979312D0,
     >          rpi=1.77245385090551588D0, pi5=34.9868366552497108D0,
     >    pi52=34.9868366552497108D0)
      parameter(pis32=5.56832799683170698E0,piss=3.14159265358979312E0,
     >          rpis=1.77245385090551588E0, pis5=34.9868366552497108E0,
     >    pis52=34.9868366552497108E0)


c Angular momenta : up to f functions ( could be easily extended if
c necessary)
c
c      common /fit/ Nang,dens1,integ1,Iexch1,igridc,igrid2c
c      common /cav/ a01,epsilon1,field1,exter1,Fx1,Fy1,Fz1
c
c      common /index/ ii,iid
c
      Data Num /1,3,6,10/
      dimension jatom(2,100),coef(100),dist(100,3),distt(100)
      integer ndis,nsteps
      real*8 kjar,xini,xfinal   

      integer, dimension(:), ALLOCATABLE :: natomc,nnps,nnpp,nnpd,nns
      integer, dimension(:), ALLOCATABLE :: nnd,nnp
      real*8, dimension (:), ALLOCATABLE :: atmin
      integer, dimension(:,:), ALLOCATABLE :: jatc
      integer kknums,kknumd
      integer, dimension (:), ALLOCATABLE :: kkind,kkinds
      real*8     rmax, rmaxs
      real*8, dimension (:), ALLOCATABLE :: cool
      real*4, dimension (:), ALLOCATABLE :: cools
      parameter rmintsol=16.0D0

      end module garcha_mod
