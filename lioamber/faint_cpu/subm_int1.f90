!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
      module subm_int1; contains
       subroutine int1(En, RMM, Smat, Nuc, a, c, d, r, Iz, ncont, NORM,&
                       & natom, M, Md )
!------------------------------------------------------------------------------!
!
!      Integrals subroutine
!      1 e integrals
!      using the Obara-Saika recursive method.
!      (See Helgaker, "Molecular Electronic Structure Theory" (2000). pg 339)
!
!      loop over all basis functions
!      now the basis is supposed to be ordered according to the type,
!      all s, then all p, then all d, .....
!      inside each type, are ordered in shells
!      px,py,pz , dx2,dxy,dyy,dzx,dzy,dzz, .....
!
!      ns ... marker for end of s
!      np ... marker for end of p
!      nd ... marker for end of d
!
!      r(Nuc(i),j) j component of position of nucleus i , j=1,3
!      Input : basis function information
!      Output: F matrix, and S matrix
!
!------------------------------------------------------------------------------!
       use liotemp   , only: FUNCT
       use garcha_mod, only: nshell
       use constants_mod, only: pi, pi32
       implicit none

!      Input quantities (ex-garchamod variables)
        double precision, allocatable, intent(inout) :: RMM(:)
        double precision, allocatable, intent(inout) :: Smat(:,:)
        double precision,              intent(inout) :: En
        integer,                       intent(inout) :: natom

        double precision, allocatable, intent(in) :: d(:,:)
        double precision, intent(in) :: r(natom,3)
        double precision, allocatable, intent(in) :: a(:,:)
        double precision, allocatable, intent(in) :: c(:,:)
        integer, intent(in) :: Nuc(M)
        integer, intent(in) :: Iz(natom)
        integer, allocatable, intent(in) :: ncont(:)
        integer, intent(in) :: M
        integer, intent(in) :: Md
        logical, intent(in) :: NORM


       integer :: natomold, igpu
       integer :: n, i, j, k, ii, jj, ni, nj, iin
       integer :: l1, l2, l3, l4, l12, l34
       integer :: MM, MMd, ns, np, nd
       integer :: M1, M2, M3, M5, M7, M9, M11

       double precision  :: E1, ovlap
       double precision  :: Q(3), term, temp, sq3, alf, alf2, cc, ccoef
       double precision  :: f1, f2, tn, tna, u, z2, zij
       double precision  :: ss, ps, dd, p0s, p1s, p2s, p3s
       double precision  :: pi0p, pi1p, piks, pikpk, pipk, pis
       double precision  :: pj0s, pj1s, pj2s, pj0p, pj1p, pjkpk
       double precision  :: pjks, pjpk, pjs, pks, sks
       double precision  :: dijs, dijpk, dijks, dijkpk
       double precision  :: d0s, d0p, d1p, d1s, d2s
       double precision  :: t0, t1, t2

       integer, allocatable, dimension(:) :: Iaux
       double precision , allocatable, dimension(:) :: s0s, s1s, s2s, s3s, s4s
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
       allocate(s0s(natom),s1s(natom),s2s(natom))
!       allocate(s3s(natom),s4s(natom),Iaux(natom))
       allocate(s3s(natom))
       allocate(s4s(natom))
!       allocate(Iaux(natom))
       if (.not.allocated(Smat)) allocate(Smat(M,M))

!-------------------------------------------------------------------------------
!      Distance between pairs of centers
!      BSSE
!-------------------------------------------------------------------------------
!      Sets to 0 atomic charges, but since they are used to
!      construct the grids, they are stored in an auxiliary array
!
!      if (BSSE) then
!      do i=1,natom
!       Iaux(i)=Iz(i)
!       Iz(i)=Iz(i)*ighost(i)
!      enddo
!      endif
!-------------------------------------------------------------------------------

      if (NORM) then
        sq3=sqrt(3.D0)
      else
        sq3=1.D0
      endif

      ns=nshell(0)
      np=nshell(1)
      nd=nshell(2)
      MM=M*(M+1)/2
      MMd=Md*(Md+1)/2
      M2=2*M

! Pointers
!---------
! P
      M1=1
! Pnew
      M3=M1+MM
! S, F also uses the same position after S was used
      M5=M3+MM
! G
      M7=M5+MM
! Gm
      M9=M7+MMd
! H
      M11=M9+MMd

!-------------------------------------------------------------------------------
! Overlap ,Kinetic energy and Nuclear Attraction matrix elements evaluation
!-------------------------------------------------------------------------------
! Overlap matrix will be kept, kinetic energy and nuclear attraction matrix
! elements wont, they're stored in Fock matrix and in the Energy directly
! in order to reduce the memory requirements
!-------------------------------------------------------------------------------

       Smat=0.0d0
       do i=1,MM
         RMM(M5+i-1)=0.D0
         RMM(M11+i-1)=0.D0
       enddo
!
!     do 50 i=1,natom
!       do 50 j=1,natom
!         d(i,j)=(r(i,1)-r(j,1))**2+(r(i,2)-r(j,2))**2+(r(i,3)-r(j,3))**2
!       end do
!     end do
!
! Nuclear Repulsion part ------------------------------------------
      En=0.D0
      do i=1,natom
      do j=1,i-1
           En=En+Iz(i)*Iz(j)/sqrt(d(i,j))
      enddo
      enddo

      call aint_query_gpu_level(igpu)
      ! Doing nuclear attraction part on GPU - KE part still is done here.
      if (igpu.gt.3) then
        natomold = natom
        natom = 0
      endif

! First loop (s|s) case -------------------------------------------

      do i=1,ns                                                                 ! Over contractions
      do j=1,i                                                                ! Over contractions
        dd=d(Nuc(i),Nuc(j))                                                     ! Over primitives
        do ni=1,ncont(i)                                                      ! Over primitives
        do nj=1,ncont(j)

!          (0|0) calculation
           zij=a(i,ni)+a(j,nj)
           alf=a(i,ni)*a(j,nj)/zij
           alf2=alf*2.D0
           Q(1)=(a(i,ni)*r(Nuc(i),1)+a(j,nj)*r(Nuc(j),1))/zij
           Q(2)=(a(i,ni)*r(Nuc(i),2)+a(j,nj)*r(Nuc(j),2))/zij
           Q(3)=(a(i,ni)*r(Nuc(i),3)+a(j,nj)*r(Nuc(j),3))/zij
           ccoef=c(i,ni)*c(j,nj)

           ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
           ovlap=ss
           tn=alf*(3.D0-alf2*dd)*ovlap
!
           k=i+((M2-j)*(j-1))/2
           RMM(M5+k-1)=RMM(M5+k-1)+ccoef*ovlap
           Smat(i,j)=Smat(i,j)+ccoef*ovlap
           Smat(j,i)=Smat(j,i)+ccoef*ovlap

!          loop over nuclei, nuclear attraction matrix elements
!          tna: accumulates nuc. attraction over all nuclei

           tna=0.D0

           do n=1,natom
              u=(Q(1)-r(n,1))**2+(Q(2)-r(n,2))**2+(Q(3)-r(n,3))**2
              u=u*zij

              s0s(n)=Iz(n)*2.*sqrt(zij/pi)*ss*FUNCT(0,u)
              tna=tna-s0s(n)
           enddo ! n
!
           term=ccoef*(tn+tna)
           RMM(M11+k-1)=RMM(M11+k-1)+ term
         end do  ! nj
         end do  ! ni
      end do  ! j
      end do  ! i

!
!-------------------------------------------------------------------------------
! second loop  (p|s) case
!-------------------------------------------------------------------------------

     do i=ns+1,ns+np,3                                                         ! Over contractions
     do j=1,ns                                                               ! Over contractions
        dd=d(Nuc(i),Nuc(j))
        do ni=1,ncont(i)                                                      ! Over primitives
        do nj=1,ncont(j)                                                    ! Over primitives

           zij=a(i,ni)+a(j,nj)
           Q(1)=(a(i,ni)*r(Nuc(i),1)+a(j,nj)*r(Nuc(j),1))/zij
           Q(2)=(a(i,ni)*r(Nuc(i),2)+a(j,nj)*r(Nuc(j),2))/zij
           Q(3)=(a(i,ni)*r(Nuc(i),3)+a(j,nj)*r(Nuc(j),3))/zij
           alf=a(i,ni)*a(j,nj)/zij
           alf2=2.D0*alf
           ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
           sks=alf*(3.D0-alf2*dd)*ss

  !        Loop over nuclei, part common for all shell
           do n=1,natom
             u=(Q(1)-r(n,1))**2+(Q(2)-r(n,2))**2+(Q(3)-r(n,3))**2
             u=u*zij
             temp=2.D0*sqrt(zij/pi)*ss

             s0s(n)=temp*FUNCT(0,u)
             s1s(n)=temp*FUNCT(1,u)
           end do ! n
!
           ccoef=c(i,ni)*c(j,nj)
!
!          l2: different p in the p shell ( x,y,z respectively)
!
           do l2=1,3
              t1=Q(l2)-r(Nuc(i),l2)
              ovlap=t1*ss
              tn=t1*sks+alf2*ovlap
              iin=i+l2-1
!            ii index , taking into account different components of the shell
!
              k=iin+((M2-j)*(j-1))/2
              RMM(M5+k-1)=RMM(M5+k-1)+ovlap*ccoef
              Smat(iin,j)=Smat(iin,j)+ovlap*ccoef
              Smat(j,iin)=Smat(j,iin)+ovlap*ccoef
!             loop over nuclei, specific part

              tna=0.D0
              do n=1,natom
                term=t1*s0s(n)-(Q(l2)-r(n,l2))*s1s(n)
                tna=tna-Iz(n)*term
              end do ! n

              term=ccoef*(tn+tna)
              RMM(M11+k-1)=RMM(M11+k-1)+term
           end do ! l2
         end do ! nj
         end do  ! ni
      end do ! j
      end do ! i

!-------------------------------------------------------------------
!
! (p|p) case
!
      do i=ns+1,ns+np,3
      do j=ns+1,i,3
         dd=d(Nuc(i),Nuc(j))
         do ni=1,ncont(i)
         do nj=1,ncont(j)

            zij=a(i,ni)+a(j,nj)
            z2=2.D0*zij
            Q(1)=(a(i,ni)*r(Nuc(i),1)+a(j,nj)*r(Nuc(j),1))/zij
           Q(2)=(a(i,ni)*r(Nuc(i),2)+a(j,nj)*r(Nuc(j),2))/zij
           Q(3)=(a(i,ni)*r(Nuc(i),3)+a(j,nj)*r(Nuc(j),3))/zij
           alf=a(i,ni)*a(j,nj)/zij
           alf2=2.D0*alf
           ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
           sks=alf*(3.D0-alf2*dd)*ss

!          loop over nuclei, part common for all shell
           do n=1,natom
             u=(Q(1)-r(n,1))**2+(Q(2)-r(n,2))**2+(Q(3)-r(n,3))**2
             u=u*zij
             temp=2.D0*sqrt(zij/pi)*ss

             s0s(n)=temp*FUNCT(0,u)
             s1s(n)=temp*FUNCT(1,u)
             s2s(n)=temp*FUNCT(2,u)
           end do

           ccoef=c(i,ni)*c(j,nj)

           do l1=1,3
             t1=Q(l1)-r(Nuc(i),l1)
             ps=ss*t1
             pks=sks*t1+alf2*ps

             do l2=1,3

               t1=Q(l2)-r(Nuc(j),l2)
               ovlap=t1*ps
               tn=t1*pks+alf2*ovlap

               if (l1.eq.l2) then
                 ovlap=ovlap+ss/z2
                 tn=tn+(sks+alf2*ss)/z2
               endif

               iin=i+l1-1
               jj=j+l2-1
               term=tn*ccoef

               if (iin.ge.jj) then
                 k=iin+((M2-jj)*(jj-1))/2
                 RMM(M5+k-1)=RMM(M5+k-1)+ovlap*ccoef
                 Smat(iin,jj)=Smat(iin,jj)+ovlap*ccoef
                 Smat(jj,iin)=Smat(jj,iin)+ovlap*ccoef
                 RMM(M11+k-1)=RMM(M11+k-1)+term
               endif

              end do ! l2
            end do ! l1


!           loop over nuclei ( specific part)
            do n=1,natom
            do l1=1,3
               t1=Q(l1)-r(Nuc(i),l1)
               t2=Q(l1)-r(n,l1)
               p0s=t1*s0s(n)-t2*s1s(n)
               p1s=t1*s1s(n)-t2*s2s(n)

               do l2=1,3
                 t1=Q(l2)-r(Nuc(j),l2)
                 t2=Q(l2)-r(n,l2)
                 tna=t1*p0s-t2*p1s

                 if (l1.eq.l2) then
                   tna=tna+(s0s(n)-s1s(n))/z2
                 endif

                 iin=i+l1-1
                 jj=j+l2-1

                 if (iin.ge.jj) then
                   k=iin+((M2-jj)*(jj-1))/2
                   term=-tna*ccoef*Iz(n)
                   RMM(M11+k-1)=RMM(M11+k-1)+term
                 endif

               end do ! l2
             end do ! l1
             end do ! n
!-----------------------
          end do ! nj
          end do ! ni
      end do ! j
      end do ! i
!-------------------------------------------------------------------
! (d|s) case

      do i=ns+np+1,M,6
      do j=1,ns
         dd=d(Nuc(i),Nuc(j))
         do ni=1,ncont(i)
         do nj=1,ncont(j)

           zij=a(i,ni)+a(j,nj)
           z2=2.D0*zij
           Q(1)=(a(i,ni)*r(Nuc(i),1)+a(j,nj)*r(Nuc(j),1))/zij
           Q(2)=(a(i,ni)*r(Nuc(i),2)+a(j,nj)*r(Nuc(j),2))/zij
           Q(3)=(a(i,ni)*r(Nuc(i),3)+a(j,nj)*r(Nuc(j),3))/zij
           alf=a(i,ni)*a(j,nj)/zij
           alf2=2.D0*alf
           ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
           sks=alf*(3.D0-alf2*dd)*ss

!          Loop over nuclei, part common for all shell
           do n=1,natom
             u=(Q(1)-r(n,1))**2+(Q(2)-r(n,2))**2+(Q(3)-r(n,3))**2
             u=u*zij
             temp=2.D0*sqrt(zij/pi)*ss
             s0s(n)=temp*FUNCT(0,u)
             s1s(n)=temp*FUNCT(1,u)
             s2s(n)=temp*FUNCT(2,u)
           end do

           ccoef=c(i,ni)*c(j,nj)

           do l1=1,3
             t1=Q(l1)-r(Nuc(i),l1)
             ps=ss*t1
             pks=sks*t1+alf2*ps

             do l2=1,l1

               t1=Q(l2)-r(Nuc(i),l2)
               ovlap=t1*ps
               tn=t1*pks

               f1=1.D0
               if (l1.eq.l2) then
                 ovlap=ovlap+ss/z2
                 tn=tn+sks/z2-alf2*ss/(2.*a(i,ni))
                 f1=sq3
               endif

               tn=tn+alf2*ovlap

               l12=l1*(l1-1)/2+l2
!                Ordering of d shell should be:
!                xx,yx,yy,zx,zy,zz ( 11, 21, 22, 31, 32, 33 )
               iin=i+l12-1

               cc=ccoef/f1
               term=cc*tn
               k=iin+((M2-j)*(j-1))/2
               RMM(M5+k-1)=RMM(M5+k-1)+ovlap*cc
               Smat(iin,j)=Smat(iin,j)+ovlap*cc
               Smat(j,iin)=Smat(j,iin)+ovlap*cc
               RMM(M11+k-1)=RMM(M11+k-1)+term
           end do
           end do

!          Nuclear attraction part
           do n=1,natom
           do l1=1,3
              t1=Q(l1)-r(Nuc(i),l1)
              t2=Q(l1)-R(n,l1)
              p0s=t1*s0s(n)-t2*s1s(n)
              p1s=t1*s1s(n)-t2*s2s(n)
              do l2=1,l1
                 t1=Q(l2)-r(Nuc(i),l2)
                 t2=Q(l2)-r(n,l2)
                 tna=t1*p0s-t2*p1s

                 f1=1.D0
                 if (l1.eq.l2) then
                   tna=tna+(s0s(n)-s1s(n))/z2
                   f1=sq3
                 endif

                 l12=l1*(l1-1)/2+l2
!                Ordering of d shell should be:
!                xx,yx,yy,zx,zy,zz ( 11, 21, 22, 31, 32, 33 )

                 iin=i+l12-1

                 k=iin+((M2-j)*(j-1))/2
                 cc=ccoef/f1
                 term=-cc*tna*Iz(n)
                 RMM(M11+k-1)=RMM(M11+k-1)+term
               end do ! l2
            end do  ! l1
            end do ! n
!           End nuclear attr. part
         end do ! j
         end do ! i
      end do ! nj
      end do ! ni
!-----------------------------------------------------------------
!
! (d|p) case
!
     do i=ns+np+1,M,6
     do j=ns+1,ns+np,3
        dd=d(Nuc(i),Nuc(j))
        do ni=1,ncont(i)
        do nj=1,ncont(j)

           zij=a(i,ni)+a(j,nj)
           z2=2.D0*zij
           Q(1)=(a(i,ni)*r(Nuc(i),1)+a(j,nj)*r(Nuc(j),1))/zij
           Q(2)=(a(i,ni)*r(Nuc(i),2)+a(j,nj)*r(Nuc(j),2))/zij
           Q(3)=(a(i,ni)*r(Nuc(i),3)+a(j,nj)*r(Nuc(j),3))/zij
           alf=a(i,ni)*a(j,nj)/zij
           alf2=2.D0*alf
           ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
           sks=alf*(3.D0-alf2*dd)*ss

!          Loop over nuclei, part common for all shell
           do n=1,natom
             u=(Q(1)-r(n,1))**2+(Q(2)-r(n,2))**2+(Q(3)-r(n,3))**2
             u=u*zij
             temp=2.D0*sqrt(zij/pi)*ss
             s0s(n)=temp*FUNCT(0,u)
             s1s(n)=temp*FUNCT(1,u)
             s2s(n)=temp*FUNCT(2,u)
             s3s(n)=temp*FUNCT(3,u)
           end do

           ccoef=c(i,ni)*c(j,nj)

           do l1=1,3
             t1=Q(l1)-r(Nuc(i),l1)
             pis=ss*t1
             piks=sks*t1+alf2*pis
             do l2=1,l1
               t1=Q(l2)-r(Nuc(i),l2)
               pjs=ss*t1
               pjks=sks*t1+alf2*pjs

               dijs=t1*pis
               dijks=t1*piks
               f1=1.D0

               if (l1.eq.l2) then
                 f1=sq3
                 dijs=dijs+ss/z2
                 dijks=dijks+sks/z2-alf2*ss/(2.D0*a(i,ni))
               endif

               dijks=dijks+alf2*dijs

               do l3=1,3

                 t1=Q(l3)-r(Nuc(j),l3)
                 ovlap=t1*dijs
                 tn=t1*dijks

                 if (l1.eq.l3) then
                   ovlap=ovlap+pjs/z2
                   tn=tn+pjks/z2
                 endif

                 if (l2.eq.l3) then
                   ovlap=ovlap+pis/z2
                   tn=tn+piks/z2
                 endif

                 tn=tn+alf2*ovlap

                 l12=l1*(l1-1)/2+l2
                 iin=i+l12-1
                 jj=j+l3-1

                 cc=ccoef/f1
                 term=cc*tn
                 k=iin+((M2-jj)*(jj-1))/2
                 RMM(M5+k-1)=RMM(M5+k-1)+cc*ovlap
                 Smat(iin,jj)=Smat(iin,jj)+cc*ovlap
                 Smat(jj,iin)=Smat(jj,iin)+cc*ovlap
                 RMM(M11+k-1)=RMM(M11+k-1)+term
             end do
          end do
          end do

!         -----------------------
!         Nuclear attraction part
           do n=1,natom
           do l1=1,3
              t1=Q(l1)-r(Nuc(i),l1)
              t2=Q(l1)-r(n,l1)
              p0s=t1*s0s(n)-t2*s1s(n)
              p1s=t1*s1s(n)-t2*s2s(n)
              p2s=t1*s2s(n)-t2*s3s(n)

              do l2=1,l1
                 t1=Q(l2)-r(Nuc(i),l2)
                 t2=Q(l2)-r(n,l2)
                 pj0s=t1*s0s(n)-t2*s1s(n)
                 pj1s=t1*s1s(n)-t2*s2s(n)

                 f1=1.D0
                 d0s=t1*p0s-t2*p1s
                 d1s=t1*p1s-t2*p2s

                 if (l1.eq.l2) then
                   f1=sq3
                   d0s=d0s+(s0s(n)-s1s(n))/z2
                   d1s=d1s+(s1s(n)-s2s(n))/z2
                 endif

                 do l3=1,3
                   t1=Q(l3)-r(Nuc(j),l3)
                   t2=Q(l3)-r(n,l3)
                   tna=t1*d0s-t2*d1s

                   if (l1.eq.l3) then
                     tna=tna+(pj0s-pj1s)/z2
                   endif

                   if (l2.eq.l3) then
                     tna=tna+(p0s-p1s)/z2
                   endif

                   l12=l1*(l1-1)/2+l2
                   iin=i+l12-1
                   jj=j+l3-1

                   k=iin+((M2-jj)*(jj-1))/2
                   cc=ccoef/f1
                   term=-cc*tna*Iz(n)
                   RMM(M11+k-1)=RMM(M11+k-1)+term
                 end do ! l3
               end do ! l2
               end do ! l1
              end do ! n
!             End of nuclear atraction part.
!             -----------------------
           end do ! j
           end do ! i
       end do ! nj
       end do ! ni
!
!-------------------------------------------------------------------
!
! (d|d) case
!
    do i=ns+np+1,M,6
    do j=ns+np+1,i,6
       dd=d(Nuc(i),Nuc(j))
       do ni=1,ncont(i)
       do nj=1,ncont(j)

       zij=a(i,ni)+a(j,nj)
       z2=2.D0*zij
       Q(1)=(a(i,ni)*r(Nuc(i),1)+a(j,nj)*r(Nuc(j),1))/zij
       Q(2)=(a(i,ni)*r(Nuc(i),2)+a(j,nj)*r(Nuc(j),2))/zij
       Q(3)=(a(i,ni)*r(Nuc(i),3)+a(j,nj)*r(Nuc(j),3))/zij
       alf=a(i,ni)*a(j,nj)/zij
       alf2=2.D0*alf
       ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
       sks=alf*(3.D0-alf2*dd)*ss

!      Loop over nuclei, part common for all shell
       do n=1,natom
          u=(Q(1)-r(n,1))**2+(Q(2)-r(n,2))**2+(Q(3)-r(n,3))**2
          u=u*zij
          temp=2.D0*sqrt(zij/pi)*ss
          s0s(n)=temp*FUNCT(0,u)
          s1s(n)=temp*FUNCT(1,u)
          s2s(n)=temp*FUNCT(2,u)
          s3s(n)=temp*FUNCT(3,u)
          s4s(n)=temp*FUNCT(4,u)
       end do ! n

       ccoef=c(i,ni)*c(j,nj)
       t0=ss/z2

       do l1=1,3

       t1=Q(l1)-r(Nuc(i),l1)
       pis=ss*t1
       piks=sks*t1+alf2*pis
       do l2=1,l1
         t1=Q(l2)-r(Nuc(i),l2)
         pjs=ss*t1
         pjks=sks*t1+alf2*pjs

         f1=1.D0
         t1=Q(l2)-r(Nuc(i),l2)
         dijs=t1*pis
         dijks=t1*piks

         if (l1.eq.l2) then
           f1=sq3
           dijs=dijs+t0
           dijks=dijks+sks/z2-alf2*ss/(2.D0*a(i,ni))
         endif

         dijks=dijks+alf2*dijs

         do l3=1,3

           t2=Q(l3)-r(Nuc(j),l3)
           pipk=t2*pis
           pjpk=t2*pjs
           pikpk=t2*piks
           pjkpk=t2*pjks
           dijpk=t2*dijs
           dijkpk=t2*dijks

           if (l1.eq.l3) then
             pipk=pipk+t0
             dijpk=dijpk+pjs/z2
             pikpk=pikpk+sks/z2
             dijkpk=dijkpk+pjks/z2
           endif

           if (l2.eq.l3) then
             pjpk=pjpk+t0
             dijpk=dijpk+pis/z2
             pjkpk=pjkpk+sks/z2
             dijkpk=dijkpk+piks/z2
           endif

           pikpk=pikpk+alf2*pipk
           pjkpk=pjkpk+alf2*pjpk
           dijkpk=dijkpk+alf2*dijpk

           do l4=1,l3

             f2=1.D0
             t1=Q(l4)-r(Nuc(j),l4)
             ovlap=t1*dijpk
             tn=t1*dijkpk

             if (l1.eq.l4) then
               ovlap=ovlap+pjpk/z2
               tn=tn+pjkpk/z2
             endif

             if (l2.eq.l4) then
               ovlap=ovlap+pipk/z2
               tn=tn+pikpk/z2
             endif

             if (l3.eq.l4) then
               ovlap=ovlap+dijs/z2
               tn=tn+dijks/z2-alf2*dijs/(2.D0*a(j,nj))
               f2=sq3
             endif

!            l12 and l34 go from 1 to 6, spanning the d shell in
!            the order xx, xy, yy, zx, zy, zz.
!            The same order should be used in ordering the basis set,
!            before calculating the matrix elements.

             l12=l1*(l1-1)/2+l2
             l34=l3*(l3-1)/2+l4
             iin=i+l12-1
             jj=j+l34-1

             if (iin.ge.jj) then
               k=iin+((M2-jj)*(jj-1))/2
               tn= tn+alf2*ovlap

               cc=ccoef/(f1*f2)
               term=cc*tn
               RMM(M5+k-1)=RMM(M5+k-1)+ovlap*cc
               Smat(iin,jj)=Smat(iin,jj)+ovlap*cc
               Smat(jj,iin)=Smat(jj,iin)+ovlap*cc
               RMM(M11+k-1)=RMM(M11+k-1)+term
             endif
           end do ! l4
           end do ! l3
        end do ! l2
        end do ! l1

!       -----------------------------------------------
!       Loop over nuclei - Nuclear attraction part.
       do n=1,natom
       do l1=1,3
          t1=Q(l1)-r(Nuc(i),l1)
          t2=Q(l1)-r(n,l1)
          p0s=t1*s0s(n)-t2*s1s(n)
          p1s=t1*s1s(n)-t2*s2s(n)
          p2s=t1*s2s(n)-t2*s3s(n)
          p3s=t1*s3s(n)-t2*s4s(n)

          do l2=1,l1
             f1=1.D0
             t1=Q(l2)-r(Nuc(i),l2)
             t2=Q(l2)-r(n,l2)
             pj0s=t1*s0s(n)-t2*s1s(n)
             pj1s=t1*s1s(n)-t2*s2s(n)
             pj2s=t1*s2s(n)-t2*s3s(n)

             d0s=t1*p0s-t2*p1s
             d1s=t1*p1s-t2*p2s
             d2s=t1*p2s-t2*p3s
             if (l1.eq.l2) then
                f1=sq3
                d0s=d0s+(s0s(n)-s1s(n))/z2
                d1s=d1s+(s1s(n)-s2s(n))/z2
                d2s=d2s+(s2s(n)-s3s(n))/z2
             endif

             do l3=1,3
                t1=Q(l3)-r(Nuc(j),l3)
                t2=Q(l3)-r(n,l3)

                d0p=t1*d0s-t2*d1s
                d1p=t1*d1s-t2*d2s

                pi0p=t1*p0s-t2*p1s
                pi1p=t1*p1s-t2*p2s
                pj0p=t1*pj0s-t2*pj1s
                pj1p=t1*pj1s-t2*pj2s

                if (l1.eq.l3) then
                  d0p=d0p+(pj0s-pj1s)/z2
                  d1p=d1p+(pj1s-pj2s)/z2
                  pi0p=pi0p+(s0s(n)-s1s(n))/z2
                  pi1p=pi1p+(s1s(n)-s2s(n))/z2
                endif

                if (l2.eq.l3) then
                  d0p=d0p+(p0s-p1s)/z2
                  d1p=d1p+(p1s-p2s)/z2
                  pj0p=pj0p+(s0s(n)-s1s(n))/z2
                  pj1p=pj1p+(s1s(n)-s2s(n))/z2
                endif


                do l4=1,l3
                   f2=1.D0
                   t1=Q(l4)-R(Nuc(j),l4)
                   t2=Q(l4)-R(n,l4)
                   tna=t1*d0p-t2*d1p
                   if (l4.eq.l1) then
                    tna=tna+(pj0p-pj1p)/z2
                   endif

                  if (l4.eq.l2) then
                    tna=tna+(pi0p-pi1p)/z2
                  endif

                  if (l4.eq.l3) then
                    f2=sq3
                    tna=tna+(d0s-d1s)/z2
                  endif

                  cc=ccoef/(f1*f2)
                  term=-cc*Iz(n)*tna

                  l12=l1*(l1-1)/2+l2
                  l34=l3*(l3-1)/2+l4
                  iin=i+l12-1
                  jj=j+l34-1

                  if (iin.ge.jj) then
                    k=iin+((M2-jj)*(jj-1))/2
                    RMM(M11+k-1)=RMM(M11+k-1)+term
                  endif
                end do ! l4
                end do ! l3
                end do ! l2
                end do ! l1
             end do ! n
!            End nuclear attraction part
!            ---------------------------
         end do ! j
         end do ! i
      end do ! nj
      end do ! ni

!*******************************************************************************
!--- Debugging and tests -----------------------------------------------
!
!     do 1 k=1,M*(M+1)/2
!      write(*,*) k,RMM(M5+k-1),RMM(M11+k-1)
! gradient debuggings
!
!*******************************************************************************

      E1=0.D0
      do k=1,MM
        E1=E1+RMM(k)*RMM(M11+k-1)
      end do

!*******************************************************************************
!
!      write(*,*) 'E1+En =',E1+En
!
! BSSE ------------
!      if (BSSE) then
!        do i=1,natom
!         Iz(i)=Iaux(i)
!        enddo
!      endif
!
!-- prueba ----------------
!      En=En+0.0D0*(d(1,2)-2.89D0)**2
!--------------------------
!
!     write(*,*) 'matriz overlap'
!     do i=1,MM
!      write(*,*) i,RMM(M5+i-1)
!     enddo
!     do i=1,natom
!      write(*,*) i,r(i,1),r(i,2),r(i,3)
!     enddo
!     pause
!*******************************************************************************

!     Avoid double-counting diagonal elements.
      do i=1,M
        Smat(i,i)=Smat(i,i)/2
      enddo

      deallocate(s0s,s1s,s2s,s3s,s4s)
!      deallocate(Iaux)

      if (igpu.gt.3) natom = natomold
      return;end subroutine
      end module subm_int1
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
