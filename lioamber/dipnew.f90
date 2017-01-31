!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! Subroutine used for the calculation of the dipole moment in NEUTRAL          !
! (non-ionic) systems. Integrals are evaluated using the Obara Saika method.   !
! Inputs the density basis and outputs the dipole moment components.           !
! Original file: 19-1-1993                                                     !
!                                                                              !
! A loop is performed over all basis functions. Basis are supposed to be       !
! ordered according to type: first all s, then all p, then all d, etc.; and    !
! inside each type, they are ordered in shells: px, py, pz, dx2, dxy, dyy, dzx,!
! dzy, dzz, and so on.                                                         !
!                                                                              !
! ns, np, nd are markers for the end of s, p and d sections respectively.      !
! r(Nuc(i),j) is j component of position of nucleus i, j = 1..3.               !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine dipnew(uDip)

    use garcha_mod, only : RMM, NCO, Nunp, nuc, Iz, a, c, r, ncont, M, pc, d,  &
                           nshell, norm, natom, nsol, pi32, pi5, sol

    implicit none
    real*8, intent(inout) :: uDip(3)

    real*8 :: aux(3), aux1(3), aux2(3), aux3(3), aux4(3), aux5(3), aux6(3),    &
              srs(3), Q(3), uDipAt(3)
    real*8 :: sq3, alf, alf2, cc, ccoef, dd, dp, dijs, f1, f2, factor,  &
              z2, zij, Qc, ps, pis, pjs, ss, t0, t1
    integer :: M2, ns, np, nd, i, j, k, ii, jj, l1, l2, l3, l4, l12, l34, nElec, &
               n, ni, nj, icrd
     
    ns  = nshell(0)    ;
    np  = nshell(1)    ; 
    nd  = nshell(2)    ;
    M2  = 2*M   ; nElec = 2*NCO + Nunp ;
    uDip = 0.0D0
    sq3  = 1.D0  ; if (NORM) sq3 = sqrt( 3.D0 )
    
    ! First Loop: <s|s> case.
    do i=1, ns
    do j=1, i
        dd = d(Nuc(i), Nuc(j))
        do ni = 1, ncont(i)
        do nj = 1, ncont(j)
            zij   = a(i,ni) + a(j,nj)
            alf   = a(i,ni) * a(j,nj) / zij
            ss    = pi32 * exp(-alf*dd) / (zij*sqrt(zij))
            k     = i + ((M2-j) * (j-1)) /2
            ccoef = c(i,ni) * c(j,nj)
            cc    = ccoef * RMM(k)
            do icrd = 1, 3
                Q(icrd)    = (a(i,ni) * r(Nuc(i),icrd)     &
                           +  a(j,nj) * r(Nuc(j),icrd)) /zij
                srs(icrd)  = Q(icrd) * ss
                uDip(icrd) = uDip(icrd) + cc * srs(icrd)
            enddo
       enddo
       enddo
   enddo
   enddo

   ! Second Loop: <p|s> case.
   do i = ns + 1, ns+np, 3
   do j = 1, ns
       dd = d(Nuc(i), Nuc(j))
       do ni = 1, ncont(i)
       do nj = 1, ncont(j)
           zij  = a(i,ni) + a(j,nj)
           alf  = a(i,ni) * a(j,nj) /zij
           ss   = pi32 * exp(-alf*dd) / (zij*sqrt(zij))

           do icrd = 1, 3
               Q(icrd)   = (a(i,ni) *r(Nuc(i),icrd) &
                         +  a(j,nj) *r(Nuc(j),icrd)) /zij
               srs(icrd) = Q(icrd)*ss
           enddo

           ccoef = c(i,ni) * c(j,nj)
           z2    = 2.D0 * zij
           ! Different types of p in the p shell (x, y, z respectively)
           ! The l1-1 term takes into account the different components of 
           ! the shell.
           do l1 = 1, 3        
               k  = i + ((M2-j) * (j-1)) /2 + (l1-1)
               cc = ccoef * RMM(k)
               aux = (Q(l1) - r(Nuc(i),l1)) * srs
               aux(l1) = aux(l1) + ss /z2
               uDip = uDip + cc * aux
           enddo
       enddo
       enddo
   enddo
   enddo

   ! Third Loop: <p|p> case.
   do i = ns+1, ns+np, 3
   do j = ns+1, i, 3
       dd=d(Nuc(i), Nuc(j))
       do ni=1, ncont(i)
       do nj=1, ncont(j)
           zij = a(i,ni) + a(j,nj)
           z2  = 2.D0 * zij
           ss  = pi32 * exp(-alf*dd) /(zij*sqrt(zij))
           do icrd = 1, 3
               Q(icrd)   = (a(i, ni) *r(Nuc(i), icrd) &
                         +  a(j, nj) *r(Nuc(j), icrd)) /zij
               srs(icrd) = Q(icrd)*ss
           enddo

           t0  = a(i,ni) * a(j,nj)
           alf = t0 /zij
           ccoef=c(i,ni)*c(j,nj)

           do l1 = 1,3
               t1  = Q(l1) - r(Nuc(i), l1)
               aux = t1 * srs
  
               aux(l1) = aux(l1) + ss /z2
               ps = ss * t1
  
               do l2 = 1,3
                   aux1 = aux * ( Q(l2)-r(Nuc(j),l2) )
                    if (l1.eq.l2) aux1 = aux1 + srs /z2
                   
                   aux1(l2) = aux1(l2) + ps /z2
                   ii = i + l1 -1
                   jj = j + l2 -1
                 
                   if(ii.ge.jj) then
                       k    = ii + ((M2 -jj)*(jj -1))/2
                       cc   = RMM(k)*ccoef
                       uDip = uDip + cc*aux1
                   endif
               enddo
           enddo
       enddo
       enddo
   enddo
   enddo

   ! Fourth Loop: <d|s> case.
   do i = ns + np + 1, M, 6
   do j = 1, ns

       dd=d(Nuc(i),Nuc(j))
       do ni = 1, ncont(i)
       do nj = 1, ncont(j)
           zij=a(i,ni)+a(j,nj)
           z2=2.D0*zij
           Q(1)=(a(i,ni)*r(Nuc(i),1)+a(j,nj)*r(Nuc(j),1))/zij
           Q(2)=(a(i,ni)*r(Nuc(i),2)+a(j,nj)*r(Nuc(j),2))/zij
           Q(3)=(a(i,ni)*r(Nuc(i),3)+a(j,nj)*r(Nuc(j),3))/zij
           alf=a(i,ni)*a(j,nj)/zij
           alf2=2.D0*alf
           ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
           srs(1)=Q(1)*ss
           srs(2)=Q(2)*ss
           srs(3)=Q(3)*ss

           ccoef=c(i,ni)*c(j,nj)

           do l1 = 1, 3
               t1=Q(l1)-r(Nuc(i),l1)
               ps=ss*t1
               aux(1)=t1*srs(1)
               aux(2)=t1*srs(2)
               aux(3)=t1*srs(3)
               aux(l1)=aux(l1)+ss/z2

               do l2 = 1, l1
                   t1=Q(l2)-r(Nuc(i),l2)
                   aux1(1)=aux(1)*t1
                   aux1(2)=aux(2)*t1
                   aux1(3)=aux(3)*t1
                   f1=1.D0
               
                   if (l1.eq.l2) then
                       aux1(1)=aux1(1)+srs(1)/z2
                       aux1(2)=aux1(2)+srs(2)/z2
                       aux1(3)=aux1(3)+srs(3)/z2
                       f1=sq3
                   endif

                   aux1(l2)=aux1(l2)+ps/z2
                   l12=l1*(l1-1)/2+l2

                   ! The order of the d-shell should be as follows:
                   ! xx, yx, yy, zx, zy, zz (11, 21, 22, 31, 32, 33)

                   ii=i+l12-1
                   k=ii+((M2-j)*(j-1))/2
                   cc = RMM(k)*ccoef/f1

                   uDip(1)=uDip(1)+cc*aux1(1)
                   uDip(2)=uDip(2)+cc*aux1(2)
                   uDip(3)=uDip(3)+cc*aux1(3)
               enddo
           enddo
       enddo
       enddo
   enddo
   enddo

   ! Fifth Loop: <d|p> case.
   do i = ns + np + 1, M, 6
   do j = ns + 1, ns + np, 3

       dd=d(Nuc(i),Nuc(j))
       do ni=1,ncont(i)
       do nj=1,ncont(j)
           zij=a(i,ni)+a(j,nj)
           z2=2.D0*zij
           Q(1)=(a(i,ni)*r(Nuc(i),1)+a(j,nj)*r(Nuc(j),1))/zij
           Q(2)=(a(i,ni)*r(Nuc(i),2)+a(j,nj)*r(Nuc(j),2))/zij
           Q(3)=(a(i,ni)*r(Nuc(i),3)+a(j,nj)*r(Nuc(j),3))/zij
           alf=a(i,ni)*a(j,nj)/zij
           ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))
           srs(1)=Q(1)*ss
           srs(2)=Q(2)*ss
           srs(3)=Q(3)*ss

           ccoef=c(i,ni)*c(j,nj)

           do l1 = 1, 3
               t1=Q(l1)-r(Nuc(i),l1)
               pis=ss*t1

               ! aux : <pi|r|s>
               aux(1)=t1*srs(1)
               aux(2)=t1*srs(2)
               aux(3)=t1*srs(3)
               aux(l1)=aux(l1)+ss/z2

               do l2 = 1, l1
                   t1=Q(l2)-r(Nuc(i),l2)
                   pjs=ss*t1 
                   dijs=t1*pis

                   ! aux1 : <pj|r|s>
                   aux1(1)=t1*srs(1)
                   aux1(2)=t1*srs(2)
                   aux1(3)=t1*srs(3)
                   aux1(l2)=aux1(l2)+ss/z2
                
                   ! aux2 : <dij|r|s>
                   aux2(1)=t1*aux(1)
                   aux2(2)=t1*aux(2)
                   aux2(3)=t1*aux(3)
                   aux2(l2)=aux2(l2)+pis/z2
                   f1=1.D0

                   if (l1.eq.l2) then
                       f1=sq3
                       aux2(1)=aux2(1)+srs(1)/z2
                       aux2(2)=aux2(2)+srs(2)/z2
                       aux2(3)=aux2(3)+srs(3)/z2
                       dijs=dijs+ss/z2
                   endif

                   do l3 = 1, 3
                       t1=Q(l3)-r(Nuc(j),l3)
                       aux3(1)=aux2(1)*t1
                       aux3(2)=aux2(2)*t1
                       aux3(3)=aux2(3)*t1

                       if (l1.eq.l3) then
                            aux3(1)=aux3(1)+aux1(1)/z2
                            aux3(2)=aux3(2)+aux1(2)/z2
                            aux3(3)=aux3(3)+aux1(3)/z2
                       endif

                       if (l2.eq.l3) then
                            aux3(1)=aux3(1)+aux(1)/z2
                            aux3(2)=aux3(2)+aux(2)/z2
                            aux3(3)=aux3(3)+aux(3)/z2
                       endif
 
                       aux3(l3)=aux3(l3)+dijs/z2
 
                       l12=l1*(l1-1)/2+l2
                       ii=i+l12-1
                       jj=j+l3-1
           
                       k=ii+((M2-jj)*(jj-1))/2
                       cc=ccoef/f1*RMM(k)

                       uDip(1)=uDip(1)+cc*aux3(1)
                       uDip(2)=uDip(2)+cc*aux3(2)
                       uDip(3)=uDip(3)+cc*aux3(3)
                   enddo ! l1 loop
               enddo     ! l2 loop
           enddo         ! l3 loop
       enddo
       enddo
   enddo
   enddo

   ! Sixth Loop: <d|d> case.
   do i=ns+np+1,M,6
   do j=ns+np+1,M,6

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
           srs(1)=Q(1)*ss
           srs(2)=Q(2)*ss
           srs(3)=Q(3)*ss

           ccoef=c(i,ni)*c(j,nj)

           do l1=1,3
               t1=Q(l1)-r(Nuc(i),l1)
               pis=ss*t1
          
               ! Aux : <pi|r|s>
               aux(1)=t1*srs(1)
               aux(2)=t1*srs(2)
               aux(3)=t1*srs(3)
               aux(l1)=aux(l1)+ss/z2

               do l2=1,l1
                   t1=Q(l2)-r(Nuc(i),l2)
                   pjs=ss*t1
                   dijs=t1*pis
                
                   ! aux1 : <pj|r|s>
                   aux1(1)=t1*srs(1)
                   aux1(2)=t1*srs(2)
                   aux1(3)=t1*srs(3)
                   aux1(l2)=aux1(l2)+ss/z2
 
                   ! Aux2 : (dij|r|s)
                   aux2(1)=t1*aux(1)
                   aux2(2)=t1*aux(2)
                   aux2(3)=t1*aux(3)
                   aux2(l2)=aux2(l2)+pis/z2
                   f1=1.D0
 
                   if (l1.eq.l2) then
                        f1=sq3
                        aux2(1)=aux2(1)+srs(1)/z2
                        aux2(2)=aux2(2)+srs(2)/z2
                        aux2(3)=aux2(3)+srs(3)/z2
                        dijs=dijs+ss/z2
                   endif

                   do l3=1,3
                       t1=Q(l3)-r(Nuc(j),l3)
                       dp=t1*dijs
              
                       ! aux3 : (dij|r|pk)
                       aux3(1)=aux2(1)*t1
                       aux3(2)=aux2(2)*t1
                       aux3(3)=aux2(3)*t1

                       ! aux4 : (pi|r|pk)
                       aux4(1)=aux(1)*t1
                       aux4(2)=aux(2)*t1
                       aux4(3)=aux(3)*t1
 
                       ! aux5 : (pj|r|pk)
                       aux5(1)=aux1(1)*t1
                       aux5(2)=aux1(2)*t1
                       aux5(3)=aux1(3)*t1

                       if (l1.eq.l3) then
                           dp=dp+pjs/z2
                           aux3(1)=aux3(1)+aux1(1)/z2
                           aux3(2)=aux3(2)+aux1(2)/z2
                           aux3(3)=aux3(3)+aux1(3)/z2
                           aux4(1)=aux4(1)+srs(1)/z2
                           aux4(2)=aux4(2)+srs(2)/z2
                           aux4(3)=aux4(3)+srs(3)/z2
                       endif
                       if (l2.eq.l3) then
                           dp=dp+pis/z2
                           aux3(1)=aux3(1)+aux(1)/z2
                           aux3(2)=aux3(2)+aux(2)/z2
                           aux3(3)=aux3(3)+aux(3)/z2
                           aux5(1)=aux5(1)+srs(1)/z2
                           aux5(2)=aux5(2)+srs(2)/z2
                           aux5(3)=aux5(3)+srs(3)/z2
                       endif

                       aux3(l3)=aux3(l3)+dijs/z2
                       aux4(l3)=aux4(l3)+pis/z2
                       aux5(l3)=aux5(l3)+pjs/z2

                       do l4=1,l3
                           t1=Q(l4)-r(Nuc(j),l4)
 
                           ! aux3 : used here for (d|r|d)
                           aux6(1)=aux3(1)*t1
                           aux6(2)=aux3(2)*t1
                           aux6(3)=aux3(3)*t1

                           if (l1.eq.l4) then
                               aux6(1)=aux6(1)+aux5(1)/z2
                               aux6(2)=aux6(2)+aux5(2)/z2
                               aux6(3)=aux6(3)+aux5(3)/z2
                           endif

                           if (l2.eq.l4) then
                               aux6(1)=aux6(1)+aux4(1)/z2
                               aux6(2)=aux6(2)+aux4(2)/z2
                               aux6(3)=aux6(3)+aux4(3)/z2
                           endif

                           f2=1.D0
                           if (l3.eq.l4) then
                               f2=sq3
                               aux6(1)=aux6(1)+aux2(1)/z2
                               aux6(2)=aux6(2)+aux2(2)/z2
                               aux6(3)=aux6(3)+aux2(3)/z2
                           endif

                           aux6(l4)=aux6(l4)+dp/z2
                           l12=l1*(l1-1)/2+l2
                           l34=l3*(l3-1)/2+l4
                           ii=i+l12-1
                           jj=j+l34-1

                           if (ii.ge.jj) then
                               k=ii+((M2-jj)*(jj-1))/2
                               cc=ccoef/(f1*f2)*RMM(k)
                               uDip(1)=uDip(1)+cc*aux6(1)
                               uDip(2)=uDip(2)+cc*aux6(2)
                               uDip(3)=uDip(3)+cc*aux6(3)
                           endif
                       enddo
                   enddo
               enddo
           enddo
       enddo
       enddo
   enddo
   enddo

   uDipAt = 0.0D0
   Qc     = 0.0D0
 
   do i=1,natom
       Qc  = Qc  + Iz(i)
       uDipAt(1) = uDipAt(1) + Iz(i)*r(i,1)
       uDipAt(2) = uDipAt(2) + Iz(i)*r(i,2)
       uDipAt(3) = uDipAt(3) + Iz(i)*r(i,3)
   enddo

   Qc=Qc-nElec
   if (sol) then
       do k=1,Nsol
           Qc = Qc + pc(k)
       enddo
   endif

! Factor : For charged species dipole moment depends on the definition of the  !
! coordinates' origin. Using this factor, it is defined with respect to the    !
! center of charge (important in Reaction Field calculations). For neutral     !
! systems this is not necessary.                                               !

   factor = (Qc + nElec)/nElec
   uDip(1) = (uDipAt(1) - uDip(1)*factor) * 2.54D0
   uDip(2) = (uDipAt(2) - uDip(2)*factor) * 2.54D0
   uDip(3) = (uDipAt(3) - uDip(3)*factor) * 2.54D0
 
   return
end subroutine dipnew
