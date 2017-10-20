!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!%% DIP.F90 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
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
subroutine dip(uDip)

    use garcha_mod, only : RMM, NCO, Nunp, nuc, Iz, a, c, r, ncont, M, pc, d,  &
                           nshell, norm, natom, nsol, pi32, pi5, sol

    implicit none
    real*8, intent(inout) :: uDip(3)

    real*8  :: P_density(M*(M+1)/2)  ! Until we get rid of RMM.
    real*8  :: aux(3), aux1(3), aux2(3), aux3(3), aux4(3), aux5(3), aux6(3),   &
               srs(3), Q(3), uDipAt(3)
    real*8  :: sq3, alf, alf2, cc, cCoef, dd, dp, dijs, f1, f2, factor, z2,    &
               zij, Qc, ps, pis, pjs, ss, t0, t1
    integer :: M2, ns, np, nd, i, j, k, ii, jj, l1, l2, l3, l4, l12, l34, n,   &
               ni, nj, iCrd, nElec
     
    ! Constants
    ns  = nshell(0)   
    np  = nshell(1)    
    nd  = nshell(2)    
    M2  = 2*M    
    uDip = 0.0D0
    sq3  = 1.D0  
    if (NORM) sq3 = sqrt( 3.D0 )
    nElec = 2*NCO + Nunp
    
    ! This is until we get rid of RMM.
    k = M*(M+1) /2
    do i=1, k
       P_density(i) = RMM(i)
    enddo 

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
            cCoef = c(i,ni) * c(j,nj)
            cc    = cCoef * P_density(k)
            do iCrd = 1, 3
                Q(iCrd)    = (a(i,ni) * r(Nuc(i),iCrd)     &
                           +  a(j,nj) * r(Nuc(j),iCrd)) /zij
                srs(iCrd)  = Q(iCrd) * ss
                uDip(iCrd) = uDip(iCrd) + cc * srs(iCrd)
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

           do iCrd = 1, 3
               Q(iCrd)   = (a(i,ni) *r(Nuc(i),iCrd)     &
                         +  a(j,nj) *r(Nuc(j),iCrd)) /zij
               srs(iCrd) = Q(iCrd)*ss
           enddo

           cCoef = c(i,ni) * c(j,nj)
           z2    = 2.D0 * zij
           ! Different types of p in the p shell (x, y, z respectively)
           ! The l1-1 term takes into account the different components of 
           ! the shell.
           do l1 = 1, 3        
               k  = i + ((M2-j) * (j-1)) /2 + (l1-1)
               cc = cCoef * P_density(k)
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
           t0  = a(i,ni) * a(j,nj)
           alf = t0 /zij
           ss  = pi32 * exp(-alf*dd) /(zij*sqrt(zij))
           do iCrd = 1, 3
               Q(iCrd)   = (a(i, ni) *r(Nuc(i), iCrd)     &
                         +  a(j, nj) *r(Nuc(j), iCrd)) /zij
               srs(iCrd) = Q(iCrd)*ss
           enddo

           z2  = 2.D0 * zij
           cCoef=c(i,ni)*c(j,nj)

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
                       cc   = P_density(k)*cCoef
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
           zij = a(i,ni) + a(j,nj)
           alf=a(i,ni)*a(j,nj)/zij
           ss=pi32*exp(-alf*dd)/(zij*sqrt(zij))

           do iCrd = 1, 3
               Q(iCrd)   = (a(i, ni) *r(Nuc(i), iCrd)     &
                         +  a(j, nj) *r(Nuc(j), iCrd)) /zij
               srs(iCrd) = Q(iCrd)*ss
           enddo

           z2 = 2.D0*zij
           alf2=2.D0*alf
           cCoef=c(i,ni)*c(j,nj)

           do l1 = 1, 3
               t1=Q(l1)-r(Nuc(i),l1)
               ps=ss*t1
               aux(1)=t1*srs(1)
               aux(2)=t1*srs(2)
               aux(3)=t1*srs(3)
               aux(l1)=aux(l1)+ss/z2

               do l2 = 1, l1
                   t1   = Q(l2) - r(Nuc(i),l2)
                   aux1 = aux*t1
                   f1   = 1.D0
               
                   if (l1.eq.l2) then
                       aux1 = aux1 + srs/z2
                       f1   = sq3
                   endif

                   aux1(l2) = aux1(l2) + ps /z2
                   l12      = l1*(l1-1) /2 + l2

                   ! The order of the d-shell should be as follows:
                   ! xx, yx, yy, zx, zy, zz (11, 21, 22, 31, 32, 33)
                   ii = i  + l12 - 1
                   k  = ii + ((M2-j)*(j-1))/2
                   cc = P_density(k) * cCoef /f1

                   uDip = uDip + cc *aux1
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
           zij = a(i,ni) + a(j,nj)
           alf = a(i,ni) *a(j,nj) /zij
           ss  = pi32*exp(-alf*dd)/(zij*sqrt(zij))

           do iCrd = 1, 3
               Q(iCrd)   = (a(i, ni) *r(Nuc(i), iCrd)     &
                         +  a(j, nj) *r(Nuc(j), iCrd)) /zij
               srs(iCrd) = Q(iCrd)*ss
           enddo

           z2  = 2.D0*zij
           cCoef=c(i,ni)*c(j,nj)

           do l1 = 1, 3
               t1  = Q(l1) - r(Nuc(i),l1)
               pis = ss * t1

               ! aux : <pi|r|s>
               aux    = t1 *srs
               aux(l1)= aux(l1) + ss /z2

               do l2 = 1, l1
                   t1   = Q(l2) - r(Nuc(i),l2)
                   pjs  = ss *t1 
                   dijs = t1 *pis

                   ! aux1 : <pj|r|s>
                   aux1     = t1 *srs
                   aux1(l2) = aux1(l2) + ss /z2
                
                   ! aux2 : <dij|r|s>
                   aux2     = t1 *aux
                   aux2(l2) = aux2(l2) + pis /z2
                   f1 = 1.D0

                   if (l1.eq.l2) then
                       f1   = sq3
                       aux2 = aux2 + srs /z2
                       dijs = dijs + ss /z2
                   endif

                   do l3 = 1, 3
                       t1   = Q(l3) - r(Nuc(j),l3)
                       aux3 = aux2  * t1

                       if (l1.eq.l3) then
                            aux3 = aux3 + aux1 /z2
                       endif

                       if (l2.eq.l3) then
                            aux3 = aux3 + aux  /z2
                       endif
 
                       aux3(l3) = aux3(l3) + dijs /z2
 
                       l12 = l1*(l1-1) /2 + l2
                       ii  = i + l12 - 1
                       jj  = j + l3  - 1
           
                       k  = ii + ((M2-jj)*(jj-1)) /2
                       cc = cCoef /f1 *P_density(k)

                       uDip = uDip + cc *aux3
                   enddo
               enddo    
           enddo        
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
           zij = a(i,ni) + a(j,nj)
           alf = a(i,ni) * a(j,nj) /zij
           ss  = pi32*exp(-alf*dd)/(zij*sqrt(zij))

           do iCrd = 1, 3
               Q(iCrd)   = (a(i, ni) *r(Nuc(i), iCrd)     &
                         +  a(j, nj) *r(Nuc(j), iCrd)) /zij
               srs(iCrd) = Q(iCrd)*ss
           enddo

           alf2  = 2.D0 *alf
           z2    = 2.D0 *zij
           cCoef = c(i,ni)*c(j,nj)

           do l1=1,3
               t1  = Q(l1) - r(Nuc(i),l1)
               pis = ss*t1
          
               ! aux : <pi|r|s>
               aux     = t1*srs
               aux(l1) = aux(l1) + ss /z2

               do l2=1,l1
                   t1   = Q(l2) - r(Nuc(i),l2)
                   pjs  = ss *t1
                   dijs = t1 *pis
                
                   ! aux1 : <pj|r|s>
                   aux1     = t1*srs
                   aux1(l2) = aux1(l2)+ss/z2
 
                   ! Aux2 : (dij|r|s)
                   aux2     = t1*aux
                   aux2(l2) = aux2(l2)+pis/z2
                   f1 = 1.D0
 
                   if (l1.eq.l2) then
                        f1 = sq3
                        aux2 = aux2 + srs /z2
                        dijs = dijs + ss  /z2
                   endif

                   do l3=1,3
                       t1 = Q(l3) - r(Nuc(j),l3)
                       dp = t1 *dijs
              
                       ! aux3 : (dij|r|pk)
                       aux3 = aux2 * t1

                       ! aux4 : (pi|r|pk)
                       aux4 = aux  * t1
 
                       ! aux5 : (pj|r|pk)
                       aux5 = aux1 * t1

                       if (l1.eq.l3) then
                           dp = dp + pjs /z2
                           aux3 = aux3 + aux1 /z2
                           aux4 = aux4 + srs  /z2
                       endif
                       if (l2.eq.l3) then
                           dp = dp + pis /z2
                           aux3 = aux3 + aux /z2
                           aux5 = aux5 + srs /z2
                       endif

                       aux3(l3) = aux3(l3) + dijs /z2
                       aux4(l3) = aux4(l3) + pis  /z2
                       aux5(l3) = aux5(l3) + pjs  /z2

                       do l4=1,l3
                           t1 = Q(l4) - r(Nuc(j),l4)
 
                           ! aux3 : used here for (d|r|d)
                           aux6 = aux3 * t1

                           if (l1.eq.l4) then
                               aux6 = aux6 + aux5 /z2
                           endif

                           if (l2.eq.l4) then
                               aux6 = aux6 + aux4 /z2
                           endif

                           f2=1.D0
                           if (l3.eq.l4) then
                               f2   = sq3
                               aux6 = aux6 + aux2 /z2
                           endif

                           aux6(l4) = aux6(l4) + dp /z2
                           l12 = l1 *(l1-1) /2 + l2
                           l34 = l3 *(l3-1) /2 + l4
                           ii  = i + l12 -1
                           jj  = j + l34 -1

                           if (ii.ge.jj) then
                               k    = ii + ((M2-jj)*(jj-1)) /2
                               cc   = cCoef / (f1*f2) * P_density(k)
                               uDip = uDip + cc*aux6
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
   uDip = (uDipAt - uDip*factor) * 2.54D0
 
   return
end subroutine dip
