subroutine obtain_sigma(wfunc,wfunc_old,coef,coef_old,sigma,&
                        ndets,all_states,M,NCO,Nvirt)
use garcha_mod, only: r, natom
use basis_data, only: a, c
use fstsh_data, only: Sovl_old, Sovl_now, a_old, c_old, r_old, tsh_file, &
                      tsh_nucStep, phases_old, tsh_time_dt
   implicit none

   integer, intent(in)  :: ndets, all_states, M, NCO, Nvirt
   LIODBLE, intent(in)  :: wfunc(ndets,all_states), wfunc_old(ndets,all_states)
   LIODBLE, intent(in)  :: coef(M,M), coef_old(M,M)
   LIODBLE, intent(out) :: sigma(all_states,all_states)

   integer, dimension(:,:), allocatable :: coupling_screening
   LIODBLE, dimension(:,:), allocatable :: diff, Smat, Sbig
   LIODBLE, dimension(:,:), allocatable :: part_cio
   LIODBLE, dimension(:), allocatable :: phases

   integer :: ii, jj

   if ( tsh_nucStep == 0 ) return
   
   ! Form distances differences
   allocate(diff(natom,natom))
   call form_diff(r, r_old, diff, natom)

   ! Obtain Overlap at differences times
   allocate(Smat(M,M)) ! S[t,t-dt]
   call form_Sprojection(Smat,diff,r,r_old,a,a_old,c,c_old)
   deallocate(diff)

   ! Form S big ( with all times )
   allocate(Sbig(M*2,M*2))
   Sbig(1:M,1:M) = Sovl_now(:,:)         ! [t,t]
   Sbig(M+1:M*2,M+1:M*2) = Sovl_old(:,:) ! [t-dt,t-dt]
   Sbig(1:M,M+1:M*2) = Smat(:,:)         ! [t,t-dt]
   Sbig(M+1:M*2,1:M) = -Smat(:,:)        ! [t-dt,t]
   deallocate(Smat)

   ! Cioverlap: obtain sigma
   allocate(coupling_screening(all_states,all_states))
   call kind_coupling(coupling_screening,all_states)

   allocate(phases(all_states),part_cio(all_states,all_states));
   call g2g_cioverlap(wfunc,wfunc_old,coef,coef_old,Sbig,part_cio,&
                      coupling_screening,phases,phases_old,    &
                      M,NCO,Nvirt,all_states,ndets)
   phases_old = phases
   deallocate(phases,Sbig)
   do ii=1,all_states
      if ( 1.0d0-abs(part_cio(ii,ii)) > 0.1d0 ) write(tsh_file,"(3X,A,I2,A,F10.5)") & 
                                                "State= ",ii," Can not Overlap",part_cio(ii,ii)
      do jj=ii,all_states
          sigma(ii,jj) = (part_cio(ii,jj) - part_cio(jj,ii)) * 1.0d0 / (2.0d0*tsh_time_dt)
          sigma(jj,ii) = -sigma(ii,jj)
      enddo
   enddo

   deallocate(coupling_screening,part_cio)
end subroutine obtain_sigma

subroutine kind_coupling( Scree, all_states )
use fstsh_data, only: type_coupling, current_state
! type_coupling:
! 0 = all vs all states
! 1 = only the current state with the upper and lower states

! current_state = this includes ground state ( 1 = GS )
implicit none

   integer, intent(in)  :: all_states
   integer, intent(out) :: Scree(all_states,all_states)
  
   integer :: ii

   Scree = 0
   if ( current_state /= 1 ) then
      select case ( type_coupling )
         case(0)
            Scree = 1

         case(1)
            Scree(current_state,current_state+1)=1
            Scree(current_state+1,current_state)=1
            Scree(current_state,current_state-1)=1
            Scree(current_state-1,current_state)=1

         case default
            print*, "ERROR!: Invalue Value Of kind_coupling=",type_coupling

      end select
   else
      Scree(1,2)=1
      Scree(2,1)=1

   endif

   do ii=1,all_states
      Scree(ii,ii) = 1
   enddo
end subroutine kind_coupling

subroutine form_diff(Rc,Ro,Diff,natom)
   implicit none

   integer, intent(in)  :: natom
   LIODBLE, intent(in)  :: Rc(natom,3), Ro(natom,3)
   LIODBLE, intent(out) :: Diff(natom,natom)

   integer :: ii, jj
   LIODBLE :: tx, ty, tz

   do ii=1,natom
   do jj=1,natom
      tx = (Rc(ii,1)-Ro(jj,1))**2
      ty = (Rc(ii,2)-Ro(jj,2))**2
      tz = (Rc(ii,3)-Ro(jj,3))**2

      Diff(ii,jj) = tx + ty + tz
   enddo
   enddo
end subroutine form_diff

subroutine form_Sprojection(Smat,dif,r,rOld,a,aOld,c,cOld)
use basis_data   , only: Nuc, ncont, NORM, M, nshell
use liosubs_math , only: FUNCT
use constants_mod, only: pi, pi32
   implicit none

   LIODBLE, intent(in)  :: dif(:,:), r(:,:), a(:,:), c(:,:)
   LIODBLE, intent(in)  :: rOld(:,:), aOld(:,:), cOld(:,:)
   LIODBLE, intent(out) :: Smat(:,:)

   integer           :: my_natom, igpu, i_ind, j_ind, k_ind, ifunct, jfunct, &
                        iatom, jatom, nci, ncj, l1, l2, l3, l4, l12, l34,    &
                        MM, ns, np, nd, M2, M5, M11
   LIODBLE  :: ovlap, uf, cc_f, Q(3), temp, sq3, alf, alf2, ccoef,    &
                        t0, t1, t2, f1, f2, tn, tna, Z2, Zij, ss, ps, dd, sks, &
                        p0s, p1s, p2s, p3s, pi0p, pi1p, piks, pikpk, pipk,     &
                        pis, pj0s, pj1s, pj2s, pj0p, pj1p, pjkpk, pjks, pjpk,  &
                        pjs, pks, dijs, dijpk, dijks, dijkpk, d0s, d0p, d1p,   &
                        d1s, d2s
   LIODBLE  :: Q_w(3), alf2_w, alf_w, ccoef_w, dd_w, ovlap_w, sks_w,  &
                        ss_w, t1_w, Zij_w, cc_f_w, pks_w, ps_w, Z2_w, dijks_w, &
                        piks_w, pis_w, pjks_w, dijs_w, pjs_w


   sq3 = 1.D0
   if (NORM) sq3 = sqrt(3.D0)

   ns  = nshell(0); np = nshell(1); nd = nshell(2)
   MM  = M*(M+1)/2
   M2  = 2*M

   Smat = 0.0D0

   ! First loop - (s|s)
   do ifunct = 1, ns ! t
   do jfunct = 1, ns ! t-dt
      dd = dif(Nuc(ifunct), Nuc(jfunct))

      do nci = 1, ncont(ifunct)
      do ncj = 1, ncont(jfunct)
         Zij  = a(ifunct,nci) + aOld(jfunct,ncj)
         alf  = a(ifunct,nci) * aOld(jfunct,ncj) / Zij
         alf2 = alf * 2.D0
         Q(1) = (a(ifunct,nci) * r(Nuc(ifunct),1) + &
                 aOld(jfunct,ncj) * rOld(Nuc(jfunct),1)) / Zij
         Q(2) = (a(ifunct,nci) * r(Nuc(ifunct),2) + &
                 aOld(jfunct,ncj) * rOld(Nuc(jfunct),2)) / Zij
         Q(3) = (a(ifunct,nci) * r(Nuc(ifunct),3) + &
                 aOld(jfunct,ncj) * rOld(Nuc(jfunct),3)) / Zij

         ccoef = c(ifunct,nci) * cOld(jfunct,ncj)
         ss    = pi32 * exp(-alf*dd) / (Zij * sqrt(Zij))
         ovlap = ss
         tn    = alf * (3.D0 - alf2*dd) * ovlap

         Smat(ifunct,jfunct) = Smat(ifunct,jfunct) + ccoef * ovlap
      enddo
      enddo
   enddo
   enddo

   ! Second loop - (p|s)
   do ifunct = ns+1, ns+np, 3
   do jfunct = 1, ns
      ! nothing = <i(t)|j(t-dt)> == <p(t)|s(t-dt)>
      ! _w      = <j(t)|i(t-dt)> == <s(t)|p(t-dt)>

      dd = dif(Nuc(ifunct), Nuc(jfunct))
      dd_w = dif(Nuc(jfunct), Nuc(ifunct))

      do nci = 1, ncont(ifunct)
      do ncj = 1, ncont(jfunct)
         Zij = a(ifunct, nci) + aOld(jfunct, ncj)
         Q(1) = (a(ifunct,nci) * r(Nuc(ifunct),1) + &
                 aOld(jfunct,ncj) * rOld(Nuc(jfunct),1)) / Zij
         Q(2) = (a(ifunct,nci) * r(Nuc(ifunct),2) + &
                 aOld(jfunct,ncj) * rOld(Nuc(jfunct),2)) / Zij
         Q(3) = (a(ifunct,nci) * r(Nuc(ifunct),3) + &
                 aOld(jfunct,ncj) * rOld(Nuc(jfunct),3)) / Zij

         Zij_w = a(jfunct, ncj) + aOld(ifunct, nci)
         Q_w(1) = (a(jfunct,ncj) * r(Nuc(jfunct),1) + &
                 aOld(ifunct,nci) * rOld(Nuc(ifunct),1)) / Zij_w
         Q_w(2) = (a(jfunct,ncj) * r(Nuc(jfunct),2) + &
                 aOld(ifunct,nci) * rOld(Nuc(ifunct),2)) / Zij_w
         Q_w(3) = (a(jfunct,ncj) * r(Nuc(jfunct),3) + &
                 aOld(ifunct,nci) * rOld(Nuc(ifunct),3)) / Zij_w

         alf   = a(ifunct, nci) * aOld(jfunct,ncj) / Zij
         alf2  = 2.D0 * alf
         ss    = pi32 * exp(-alf*dd) / (Zij * sqrt(Zij))
         sks   = alf * (3.D0 - alf2*dd) * ss
         ccoef = c(ifunct,nci) * cOld(jfunct,ncj)

         alf_w   = a(jfunct, ncj) * aOld(ifunct,nci) / Zij_w
         alf2_w  = 2.D0 * alf_w
         ss_w    = pi32 * exp(-alf_w*dd_w) / (Zij_w * sqrt(Zij_w))
         sks_w   = alf_w * (3.D0 - alf2_w*dd_w) * ss_w
         ccoef_w = c(jfunct,ncj) * cOld(ifunct,nci)

         ! l2: different p in the p shell ( x,y,z respectively)
         do l2 = 1, 3
            i_ind = ifunct + l2 -1

            t1 = Q(l2) - r(Nuc(ifunct),l2)
            ovlap = t1 * ss
            Smat(i_ind, jfunct) = Smat(i_ind ,jfunct) + ovlap * ccoef

            t1_w = Q_w(l2) - rOld(Nuc(ifunct),l2)
            ovlap_w = t1_w * ss_w
            Smat(jfunct, i_ind) = Smat(jfunct, i_ind) + ovlap_w * ccoef_w
         enddo
      enddo
      enddo
   enddo
   enddo

   ! Third loop - (p|p)
   do ifunct = ns+1, ns+np , 3 ! t
   do jfunct = ns+1, ns+np , 3 ! t-dt
      dd = dif(Nuc(ifunct), Nuc(jfunct))

      do nci = 1, ncont(ifunct)
      do ncj = 1, ncont(jfunct)
         Zij = a(ifunct, nci) + aOld(jfunct, ncj)
         Z2  = 2.D0 * Zij
         Q(1) = (a(ifunct,nci) * r(Nuc(ifunct),1) + &
                 aOld(jfunct,ncj) * rOld(Nuc(jfunct),1)) / Zij
         Q(2) = (a(ifunct,nci) * r(Nuc(ifunct),2) + &
                 aOld(jfunct,ncj) * rOld(Nuc(jfunct),2)) / Zij
         Q(3) = (a(ifunct,nci) * r(Nuc(ifunct),3) + &
                 aOld(jfunct,ncj) * rOld(Nuc(jfunct),3)) / Zij

         alf   = a(ifunct, nci) * aOld(jfunct,ncj) / Zij
         alf2  = 2.D0 * alf
         ss    = pi32 * exp(-alf*dd) / (Zij * sqrt(Zij))
         sks   = alf * (3.D0 - alf2*dd) * ss
         ccoef = c(ifunct,nci) * cOld(jfunct,ncj)

         do l1 = 1, 3
            t1  = Q(l1) - r(Nuc(ifunct),l1)
            ps  = ss * t1
            pks = sks * t1 + alf2 * ps

            do l2 = 1, 3
               t1 = Q(l2) - rOld(Nuc(jfunct),l2)

               ovlap = t1 * ps
               tn    = t1 * pks + alf2 * ovlap
               if (l1 .eq. l2) then
                  ovlap = ovlap + ss / Z2
                  tn    = tn + (sks + alf2*ss) / Z2
               endif

               i_ind = ifunct + l1 -1
               j_ind = jfunct + l2 -1
               Smat(i_ind,j_ind) = Smat(i_ind,j_ind) + ovlap * ccoef
            enddo
         enddo
      enddo
      enddo
   enddo
   enddo

   ! Fourth loop - (d|s)
  do ifunct = ns+np+1, M, 6
  do jfunct = 1, ns
     ! nothing = <i(t)|j(t-dt)> == <d(t)|s(t-dt)>
     ! _w      = <j(t)|i(t-dt)> == <s(t)|d(t-dt)>

     dd = dif(Nuc(ifunct), Nuc(jfunct))
     dd_w = dif(Nuc(jfunct), Nuc(ifunct))

     do nci = 1, ncont(ifunct)
     do ncj = 1, ncont(jfunct)
        Zij  = a(ifunct, nci) + aOld(jfunct, ncj)
        Z2   = 2.D0 * Zij
        Q(1) = (a(ifunct,nci) * r(Nuc(ifunct),1) + &
                aOld(jfunct,ncj) * rOld(Nuc(jfunct),1)) / Zij
        Q(2) = (a(ifunct,nci) * r(Nuc(ifunct),2) + &
                aOld(jfunct,ncj) * rOld(Nuc(jfunct),2)) / Zij
        Q(3) = (a(ifunct,nci) * r(Nuc(ifunct),3) + &
                aOld(jfunct,ncj) * rOld(Nuc(jfunct),3)) / Zij
        Zij_w  = a(jfunct, ncj) + aOld(ifunct, nci)
        Z2_w   = 2.D0 * Zij_w
        Q_w(1) = (a(jfunct,ncj) * r(Nuc(jfunct),1) + &
                aOld(ifunct,nci) * rOld(Nuc(ifunct),1)) / Zij_w
        Q_w(2) = (a(jfunct,ncj) * r(Nuc(jfunct),2) + &
                aOld(ifunct,nci) * rOld(Nuc(ifunct),2)) / Zij_w
        Q_w(3) = (a(jfunct,ncj) * r(Nuc(jfunct),3) + &
                aOld(ifunct,nci) * rOld(Nuc(ifunct),3)) / Zij_w

        alf   = a(ifunct, nci) * aOld(jfunct,ncj) / Zij
        alf2  = 2.D0 * alf
        ss    = pi32 * exp(-alf*dd) / (Zij * sqrt(Zij))
        sks   = alf * (3.D0 - alf2*dd) * ss
        ccoef = c(ifunct,nci) * cOld(jfunct,ncj)

        alf_w   = a(jfunct, ncj) * aOld(ifunct,nci) / Zij_w
        alf2_w  = 2.D0 * alf_w
        ss_w    = pi32 * exp(-alf_w*dd_w) / (Zij_w * sqrt(Zij_w))
        sks_w   = alf_w * (3.D0 - alf2_w*dd_w) * ss_w
        ccoef_w = c(jfunct,ncj) * cOld(ifunct,nci)

        do l1 = 1, 3
           t1  = Q(l1) - r(Nuc(ifunct),l1)
           ps  = ss  * t1
           pks = sks * t1 + alf2 * ps

           t1_w  = Q_w(l1) - rOld(Nuc(ifunct),l1)
           ps_w  = ss_w  * t1_w
           pks_w = sks_w * t1_w + alf2_w * ps_w

           do l2 = 1, l1
              t1    = Q(l2) - r(Nuc(ifunct),l2)
              ovlap = t1 * ps

              t1_w    = Q_w(l2) - rOld(Nuc(ifunct),l2)
              ovlap_w = t1_w * ps_w

              f1 = 1.D0
              if (l1 .eq. l2) then
                 ovlap = ovlap + ss / Z2
                 ovlap_w = ovlap_w + ss_w / Z2_w

                 f1    = sq3
              endif
              ! Ordering of the d shell: xx, yx, yy, zx, zy, zz
              ! (11, 21, 22, 31, 32, 33)
              l12   = l1 * (l1-1)/2 + l2
              i_ind = ifunct + l12 -1

              cc_f = ccoef / f1
              Smat(i_ind,jfunct) = Smat(i_ind,jfunct)  + ovlap * cc_f

              cc_f_w = ccoef_w / f1
              Smat(jfunct,i_ind) = Smat(jfunct,i_ind)  + ovlap_w * cc_f_w
           enddo
        enddo
     enddo
     enddo
  enddo
  enddo

  ! Fifth loop - (d|p)
  do ifunct = ns+np+1, M    , 6
  do jfunct = ns+1   , ns+np, 3
     ! nothing = <i(t)|j(t-dt)> == <d(t)|p(t-dt)>
     ! _w      = <j(t)|i(t-dt)> == <p(t)|d(t-dt)>

     dd = dif(Nuc(ifunct), Nuc(jfunct))
     dd_w = dif(Nuc(jfunct), Nuc(ifunct))

     do nci = 1, ncont(ifunct)
     do ncj = 1, ncont(jfunct)
        Zij  = a(ifunct, nci) + aOld(jfunct, ncj)
        Z2   = 2.D0 * Zij
        Q(1) = (a(ifunct,nci) * r(Nuc(ifunct),1) + &
                aOld(jfunct,ncj) * rOld(Nuc(jfunct),1)) / Zij
        Q(2) = (a(ifunct,nci) * r(Nuc(ifunct),2) + &
                aOld(jfunct,ncj) * rOld(Nuc(jfunct),2)) / Zij
        Q(3) = (a(ifunct,nci) * r(Nuc(ifunct),3) + &
                aOld(jfunct,ncj) * rOld(Nuc(jfunct),3)) / Zij

        Zij_w  = a(jfunct, ncj) + aOld(ifunct, nci)
        Z2_w   = 2.D0 * Zij_w
        Q_w(1) = (a(jfunct,ncj) * r(Nuc(jfunct),1) + &
                aOld(ifunct,nci) * rOld(Nuc(ifunct),1)) / Zij_w
        Q_w(2) = (a(jfunct,ncj) * r(Nuc(jfunct),2) + &
                aOld(ifunct,nci) * rOld(Nuc(ifunct),2)) / Zij_w
        Q_w(3) = (a(jfunct,ncj) * r(Nuc(jfunct),3) + &
                aOld(ifunct,nci) * rOld(Nuc(ifunct),3)) / Zij_w

        alf   = a(ifunct, nci) * aOld(jfunct,ncj) / Zij
        alf2  = 2.D0 * alf
        ss    = pi32 * exp(-alf*dd) / (Zij * sqrt(Zij))
        sks   = alf * (3.D0 - alf2*dd) * ss
        ccoef = c(ifunct,nci) * cOld(jfunct,ncj)

        alf_w   = a(jfunct, ncj) * aOld(ifunct,nci) / Zij_w
        alf2_w  = 2.D0 * alf_w
        ss_w    = pi32 * exp(-alf_w*dd_w) / (Zij_w * sqrt(Zij_w))
        sks_w   = alf_w * (3.D0 - alf2_w*dd_w) * ss_w
        ccoef_w = c(jfunct,ncj) * cOld(ifunct,nci)

        do l1 = 1, 3
           t1   = Q(l1) - r(Nuc(ifunct),l1)
           pis  = ss  * t1
           piks = sks * t1 + alf2 * pis

           t1_w   = Q_w(l1) - rOld(Nuc(ifunct),l1)
           pis_w  = ss_w  * t1_w
           piks_w = sks_w * t1_w + alf2_w * pis_w

           do l2 = 1, l1
              t1    = Q(l2) - r(Nuc(ifunct),l2)
              pjs   = ss  * t1
              pjks  = sks * t1 + alf2 * pjs
              dijs  = t1 * pis
              dijks = t1 * piks

              t1_w    = Q_w(l2) - rOld(Nuc(ifunct),l2)
              pjs_w   = ss_w  * t1_w
              pjks_w  = sks_w * t1_w + alf2_w * pjs_w
              dijs_w  = t1_w * pis_w
              dijks_w = t1_w * piks_w

              f1    = 1.D0
              if (l1 .eq. l2) then
                 f1    = sq3
                 dijs  = dijs  + ss / Z2
                 dijks = dijks + sks/ Z2 - alf2 * ss / (2.D0*a(ifunct,nci))

                 dijs_w  = dijs_w  + ss_w / Z2_w
                 dijks_w = dijks_w + sks_w/ Z2_w - alf2_w * ss_w / (2.D0*aOld(ifunct,nci))
              endif

              dijks = dijks + alf2 * dijs
              dijks_w = dijks_w + alf2_w * dijs_w
              do l3 = 1, 3
                 t1    = Q(l3) - rOld(Nuc(jfunct),l3)
                 ovlap = t1 * dijs

                 t1_w    = Q_w(l3) - r(Nuc(jfunct),l3)
                 ovlap_w = t1_w * dijs_w
                 if (l1 .eq. l3) then
                    ovlap = ovlap + pjs  / Z2

                    ovlap_w = ovlap_w + pjs_w  / Z2_w
                 endif
                 if (l2 .eq. l3) then
                    ovlap = ovlap + pis  / Z2

                    ovlap_w = ovlap_w + pis_w  / Z2_w
                 endif

                 cc_f = ccoef / f1
                 cc_f_w = ccoef_w / f1

                 l12 = l1 * (l1-1)/2 + l2
                 i_ind = ifunct + l12 -1
                 j_ind = jfunct + l3  -1

                 Smat(i_ind,j_ind) = Smat(i_ind,j_ind) + cc_f * ovlap
                 Smat(j_ind,i_ind) = Smat(j_ind,i_ind) + cc_f_w * ovlap_w
              enddo
           enddo
        enddo
     enddo
     enddo
  enddo
  enddo
  ! Sixth and final loop - (d|d)
  do ifunct = ns+np+1, M, 6 ! t
  do jfunct = ns+np+1, M, 6 ! t-dt
     dd = dif(Nuc(ifunct), Nuc(jfunct))

     do nci = 1, ncont(ifunct)
     do ncj = 1, ncont(jfunct)
        Zij  = a(ifunct, nci) + aOld(jfunct, ncj)
        Z2   = 2.D0 * Zij
        Q(1) = (a(ifunct,nci) * r(Nuc(ifunct),1) + &
                aOld(jfunct,ncj) * rOld(Nuc(jfunct),1)) / Zij
        Q(2) = (a(ifunct,nci) * r(Nuc(ifunct),2) + &
                aOld(jfunct,ncj) * rOld(Nuc(jfunct),2)) / Zij
        Q(3) = (a(ifunct,nci) * r(Nuc(ifunct),3) + &
                aOld(jfunct,ncj) * rOld(Nuc(jfunct),3)) / Zij

        alf   = a(ifunct, nci) * aOld(jfunct,ncj) / Zij
        alf2  = 2.D0 * alf
        ss    = pi32 * exp(-alf*dd) / (Zij * sqrt(Zij))
        sks   = alf * (3.D0 - alf2*dd) * ss
        t0    = ss / Z2
        ccoef = c(ifunct,nci) * cOld(jfunct,ncj)

        do l1 = 1, 3
           t1   = Q(l1) - r(Nuc(ifunct),l1)
           pis  = ss  * t1
           piks = sks * t1 + alf2 * pis
           do l2 = 1, l1
              t1   = Q(l2) - r(Nuc(ifunct),l2)
              pjs  = ss  * t1
              pjks = sks * t1 + alf2 * pjs

              f1    = 1.D0
              dijs  = t1 * pis
              dijks = t1 * piks
              if (l1 .eq. l2) then
                 f1    = sq3
                 dijs  = dijs  + t0
                 dijks = dijks + sks / Z2 - alf2 * ss / (2.D0*a(ifunct,nci))
              endif
              dijks = dijks + alf2 * dijs

              do l3 = 1, 3
                 t2 = Q(l3) - rOld(Nuc(jfunct),l3)

                 pipk   = t2 * pis
                 pjpk   = t2 * pjs
                 pikpk  = t2 * piks
                 pjkpk  = t2 * pjks
                 dijpk  = t2 * dijs
                 dijkpk = t2 * dijks
                 if (l1 .eq. l3) then
                    pipk   = pipk   + t0
                    dijpk  = dijpk  + pjs  / Z2
                    pikpk  = pikpk  + sks  / Z2
                    dijkpk = dijkpk + pjks / Z2
                 endif
                 if (l2 .eq. l3) then
                    pjpk   = pjpk   + t0
                    dijpk  = dijpk  + pis  / Z2
                    pjkpk  = pjkpk  + sks  / Z2
                    dijkpk = dijkpk + piks / Z2
                 endif
                 pikpk  = pikpk  + alf2 * pipk
                 pjkpk  = pjkpk  + alf2 * pjpk
                 dijkpk = dijkpk + alf2 * dijpk

                 do l4 = 1, l3
                    t1    = Q(l4) - rOld(Nuc(jfunct),l4)
                    ovlap = t1 * dijpk
                    tn    = t1 * dijkpk
                    f2    = 1.D0

                    if (l1 .eq. l4) then
                       ovlap = ovlap + pjpk  / Z2
                       tn    = tn    + pjkpk / Z2
                    endif
                    if (l2 .eq. l4) then
                       ovlap = ovlap + pipk  / Z2
                       tn    = tn    + pikpk / Z2
                    endif
                    if (l3.eq.l4) then
                       ovlap = ovlap+ dijs  / Z2
                       tn    = tn   + dijks / Z2 - &
                               alf2 * dijs / (2.D0*aOld(jfunct,ncj))
                       f2    = sq3
                    endif

                    ! l12 and l34 range from 1 to 6, spanning the d shell in
                    ! the order xx, xy, yy, zx, zy, zz.
                    l12   = l1 * (l1-1)/2 + l2
                    l34   = l3 * (l3-1)/2 + l4
                    i_ind = ifunct + l12 -1
                    j_ind = jfunct + l34 -1

                    cc_f  = ccoef / (f1 * f2)
                    Smat(i_ind,j_ind) = Smat(i_ind,j_ind) + ovlap * cc_f
                 enddo
              enddo
           enddo
        enddo
     enddo
     enddo
  enddo
  enddo
end subroutine form_Sprojection

