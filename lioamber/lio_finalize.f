!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       subroutine lio_finalize()
!--------------------------------------------------------------------!
! DEALLOCATION OF GLOBAL VARIABLES
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       use garcha_mod
       implicit none
!--------------------------------------------------------------------!
       if (allocated(Smat))    deallocate(Smat)
       if (allocated(RealRho)) deallocate(RealRho)
!--------------------------------------------------------------------!
       deallocate(r,v,rqm, Em, Rm, pc,Iz, nnat,
     > af,c,a,cx,ax,cd,ad,B,Nuc,ncont,Nucx,ncontx,Nucd
     > ,ncontd, indexii, indexiid, RMM, X, XX)
       deallocate(old1,old2,old3)
       deallocate(natomc,nnps,nnpp,nnpd,nns)
       deallocate(nnd,nnp,atmin,jatc,d)
       call g2g_deinit()
       call aint_deinit()
!--------------------------------------------------------------------!
       return;end subroutine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
