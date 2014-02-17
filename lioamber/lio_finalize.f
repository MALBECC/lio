         subroutine lio_finalize() 
         use garcha_mod
      deallocate (r,v,rqm, Em, Rm, pc,Iz, nnat,
     > af,c,a,cx,ax,cd,ad,B,Nuc,ncont,Nucx,ncontx,Nucd
     >  ,ncontd, indexii, indexiid, RMM, X, XX, old1, old2, old3)
         deallocate (natomc,nnps,nnpp,nnpd,nns)
         deallocate (nnd,nnp,atmin,jatc,d)
      call g2g_deinit () 
      end subroutine lio_finalize 
 
