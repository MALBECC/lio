c subroutine that calculates the rcut and block QM-MM only in the first step

	subroutine qmmm_lst_blk(na_u,nac,natot,nroaa,atxres,rclas,
     .  rcorteqmmm,radiobloqmmm,blockqmmm,listqmmm,rcorteqm,slabel)

	use precision 
	use ionew
	use sys
	implicit none
        integer na_u,nac,natot,iu,nroaa,atxres(20000)
        double precision cm(3,20000),dist,dist2,rcorteqm,
     .  rcorteqmmm,radiobloqmmm,rclas(3,natot),Ang
        integer i,j,k,l,count,blockqmmm(nac),listqmmm(nac)
        logical bloqmmm,liqmmm,found
        character slabel*20,fname*24,paste*24
        external paste
        Ang    = 1._dp / 0.529177_dp

c change units
        rclas=rclas/Ang

c calculate center of masses of all residues
	k=na_u+1
	do i=1,nroaa
	 do j=1,atxres(i)
	cm(1:3,i)=cm(1:3,i)+rclas(1:3,k)
	k=k+1
	 enddo
	cm(1:3,i)=cm(1:3,i)/atxres(i)
	enddo

c rcut QM-MM
        found=.false.
        liqmmm=.false.
        if(rcorteqmmm.le.99.0) liqmmm=.true.
        if(liqmmm) then
c find file name
        fname = paste(slabel,'.lst')
c check if the input file exists
        inquire( file=fname, exist=found )
       if (found) then
c Open file
        call io_assign( iu )
        open( iu, file=fname, status='old' )
c read listqmmm
        do i=1,nac
        read(iu,*,err=1,end=1) listqmmm(i)
        enddo
c Close file
        call io_close( iu )
        write(6,'(/a)') 'qm-mm: Reading QM-MM neighbours list from file'
       else
c QM-MM neigh list
      write(6,'(/a,f12.6)') 
     .'qm-mm: cut off radius QM-MM (Ang):',rcorteqmmm
      if(rcorteqm.ge.rcorteqmmm) then
      call die("qm-mm: cut off QM-MM must be greater than cut off QM")
      endif
        dist=0.0
        dist2=rcorteqmmm**2
	listqmmm=1
	do l=1,na_u
	k=1
	 do i=1,nroaa
	  dist=(rclas(1,l)-cm(1,i))**2+
     .         (rclas(2,l)-cm(2,i))**2+
     .         (rclas(3,l)-cm(3,i))**2
	  do j=1,atxres(i)
	   if(dist.le.dist2) listqmmm(k)=0
	   k=k+1
	  enddo
	 enddo
	enddo
	count=0
	do i=1,nac
	if(listqmmm(i).eq.0) count=count+1
	enddo
	write(6,'(/a,2x,i5)') 
     .  'qm-mm: Number of QM-MM interacting atoms:',count
c Open file
        call io_assign( iu )
        open( iu, file=fname, form='formatted', status='unknown' )
c write listqmmm
        do i=1,nac
        write(iu,*) listqmmm(i)
        enddo
c Close file
        call io_close( iu )
c        write(6,'(/a)') 'qm-mm: Writing QM-MM neighbours list in file'
       endif
        endif

c block QM-MM
        found=.false.
        bloqmmm=.false.
        if(radiobloqmmm.le.99.0) bloqmmm=.true.
        if(bloqmmm) then
c find file name
        fname = paste(slabel,'.blk')
c check if the input file exists
        inquire( file=fname, exist=found )
       if (found) then
c Open file
        call io_assign( iu )
        open( iu, file=fname, status='old' )
c read blockqmmm
        do i=1,nac
        read(iu,*,err=2,end=2) blockqmmm(i)
        enddo
c Close file
        call io_close( iu )
        write(6,'(/a)') 'qm-mm: Reading blocked QM-MM atoms from file'
       else
c fixing MM atoms beyond block cut off  
        write(6,'(/a,f12.6)') 
     .  'qm-mm: cut off radius Block (Ang):',radiobloqmmm
        dist=0.0
        dist2=radiobloqmmm**2
        blockqmmm=1
        do l=1,na_u
        k=1
         do i=1,nroaa
          dist=(rclas(1,l)-cm(1,i))**2+
     .         (rclas(2,l)-cm(2,i))**2+
     .         (rclas(3,l)-cm(3,i))**2
          do j=1,atxres(i)
           if(dist.le.dist2) blockqmmm(k)=0
           k=k+1
          enddo
         enddo
        enddo
        count=0
        do i=1,nac
        if(radiobloqmmm.eq.0.0) blockqmmm(i)=1
        if(blockqmmm(i).eq.0) count=count+1
        enddo
        write(6,'(/a,2x,i5)') 
     .  'qm-mm: Number of QM-MM free atoms:',count 
c Open file
        call io_assign( iu )
        open( iu, file=fname, form='formatted', status='unknown' )
c write blockqmmm
        do i=1,nac
        write(iu,*) blockqmmm(i)
        enddo
c Close file
        call io_close( iu )
c        write(6,'(/a)') 'qm-mm: Writing blocked QM-MM atoms in file'
       endif
        endif

c change units
      rclas=rclas*Ang
      rcorteqm=rcorteqm*Ang       
      call flush(6)

      return
 1    stop 'qm-mm: Problem reading from .lst file'
 2    stop 'qm-mm: Problem reading from .blk file'
      end

