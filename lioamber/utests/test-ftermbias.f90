!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       include "../vvfterm.f"
       program vvfterm_test01
       implicit none
       real*8   :: sqsmat(3,3),outmat(3,3)
       integer  :: vector(3)
       real*8   :: weight

       write(*,*) ''
       write(*,*) ''


!------------------------------------------------------------------------------!
       sqsmat(1,1)=0.0d0
       sqsmat(1,2)=0.0d0
       sqsmat(1,3)=0.0d0
       sqsmat(2,1)=0.0d0
       sqsmat(2,2)=0.0d0
       sqsmat(2,3)=0.0d0
       sqsmat(3,1)=0.0d0
       sqsmat(3,2)=0.0d0
       sqsmat(3,3)=0.0d0

       vector(1)=1
       vector(2)=1
       vector(3)=1
       weight=1.0d0

       outmat(:,:)=0.0d0
       call vvfterm(3,sqsmat,vector,weight,outmat)
       write(*,*) 'null matrix from null mat:'
       write(*,*) outmat(1,1),outmat(1,2),outmat(1,3)
       write(*,*) outmat(2,1),outmat(2,2),outmat(2,3)
       write(*,*) outmat(3,1),outmat(3,2),outmat(3,3)
       write(*,*) ''


!------------------------------------------------------------------------------!
       sqsmat(1,1)=1.0d0
       sqsmat(1,2)=3.0d0
       sqsmat(1,3)=2.0d0
       sqsmat(2,1)=5.0d0
       sqsmat(2,2)=6.0d0
       sqsmat(2,3)=3.0d0
       sqsmat(3,1)=4.0d0
       sqsmat(3,2)=2.0d0
       sqsmat(3,3)=1.0d0

       vector(1)=0
       vector(2)=0
       vector(3)=0
       weight=1.0d0

       outmat(:,:)=0.0d0
       call vvfterm(3,sqsmat,vector,weight,outmat)
       write(*,*) 'null matrix from null vec:'
       write(*,*) outmat(1,1),outmat(1,2),outmat(1,3)
       write(*,*) outmat(2,1),outmat(2,2),outmat(2,3)
       write(*,*) outmat(3,1),outmat(3,2),outmat(3,3)
       write(*,*) ''


!------------------------------------------------------------------------------!
       sqsmat(1,1)=1.0d0
       sqsmat(1,2)=3.0d0
       sqsmat(1,3)=2.0d0
       sqsmat(2,1)=5.0d0
       sqsmat(2,2)=6.0d0
       sqsmat(2,3)=3.0d0
       sqsmat(3,1)=4.0d0
       sqsmat(3,2)=2.0d0
       sqsmat(3,3)=1.0d0

       vector(1)=1
       vector(2)=1
       vector(3)=1
       weight=0.0d0

       outmat(:,:)=0.0d0
       call vvfterm(3,sqsmat,vector,weight,outmat)
       write(*,*) 'null matrix from null weight:'
       write(*,*) outmat(1,1),outmat(1,2),outmat(1,3)
       write(*,*) outmat(2,1),outmat(2,2),outmat(2,3)
       write(*,*) outmat(3,1),outmat(3,2),outmat(3,3)
       write(*,*) ''


!------------------------------------------------------------------------------!
       sqsmat(1,1)=2.0d0
       sqsmat(1,2)=0.0d0
       sqsmat(1,3)=0.0d0
       sqsmat(2,1)=0.0d0
       sqsmat(2,2)=3.0d0
       sqsmat(2,3)=0.0d0
       sqsmat(3,1)=0.0d0
       sqsmat(3,2)=0.0d0
       sqsmat(3,3)=4.0d0

       vector(1)=1
       vector(2)=1
       vector(3)=1
       weight=1.0d0

       outmat(:,:)=0.0d0
       call vvfterm(3,sqsmat,vector,weight,outmat)
       write(*,*) '4,9,16 id matrix:'
       write(*,*) outmat(1,1),outmat(1,2),outmat(1,3)
       write(*,*) outmat(2,1),outmat(2,2),outmat(2,3)
       write(*,*) outmat(3,1),outmat(3,2),outmat(3,3)
       write(*,*) ''


!------------------------------------------------------------------------------!
       sqsmat(1,1)=2.0d0
       sqsmat(1,2)=0.0d0
       sqsmat(1,3)=0.0d0
       sqsmat(2,1)=0.0d0
       sqsmat(2,2)=3.0d0
       sqsmat(2,3)=0.0d0
       sqsmat(3,1)=0.0d0
       sqsmat(3,2)=0.0d0
       sqsmat(3,3)=4.0d0

       vector(1)=1
       vector(2)=1
       vector(3)=0
       weight=2.0d0

       outmat(:,:)=0.0d0
       call vvfterm(3,sqsmat,vector,weight,outmat)
       write(*,*) '2*4,2*9,0 id matrix:'
       write(*,*) outmat(1,1),outmat(1,2),outmat(1,3)
       write(*,*) outmat(2,1),outmat(2,2),outmat(2,3)
       write(*,*) outmat(3,1),outmat(3,2),outmat(3,3)
       write(*,*) ''


!------------------------------------------------------------------------------!
       sqsmat(1,1)=1.0d0
       sqsmat(1,2)=1.0d0
       sqsmat(1,3)=1.0d0
       sqsmat(2,1)=1.0d0
       sqsmat(2,2)=1.0d0
       sqsmat(2,3)=1.0d0
       sqsmat(3,1)=1.0d0
       sqsmat(3,2)=1.0d0
       sqsmat(3,3)=1.0d0

       vector(1)=1
       vector(2)=1
       vector(3)=1
       weight=1.0d0

       call vvfterm(3,sqsmat,vector,weight,outmat)
       write(*,*) 'previows + 3 matrix:'
       write(*,*) outmat(1,1),outmat(1,2),outmat(1,3)
       write(*,*) outmat(2,1),outmat(2,2),outmat(2,3)
       write(*,*) outmat(3,1),outmat(3,2),outmat(3,3)
       write(*,*) ''


!------------------------------------------------------------------------------!
       write(*,*) ''
       write(*,*) ''
       return; end program
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
