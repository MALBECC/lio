!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       include "../lowdin_calc.f"
       program lowdin_calc_test01
       implicit none
       real*8   :: rhomat(2,2),sqsmat(2,2)
       integer  :: atomorb(2)
       real*8   :: outvec(2)

       write(*,*) ''
       write(*,*) ''


!------------------------------------------------------------------------------!
       rhomat(1,1)=1.0d0
       rhomat(1,2)=1.0d0
       rhomat(2,1)=1.0d0
       rhomat(2,2)=1.0d0

       sqsmat(1,1)=0.0d0
       sqsmat(1,2)=0.0d0
       sqsmat(2,1)=0.0d0
       sqsmat(2,2)=0.0d0

       atomorb(1)=1
       atomorb(2)=2

       outvec(:)=0.0d0
       call lowdin_calc(2,2,rhomat,sqsmat,atomorb,outvec)
       write(*,*) 'null vector from null sqs:'
       write(*,*) outvec(1),outvec(2)
       write(*,*) ''


!------------------------------------------------------------------------------!
       rhomat(1,1)=0.0d0
       rhomat(1,2)=0.0d0
       rhomat(2,1)=0.0d0
       rhomat(2,2)=0.0d0

       sqsmat(1,1)=1.0d0
       sqsmat(1,2)=1.0d0
       sqsmat(2,1)=1.0d0
       sqsmat(2,2)=1.0d0

       atomorb(1)=1
       atomorb(2)=2

       outvec(:)=0.0d0
       call lowdin_calc(2,2,rhomat,sqsmat,atomorb,outvec)
       write(*,*) 'null vector from null rho:'
       write(*,*) outvec(1),outvec(2)
       write(*,*) ''


!------------------------------------------------------------------------------!
       rhomat(1,1)=1.0d0
       rhomat(1,2)=0.0d0
       rhomat(2,1)=0.0d0
       rhomat(2,2)=1.0d0

       sqsmat(1,1)=2.0d0
       sqsmat(1,2)=2.0d0
       sqsmat(2,1)=1.0d0
       sqsmat(2,2)=1.0d0

       atomorb(1)=1
       atomorb(2)=2

       outvec(:)=0.0d0
       call lowdin_calc(2,2,rhomat,sqsmat,atomorb,outvec)
       write(*,*) 'expected result is (6,3) :'
       write(*,*) outvec(1),outvec(2)
       write(*,*) ''


!------------------------------------------------------------------------------!
       rhomat(1,1)= 2.0d0
       rhomat(1,2)= 1.0d0
       rhomat(2,1)=-1.0d0
       rhomat(2,2)= 3.0d0

       sqsmat(1,1)=2.0d0
       sqsmat(1,2)=2.0d0
       sqsmat(2,1)=1.0d0
       sqsmat(2,2)=1.0d0

       atomorb(1)=1
       atomorb(2)=2

       outvec(:)=0.0d0
       call lowdin_calc(2,2,rhomat,sqsmat,atomorb,outvec)
       write(*,*) 'expected result is (12,6) :'
       write(*,*) outvec(1),outvec(2)
       write(*,*) ''


!------------------------------------------------------------------------------!
       rhomat(1,1)= 2.0d0
       rhomat(1,2)= 1.0d0
       rhomat(2,1)=-1.0d0
       rhomat(2,2)= 3.0d0

       sqsmat(1,1)=2.0d0
       sqsmat(1,2)=2.0d0
       sqsmat(2,1)=1.0d0
       sqsmat(2,2)=1.0d0

       atomorb(1)=1
       atomorb(2)=1

       outvec(:)=0.0d0
       call lowdin_calc(2,1,rhomat,sqsmat,atomorb,outvec)
       write(*,*) 'expected result is (18,0) :'
       write(*,*) outvec(1),outvec(2)
       write(*,*) ''


!------------------------------------------------------------------------------!
       write(*,*) ''
       write(*,*) ''
       return; end program
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
