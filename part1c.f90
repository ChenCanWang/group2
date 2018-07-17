      program part1c

      implicit none
      integer i,j,ntmp
      real*8 H(6,6),eigenv(6),work(18)
      real*8 g,exact

      !g=-1.
      write(*,*)' Enter interaction strength'
      read(*,*),g
      write(*,"('g=',f8.3)")g

      do i=1,6
        do j=1,6
           if(i==j)then
             if(i<=3)then
               H(i,j)=2.*i-g
             else
               H(i,j)=2.*(i-1)-g
             endif
           else
             H(i,j)=-g/2.
           endif
           H(1,6)=0.; H(2,5)=0.; H(3,4)=0.
           H(6,1)=0.; H(5,2)=0.; H(4,3)=0.
        enddo
      enddo
      
      write(6,*)'Matrix'
      do i=1,6
        do j=1,6 
          write(6,'(1x,f8.3\)') H(i,j)
        enddo
        write(6,*)' '
      enddo

      call dsyev('V','U',6,H,6,eigenv,work,18,ntmp)
      write(6,*)'eigenvalues'
      write(6,*)eigenv
      
      exact=minval(eigenv)-(2.-g)
      write(6,*)'Correlation energy'
      write(6,'(1x,f8.3)')exact

      end
           
             
  
