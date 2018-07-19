      program configuration


      implicit none
      integer :: i,j,k
      integer :: np,nl,pair,numcon
      real*8 :: factorial
      integer,allocatable,dimension(:,:) :: states
      integer,allocatable,dimension(:) :: sp

      write(*,*)' Enter number of particles'
      read(*,*),np
      !write(*,"('g=',f8.3)")g
      
      write(*,*)' Enter number of lines'
      read(*,*),nl 
      
! the number of configurations.
      numcon=int(factorial(nl)/factorial(np)/factorial(nl-np))

      write(6,"('number of configurations = 'I4)")numcon

      allocate(sp(np),states(numcon,np))    
      sp=0; states=0
      
      do i=1,np
         sp(i)=i
      enddo
      
      states(1,:)=sp

      k=1
     
10    k=k+1
      i=np
      if(sp(i)==nl)then
        goto 11
      else
        goto 12
      endif
      
11    if((i-1)<1)then
        goto 13
      else
        goto 14
      endif

12    sp(i)=sp(i)+1
      do j=1,np-i
         sp(i+j)=sp(i)+j
      enddo
      states(k,:)=sp
      goto 10

14    i=i-1
      if(sp(i)+1==sp(i+1))then
        goto 11
      else
        goto 12
      endif

13    continue
     
      do k=1,numcon
         write(6,*),states(k,:)
      enddo
       

      end
       
      recursive function factorial(n) result(res)
      implicit none
      integer :: n
      real*8 :: res

      if (n==0 .or. n==1) then
         res=1.0
      else
         res=n*factorial(n-1)
      end if
      
      end function factorial
