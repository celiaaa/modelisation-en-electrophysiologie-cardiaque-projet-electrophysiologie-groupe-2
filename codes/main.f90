program main
  use donnees
  use stockage

  implicit none

  character(25) :: output = 'output.dat'
  
  !type(matrice) :: A
  real*8, dimension(:), allocatable :: u , unext , uprev
  integer :: i,j,k,time_to_save,n,nb_iter
  real*8 :: dx,dt,tmax,xmin,xmax,a

  open(11,file='parametres.dat',action='read',status='old')

  read(11,*) xmin
  read(11,*) xmax
  read(11,*) n
  read(11,*) 
  read(11,*) tmax
  read(11,*) time_to_save

  close(11)
  
  dx = (xmax-xmin)/n
  dt = dx*dx/2.
  nb_iter = ceiling(tmax/dt)
  print*, 'dx / dt', dx , '   /   ', dt
  print*, 'T =',tmax, ' / ', nb_iter , ' iterations'
  allocate(uprev(0:N),unext(0:N),u(0:N))
  a = 2*dt/dx/dx

  open(11,file=trim(output),action='write',status='unknown')
  
  uprev = 0.
  do i = ceiling(N/2.)-10, ceiling(N/2.)+10,1
     uprev(i) = 100.
  enddo

!  write(11,'(3E22.15)') 0. , xmin , uprev(0)
  do i = 0,N
     write(11,'(3E22.15)') 0., real(i*dx) , uprev(i)
  enddo
  
  write(11,*)
  
  u(0) = uprev(0) + a*(uprev(1)-uprev(0))
  u(n) = uprev(n) + a*(uprev(N)-uprev(N-1))
  do i = 1,N-1
     u(i) = uprev(i) + (dt/dx/dx)*(uprev(i+1)-2*uprev(i)+uprev(i-1))
  enddo

  ! do i = 0,N
  !    write(11,'(2E22.15)')  real(i*dx) , u(i)
  ! enddo
  ! write(11,'(3E22.15)') dt, xmax , u(N)
  ! write(11,*)
  
  do k = 1, nb_iter-1

     unext(0) = uprev(0) + 2*a*(u(1)-u(0))
     unext(N) = uprev(N) + 2*a*(u(N-1)-u(N))
     do i = 1,N-1

        unext(i) = uprev(i)*(1-a)/(1+a) + (a/(1+a))*(u(i-1)+u(i+1))

     enddo
     
     if (mod(k,time_to_save) == 0) then

        !write(11,*)  real(k*dt) , xmin , unext(0)
        do i = 0,N
           write(11,'(3E22.15)') real(k*dt) , real(i*dx) , unext(i)
        enddo
        
        write(11,*)
        
     endif

  end do

  close(11)
end program main
