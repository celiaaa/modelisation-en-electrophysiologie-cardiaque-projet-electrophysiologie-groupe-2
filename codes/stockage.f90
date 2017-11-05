module stockage
  implicit none

  abstract interface
     function rhs_of_ode(x)
       use donnees
       real*8,intent(in) :: x
       real*8 :: rhs_of_ode
     end function rhs_of_ode
  end interface

  type matrice
     real*8,dimension(:,:),allocatable :: data
     integer,dimension(:),allocatable :: dist
  end type matrice

  real*8,parameter :: pi=3.14159265359

contains
  
  subroutine affiche(mat)
    type(matrice),intent(in) :: mat 
    integer :: i,j,k,m,nbdi

    m = size(mat%data(:,1)) !Nombre de lignes                                   
    nbdi = size(mat%data(1,:)) !Nombres de diagonales                           

    print "('Matrice de taille',I4)", m

    do i = 1,m
       do k = 1,nbdi
          j = i+mat%dist(k)
          print*, 'A(' , i,j,') = ',mat%data(i,k)
       enddo
    enddo
  end subroutine affiche

  function produit(A,x) !!Produit Ax = b                                                                                                            
    type(matrice),intent(in) :: A 
    real*8,dimension(:),intent(in) :: x 
    real*8,dimension(:),allocatable :: produit

    integer :: i,j,m,nbdi,k
    real*8,dimension(:),allocatable :: x_etendu

    m = size(A%data(:,1)) !Nombre de lignes                                                                                                           
    nbdi = size(A%data(1,:)) !Nombre de diago non nulles                                                                                              

    allocate(produit(m))

    !V‚rification des tailles                                                                                                                         
    if (m .ne. size(x) .or. m  .ne. size(produit)) then
       print*,"Les tailles ne sont pas compatibles"
       !   call exit(1)                                                                                                                                  
    endif

    allocate(x_etendu(A%dist(1)+1:m+A%dist(nbdi)))
    x_etendu = 0
    x_etendu(1:m) = x


    produit=0
    do i = 1,m
       do k = 1,nbdi
          produit(i) = produit(i) + A%data(i,k)*x_etendu(i+A%dist(k))
       enddo
    enddo

    deallocate(x_etendu)
  end function produit

  subroutine gradconj(A,b,x)
    implicit none
    type(matrice),intent(in) :: a
    type(matrice) :: C,D !Pour préconditionnement                                                                                                     
    real*8, dimension(:),intent(in) :: b
    real*8,dimension(:),intent(inout) :: x
    real*8,dimension(:), allocatable :: y,y1,h0,h,g,g0,x0
    real :: p, gamma
    integer :: k

    allocate(x0(size(b)),g0(size(b)),h0(size(b)),g(size(b)))
    allocate(h(size(b)),y1(size(b)),y(size(b)))

    x0 = 0.

    allocate(c%dist(1),c%data(size(b),1))
    allocate(d%dist(1),d%data(size(b),1))

    c%data(:,1) = 1./a%data(:,3)
    c%dist(1) = 0. ; d%dist(1) = 0.

    d%data(:,1) = -c%data(:,1)

    g0 = -b
    h0 = produit(d,g0)

    do k=0,size(b)
       y=produit(a,h0)

       p = -(dot_product(g0,h0)/dot_product(h0,y))

       x=x+p*h0

       g=g0+p*y
       y1 = produit(c,g) ; y=produit(c,g0)

       gamma=dot_product(y1,g)/dot_product(y,g0)
       y=produit(d,g)

       h=y+gamma*h0

       if (dot_product(y1,g) .lt. 0.0000001) then
          EXIT
       endif

       h0=h ; g0=g
    enddo

  end subroutine gradconj

  subroutine init_A(N,dt,dx,alpha,A)
    procedure(rhs_of_ode) :: alpha
    integer, intent(in) :: N
    real*8, intent(in) :: dt, dx
    type(matrice), intent(out) :: A

    integer :: i

    allocate(A%dist(3)) ; allocate(A%data(0:N,3))

    A%dist(:) = (/-1,0,1/)
    print*, 'test1'
    A%data = 0.
    print*, 'test2'
    do i = 1,N-1
       A%data(i,2) = 2*alpha((i+0.5)*dx) + 2*alpha((i-0.5)*dx)
       A%data(i,1) = -alpha((i-0.5)*dx)
       A%data(i,3) = -alpha((i+0.5)*dx)
    enddo
    print*, 'test3'
    i = 0
    A%data(i,1) = 0.
    A%data(i,2) = 2*alpha(0.5*dx)
    A%data(i,3) = -2*alpha(0.5*dx)

    i = N
    A%data(i,1) = -2*alpha((N-0.5)*dx)
    A%data(i,2) = 2*alpha((N-0.5)*dx)
    A%data(i,3) = 0.

  end subroutine init_A

end module stockage
