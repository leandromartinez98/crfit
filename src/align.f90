!
! Subroutine align: aligns to structures
!
! Returns x aligned to y
!
! Requires subroutine dsyev from Lapack (compile with -llapack)
!
! L. Martinez, Institute of Chemistry - University of Campinas
! April 09, 2015
!

subroutine align(n,x,y)

  implicit none 
  integer :: n, i, j, k
  double precision :: x(n,3), y(n,3), cmx(3), cmy(3) , q(4,4), u(3,3), a(4)
  double precision, allocatable :: xm(:), ym(:), zm(:), xp(:), yp(:), zp(:),&
                                   xnew(:,:)

  ! For dsyev
  double precision :: work(12)
  integer :: info

  ! Test if xm and xp are already allocated, if not, allocate them

  if ( .not. allocated(xm) ) then 
    allocate(xm(n), ym(n), zm(n), xp(n), yp(n), zp(n), xnew(n,3))
  else if ( size(xm) /= n ) then
    deallocate(xm, ym, zm, xp, yp, zp, xnew)
    allocate(xm(n), ym(n), zm(n), xp(n), yp(n), zp(n), xnew(n,3))
  end if

  ! Computing centroid

  do i = 1, 3
    cmx(i) = 0.d0
    cmy(i) = 0.d0
  end do

  do i = 1, n
    do j = 1, 3
      cmx(j) = cmx(j) + x(i,j)
      cmy(j) = cmy(j) + y(i,j)
    end do
  end do
      
  do i = 1, 3
    cmx(i) = cmx(i) / dble(n)
    cmy(i) = cmy(i) / dble(n)
  end do

  ! Translating both sets to the origin

  do i = 1, n
    do j = 1, 3
      x(i,j) = x(i,j) - cmx(j)
      y(i,j) = y(i,j) - cmy(j)
    end do
  end do

 ! Computing the quaternion matrix

  do i = 1, n
    xm(i) = y(i,1) - x(i,1)
    ym(i) = y(i,2) - x(i,2)
    zm(i) = y(i,3) - x(i,3)
    xp(i) = y(i,1) + x(i,1)
    yp(i) = y(i,2) + x(i,2)
    zp(i) = y(i,3) + x(i,3)
  end do

  do i = 1, 4
    do j = 1, 4
      q(i,j) = 0.
    end do
  end do

  do i = 1, n
    q(1,1) = q(1,1) + xm(i)**2 + ym(i)**2 + zm(i)**2
    q(1,2) = q(1,2) + yp(i)*zm(i) - ym(i)*zp(i)
    q(1,3) = q(1,3) + xm(i)*zp(i) - xp(i)*zm(i)
    q(1,4) = q(1,4) + xp(i)*ym(i) - xm(i)*yp(i)
    q(2,2) = q(2,2) + yp(i)**2 + zp(i)**2 + xm(i)**2
    q(2,3) = q(2,3) + xm(i)*ym(i) - xp(i)*yp(i)
    q(2,4) = q(2,4) + xm(i)*zm(i) - xp(i)*zp(i)
    q(3,3) = q(3,3) + xp(i)**2 + zp(i)**2 + ym(i)**2
    q(3,4) = q(3,4) + ym(i)*zm(i) - yp(i)*zp(i)
    q(4,4) = q(4,4) + xp(i)**2 + yp(i)**2 + zm(i)**2
  end do
  q(2,1) = q(1,2)
  q(3,1) = q(1,3)
  q(3,2) = q(2,3)
  q(4,1) = q(1,4)
  q(4,2) = q(2,4)
  q(4,3) = q(3,4)          

  ! Computing the eigenvectors 'a' and eigenvalues 'q' of the q matrix
  ! 'q' contains the eigenvectors associates with eigenvalues in ascending order

  call dsyev('V','U',4,q,4,a,work,12,info)

  ! Compute rotation matrix
  
  u(1,1) = q(1,1)**2 + q(2,1)**2 - q(3,1)**2 - q(4,1)**2
  u(1,2) = 2. * ( q(2,1)*q(3,1) + q(1,1)*q(4,1) )
  u(1,3) = 2. * ( q(2,1)*q(4,1) - q(1,1)*q(3,1) )
  u(2,1) = 2. * ( q(2,1)*q(3,1) - q(1,1)*q(4,1) )
  u(2,2) = q(1,1)**2 + q(3,1)**2 - q(2,1)**2 - q(4,1)**2
  u(2,3) = 2. * ( q(3,1)*q(4,1) + q(1,1)*q(2,1) )
  u(3,1) = 2. * ( q(2,1)*q(4,1) + q(1,1)*q(3,1) )
  u(3,2) = 2. * ( q(3,1)*q(4,1) - q(1,1)*q(2,1) )
  u(3,3) = q(1,1)**2 + q(4,1)**2 - q(2,1)**2 - q(3,1)**2      

  ! Rotate vector x (will be stored in xnew)

  do i = 1, n
    do j = 1, 3
      xnew(i,j) = 0.d0
      do k = 1, 3
        xnew(i,j) = xnew(i,j) + u(j,k) * x(i,k)
      end do
    end do
  end do

  ! Translate vector to the centroid of y (and restore y)

  do i = 1, n
    do j = 1, 3
      x(i,j) = xnew(i,j) + cmy(j)
      y(i,j) = y(i,j) + cmy(j)
    end do
  end do

  return 
end subroutine align
