!
! Subroutine compute_gbond: Compute the gradient of the potential relative
!            to bonds.
! 
! L. Martinez, Nov 03, 2013
! Institute of Chemistry - State University of Campinas - Brazil
!

subroutine compute_gbond(x,g)
  
  use force_field
  implicit none
  integer :: i, ix, iy, iz, jx, jy, jz
  double precision :: x(natoms*3), d, g(natoms*3), gx, gy, gz, gc, dx, dy, dz

  do i = 1, nbonds

    ix = (ibond(i,1)-1)*3 + 1
    iy = ix + 1
    iz = ix + 2

    jx = (ibond(i,2)-1)*3 + 1
    jy = jx + 1
    jz = jx + 2

    dx = x(ix) - x(jx)
    dy = x(iy) - x(jy)
    dz = x(iz) - x(jz)

    d = dsqrt( dx**2 + dy**2 + dz**2 )

    gc = 2.d0*fbond(i,1)*(d-fbond(i,2)) / d

    gx = gc*dx
    gy = gc*dy
    gz = gc*dz

    g(ix) = g(ix) + gx
    g(iy) = g(iy) + gy
    g(iz) = g(iz) + gz

    g(jx) = g(jx) - gx
    g(jy) = g(jy) - gy
    g(jz) = g(jz) - gz

  end do

  return
end subroutine compute_gbond


