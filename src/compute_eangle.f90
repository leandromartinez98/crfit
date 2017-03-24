!
! Subroutine compute_eangle: Computes the angle energy
! 
! L. Martinez, Nov 28, 2013
! Institute of Chemistry - State University of Campinas - Brazil
!

subroutine compute_eangle(x,eangle)
  
  use force_field
  implicit none
  integer :: i, ix, iy, iz, jx, jy, jz, kx, ky, kz
  double precision :: eangle, x(*), theta, costheta, &
                      v1(3), v2(3), v1norm, v2norm, d
  
  eangle = 0.d0
  do i = 1, nangles

    ix = (iangle(i,1)-1)*3 + 1
    iy = ix + 1
    iz = ix + 2

    jx = (iangle(i,2)-1)*3 + 1
    jy = jx + 1
    jz = jx + 2
 
    kx = (iangle(i,3)-1)*3 + 1
    ky = kx + 1
    kz = kx + 2

    v1(1) = x(ix) - x(jx)
    v1(2) = x(iy) - x(jy)
    v1(3) = x(iz) - x(jz)
    v2(1) = x(kx) - x(jx)
    v2(2) = x(ky) - x(jy)
    v2(3) = x(kz) - x(jz)

    v1norm = v1(1)**2 + v1(2)**2 + v1(3)**2 
    v2norm = v2(1)**2 + v2(2)**2 + v2(3)**2 
    costheta = (v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)) / dsqrt(v1norm*v2norm)
    if ( costheta > 1.d0 ) then 
      costheta = 1.d0
    else if ( costheta < -1.d0 ) then
      costheta = -1.d0
    end if

    theta = dacos(costheta)
    eangle = eangle + fangle(i,1)*(theta-fangle(i,2) )**2

  end do

  ! Urey-Bradley therms

  do i = 1, nureybradley

    ix = (iureybradley(i,1)-1)*3 + 1
    iy = ix + 1
    iz = ix + 2

    kx = (iureybradley(i,2)-1)*3 + 1
    ky = kx + 1
    kz = kx + 2

    d = dsqrt((x(ix) - x(kx))**2 + (x(iy) - x(ky))**2 + (x(iz) - x(kz))**2)

    eangle = eangle + fureybradley(i,1)*(d-fureybradley(i,2))**2

  end do

  return
end subroutine compute_eangle

