!
! Subroutine compute_gangle: Computes the gradient relative to angles
! 
! L. Martinez, Nov 28, 2013
! Institute of Chemistry - State University of Campinas - Brazil
!

subroutine compute_gangle(x,g)
  
  use force_field
  implicit none
  integer :: i, ix, iy, iz, jx, jy, jz, kx, ky, kz
  double precision :: g(natoms*3), x(natoms*3), theta, costheta, dot_product, &
                      v1(3), v2(3), v1norm_inv, v2norm_inv, v1norm2, v2norm2, &
                      dedtheta, dthetadcostheta, a, b, d, &
                      v1v2norm_inv, dd, ddv1, ddv2, dotv1, dotv2, cv1, cv2

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

    v1norm2 = v1(1)**2 + v1(2)**2 + v1(3)**2
    v2norm2 = v2(1)**2 + v2(2)**2 + v2(3)**2
    v1norm_inv = 1.d0 / dsqrt(v1norm2)
    v2norm_inv = 1.d0 / dsqrt(v2norm2)
    v1v2norm_inv = v1norm_inv * v2norm_inv

    dot_product = v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)
    costheta = dot_product * v1norm_inv * v2norm_inv
    if ( costheta > 1.d0 ) then 
      costheta = 1.d0
    else if ( costheta < -1.d0 ) then
      costheta = -1.d0
    end if
    theta = dacos(costheta)

    dedtheta = 2.d0*fangle(i,1)*(theta-fangle(i,2))
    dthetadcostheta = -1.d0 / ( dsqrt( 1.d0 - costheta**2 ) )
    dd = dedtheta * dthetadcostheta
    ddv1 = dd*v1norm_inv
    ddv2 = dd*v2norm_inv
    dotv1 = dot_product * v1norm_inv**3
    dotv2 = dot_product * v2norm_inv**3
    cv1 = costheta / v1norm2
    cv2 = costheta / v2norm2

    g(ix) = g(ix) + ddv2 * ( v2(1)*v1norm_inv - dotv1 * v1(1) )
    g(iy) = g(iy) + ddv2 * ( v2(2)*v1norm_inv - dotv1 * v1(2) )
    g(iz) = g(iz) + ddv2 * ( v2(3)*v1norm_inv - dotv1 * v1(3) )

    g(kx) = g(kx) + ddv1 * ( v1(1)*v2norm_inv - dotv2 * v2(1) )
    g(ky) = g(ky) + ddv1 * ( v1(2)*v2norm_inv - dotv2 * v2(2) )
    g(kz) = g(kz) + ddv1 * ( v1(3)*v2norm_inv - dotv2 * v2(3) )

    g(jx) = g(jx) + dd * ( v1(1)*cv1 + v2(1)*cv2 - (v2(1)+v1(1))*v1v2norm_inv )
    g(jy) = g(jy) + dd * ( v1(2)*cv1 + v2(2)*cv2 - (v2(2)+v1(2))*v1v2norm_inv )
    g(jz) = g(jz) + dd * ( v1(3)*cv1 + v2(3)*cv2 - (v2(3)+v1(3))*v1v2norm_inv )

  end do

  ! Urey-Bradley terms

  do i = 1, nureybradley

    ix = (iureybradley(i,1)-1)*3 + 1
    iy = ix + 1
    iz = ix + 2

    kx = (iureybradley(i,2)-1)*3 + 1
    ky = kx + 1
    kz = kx + 2

    v1(1) = x(ix) - x(kx)
    v1(2) = x(iy) - x(ky)
    v1(3) = x(iz) - x(kz)

    d = dsqrt(v1(1)**2 + v1(2)**2 + v1(3)**2)

    a = 2.d0*fureybradley(i,1)*(d - fureybradley(i,2))/d

    b = a*v1(1)
    g(ix) = g(ix) + b
    g(kx) = g(kx) - b

    b = a*v1(2)
    g(iy) = g(iy) + b
    g(ky) = g(ky) - b

    b = a*v1(3)
    g(iz) = g(iz) + b
    g(kz) = g(kz) - b

  end do

  return
end subroutine compute_gangle

