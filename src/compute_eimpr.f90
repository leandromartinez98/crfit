!
! Subroutine compute_eimpr: Computes the improper angle energy
! 
! L. Martinez, Dez 02, 2013
! Institute of Chemistry - State University of Campinas - Brazil
!

subroutine compute_eimpr(x,eimpr)
  
  use force_field
  implicit none
  integer :: i, j, k, ix, iy, iz, jx, jy, jz, kx, ky, kz, lx, ly, lz
  double precision :: eimpr, x(*), phi, cosphi, &
                      v1(3), v2(3), v3(3), c1(3), c2(3), norm2_c1, norm2_c2
  
  eimpr = 0.d0
  do i = 1, nimpr

    ix = (iimpr(i,1)-1)*3 + 1
    iy = ix + 1
    iz = ix + 2

    jx = (iimpr(i,2)-1)*3 + 1
    jy = jx + 1  
    jz = jx + 2

    kx = (iimpr(i,3)-1)*3 + 1
    ky = kx + 1
    kz = kx + 2

    lx = (iimpr(i,4)-1)*3 + 1
    ly = lx + 1
    lz = lx + 2

    v1(1) = x(ix) - x(jx)
    v1(2) = x(iy) - x(jy)
    v1(3) = x(iz) - x(jz)
  
    v2(1) = x(jx) - x(kx)
    v2(2) = x(jy) - x(ky)
    v2(3) = x(jz) - x(kz)

    v3(1) = x(kx) - x(lx)
    v3(2) = x(ky) - x(ly)
    v3(3) = x(kz) - x(lz)
    
    c1(1) = v1(2)*v2(3) - v1(3)*v2(2) 
    c1(2) = v1(3)*v2(1) - v1(1)*v2(3)
    c1(3) = v1(1)*v2(2) - v1(2)*v2(1)
    norm2_c1 =  c1(1)**2 + c1(2)**2 + c1(3)**2

    c2(1) = v2(2)*v3(3) - v2(3)*v3(2) 
    c2(2) = v2(3)*v3(1) - v2(1)*v3(3)
    c2(3) = v2(1)*v3(2) - v2(2)*v3(1)
    norm2_c2 =  c2(1)**2 + c2(2)**2 + c2(3)**2

    cosphi = ( c1(1)*c2(1) + c1(2)*c2(2) + c1(3)*c2(3) ) / &
               dsqrt( norm2_c1 * norm2_c2 )
                  
    if ( cosphi > 0.d0 ) cosphi = dmin1(cosphi,1.d0-1.d-10)
    if ( cosphi < 0.d0 ) cosphi = dmax1(cosphi,-1.d0+1.d-10)
    phi = dacos(cosphi)

    do j = 0, iimpr(i,6) - 1
      k = iimpr(i,5) + j
      eimpr = eimpr + fimpr(k,1)*(phi - fimpr(k,2))**2
    end do

  end do

  return
end subroutine compute_eimpr






