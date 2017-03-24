!
! Subroutine compute_gimpr: Computes the improper dihedral angle energy gradient
! 
! L. Martinez, Dez 03, 2013
! Institute of Chemistry - State University of Campinas - Brazil
!

subroutine compute_gimpr(x,g)
  
  use force_field
  implicit none
  integer :: i, j, k, ix, iy, iz, jx, jy, jz, kx, ky, kz, lx, ly, lz
  double precision :: g(natoms*3), x(natoms*3), phi, dot_product, norm2_c1, norm2_c2,&
                      v1(3), v2(3), v3(3), c1(3), c2(3), &
                      de_dphi, dphi_dcosphi, dout, pnorm_inv, cosphi, dotn1, dotn2
  double precision :: c11v12, c11v13, c11v22, c11v23, c11v32, c11v33, &
                      c12v11, c12v13, c12v21, c12v23, c12v31, c12v33, &
                      c13v11, c13v12, c13v21, c13v22, c13v31, c13v32, & 
                      c21v12, c21v13, c21v22, c21v23, c21v32, c21v33, &
                      c22v11, c22v13, c22v21, c22v23, c22v31, c22v33, &
                      c23v11, c23v12, c23v21, c23v22, c23v31, c23v32

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
    pnorm_inv = 1.d0 / dsqrt(norm2_c1*norm2_c2)

    dot_product = ( c1(1)*c2(1) + c1(2)*c2(2) + c1(3)*c2(3) ) 
    cosphi = dot_product * pnorm_inv

    if ( cosphi > 0.d0 ) cosphi = dmin1(cosphi,1.d0-1.d-10)
    if ( cosphi < 0.d0 ) cosphi = dmax1(cosphi,-1.d0+1.d-10)
    phi = dacos(cosphi)
    dphi_dcosphi = -1.d0*pnorm_inv / dsqrt( 1.d0 - cosphi**2 )
    dotn1 = dot_product/norm2_c1
    dotn2 = dot_product/norm2_c2

    c11v12 = c1(1)*v1(2); c11v13 = c1(1)*v1(3); c11v22 = c1(1)*v2(2)
    c11v23 = c1(1)*v2(3); c11v32 = c1(1)*v3(2); c11v33 = c1(1)*v3(3)
    c12v11 = c1(2)*v1(1); c12v13 = c1(2)*v1(3); c12v21 = c1(2)*v2(1)
    c12v23 = c1(2)*v2(3); c12v31 = c1(2)*v3(1); c12v33 = c1(2)*v3(3)
    c13v11 = c1(3)*v1(1); c13v12 = c1(3)*v1(2); c13v21 = c1(3)*v2(1)
    c13v22 = c1(3)*v2(2); c13v31 = c1(3)*v3(1); c13v32 = c1(3)*v3(2)
    c21v12 = c2(1)*v1(2); c21v13 = c2(1)*v1(3); c21v22 = c2(1)*v2(2)
    c21v23 = c2(1)*v2(3); c21v32 = c2(1)*v3(2); c21v33 = c2(1)*v3(3)
    c22v11 = c2(2)*v1(1); c22v13 = c2(2)*v1(3); c22v21 = c2(2)*v2(1)
    c22v23 = c2(2)*v2(3); c22v31 = c2(2)*v3(1); c22v33 = c2(2)*v3(3)
    c23v11 = c2(3)*v1(1); c23v12 = c2(3)*v1(2); c23v21 = c2(3)*v2(1)
    c23v22 = c2(3)*v2(2); c23v31 = c2(3)*v3(1); c23v32 = c2(3)*v3(2)

    do j = 0, iimpr(i,6) - 1
      k = iimpr(i,5) + j

      de_dphi = 2.d0*fimpr(k,1)*(phi-fimpr(k,2))
      dout = de_dphi*dphi_dcosphi

      g(ix) = g(ix) + dout*(c23v22-c22v23 + dotn1*(c12v23-c13v22))
      g(iy) = g(iy) + dout*(c21v23-c23v21 + dotn1*(c13v21-c11v23))
      g(iz) = g(iz) + dout*(c22v21-c21v22 + dotn1*(c11v22-c12v21))

      g(jx) = g(jx) + dout*(c22v13+c22v23+c13v32-c12v33-c23v22-c23v12 + &
                   dotn1*(c13v22+c13v12-c12v13-c12v23) + dotn2*(c22v33-c23v32))
      g(jy) = g(jy) + dout*(c23v11+c23v21+c11v33-c13v31-c21v23-c21v13 + &  
                   dotn1*(c11v23+c11v13-c13v11-c13v21) + dotn2*(c23v31-c21v33))
      g(jz) = g(jz) + dout*(c21v12+c21v22+c12v31-c11v32-c22v21-c22v11 + &  
                   dotn1*(c12v21+c12v11-c11v12-c11v22) + dotn2*(c21v32-c22v31))

      g(kx) = g(kx) + dout*(c12v23+c12v33+c23v12-c22v13-c13v32-c13v22 + &
                   dotn2*(c23v32+c23v22-c22v23-c22v33) + dotn1*(c12v13-c13v12))
      g(ky) = g(ky) + dout*(c13v21+c13v31+c21v13-c23v11-c11v33-c11v23 + &
                   dotn2*(c21v33+c21v23-c23v21-c23v31) + dotn1*(c13v11-c11v13))
      g(kz) = g(kz) + dout*(c11v22+c11v32+c22v11-c21v12-c12v31-c12v21 + &
                   dotn2*(c22v31+c22v21-c21v22-c21v32) + dotn1*(c11v12-c12v11))

      g(lx) = g(lx) + dout*(c13v22-c12v23 + dotn2*(c22v23-c23v22))
      g(ly) = g(ly) + dout*(c11v23-c13v21 + dotn2*(c23v21-c21v23))
      g(lz) = g(lz) + dout*(c12v21-c11v22 + dotn2*(c21v22-c22v21))

    end do
  end do

  return
end subroutine compute_gimpr

