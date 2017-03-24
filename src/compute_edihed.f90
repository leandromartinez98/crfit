!
! Subroutine compute_edihed: Computes the dihedral angle energy
! 
! L. Martinez, Dez 02, 2013
! Institute of Chemistry - State University of Campinas - Brazil
!

subroutine compute_edihed(x,edihed)
  
  use force_field
  implicit none
  integer :: i, j, k
  double precision :: edihed, x(natoms*3), phi, dihedangle
  
  edihed = 0.d0
  do i = 1, ndihed
    phi = dihedangle(idihed(i,1),idihed(i,2),idihed(i,3),idihed(i,4),x)
    do j = 0, idihed(i,6) - 1
      k = idihed(i,5) + j
      edihed = edihed + fdihed(k,1)*(1+dcos(fdihed(k,2)*phi-fdihed(k,3)))
    end do
  end do

  return
end subroutine compute_edihed

!
! Function that computes a dihedral angle given the indices of the atoms
!

double precision function dihedangle(i,j,k,l,x)

  implicit none
  integer :: i, j, k, l, ix, iy, iz, jx, jy, jz, kx, ky, kz, lx, ly, lz
  double precision :: x(*), v1(3), v2(3), v3(3), c1(3), c2(3), c3(3), &
                      norm2_c1, norm2_c2, norm2_c3, cosphi, sinphi
                      
  ix = (i-1)*3 + 1 
  iy = ix + 1
  iz = ix + 2

  jx = (j-1)*3 + 1
  jy = jx + 1
  jz = jx + 2

  kx = (k-1)*3 + 1
  ky = kx + 1
  kz = kx + 2

  lx = (l-1)*3 + 1
  ly = lx + 1
  lz = lx + 2

  v1(1) = x(jx) - x(ix)
  v1(2) = x(jy) - x(iy)
  v1(3) = x(jz) - x(iz)
  
  v2(1) = x(kx) - x(jx)
  v2(2) = x(ky) - x(jy)
  v2(3) = x(kz) - x(jz)

  v3(1) = x(lx) - x(kx)
  v3(2) = x(ly) - x(ky)
  v3(3) = x(lz) - x(kz)

  c1(1) = v1(2)*v2(3) - v1(3)*v2(2) 
  c1(2) = v1(3)*v2(1) - v1(1)*v2(3)
  c1(3) = v1(1)*v2(2) - v1(2)*v2(1)
  norm2_c1 =  c1(1)**2 + c1(2)**2 + c1(3)**2

  c2(1) = v2(2)*v3(3) - v2(3)*v3(2) 
  c2(2) = v2(3)*v3(1) - v2(1)*v3(3)
  c2(3) = v2(1)*v3(2) - v2(2)*v3(1)
  norm2_c2 =  c2(1)**2 + c2(2)**2 + c2(3)**2

  c3(1) = v2(2)*c1(3) - v2(3)*c1(2) 
  c3(2) = v2(3)*c1(1) - v2(1)*c1(3)
  c3(3) = v2(1)*c1(2) - v2(2)*c1(1)
  norm2_c3 =  c3(1)**2 + c3(2)**2 + c3(3)**2

  cosphi = ( c1(1)*c2(1) + c1(2)*c2(2) + c1(3)*c2(3) ) / &
             dsqrt( norm2_c1 * norm2_c2 )
  sinphi = ( c3(1)*c2(1) + c3(2)*c2(2) + c3(3)*c2(3) ) / &
             dsqrt( norm2_c3 * norm2_c2 )

  dihedangle = atan2(sinphi,cosphi)

end function dihedangle




