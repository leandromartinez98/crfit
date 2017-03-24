!
! Subroutine compute_gnb: Computes the gradient relative to non-bonded 
!                         interactions
! 
! L. Martinez, Dez 4, 2013
! Institute of Chemistry - State University of Campinas - Brazil
!

subroutine compute_gnb(x,g)
  
  use force_field
  implicit none
  integer :: inb, kbb, ix, iy, iz, jx, jy, jz
  double precision :: d2, d, c, dinv6, &
                      d3, dx, dy, dz, qqd3, qqd3x, qqd3y, qqd3z, &
                      g(natoms*3), x(natoms*3)

  ! Gradient for non-bonded interactions for all atoms

  if ( ffcomp(7) ) then

    do inb = 1, n_nonbonded

      ix = (ijnonbonded(inb,1)-1)*3 + 1
      iy = ix + 1
      iz = ix + 2

      jx = (ijnonbonded(inb,2)-1)*3 + 1
      jy = jx + 1
      jz = jx + 2

      dx = x(ix) - x(jx)
      dy = x(iy) - x(jy)
      dz = x(iz) - x(jz)

      d2 = dx**2 + dy**2 + dz**2

      if ( d2 > cutoff2 ) cycle

      d = dsqrt(d2)
      
      d3 = d**3
      qqd3 = qq(inb)/d3
      qqd3x = dx*qqd3
      qqd3y = dy*qqd3
      qqd3z = dz*qqd3

      g(ix) = g(ix) - qqd3x
      g(iy) = g(iy) - qqd3y
      g(iz) = g(iz) - qqd3z
      g(jx) = g(jx) + qqd3x
      g(jy) = g(jy) + qqd3y
      g(jz) = g(jz) + qqd3z

      dinv6 = (1.d0/d3)**2
      c = -12.d0*epseps(inb)*ss6(inb)*dinv6/d2*( ss6(inb)*dinv6 - 1.d0 )

      g(ix) = g(ix) + c*dx
      g(iy) = g(iy) + c*dy
      g(iz) = g(iz) + c*dz
      g(jx) = g(jx) - c*dx
      g(jy) = g(jy) - c*dy
      g(jz) = g(jz) - c*dz

    end do
  
   ! Gradient for backbone non-bonded interactions only

   else

     do kbb = 1, n_nbbackbone

      inb = inbbackbone(kbb)

      ix = (ijnonbonded(inb,1)-1)*3 + 1
      iy = ix + 1
      iz = ix + 2

      jx = (ijnonbonded(inb,2)-1)*3 + 1
      jy = jx + 1
      jz = jx + 2

      dx = x(ix) - x(jx)
      dy = x(iy) - x(jy)
      dz = x(iz) - x(jz)

      d2 = dx**2 + dy**2 + dz**2

      if ( d2 > cutoff2 ) cycle

      d = dsqrt(d2)
      
      d3 = d**3
      qqd3 = qq(inb)/d3
      qqd3x = dx*qqd3
      qqd3y = dy*qqd3
      qqd3z = dz*qqd3

      g(ix) = g(ix) - qqd3x
      g(iy) = g(iy) - qqd3y
      g(iz) = g(iz) - qqd3z
      g(jx) = g(jx) + qqd3x
      g(jy) = g(jy) + qqd3y
      g(jz) = g(jz) + qqd3z

      dinv6 = (1.d0/d3)**2
      c = -12.d0*epseps(inb)*ss6(inb)*dinv6/d2*( ss6(inb)*dinv6 - 1.d0 )

      g(ix) = g(ix) + c*dx
      g(iy) = g(iy) + c*dy
      g(iz) = g(iz) + c*dz
      g(jx) = g(jx) - c*dx
      g(jy) = g(jy) - c*dy
      g(jz) = g(jz) - c*dz

    end do

   end if

  return
end subroutine compute_gnb














