!
! Subroutine compute_gnbcell: Computes the non-bonded gradient with short-range modification,
!                             and a cutoff, using the linked-cell method to avoid the evaluation
!                             of distances greater than the cutoff
! 
! L. Martinez, Aug 18, 2014
! Institute of Chemistry - State University of Campinas - Brazil
!

subroutine compute_gnbcell(x,g)
  
  use force_field
  use linkedcells
  implicit none
  integer :: i, j, k, iatom
  double precision :: x(natoms*3), g(natoms*3)

  ! Initialize the linked cells with the atoms within

  call cells(x)

  ! Compute the gradient

  do i = 1, nboxes(1)
    do j = 1, nboxes(2)
      do k = 1, nboxes(3)

        iatom = iatomfirst(i,j,k)
        do while( iatom /= 0 ) 

          ! If backbone only computations are being performed, skip this
          ! atom if it is not a backbone atom

          if ( .not. ffcomp(7) .and. .not. isbackbone(iatom) ) then 
            iatom = iatomnext(iatom)
            cycle
          end if

          ! Interactions inside box

          call gnbcell(iatom,iatomnext(iatom),x,g)

          ! Interactions of boxes that share faces

          call gnbcell(iatom,iatomfirst(i+1,j,k),x,g)
          call gnbcell(iatom,iatomfirst(i,j+1,k),x,g)
          call gnbcell(iatom,iatomfirst(i,j,k+1),x,g)

          ! Interactions of boxes that share axes

          call gnbcell(iatom,iatomfirst(i+1,j+1,k),x,g)
          call gnbcell(iatom,iatomfirst(i+1,j,k+1),x,g)
          call gnbcell(iatom,iatomfirst(i+1,j-1,k),x,g)
          call gnbcell(iatom,iatomfirst(i+1,j,k-1),x,g)
          call gnbcell(iatom,iatomfirst(i,j+1,k+1),x,g)
          call gnbcell(iatom,iatomfirst(i,j+1,k-1),x,g)

          ! Interactions of boxes that share vertices

          call gnbcell(iatom,iatomfirst(i+1,j+1,k+1),x,g)
          call gnbcell(iatom,iatomfirst(i+1,j+1,k-1),x,g)
          call gnbcell(iatom,iatomfirst(i+1,j-1,k+1),x,g)
          call gnbcell(iatom,iatomfirst(i+1,j-1,k-1),x,g)

          iatom = iatomnext(iatom)
        end do

      end do
    end do
  end do

  return
end subroutine compute_gnbcell

!
! Function that computes the interactions of atom iatom with all
! atoms jatom of the box with first atom jatom entered
!

subroutine gnbcell(iatom,jatom,x,g)
  
  use force_field, only : natoms, ffcomp, isbackbone
  use linkedcells, only : iatomnext
  implicit none
  integer :: iatom, jatom, j
  double precision :: x(natoms*3), g(natoms*3)

  j = jatom
  do while( j /= 0 ) 

    if ( .not. ffcomp(7) .and. .not. isbackbone(j) ) then
      j = iatomnext(j)
      cycle
    end if

    call gnbpair(iatom,j,x,g)
    j = iatomnext(j)

  end do

end subroutine gnbcell

!
! Function that computes, for a pair of atoms, the non-bonded gradient
!

subroutine gnbpair(iatom,j,x,g)
 
  use force_field
  use linkedcells, only : knonbonded
  implicit none
  integer :: inb, iatom, j, ix, iy, iz, jx, jy, jz
  double precision :: x(natoms*3), g(natoms*3)
  double precision :: d2, d, d3, qqd3, qqd3x, qqd3y, qqd3z, dx, dy, dz,&
                      c, dinv6  

  inb = knonbonded(iatom,j) 

  if ( inb == 0 ) return

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

  if ( d2 > cutoff2 ) return

  d = dsqrt(d2)

  if ( d > cutnb(inb,1) ) then

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

  else

    g(ix) = g(ix) + cutnb(inb,2)*dx/d
    g(iy) = g(iy) + cutnb(inb,2)*dy/d
    g(iz) = g(iz) + cutnb(inb,2)*dz/d
    g(jx) = g(jx) - cutnb(inb,2)*dx/d
    g(jy) = g(jy) - cutnb(inb,2)*dy/d
    g(jz) = g(jz) - cutnb(inb,2)*dz/d

  end if

end subroutine gnbpair 
