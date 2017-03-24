!
! Subroutine compute_enbcell: Computes the non-bonded energy with short-range modification,
!                             and a cutoff, using the linked-cell method to avoid the evaluation
!                             of distances greater than the cutoff
! 
! L. Martinez, Aug 15, 2014
! Institute of Chemistry - State University of Campinas - Brazil
!
! San Lorenzo Campeón! (13/08/2014)
!

subroutine compute_enbcell(x,etot)
  
  use force_field
  use linkedcells
  implicit none
  integer :: i, j, k, iatom 
  double precision :: x(natoms*3), etot, enbcell

  ! Initialize the linked cells with the atoms within

  call cells(x)

  ! Compute the energy

  etot = 0.d0
  do i = 1, nboxes(1)
    do j = 1, nboxes(2)
      do k = 1, nboxes(3)

        iatom = iatomfirst(i,j,k)

        do while( iatom /= 0 ) 

          ! If only computing for backbone atoms, cycle (not the best way
          ! to do this, but ok for the moment. 

          if ( .not. ffcomp(7) .and. .not. isbackbone(iatom) ) then
            iatom = iatomnext(iatom)
            cycle
          end if

          ! Interactions inside box

          etot = etot + enbcell(iatom,iatomnext(iatom),x)

          ! Interactions of boxes that share faces

          etot = etot + enbcell(iatom,iatomfirst(i+1,j,k),x)
          etot = etot + enbcell(iatom,iatomfirst(i,j+1,k),x)
          etot = etot + enbcell(iatom,iatomfirst(i,j,k+1),x)

          ! Interactions of boxes that share axes

          etot = etot + enbcell(iatom,iatomfirst(i+1,j+1,k),x)
          etot = etot + enbcell(iatom,iatomfirst(i+1,j,k+1),x)
          etot = etot + enbcell(iatom,iatomfirst(i+1,j-1,k),x)
          etot = etot + enbcell(iatom,iatomfirst(i+1,j,k-1),x)
          etot = etot + enbcell(iatom,iatomfirst(i,j+1,k+1),x)
          etot = etot + enbcell(iatom,iatomfirst(i,j+1,k-1),x)

          ! Interactions of boxes that share vertices

          etot = etot + enbcell(iatom,iatomfirst(i+1,j+1,k+1),x)
          etot = etot + enbcell(iatom,iatomfirst(i+1,j+1,k-1),x)
          etot = etot + enbcell(iatom,iatomfirst(i+1,j-1,k+1),x)
          etot = etot + enbcell(iatom,iatomfirst(i+1,j-1,k-1),x)

          iatom = iatomnext(iatom)
        end do

      end do
    end do
  end do

  return
end subroutine compute_enbcell

!
! Function that computes the interactions of atom iatom with all
! atoms jatom of the box with first atom jatom entered
!

double precision function enbcell(iatom,jatom,x)
  
  use force_field, only : natoms, ffcomp, isbackbone
  use linkedcells, only : iatomnext
  implicit none
  integer :: iatom, jatom, j
  double precision :: x(natoms*3), enbpair

  j = jatom
  enbcell = 0.d0
  do while( j /= 0 ) 

    ! Compute only for backbone atoms if that is the case
    if ( .not. ffcomp(7) .and. .not. isbackbone(j) ) then
      j = iatomnext(j)
      cycle
    end if

    enbcell = enbcell + enbpair(iatom,j,x)
    j = iatomnext(j)

  end do

end function enbcell

!
! Function that computes, for a pair of atoms, the non-bonded interaction
! This function contains the decision based on the exclusion list
!

double precision function enbpair(iatom,j,x)
 
  use force_field
  use linkedcells, only : knonbonded
  implicit none
  integer :: inb, iatom, j, ix, iy, iz, jx, jy, jz
  double precision :: x(natoms*3), d2, d, elec_energy, vdw_energy, p6, p12

  inb = knonbonded(iatom,j) 

  if ( inb == 0 ) then
    enbpair = 0.d0
    return
  end if

  ix = (ijnonbonded(inb,1)-1)*3 + 1
  iy = ix + 1
  iz = ix + 2

  jx = (ijnonbonded(inb,2)-1)*3 + 1
  jy = jx + 1
  jz = jx + 2

  d2 = (x(ix) - x(jx))**2 + & 
       (x(iy) - x(jy))**2 + &
       (x(iz) - x(jz))**2

  if ( d2 > cutoff2 ) then
    enbpair = 0.d0 
    return
  end if

  d = dsqrt(d2)

  if ( d > cutnb(inb,1) ) then

    elec_energy = qq(inb)/d

    p6 = ss6(inb) / d2**3
    p12 = p6*p6
    vdw_energy = epseps(inb)*( p12 - 2.d0*p6 )

    enbpair = elec_energy + vdw_energy

  else

    enbpair = cutnb(inb,2)*d + cutnb(inb,3)

  end if

end function enbpair 

!
! Subroutine that sets up the cells for the linked cell method
!

subroutine cells(x)
  
  use force_field, only : natoms, cutoff
  use linkedcells
  integer :: i, j, k, ix, iy, iz, iboxx, iboxy, iboxz
  double precision :: x(natoms*3)
  double precision :: xmin(3), xmax(3)

  ! Compute maximum and minimum coordinates for all atoms, to build
  ! the bounding box

  i = 1
  ix = (i-1)*3 + 1
  iy = ix + 1
  iz = ix + 2
  xmin(1) = x(ix)
  xmin(2) = x(iy)
  xmin(3) = x(iz)
  xmax(1) = x(ix)
  xmax(2) = x(iy)
  xmax(3) = x(iz)
  do i = 2, natoms
    ix = (i-1)*3 + 1
    iy = ix + 1
    iz = ix + 2
    xmin(1) = dmin1(x(ix),xmin(1))
    xmin(2) = dmin1(x(iy),xmin(2))
    xmin(3) = dmin1(x(iz),xmin(3))
    xmax(1) = dmax1(x(ix),xmax(1))
    xmax(2) = dmax1(x(iy),xmax(2))
    xmax(3) = dmax1(x(iz),xmax(3))
  end do

  ! To avoid problems if x == xmax or x == xmin in setting the boxes:

  xmin(1) = xmin(1) - 1.d-5*(xmax(1)-xmin(1))
  xmin(2) = xmin(2) - 1.d-5*(xmax(2)-xmin(2))
  xmin(3) = xmin(3) - 1.d-5*(xmax(3)-xmin(3))
  xmax(1) = xmax(1) + 1.d-5*(xmax(1)-xmin(1))
  xmax(2) = xmax(2) + 1.d-5*(xmax(2)-xmin(2))
  xmax(3) = xmax(3) + 1.d-5*(xmax(3)-xmin(3))

  ! Setting the number of cells in each direction

  nboxes(1) = int( (xmax(1)-xmin(1)) / cutoff ) + 1
  if ( nboxes(1) <= maxboxes ) then
    boxcut(1) = cutoff
  else
    nboxes(1) = maxboxes
    boxcut(1) = (xmax(1) - xmin(1)) / maxboxes
  end if
    
  nboxes(2) = int( (xmax(2)-xmin(2)) / cutoff ) + 1
  if ( nboxes(2) <= maxboxes ) then
    boxcut(2) = cutoff
  else
    nboxes(2) = maxboxes
    boxcut(2) = (xmax(2) - xmin(2)) / maxboxes
  end if

  nboxes(3) = int( (xmax(3)-xmin(3)) / cutoff ) + 1
  if ( nboxes(3) <= maxboxes ) then
    boxcut(3) = cutoff
  else
    nboxes(3) = maxboxes
    boxcut(3) = (xmax(3) - xmin(3)) / maxboxes
  end if

  ! Initializatin iatomfix array

  do i = 0, nboxes(1) + 1
    do j = 0, nboxes(2) + 1
      do k = 0, nboxes(3) + 1
        iatomfirst(i,j,k) = 0
      end do
    end do
  end do

  ! Now set in which box each atom is 

  do i = 1, natoms
    ix = (i-1)*3 + 1
    iy = ix + 1
    iz = ix + 2
    iboxx = int( (x(ix)-xmin(1)) / boxcut(1) ) + 1
    iboxy = int( (x(iy)-xmin(2)) / boxcut(2) ) + 1
    iboxz = int( (x(iz)-xmin(3)) / boxcut(3) ) + 1
    iatomnext(i) = iatomfirst(iboxx,iboxy,iboxz)
    iatomfirst(iboxx,iboxy,iboxz) = i
  end do

  return
end subroutine cells

