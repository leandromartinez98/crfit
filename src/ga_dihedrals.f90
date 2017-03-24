!
! subroutine nc_dihedrals: Reads the protein structure and annotates
!                          which are the atoms N and CA of a Phi dihedral
!                          angle. This discriminates which atoms are
!                          to be rotated when a dihedral is rigid-body
!                          changed on the genetic algorithm code. 
!
! L. Martinez
! Institute of Chemistry - State University of Campinas, UNICAMP - Brazil
! Jun 24, 2014 (Italia 0 x Uruguai 1)
!
!

subroutine nc_dihedrals()

  use force_field, only : natoms, nresidues, name, residue_number
  use dihedrals
  implicit none 
  integer :: i, j

  allocate( backbone(nresidues,3) )

  ! All atoms with indexes greater than that of the ca_atom are rotated, the
  ! axis of rotation been the N-CA bond. This strategy depends on the PSF
  ! file being ordered such that every residue starts with the N atom.
  ! Apparently this is standard, at least of PSF files generated with the
  ! CHARMM psfgen tool. The last residue of the sequence is in a different
  ! order, but that will be ignored because there is no relevant dihedral
  ! to be rotated.  

  do i = 1, natoms
    j = residue_number(i) - residue_number(1) + 1
    if ( name(i) == "N" ) backbone(j,1) = i
    if ( name(i) == "CA" ) backbone(j,2) = i
    if ( name(i) == "C" ) backbone(j,3) = i
  end do

  return
end subroutine nc_dihedrals

!
! Function that computes the PHI dihedral angle of a residue
!

double precision function compute_phi(iresidue,x)

  use dihedrals
  implicit none
  integer :: iresidue, i, j, k, l
  double precision :: x(*), dihedangle

  i = backbone(iresidue-1,3)
  j = backbone(iresidue,1)
  k = backbone(iresidue,2)
  l = backbone(iresidue,3)
  compute_phi = dihedangle(i,j,k,l,x)

  return
end function compute_phi

!
! Function that computes the PSI dihedral angle of a residue
!

double precision function compute_psi(iresidue,x)

  use dihedrals
  implicit none
  integer :: iresidue, i, j, k, l
  double precision :: x(*), dihedangle

  i = backbone(iresidue,1)
  j = backbone(iresidue,2)
  k = backbone(iresidue,3)
  l = backbone(iresidue+1,1)
  compute_psi = dihedangle(i,j,k,l,x)

  return
end function compute_psi

!
! Rotate PHI dihedral angle of a residue
!

subroutine rotate_phi(iresidue,x,phi)

  use force_field, only : natoms
  use dihedrals
  integer :: iresidue, i1, i2
  double precision :: phi, x(natoms*3)  

  i1 = backbone(iresidue,1)
  i2 = backbone(iresidue,2)
  call rotate_around_axis(i1,i2,x,phi)

  return
end subroutine rotate_phi

!
! Rotate PSI dihedral angle of a residue
!

subroutine rotate_psi(iresidue,x,psi)

  use force_field, only : natoms
  use dihedrals
  integer :: iresidue, i1, i2
  double precision :: psi, x(natoms*3)  

  i1 = backbone(iresidue,2)
  i2 = backbone(iresidue,3)
  call rotate_around_axis(i1,i2,x,psi)

  return
end subroutine rotate_psi

!
!
! subroutine rotate_around_axis: Given the indexes of two atoms, rotates rigidly
!                                the rest of the structure according around the axis
!                                of this dihedral 
!

subroutine rotate_around_axis(i1,i2,x,theta)

  use force_field, only : natoms, iorder
  implicit none
  integer :: i1, i2, j1, j2, ii, i, j, k
  double precision :: x(natoms*3), theta, ux, uy, uz, axis_norm, rotmat(3,3), &
                      cost, sint, onemcost, ux2, uy2, uz2, xtmp(3)

  ! Compute rotation matrix. Axis vector is the unitary vector connecting
  ! atoms N and C with indices in and ic provided. 

  j1 = (i1-1)*3
  j2 = (i2-1)*3
  ux = x(j2+1) - x(j1+1)
  uy = x(j2+2) - x(j1+2)
  uz = x(j2+3) - x(j1+3)
  ux2 = ux**2
  uy2 = uy**2
  uz2 = uz**2
  axis_norm = dsqrt( ux2 + uy2 + uz2 )
  ux = ux / axis_norm
  uy = uy / axis_norm
  uz = uz / axis_norm
  ux2 = ux**2
  uy2 = uy**2
  uz2 = uz**2
  
  cost = dcos(theta)
  sint = dsin(theta)
  onemcost = 1 - cost

  rotmat(1,1) = cost + ux2 * onemcost
  rotmat(1,2) = ux*uy*onemcost - uz*sint
  rotmat(1,3) = ux*uz*onemcost + uy*sint
  rotmat(2,1) = uy*ux*onemcost + uz*sint
  rotmat(2,2) = cost + uy2*onemcost
  rotmat(2,3) = uy*uz*onemcost - ux*sint
  rotmat(3,1) = uz*ux*onemcost - uy*sint
  rotmat(3,2) = uz*uy*onemcost + ux*sint
  rotmat(3,3) = cost + uz2*onemcost

  ! Apply rotaion to all atoms with indices greater than i1

  do ii = 1, natoms
    if ( iorder(ii) > iorder(i1) ) then
      i = (ii-1)*3 + 1
      j = i + 1
      k = i + 2
      xtmp(1) = x(i) - x(j1+1)
      xtmp(2) = x(j) - x(j1+2)
      xtmp(3) = x(k) - x(j1+3)
      xtmp = matmul(rotmat,xtmp)
      x(i) = xtmp(1) + x(j1+1) 
      x(j) = xtmp(2) + x(j1+2)
      x(k) = xtmp(3) + x(j1+3)
    end if
  end do

  return
end subroutine rotate_around_axis

