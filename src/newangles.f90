!
! Create new phi and psi angles for a residue, using more or less the
! knowlege of a Ramachandran plot
!

subroutine newangles(phi,psi)

  implicit none
  double precision, parameter :: torad = 3.141592d0/180.d0
  double precision :: structure, random, phi, psi

  call random_number(structure)

  ! Beta-sheet

  if ( structure <= 0.4d0 ) then
    call random_number(random)
    phi = -130.d0 + 70.d0*random
    call random_number(random)
    psi = 90.d0 + 80.d0*random

  ! Right-handed alpha-helix

  else if ( structure > 0.4d0 .and. structure <= 0.8d0 ) then
    call random_number(random)
    phi = -145.d0 + 100.d0*random
    call random_number(random)
    psi = -67.d0 + 34.d0*random

  ! Left handed alpha-helix

  else 
    call random_number(random)
    phi = 45.d0 + 14.d0*random
    call random_number(random)
    psi = 28.d0 + 65.d0*random
  end if

  phi = phi*torad
  psi = psi*torad
  
  return
end subroutine newangles

