!
! Subroutine compute_bonds: Computes the covalent bond interaction energy
! 
! L. Martinez, Nov 28, 2013
! Institute of Chemistry - State University of Campinas - Brazil
!

subroutine compute_ebond(x,ebond)
  
  use force_field
  implicit none
  integer :: i, ix, jx
  double precision :: ebond, x(*), d
  
  ebond = 0.d0
  do i = 1, nbonds
    ix = (ibond(i,1) - 1)*3
    jx = (ibond(i,2) - 1)*3
    d = dsqrt( (x(ix+1) - x(jx+1))**2 + &
               (x(ix+2) - x(jx+2))**2 + &
               (x(ix+3) - x(jx+3))**2 )
    ebond = ebond + fbond(i,1)*(d - fbond(i,2))**2
  end do

  return
end subroutine compute_ebond

