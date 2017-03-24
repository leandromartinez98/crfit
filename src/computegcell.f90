!
! Subroutine computegcell: Computes the gradient using the linked cell
!                          method
!
! L. Martínez, Aug 18, 2014
! Institute of Chemistry - State University of Campinas - Brazil
!

subroutine computegcell(n,x,g)

  use force_field, only : ffcomp
  implicit none
  integer :: i, n
  double precision :: x(*), g(*)

  ! Important: in the sub-routines, the gradient is not zeroed, so that 
  ! it accumulates in g

  do i = 1, n
    g(i) = 0.d0
  end do

  if ( ffcomp(1) ) call compute_gbond(x,g)
  if ( ffcomp(2) ) call compute_gangle(x,g)
  if ( ffcomp(3) ) call compute_gdihed(x,g)
  if ( ffcomp(4) ) call compute_gimpr(x,g)
  if ( ffcomp(5) ) call compute_gnbcell(x,g)  
  if ( ffcomp(6) ) call compute_gconst(x,g)
  
  return
end subroutine computegcell

