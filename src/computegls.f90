!
! Subroutine computegls: Computes the total gradient of the objective function,
!                        for the short-range modified function
!
! L. Martínez, Dec 9, 2013
! Institute of Chemistry - State University of Campinas - Brazil
!

subroutine computegls(n,x,g)

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
  if ( ffcomp(5) ) call compute_gnbls(x,g)  
  if ( ffcomp(6) ) call compute_gconst(x,g)
  
  return
end subroutine computegls

