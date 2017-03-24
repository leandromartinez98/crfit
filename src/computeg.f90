!
! Subroutine computeg: Computes the total gradient of the objective function
!
! L. Martínez, Dec 9, 2013
! Institute of Chemistry - State University of Campinas - Brazil
!

subroutine computeg(n,x,g)

  use force_field, only : ffcomp, useroutine
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
  if ( ffcomp(5) ) then 
    if ( useroutine == 1 ) then
      call compute_gnb(x,g)  
    else if ( useroutine == 2 ) then
      call compute_gnbls(x,g)  
    else if ( useroutine == 3 ) then
      call compute_gnbcell(x,g)  
    end if
  end if
  if ( ffcomp(6) ) call compute_gconst(x,g)
  
  return
end subroutine computeg

