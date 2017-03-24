!
! Subroutine computefls: Computes the total objective function value
! using the short-range modification for non-bonded interactions
!
! L. Martínez, Dec 9, 2013
! Institute of Chemistry - State University of Campinas - Brazil
!

subroutine computefls(x,f)

  use force_field, only : ffcomp
  implicit none
  double precision :: x(*), f, f1

  f = 0.d0

  if ( ffcomp(1) ) then 
    call compute_ebond(x,f1)
    f = f + f1
  end if

  if ( ffcomp(2) ) then
    call compute_eangle(x,f1)
    f = f + f1
  end if
 
  if ( ffcomp(3) ) then
    call compute_edihed(x,f1)
    f = f + f1
  end if
  
  if ( ffcomp(4) ) then
    call compute_eimpr(x,f1)
    f = f + f1
  end if

  if ( ffcomp(5) ) then
    call compute_enbls(x,f1)  
    f = f + f1
  end if

  if ( ffcomp(6) ) then
    call compute_fconst(x,f1)
    f = f + f1
  end if
  
  return
end subroutine computefls

