!
! Subroutine computef: Computes the total objective function value
!
! L. Martínez, Dec 9, 2013
! Institute of Chemistry - State University of Campinas - Brazil
!

subroutine computef(x,f)

  use force_field, only : ffcomp, useroutine, ecomp
  implicit none
  double precision :: x(*), f

  ecomp(1) = 0.d0
  ecomp(2) = 0.d0
  ecomp(3) = 0.d0
  ecomp(4) = 0.d0
  ecomp(5) = 0.d0
  ecomp(6) = 0.d0
  ecomp(7) = 0.d0
  ecomp(8) = 0.d0

  if ( ffcomp(1) ) call compute_ebond(x,ecomp(1))
  if ( ffcomp(2) ) call compute_eangle(x,ecomp(2))
  if ( ffcomp(3) ) call compute_edihed(x,ecomp(3))
  if ( ffcomp(4) ) call compute_eimpr(x,ecomp(4))
  if ( ffcomp(5) ) then
    if ( useroutine == 1 ) then
      call compute_enb(x,ecomp(5),ecomp(6))  
      ecomp(8) = ecomp(5) + ecomp(6)
    else if ( useroutine == 2 ) then
      call compute_enbls(x,ecomp(8))
    else if ( useroutine == 3 ) then
      call compute_enbcell(x,ecomp(8))
    end if
  end if
  if ( ffcomp(6) ) call compute_fconst(x,ecomp(7))

  f = ecomp(1) + ecomp(2) + ecomp(3) + ecomp(4) + ecomp(7) + ecomp(8)
  
  return
end subroutine computef

