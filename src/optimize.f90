!
! Subroutine optimize: Chooses which optimization code will be used
!

subroutine optimize(n,x,f,outputunit,optpars,minmethod)

  implicit none
  integer :: optpars(10)
  integer :: n, outputunit, minmethod
  double precision :: x(n), f
  external :: computef, computeg

  ! Use spg method

  if ( minmethod == 1 ) then
    call callspg(n,x,f,outputunit,computef,computeg,optpars)
  end if

  ! Use cgnewton method

  if ( minmethod == 2 ) then
    call callcgnewton(n,x,f,outputunit,computef,computeg,optpars)
  end if

return
end subroutine optimize
