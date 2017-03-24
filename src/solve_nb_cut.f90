!
! subroutine solve_nb_cut: finds the point in which the non-bonded potential
!                          will be cut to a straight line
!
! L. Martinez, Dec 10, 2013
! Institute of Chemistry, State University of Campinas
!

subroutine solve_nb_cut(eps,sig6,q,xcut,dydx,ycut)

  implicit none
  double precision :: eps, q, sig6, sig12, xcut, ycut, dydx, &
                      lbound, ubound, precision, y

  sig12 = sig6**2
  precision = 1.d-12
  lbound = 0.d0
  ubound = 50.d0
  xcut = lbound + ( ubound - lbound ) / 2.d0
  if ( eps > 0.d0 .and. sig6 > 0.d0 ) then
    do
      y = 13.d0*eps*sig12/(xcut**12) - 14.d0*eps*sig6/(xcut**6) + 2.d0*q/xcut - ycut
      if ( dabs(y) < precision ) exit
      if ( y > 0.d0 ) then
        lbound = xcut
      else
        ubound = xcut
      end if
      xcut = lbound + ( ubound - lbound ) / 2.d0
      if ( dabs(xcut) < precision ) then
        write(*,*) ' ERROR: xcut < precision in solve_nb_cut. Contact the developer. '
        call stop_all()
      end if
    end do
    dydx = eps*(-12.d0*sig12/xcut**13 + 12.d0*sig6/xcut**7) - q/xcut**2
    if ( dydx > 0.d0 ) then
      write(*,*) ' ERROR: The linear cut approximation for short distances '
      write(*,*) '        has a positive slope, and there is a non-null vdW '
      write(*,*) '        interaction. Contact the developer. '
      call stop_all()
    end if
  else
    do
      if ( q > 0.d0 ) then
        y = 2.d0*q/xcut - ycut
        if ( dabs(y) < precision ) exit
        if ( y > 0.d0 ) then
          lbound = xcut
        else
          ubound = xcut
        end if
      end if
      if ( q < 0.d0 ) then
        y = 2.d0*q/xcut + ycut
        if ( dabs(y) < precision ) exit
        if ( y > 0.d0 ) then
          ubound = xcut
        else
          lbound = xcut
        end if
      end if
      xcut = lbound + ( ubound - lbound ) / 2.d0
    end do
    dydx = -q/xcut**2
  end if

  return
end subroutine solve_nb_cut


