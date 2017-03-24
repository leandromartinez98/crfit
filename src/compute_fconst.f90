!
! Subroutine compute_fconst: Computes the function value relative to constraints 
! 
! L. Martinez, Jun 2, 2014
! Institute of Chemistry - State University of Campinas - Brazil
!

subroutine compute_fconst(x,fconst)
  
  use force_field, only : natoms
  use constraints
  use flashsortvars
  implicit none
  integer :: i, ix, iy, iz, jx, jy, jz, iconst, igroup
  double precision :: x(natoms*3), fconst, d, vmin

  ! Compute the individual constraint violations

  do i = 1, ntotconst
    violation(i) = 0.d0
  end do

  iconst = 0
  do igroup = 1, ncgroups
    do i = 1, ncpergroup(igroup)
      iconst = iconst + 1

      ix = (iconstraint(iconst,1)-1)*3+1
      iy = ix + 1
      iz = ix + 2

      jx = (iconstraint(iconst,2)-1)*3+1
      jy = jx + 1
      jz = jx + 2

      d = (x(jx)-x(ix))**2 + (x(jy)-x(iy))**2 + (x(jz)-x(iz))**2

      if ( d < dconstraint(iconst,1) ) then
        violation(iconst) = kconstraint(iconst)*( dconstraint(iconst,1) - d )
      end if

      if ( d > dconstraint(iconst,2) ) then
        violation(iconst) = kconstraint(iconst)*( d - dconstraint(iconst,2) )
      end if 
    end do
  end do

  ! Compute the function value

  fconst = 0.d0
  iconst = 1
  do igroup = 1, ncgroups

    if ( ncpergroup(igroup) == ncconsider(igroup) ) then

      do i = 1, ncpergroup(igroup)
        fconst = fconst + violation(iconst)
        iconst = iconst + 1
      end do

    else if ( ncconsider(igroup) == 1 ) then

      vmin = violation(iconst)
      do i = 2, ncpergroup(igroup)
        iconst = iconst + 1
        vmin = dmin1(vmin,violation(iconst))
      end do
      fconst = fconst + vmin

    else
      
      do i = 1, ncpergroup(igroup) 
        viol(i) = violation(iconst)
        iconst = iconst + 1
      end do
      call flashsort(viol, ncpergroup(igroup), lflash, mflash, indflash, indflashback)
      do i = 1, ncconsider(igroup)
        fconst = fconst + viol(i)
      end do

    end if

  end do

  return
end subroutine compute_fconst

