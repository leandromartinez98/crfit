!
! Subroutine compute_gconst: Computes the gradient relative constraints
! 
! L. Martinez, Jul, 2013
! Institute of Chemistry - State University of Campinas - Brazil
!

subroutine compute_gconst(x,g)
  
  use force_field, only : natoms
  use constraints
  use flashsortvars
  implicit none
  integer :: i, ix, iy, iz, jx, jy, jz, igroup, iconst, iconsider, ifirst
  double precision :: d, g(natoms*3), x(natoms*3), vmin

  ! Compute individual constraint violations

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

  ! Add contribution to the gradient

  iconst = 1 
  do igroup = 1, ncgroups

    ! If all constraints of this group are to be considered

    if ( ncpergroup(igroup) == ncconsider(igroup) ) then

      do i = 1, ncconsider(igroup)
        call gpair(iconst,x,g)
        iconst = iconst + 1
      end do

    ! If only one constraint of this group is to be considered

    else if ( ncpergroup(igroup) == 1 ) then

      iconsider = iconst
      vmin = violation(iconst)
      do i = 2, ncpergroup(igroup) 
        if ( violation(iconst) < vmin ) then
          vmin = violation(iconst)
          iconsider = iconst
        end if
        iconst = iconst + 1
      end do
      call gpair(iconsider,x,g) 

    ! Otherwise, sort the violations and add the lowest ones 

    else

      ! Sort the constraint violations from lowest to highest

      ifirst = iconst
      do i = 1, ncpergroup(igroup) 
        viol(i) = violation(iconst)
        iconst = iconst + 1
      end do
      call flashsort(viol, ncpergroup(igroup), lflash, mflash, indflash, indflashback)

      ! Now, add to the gradient the contribution of the ncconsider smallest violations

      do i = 1, ncconsider(igroup)
        iconsider = ifirst + indflash(i) - 1
        call gpair(iconsider,x,g)
      end do

    end if

  end do

  return
end subroutine compute_gconst

! Subroutine that given the index of the constraint, adds its contribution
! to the gradient

subroutine gpair(iconst,x,g)

  use constraints
  implicit none
  integer :: iconst, ix, iy, iz, jx, jy, jz
  double precision :: x(*), g(*), d

  ix = (iconstraint(iconst,1)-1)*3+1
  iy = ix + 1
  iz = ix + 2

  jx = (iconstraint(iconst,2)-1)*3+1
  jy = jx + 1
  jz = jx + 2

  d = (x(jx)-x(ix))**2 + (x(jy)-x(iy))**2 + (x(jz)-x(iz))**2

  if ( d < dconstraint(iconst,1) ) then

    g(ix) = g(ix) + 2.d0*kconstraint(iconst)*( x(jx) - x(ix) ) 
    g(iy) = g(iy) + 2.d0*kconstraint(iconst)*( x(jy) - x(iy) ) 
    g(iz) = g(iz) + 2.d0*kconstraint(iconst)*( x(jz) - x(iz) ) 

    g(jx) = g(jx) - 2.d0*kconstraint(iconst)*( x(jx) - x(ix) ) 
    g(jy) = g(jy) - 2.d0*kconstraint(iconst)*( x(jy) - x(iy) ) 
    g(jz) = g(jz) - 2.d0*kconstraint(iconst)*( x(jz) - x(iz) ) 

  end if

  if ( d > dconstraint(iconst,2) ) then

    g(ix) = g(ix) - 2.d0*kconstraint(iconst)*( x(jx) - x(ix) ) 
    g(iy) = g(iy) - 2.d0*kconstraint(iconst)*( x(jy) - x(iy) ) 
    g(iz) = g(iz) - 2.d0*kconstraint(iconst)*( x(jz) - x(iz) ) 

    g(jx) = g(jx) + 2.d0*kconstraint(iconst)*( x(jx) - x(ix) ) 
    g(jy) = g(jy) + 2.d0*kconstraint(iconst)*( x(jy) - x(iy) ) 
    g(jz) = g(jz) + 2.d0*kconstraint(iconst)*( x(jz) - x(iz) ) 

  end if

  return
end subroutine gpair














