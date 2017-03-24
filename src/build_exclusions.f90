!
! Subroutine build_exclusions: Build the lists of atoms excluded from the computation
!     of non-bonded interactions, from the bond data (excluded 1-4).
!
! L. Martinez, Nov 27, 2013
! Institute of Chemistry - State University of Campinas - Brazil
! 

subroutine build_exclusions()

  use force_field
  use linkedcells
  implicit none

  ! To build exclusion lists, we will create temporary vectors that have, for
  ! each atom, which atoms are bound to it. 

  integer :: i, j, k, ii, jj, kk, kbb
  integer, allocatable :: bond_list(:,:), bond_count(:)
  logical :: scale14

  ! The electrostatic constant
  double precision :: xcut, ycut, dydx
  double precision, parameter :: qelectrostatic = 332.05382d0

  allocate ( bond_list(natoms,8), bond_count(natoms) )

  write(*,*) ' Computing non-bonded interaction list. '

  ! Reseting bond count

  do i = 1, natoms
    bond_count(i) = 0
  end do

  ! Filling up the bond_list vector

  do i = 1, nbonds
    ii = ibond(i,1)
    jj = ibond(i,2)
    bond_count(ii) = bond_count(ii) + 1
    if ( bond_count(ii) > 8 ) then
      write(*,*) ' ERROR: Found more than 8 bonds for atom: ', ii
      call stop_all()
    end if
    bond_list(ii,bond_count(ii)) = ibond(i,2)
    bond_count(jj) = bond_count(jj) + 1
    if ( bond_count(ii) > 8 ) then
      write(*,*) ' ERROR: Found more than 8 bonds for atom: ', jj
      call stop_all()
    end if
    bond_list(jj,bond_count(jj)) = ibond(i,1)
  end do

  !
  ! Building inclusion lists
  ! 

  ! Counting the total number of interactions to allocate include vector

  n_nonbonded = 0
  do i = 1, natoms - 1
    jloop1: do j = i + 1, natoms

      do ii = 1, bond_count(i)

        ! 1-2 exclusions

        if ( j == bond_list(i,ii) ) cycle jloop1
 
        ! 1-3 exclusions

        do jj = 1, bond_count(bond_list(i,ii))
          if ( j == bond_list(bond_list(i,ii),jj) ) cycle jloop1
        end do

      end do

      ! If arrived here, j is not listed in 1-2 nor 1-3 exclusions, so include it

      n_nonbonded = n_nonbonded + 1 
   
      ! If i and j are backbone atoms, include them into the backbone
      ! interaction list

      if ( isbackbone(i) .and. isbackbone(j) ) then
        n_nbbackbone = n_nbbackbone + 1
      end if

    end do jloop1
  end do 
  write(*,*) ' Total number of non-bonded interactions: ', n_nonbonded
  write(*,*) ' Number of backbone non-bonded interactions: ', n_nbbackbone

  ! This is to initialize the vectors of the linked cell method

  maxboxes = min0(int(float(nresidues)/cutoff) + 1, 10) 

  ! Allocating include vector and filling it

  allocate( ijnonbonded(n_nonbonded,2), qq(n_nonbonded), qqreal(n_nonbonded), &
            epseps(n_nonbonded), ss6(n_nonbonded), cutnb(n_nonbonded,3), &
            knonbonded(natoms,natoms), inbbackbone(n_nbbackbone), &
            iatomfirst(0:maxboxes+1,0:maxboxes+1,0:maxboxes+1), iatomnext(natoms) )

  ! Build exclusion lists

  write(*,*) ' Building exclusion lists ... '

  k = 0
  kbb = 0
  do i = 1, natoms - 1
    jloop2: do j = i + 1, natoms

      knonbonded(i,j) = 0
      knonbonded(j,i) = 0

      scale14 = .false.
      do ii = 1, bond_count(i)

        ! 1-2 exclusions

        if ( j == bond_list(i,ii) ) cycle jloop2

        ! 1-3 exclusions

        do jj = 1, bond_count(bond_list(i,ii))
          if ( j == bond_list(bond_list(i,ii),jj) ) cycle jloop2

          ! 1-4 scaling
          do kk = 1, bond_count(bond_list(bond_list(i,ii),jj))
            if ( j == bond_list(bond_list(bond_list(i,ii),jj),kk) ) scale14 = .true.
          end do

        end do

      end do

      ! If arrived here, j is not listed in 1-2 nor 1-3 exclusions, so include it

      k = k + 1 
      ijnonbonded(k,1) = i
      ijnonbonded(k,2) = j

      ! And lets save time, and compute the combinations of the interactions parameters
      ! between pairs

      if ( .not. scale14 ) then
        if ( isbackbone(i) .and. isbackbone(j) ) then
          qqreal(k) = qelectrostatic*charge(i)*charge(j)
        else
          qqreal(k) = (1.d0/dielectric)*qelectrostatic*charge(i)*charge(j)
        end if
        qq(k) = qqreal(k)
        epseps(k) = dsqrt( eps(i) * eps(j) )
        ss6(k) = (sig(i) + sig(j))**6
      else 
        if ( isbackbone(i) .and. isbackbone(j) ) then
          qqreal(k) = qelectrostatic*charge(i)*charge(j)
        else
          qqreal(k) = (1.d0/dielectric)*qelectrostatic*charge(i)*charge(j)
        end if
        qq(k) = qqreal(k)
        epseps(k) = dsqrt( eps14(i) * eps14(j) )
        ss6(k) = (sig14(i) + sig14(j))**6
      end if

      ! And the parameters of the short distance linear extrapolation

      ycut = 100.d0
      call solve_nb_cut(epseps(k),ss6(k),qq(k),xcut,dydx,ycut)
      cutnb(k,1) = xcut
      cutnb(k,2) = dydx
      cutnb(k,3) = ycut

      ! Annotate the index of this interaction for linked cell computation

      knonbonded(i,j) = k
      knonbonded(j,i) = k

      if ( isbackbone(i) .and. isbackbone(j) ) then
        kbb = kbb + 1
        inbbackbone(kbb) = k
      end if

    end do jloop2
  end do 

  ! We don't need the bond lists anymore, so get rid of them

  deallocate( bond_count, bond_list )

  return
end subroutine build_exclusions

!
! Subroutine updatels: Updates the linear extraploation at short
! distances if the charges are changed (for turning sidechains neutral
! or polar, or vice-versa)
!
! L. Martinez, Jan 15, 2015 (Est-ce que je suis Charlie?)
!

subroutine updatels()

  use force_field
  implicit none
  integer :: inb, i, j
  double precision :: xcut, ycut, dydx

  do inb = 1, n_nonbonded

    i = ijnonbonded(inb,1)
    j = ijnonbonded(inb,2)
    
    if ( ffcomp(8) ) then 
      qq(inb) = qqreal(inb)
    else
      if ( isbackbone(i) .and. isbackbone(j) ) then
        qq(inb) = qqreal(inb)
      else
        qq(inb) = 0.d0
      end if
    end if

    ycut = 100.d0
    call solve_nb_cut(epseps(inb),ss6(inb),qq(inb),xcut,dydx,ycut)
    cutnb(inb,1) = xcut
    cutnb(inb,2) = dydx
    cutnb(inb,3) = ycut

  end do

return
end





