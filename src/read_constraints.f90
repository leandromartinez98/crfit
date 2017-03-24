!
! Subroutine read_constraints: Reads the file containing the user-defined
!                              distance constraints
!
! L. Martinez, Dec 17, 2013
! Institute of Chemistry, State University of Campinas
!

subroutine read_constraints(constraint_file)

  use force_field
  use constraints
  implicit none
  integer :: iconst, i, j, ioer, ires1, ires2, ioerr, ifile
  double precision :: dmin, dmax, kread
  character(len=200) :: record, constraint_file(*), keyword
  character(len=5) :: res1, name1, chain1, res2, name2, chain2

  ! Counting the number of constraints and the number of groups
  ! of constraints

  ncgroups = 0
  ntotconst = 0
  do ifile = 1, nconstfiles

    record = constraint_file(ifile)
    open(10,file=record,status="old",action="read",iostat=ioerr)
    if ( ioerr /= 0 ) then
      write(*,*) ' ERROR: Could not find or open constraint file: ', trim(constraint_file(ifile))
      call stop_all()
    end if
    
    do
      read(10,"( a200 )",iostat=ioer) record
      if ( ioer /= 0 ) exit
      if ( record(1:1) == "#" .or. len(trim(record)) < 1 ) cycle
      if ( keyword(record) == "set" ) then
        read(record,*,iostat=ioer) res1, i, res1, i
        if ( ioer /= 0 ) then
          write(*,*) ' ERROR: Fail to read constraint set definition in file: ', &
                     trim(constraint_file(ifile))
          call stop_all()
        end if
        ncgroups = ncgroups + 1
        cycle
      end if
      read(record,*,iostat=ioer) i, j, dmin, dmax, kread
      if ( ioer == 0 ) then
        ntotconst = ntotconst + 1
        cycle
      end if
      read(record,*,iostat=ioer) res1, ires1, chain1, name1, res2, &
                                       ires2, chain2, name2, dmin, dmax, kread
      if ( ioer == 0 ) then
        ntotconst = ntotconst + 1
        cycle
      end if
    end do
    close(10)

  end do

  allocate( iconstraint(ntotconst,2),dconstraint(ntotconst,2), ncpergroup(ncgroups), &
            violation(ntotconst), kconstraint(ntotconst), ncconsider(ncgroups), &
            viol(ntotconst) )

  ! Now reading the constraints

  ntotconst = 0
  ncgroups = 0
  do ifile = 1, nconstfiles

    record = constraint_file(ifile)
    open(10,file=record,status="old",action="read",iostat=ioerr)
    do
      read(10,"( a200 )",iostat=ioer) record
      if ( ioer /= 0 ) exit
      if ( record(1:1) == "#" .or. len(trim(record)) < 1 ) cycle
      if ( keyword(record) == "set" ) then
        ncgroups = ncgroups + 1
        read(record,*) res1, i, res1, ncconsider(ncgroups)
        ncpergroup(ncgroups) = 0
        cycle
      end if
      read(record,*,iostat=ioer) i, j, dmin, dmax, kread
      if ( ioer == 0 ) then
        ntotconst = ntotconst + 1
        iconstraint(ntotconst,1) = i
        iconstraint(ntotconst,2) = j
        dconstraint(ntotconst,1) = dmin**2
        dconstraint(ntotconst,2) = dmax**2
        kconstraint(ntotconst) = kread
        ncpergroup(ncgroups) = ncpergroup(ncgroups) + 1
        cycle
      end if
      read(record,*,iostat=ioer) res1, ires1, chain1, name1, res2, &
                                       ires2, chain2, name2, dmin, dmax, kread
      if ( ioer == 0 ) then
        ntotconst = ntotconst + 1
        i = 1
        do while( i <= natoms )  
          if ( residue(i) == res1 .and. &
               residue_number(i) == ires1 .and. & 
               !chain(i) = chain1 .and. &
               name(i) == name1 ) exit
          i = i + 1
          if ( i > natoms ) then
            write(*,*) ' ERROR: Could not find atom: ', res1, ires1, chain1, name1
            call stop_all()
          end if
        end do
        j = 1
        do while( i <= natoms )
          if ( residue(j) == res2 .and. &
               residue_number(j) == ires2 .and. & 
               !chain(j) = chain2 .and. &
               name(j) == name2 ) exit
          j = j + 1
          if ( j > natoms ) then
            write(*,*) ' ERROR: Could not find atom: ', res2, ires2, chain2, name2
            call stop_all()
          end if
        end do
        iconst = iconst + 1 
        iconstraint(ntotconst,1) = i
        iconstraint(ntotconst,2) = j
        dconstraint(ntotconst,1) = dmin**2
        dconstraint(ntotconst,2) = dmax**2
        kconstraint(ntotconst) = kread 
        ncpergroup(ncgroups) = ncpergroup(ncgroups) + 1
        cycle
      end if
    end do
    close(10)

  end do

  ! Print error if some group of constraints was not defined correctly

  do i = 1, ncgroups
    if ( ncpergroup(i) == 0 ) then 
      write(*,*) ' ERROR: Number of constraints of set ', i, ' is zero. '
      call stop_all()
    end if
    if ( ncconsider(i) == 0 ) then 
      write(*,*) ' WARNING: Not considering any constraint of set ', i
    end if
    if ( ncconsider(i) > ncpergroup(i) ) then
      write(*,*) ' ERROR: Number os constraints to be considered in set ', i,&
                 '        is greater than total number of constraints of this set. '
      call stop_all()
    end if
  end do

  write(*,*) ' Total number of constraints: ', ntotconst
  write(*,*) ' Total number of groups of constraints: ', ncgroups

  if ( ntotconst == 0 ) then
    write(*,*) ' ERROR: Set to use constraints, but no constraint could be read. '
    call stop_all()
  end if

end subroutine read_constraints







