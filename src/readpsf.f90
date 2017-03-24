!
! Subroutine readpsf: Read all that from the psf file and allocates
!                     some arrays of the force_field module.
!
! L. Martinez, Nov 11, 2013
! Institute of Chemistry - State University of Campinas - Brazil
!

subroutine readpsf(psf_file)

  use force_field
  implicit none
  integer :: itemp, i, j, ii, iline, icount, ioerr, iiread(50), &
             iicount
  character(len=200) :: psf_file, record, keyword
  
  open(10,file=psf_file,action='read',status='old',iostat=ioerr)
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Could not open file: ', trim(psf_file)
    call stop_all()
  end if
  
  ! Get the number of atoms, bonds, etc, from the PSF file headers
  
  natoms = 0
  do 
    read(10,"( a200 )",iostat=ioerr) record 
    if ( ioerr /= 0 ) exit
  
    read(record,*,iostat=ioerr) itemp, keyword
    if ( ioerr /= 0 ) cycle
  
    if ( keyword == "!NATOM" ) natoms = itemp 
    if ( keyword == "!NBOND:" ) nbonds = itemp
    if ( keyword == "!NTHETA:" ) nangles = itemp
    if ( keyword == "!NPHI:" ) ndihed = itemp
    if ( keyword == "!NIMPHI:" ) nimpr = itemp
    if ( keyword == "!NCRTERM:" ) ncrterm = itemp
  
  end do

  if ( natoms == 0 ) then
    write(*,*) ' ERROR: Flag !NATOM not found in PSF file. '
    write(*,*) '        Could not read number of atoms. '
    call stop_all()
  end if
  
  ! Allocate all arrays
  
  allocate( & 
            segment(natoms), residue(natoms), name(natoms), type(natoms), &
            residue_number(natoms), charge(natoms), mass(natoms), &
            eps(natoms), sig(natoms), eps14(natoms), sig14(natoms), &
            ibond(nbonds,2), fbond(nbonds,2), &
            iangle(nangles,3), fangle(nangles,2), &  
            idihed(ndihed,6), &
            iimpr(nimpr,6), &
            icrterm(ncrterm,6), &
            iorder(natoms), &
            isbackbone(natoms) &
          )
  
  ! Read PSF file data
  
  rewind(10)
  do 
  
    read(10,"( a200 )",iostat=ioerr) record 
    iline = iline + 1
    if ( ioerr /= 0 ) exit
  
    read(record,*,iostat=ioerr) itemp, keyword
    if ( ioerr /= 0 ) cycle
  
    ! Read atom parameters
  
    if ( keyword == "!NATOM" ) then 
      do i = 1, natoms
        read(10,"( a200 )",iostat=ioerr) record
        iline = iline + 1
        if ( ioerr /= 0 ) call ioerror(1,record,iline)
        read(record,*,iostat=ioerr) itemp, segment(i), residue_number(i), residue(i), &
                                    name(i), type(i), charge(i), mass(i) 
        if ( ioerr /= 0 ) call ioerror(1,record,iline)
        if ( name(i) == 'HN' .or. &
             !name(i) == 'HT1' .or. &
             !name(i) == 'HT2' .or. &
             !name(i) == 'HT3' .or. &
             !( name(i) == 'HN1' .and. residue(i) == 'PRO' ) .or. &
             !( name(i) == 'HN2' .and. residue(i) == 'PRO' ) .or. &
             name(i) == 'N' .or. &
             name(i) == 'HA' .or. &
             ( name(i) == 'HA1' .and. residue(i) == 'GLY' ) .or. &
             ( name(i) == 'HA2' .and. residue(i) == 'GLY' ) .or. &
             name(i) == 'CA' .or. &
             name(i) == 'C' .or. & 
             !name(i) == 'OT1' .or. &
             !name(i) == 'OT2' .or. &
             name(i) == 'O' ) then
          isbackbone(i) = .true.
        else
          isbackbone(i) = .false.
        end if
      end do
    end if
  
    ! Read bonds
  
    if ( keyword == "!NBOND:" ) then
      icount = 0
      iicount = 0
      do while( (icount / 2) < nbonds ) 
        read(10,"( a200 )",iostat=ioerr) record
        iline = iline + 1
        if ( ioerr /= 0 ) call ioerror(1,record,iline)
  
        ioerr = 0
        j = 0
        do while( ioerr == 0 )
          j = j + 1
          read(record,*,iostat=ioerr) (itemp, i = 1, j)
        end do
        j = j - 1
        read(record,*) (iiread(i), i = 1, j)
  
        do i = 1, j
          icount = icount + 1
          iicount = iicount + 1
          ibond((icount-iicount)/2+1,iicount) = iiread(i)
          if ( iicount == 2 ) iicount = 0
        end do
  
      end do
      if ( icount/2 /= nbonds ) then
        write(*,*) ' Error reading the index of atoms in bonds, the number '
        write(*,*) ' of bonds found is not equal to the number of bonds '
        write(*,*) ' specified in the PSF header. '
        call stop_all()
      end if
    end if
  
    ! Read angles
  
    if ( keyword == "!NTHETA:" ) then
      icount = 0     
      iicount = 0
      do while( (icount / 3) < nangles ) 
        read(10,"( a200 )",iostat=ioerr) record
        iline = iline + 1
        if ( ioerr /= 0 ) call ioerror(1,record,iline)
  
        ioerr = 0
        j = 0
        do while( ioerr == 0 )
          j = j + 1
          read(record,*,iostat=ioerr) (itemp, i = 1, j)
        end do
        j = j - 1
        read(record,*) (iiread(i), i = 1, j)
  
        do i = 1, j
          icount = icount + 1
          iicount = iicount + 1
          iangle((icount-iicount)/3+1,iicount) = iiread(i)
          if( iicount == 3 ) iicount = 0
        end do
  
      end do
      if ( icount/3 /= nangles ) then
        write(*,*) ' Error reading the index of atoms in angles, the number '
        write(*,*) ' of angles found is not equal to the number of angles '
        write(*,*) ' specified in the PSF header. '
        call stop_all()
      end if
    end if
  
    ! Read dihedrals
  
    if ( keyword == "!NPHI:" ) then
      icount = 0     
      iicount = 0
      do while( (icount / 4) < ndihed ) 
        read(10,"( a200 )",iostat=ioerr) record
        iline = iline + 1
        if ( ioerr /= 0 ) call ioerror(1,record,iline)
  
        ioerr = 0
        j = 0
        do while( ioerr == 0 )
          j = j + 1
          read(record,*,iostat=ioerr) (itemp, i = 1, j)
        end do
        j = j - 1
        read(record,*) (iiread(i), i = 1, j)
  
        do i = 1, j
          icount = icount + 1
          iicount = iicount + 1
          idihed((icount-iicount)/4+1,iicount) = iiread(i)
          if( iicount == 4 ) iicount = 0
        end do
  
      end do
      if ( icount/4 /= ndihed ) then
        write(*,*) ' Error reading the index of atoms in dihedrals, the number '
        write(*,*) ' of dihedrals found is not equal to the number of dihedrals '
        write(*,*) ' specified in the PSF header. '
        call stop_all()
      end if
    end if
  
    ! Read impropers
  
    if ( keyword == "!NIMPHI:" ) then
      icount = 0     
      iicount = 0
      do while( (icount / 4) < nimpr ) 
        read(10,"( a200 )",iostat=ioerr) record
        iline = iline + 1
        if ( ioerr /= 0 ) call ioerror(1,record,iline)
  
        ioerr = 0
        j = 0
        do while( ioerr == 0 )
          j = j + 1
          read(record,*,iostat=ioerr) (itemp, i = 1, j)
        end do
        j = j - 1
        read(record,*) (iiread(i), i = 1, j)
  
        do i = 1, j
          icount = icount + 1
          iicount = iicount + 1
          iimpr((icount-iicount)/4+1,iicount) = iiread(i)
          if( iicount == 4 ) iicount = 0
        end do
  
      end do
      if ( icount/4 /= nimpr ) then
        write(*,*) ' Error reading the index of atoms in improper dihedrals, the number '
        write(*,*) ' of improper dihedrals found is not equal to the number of improper dihedrals '
        write(*,*) ' specified in the PSF header. '
        call stop_all()
      end if
    end if
  
    ! Read cross-terms
  
!    if ( keyword == "!NCRTERM:" ) then
!      icount = 0     
!      iicount = 0
!      do while( (icount / 4) < ncrterm ) 
!        read(10,"( a200 )",iostat=ioerr) record
!        iline = iline + 1
!        if ( ioerr /= 0 ) call ioerror(1,record,iline)
!  
!        ioerr = 0
!        j = 0
!        do while( ioerr == 0 )
!          j = j + 1
!          read(record,*,iostat=ioerr) (itemp, i = 1, j)
!        end do
!        j = j - 1
!        read(record,*) (iiread(i), i = 1, j)
!  
!        do i = 1, j
!          icount = icount + 1
!          iicount = iicount + 1
!          icrterm((icount-iicount)/4+1,iicount) = iiread(i)
!          if( iicount == 4 ) iicount = 0
!        end do
!  
!      end do
!      if ( icount/4 /= ncrterm ) then
!        write(*,*) ' Error reading the index of atoms in cross-terms, the number '
!        write(*,*) ' of cross-terms found is not equal to the number of cross-terms'
!        write(*,*) ' specified in the PSF header. '
!        call stop_all()
!      end if
!    end if
  
  end do
  close(10)

  nresidues = residue_number(natoms)-residue_number(1) + 1

  ! Count the number of atoms per residue

  allocate( natres(nresidues) )
  ii = 0
  do i = 1, natoms
    if ( residue_number(i)-residue_number(1)+1 /= ii ) then
      ii = ii + 1
      natres(ii) = 0
    end if
    natres(ii) = natres(ii) + 1
  end do

  ! Compute order of atoms which is practical for the rotation of phi
  ! psi dihedrals

  icount = 0
  iicount = 0

  ! First residue

  ii = 0
  do j = 1, natres(1)
    icount = icount + 1
    if ( name(icount) == 'CA' ) then
      iorder(icount) = iicount + natres(1) - 2
    else if ( name(icount) == 'C' ) then
      iorder(icount) = iicount + natres(1) - 1 
    else if ( name(icount) == 'O' ) then
      iorder(icount) = iicount + natres(1)
    else                                    
      ii = ii + 1
      iorder(icount) = iicount + ii
    end if
  end do
  iicount = iicount + natres(1)

  ! Residues in the midle

  do i = 2, nresidues - 1
    ii = 2           
    do j = 1, natres(i)
      icount = icount + 1
      if ( name(icount) == 'HN' ) then 
        iorder(icount) = iicount + 1
      else if ( name(icount) == 'N' ) then
        iorder(icount) = iicount + 2
      else if ( name(icount) == 'CA' ) then
        iorder(icount) = iicount + natres(i) - 2
      else if ( name(icount) == 'C' ) then
        iorder(icount) = iicount + natres(i) - 1
      else if ( name(icount) == 'O' ) then
        iorder(icount) = iicount + natres(i) 
      else                                    
        ii = ii + 1
        iorder(icount) = iicount + ii
      end if
    end do
    iicount = iicount + natres(i)
  end do

  ! Last residue

  ii = 2
  do j = 1, natres(nresidues)
    icount = icount + 1
    if ( name(icount) == 'HN' ) then
      iorder(icount) = iicount + 1
    else if ( name(icount) == 'N' ) then
      iorder(icount) = iicount + 2
    else                                    
      ii = ii + 1
      iorder(icount) = iicount + ii
    end if
  end do

  return  
end subroutine readpsf

!
! Suborutine that prints error messages on input or output
!

subroutine ioerror(ierror,record,iline)

  implicit none
  integer :: iline, ierror
  character(len=200) :: record
  
  if ( ierror == 1 ) then
    write(*,*) ' Error reading PSF file at line ', iline
    write(*,*) ' Content: ', trim(record(1:len(record)))
  end if
  call stop_all()

end subroutine ioerror

