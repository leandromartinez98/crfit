!
! Subroutine read_coor_pdb : Read coordinates from PDB file
!
! L. Martinez, Nov 28, 2013
! Institute of Chemistry - State University of Campinas - Brazil
!

subroutine read_coor_pdb(pdb_file,natoms,x)

  implicit none
  integer :: natoms, icount, ioerr, i
  character(len=200) :: pdb_file, line
  double precision :: x(*)
  logical :: file_exist

  inquire(file=pdb_file,exist=file_exist)
  if ( .not. file_exist ) then
    write(*,*) ' ERROR: Could not open (find?) pdb file: ', trim(pdb_file)
    call stop_all()
  end if
  open(10,file=pdb_file,status='old',action='read',iostat=ioerr)
  icount = 0
  do
    read(10,"( a200 )",iostat=ioerr) line
    if ( ioerr /= 0 ) exit
    if ( line(1:4) == "ATOM" .or. line(1:6) == "HETATM" ) then
      icount = icount + 1
      read(line(31:38),*,iostat=ioerr) x((icount-1)*3+1)
      if ( ioerr /= 0 ) call pdb_read_error(i)
      read(line(39:46),*,iostat=ioerr) x((icount-1)*3+2)
      if ( ioerr /= 0 ) call pdb_read_error(i)
      read(line(47:54),*,iostat=ioerr) x((icount-1)*3+3)
      if ( ioerr /= 0 ) call pdb_read_error(i)
    end if
  end do
  close(10)

  if ( icount /= natoms ) then 
    write(*,*) ' ERROR: Number of atoms found in PDB file is different from '
    write(*,*) '        the number of atoms set in PSF file. '
    call stop_all()
  end if
  
  return
end subroutine read_coor_pdb

!
! Prints error message for wrong read and stop
!

subroutine pdb_read_error(icount)
  integer :: icount
  write(*,*) ' ERROR reading coordinates from PDB file, atom: ', icount
  call stop_all()
end subroutine pdb_read_error
  
