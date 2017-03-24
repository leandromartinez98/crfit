!
! Program createrandom: This is a extract from crfitone only to
! create a initial structure with random coordinates
!
! L. Martinez
!
! Institute of Chemistry - State University of Campinas - Brazil
!

program createrandom

  use force_field
  use flashsortvars

  implicit none

  double precision, parameter :: torad = 3.141592d0/180.d0 
  integer :: narg, i, seed 
  character(len=200) :: pdb_file, pdb_file_out, psf_file

  double precision :: random, phi, psi
  double precision, allocatable :: x(:) 

  narg = iargc()
  if ( narg /= 3 ) then
    write(*,*) ' ERROR: Run with ./createrandom psffile.psf pdbfile.pdb output.pdb '
    stop
  end if
  call getarg(1,psf_file)
  call getarg(2,pdb_file)
  call getarg(3,pdb_file_out)

  ! Readint PSF file and dihedrals

  call readpsf(psf_file)
  call nc_dihedrals()

  allocate( x(natoms*3) ) 

  ! Reading initial coordinates

  call read_coor_pdb(pdb_file,natoms,x)

  ! Initialize random number generator

  call seed_from_time(seed)
  call init_random_number(seed)

  ! Initial guess
  
  do i = 1, nresidues-1
    call random_number(random)
    phi = (-180.d0 + 360.d0*random)*torad
    call rotate_phi(i+1,x,phi)
    call random_number(random)
    psi = (-180.d0 + 360.d0*random)*torad
    call rotate_psi(i,x,psi)
  end do

  call write_pdb(pdb_file,pdb_file_out,x)

end program createrandom
