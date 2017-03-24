!
! Program crfitONE
!
! L. Martinez
!
! Institute of Chemistry - State University of Campinas - Brazil
!

program crfitone

  use force_field
  use constraints
  use linkedcells
  use flashsortvars

  implicit none

  double precision, parameter :: torad = 3.141592d0/180.d0 
  integer :: n_parameter_file, narg, i, j, iresidue, nind, &
             ii, seed, nsurvive, nca, ica, ix, iy, iz, &
             ibest, minoutunit, optpars(10), ioerr, &
             minmethod, print, maxmut, nflash, nshake, maxpop, muttype, ipop
  integer, allocatable :: nctemp(:)
  character(len=200) :: pdb_file, pdb_file_out, parameter_file(10), psf_file, &
                        addbest, record, record2, value, keyword,&
                        inputfile, minoutputfile, read_dihed_file, print_dihed_file
  character(len=200), allocatable :: constraint_file(:)

  double precision :: f, random, phi, psi, fbest, fconstbest1, fconstbest2, ftemp, &
                      probmut, mutbin, ftrue, rmsdmax
  double precision, allocatable :: x(:), g(:), xbest(:), xtmp(:), &
                                   xcatrue(:,:), xca(:,:)
  double precision :: compute_psi, compute_phi 

  logical :: linearshort, error, useinput,&
             combine, mirror, mutations, minimize, pause, testgradient, restart,&
             read_dihed, print_dihed, print_pdb, eval_input_only

  ! These external declarations are required for the gradient checking routine

  external :: computef, computeg 

  write(*,*)
  write(*,"( tr1,80('#') )")
  write(*,*)
  write(*,*) " CRFit - Modelling by LOVO constraint fitting. "
  write(*,*)
  write(*,*) " Institute of Chemistry, State University of Campinas "
  write(*,*)
  write(*,*) " THIS IS THE TRIAL --ONE-- VERSION "
  write(*,*)
  write(*,"( tr1,80('#') )")
  write(*,*)

  !
  ! Default parameters
  !

  write(*,*) ' Setting up default parameters ... '
  print = 0
  nind = 500
  maxmut = 10
  nsurvive = max0(1,int(0.3*nind)) 
  cutoff = 100.d0
  n_parameter_file = 0
  linearshort = .true.
  useinput = .true.
  combine = .true.
  mirror = .true.
  mutations = .true.
  minimize = .true.
  minmethod = 2    ! cgnewton method
  optpars(1) = 500 ! Maximum number of functional evaluations
  optpars(2) = 50  ! Maximum number of CG iterations 
  nshake = 0
  maxpop = -1
  muttype = 1
  probmut = 0.9d0
  minoutunit = 0
  do i = 1, 8
    ffcomp(i) = .true.
  end do
  useroutine = 2
  pause = .false.
  testgradient = .false.
  dielectric = 78.2d0
  restart = .false.
  nconstfiles = 0
  read_dihed = .false.
  print_dihed = .false.
  print_pdb = .false.
  eval_input_only = .false.
  seed = 1431577

  ! 
  ! Read input file 
  ! 

  write(*,*) ' Reading input file name ... '
  error = .false.
  narg = iargc()
  if ( narg == 0 ) then
    write(*,*) ' ERROR: Input file not specified. '
    stop
  end if
  call getarg(1,inputfile)

  write(*,*) ' Reading the input file ... '
  open(10,file=inputfile,action="read",status='old',iostat=ioerr)
  if ( ioerr /= 0 ) then
    write(*,*) ' ERROR: Could not open input file: ', trim(inputfile)
    call stop_all()
  end if
  do
    read(10,"( a200 )",iostat=ioerr) record
    if( ioerr /= 0 ) exit
    if(record(1:1) == "#" .or. len_trim(record) < 1 ) cycle
    select case ( keyword(record) )
      case ("psf")
        psf_file = value(record,1)
      case ("pdb")
        pdb_file = value(record,1)
      case ("output")
        print_pdb = .true.
        pdb_file_out = value(record,1)
      case ("parameters") 
        n_parameter_file = n_parameter_file + 1
        parameter_file(n_parameter_file) = value(record,1)
      case ("bonds") 
        if ( value(record,1) == "yes" ) ffcomp(1) = .true.
        if ( value(record,1) == "no" ) ffcomp(1) = .false.
      case ("angles")
        if ( value(record,1) == "yes" ) ffcomp(2) = .true.
        if ( value(record,1) == "no" ) ffcomp(2) = .false.
      case ("dihedrals")
        if ( value(record,1) == "yes" ) ffcomp(3) = .true.
        if ( value(record,1) == "no" ) ffcomp(3) = .false.
      case ("impropers")
        if ( value(record,1) == "yes" ) ffcomp(4) = .true.
        if ( value(record,1) == "no" ) ffcomp(4) = .false.
      case ("nonbonded")
        if ( value(record,1) == "yes" ) ffcomp(5) = .true.
        if ( value(record,1) == "no" ) ffcomp(5) = .false.
        if ( value(record,1) == "backbone" ) then
          ffcomp(5) = .true.
          ffcomp(7) = .false.
        end if
      case ("sidechains")
        if ( value(record,1) == "neutral" ) ffcomp(8) = .false.
        if ( value(record,1) == "polar" ) ffcomp(8) = .true.
      case ("useconstraints") 
        if ( value(record,1) == "yes" ) ffcomp(6) = .true.
        if ( value(record,1) == "no" ) ffcomp(6) = .false.
      case ("useinput")
        if ( value(record,1) == "yes" ) useinput = .true.
        if ( value(record,1) == "no" ) useinput = .false.
      case ("combine")
        if ( value(record,1) == "yes" ) combine = .true.
        if ( value(record,1) == "no" ) combine = .false.
      case ("mirror")
        if ( value(record,1) == "yes" ) mirror = .true.
        if ( value(record,1) == "no" ) mirror = .false.
      case ("mutations")
        if ( value(record,1) == "yes" ) mutations = .true.
        if ( value(record,1) == "no" ) mutations = .false.
      case ("minimize")
        if ( value(record,1) == "yes" ) minimize = .true.
        if ( value(record,1) == "no" ) minimize = .false.
      case ("pause")
        if ( value(record,1) == "yes" ) pause = .true.
        if ( value(record,1) == "no" ) pause = .false.
      case ("testgradient")
        testgradient = .true.
      case ("eval_input_only")
        eval_input_only = .true.
      case ("restart")
        restart = .true.
      case ("useroutine")
        if ( value(record,1) == "simple" ) useroutine = 1
        if ( value(record,1) == "linearshort" ) useroutine = 2
        if ( value(record,1) == "linkedcell" ) useroutine = 3
      case ("minmethod")
        if ( value(record,1) == "spg" ) minmethod = 1
        if ( value(record,1) == "cgnewton" ) minmethod = 2
      case ("minprint")
        if ( value(record,1) == "no" ) then 
          minoutunit = 0
        else if ( value(record,1) == "screen" ) then
          minoutunit = -1
        else
          minoutunit = 1
          minoutputfile = value(record,1)
        end if
      case ("constraints") 
        nconstfiles = nconstfiles + 1
      case ("print")
        record2 = value(record,1)
        read(record2,*) print
      case ("cutoff")
        record2 = value(record,1)
        read(record2,*) cutoff
      case ("nind")
        record2 = value(record,1)
        read(record2,*) nind
      case ("nshake")
        record2 = value(record,1)
        read(record2,*) nshake
      case ("maxpop")
        record2 = value(record,1)
        read(record2,*) maxpop
      case ("muttype")
        record2 = value(record,1)
        read(record2,*) muttype
      case ("mutbin")
        record2 = value(record,1)
        read(record2,*) mutbin
      case ("seed")
        if ( value(record,1) == 'random' ) then
          call seed_from_time(seed)
        else
          record2 = value(record,1)
          read(record2,*) seed
        end if
      case ("maxmut")
        record2 = value(record,1)
        read(record2,*) maxmut
      case ("probmut")
        record2 = value(record,1)
        read(record2,*) probmut
      case ("nsurvive")
        record2 = value(record,1)
        read(record2,*) nsurvive
      case ("dielectric")
        record2 = value(record,1)
        read(record2,*) dielectric
      case ("maxfeval")
        record2 = value(record,1)
        read(record2,*) optpars(1)
      case ("maxcg")
        record2 = value(record,1)
        read(record2,*) optpars(2)
      case ("read_dihed")
        read_dihed = .true.
        read_dihed_file = value(record,1)
      case ("print_dihed")
        print_dihed = .true.
        print_dihed_file = value(record,1)
      case default
        write(*,*) ' ERROR: Unrecognized keyword: ', keyword(record)
        error = .true.
     end select
  end do
  close(10)

  if ( nind /= 1 ) then
    write(*,*) ' ERROR: This is CRFitONE, and must be run with nind = 1 '
    stop
  end if

  ! Stop if errors were found in the input file

  if ( error ) stop

  ! Print input parameters

  write(*,*)
  write(*,*) ' Force-field options: '
  write(*,*) 
  write(*,*) ' Consider bonds: ', ffcomp(1)
  write(*,*) ' Consider angles: ', ffcomp(2)
  write(*,*) ' Consider dihedrals: ', ffcomp(3)
  write(*,*) ' Consider impropers: ', ffcomp(4)
  write(*,*) ' Consider non-bonded interactions: ', ffcomp(5)
  write(*,*) ' -- for all atoms: ', ffcomp(7)
  write(*,*) ' Use constraints: ', ffcomp(6)
  write(*,*) ' Side chains are polar: ', ffcomp(8)
  write(*,*) ' Linear short-range modification: ', linearshort
  write(*,*) ' Long-range cutoff: ', cutoff
  write(*,*) ' Number of individuals of genetic algorithm: ', nind
  cutoff2 = cutoff**2

  if ( useroutine == 1 ) then
    write(*,*) ' Using standard non-bonded force field computation. '
  else if ( useroutine == 2 ) then
    write(*,*) ' Using non-bonded force field with linearization at short distances. '
  else if ( useroutine == 3 ) then
    write(*,*) ' Using linked cell method with linearization at short distances. '
  end if

  write(*,*) ' Reading PSF file: ', trim(psf_file)
  call readpsf(psf_file)
  write(*,*)
  write(*,*) ' Structure summary: '
  write(*,*)  natoms, ' atoms. '
  write(*,*)  nresidues, ' residues. '
  write(*,*)  nbonds, ' bonds. '
  write(*,*)  nangles, ' angles. '
  write(*,*)  ndihed, ' dihedral angles. '
  write(*,*)  nimpr, ' improper dihedral angles. '

  ! Compute order of atoms for dihedral rotations

  call nc_dihedrals()

  allocate( x(natoms*3), g(natoms*3) ) 

  ! Reading parameters, initial coordinates, and constraints

  call readprm(parameter_file,n_parameter_file)

  write(*,*) 
  write(*,*) ' Reading coordinates from PDB file: ', trim(pdb_file)
  call read_coor_pdb(pdb_file,natoms,x)

  ! Read constraint data

  if ( ffcomp(6) ) then

    if ( nconstfiles > 0 ) then
      allocate( constraint_file(nconstfiles) )
      open(10,file=inputfile,action='read') 
      nconstfiles = 0
      do 
        read(10,"( a200 )",iostat=ioerr) record
        if ( ioerr /= 0 ) exit
        if ( keyword(record) == "constraints" ) then
          nconstfiles = nconstfiles + 1
          constraint_file(nconstfiles) = value(record,1)
        end if
      end do
      close(10)
      write(*,"( /, '  Reading constraint files ... ' )")
      write(*,*) ' -- Number of files: ', nconstfiles
      call read_constraints(constraint_file)
    else
      write(*,*) ' ERROR: Set to use constraints, but no constraint file is defined. '
      stop
    end if

  else

     ncgroups = 0
     ntotconst = 0

  end if

  ! Computing the number of CA atoms

  nca = 0
  do i = 1, natoms
    if ( name(i) == 'CA' ) nca = nca + 1
  end do
  allocate( xcatrue(nca,3), xca(nca,3) )

  ! For flashsort (used for sorting individuals here, and constraints in compute_fconst)

  nflash = max(ntotconst,nind)
  allocate( indflash(nflash), lflash(nflash), indflashback(nflash) )
  mflash = 1 + nflash / 10

  ! Print constraint parameters

  do i = 1, ncgroups
    write(*,*) ' Constraints of group: ', i
    write(*,*) ' -- Number of constraints in this group: ', ncpergroup(i)
    write(*,*) ' -- Number of constraints for LOVO evaluation: ', ncconsider(i)
  end do

  ! Build non-bonded pair list
                             
  call build_exclusions(error)
  if ( error ) stop

  ! If the side chains are not polar, update the linear extrapolation of 
  ! short distance interactions
  
  if ( .not. ffcomp(8) ) call updatels(error)

  ! 
  ! Write properties of the setup of the genetic algorithm
  ! 

  write(*,*)
  write(*,*) ' Genetic algorithm options: '
  write(*,*)
  write(*,*) ' Use parents to build new individuals (combine): ', combine
  write(*,*) ' Mirror individuals sometimes: ', mirror
  write(*,*) ' Perform mutations: ', mutations
  write(*,*) ' Stop at maximum number of populations: ', maxpop
  write(*,*) ' Maximum number of mutations per individual: ', maxmut
  write(*,*) ' Probability of mutating an individual: ', probmut
  write(*,*) ' Type of mutation (0: random; 1: Ramachandran): ', muttype
  write(*,*) ' Maximum size of the perturbation (degrees): ', mutbin*360.

  ! Allocate local vectors that are required for the GA evaluation

  allocate( xbest(natoms*3), nctemp(ncgroups), xtmp(natoms*3) )

  ! Define how the optimization method will output information

  if ( minoutunit > 0 ) open(minoutunit,file=minoutputfile) 

  ! Compute the function value for input structure

  ii = useroutine
  useroutine = 1
  call computef(x,f)
  useroutine = ii
  write(*,*) 
  write(*,*) ' Energy of input structure = ', f
  write(*,*) ' ---- Bond energy: ', ecomp(1)
  write(*,*) ' ---- Angle energy: ', ecomp(2)
  write(*,*) ' ---- Dihedral angle energy: ', ecomp(3)
  write(*,*) ' ---- Improprer dihedral angle energy: ', ecomp(4)
  write(*,*) ' ---- VDW non-bonded energy = ', ecomp(5)
  write(*,*) ' ---- Electrostatic non-bonded energy = ', ecomp(6)
  write(*,*) ' ---- Constraint violation energy: ', ecomp(7)
  if ( useroutine /= 1 ) then 
    call computef(x,f)
    write(*,*) ' ---- Function value: ', f
    write(*,*) ' ---- Non-bonded function value: ', ecomp(8)
  end if

  ! Just evaluate the input structure and exit

  if ( eval_input_only ) then
    if ( print_dihed ) then
      open(10,file=print_dihed_file)
      do i = 1, nresidues-1
        write(10,*) compute_phi(i+1,x), compute_psi(i,x)
      end do
      close(10)
    end if
    call stop_all()
  end if

  ! Write the input structure to the the trajectory file (first write)

  if ( print_pdb .and. .not. restart ) then
    call system("rm -f ./traj.pdb")
    call write_pdb(pdb_file,addbest(pdb_file_out),x)
    call system("cat ./testeBEST.pdb >> ./traj.pdb")
  end if

  ! Save the initial structure as the reference structure

  ftrue = f
  ica = 0 
  do i = 1, natoms
    ix = (i-1)*3 + 1
    iy = ix + 1
    iz = ix + 2
    if ( name(i) == 'CA' ) then
      ica = ica + 1
      xcatrue(ica,1) = x(ix)
      xcatrue(ica,2) = x(iy)
      xcatrue(ica,3) = x(iz)
    end if
  end do

  ! Initialize random number generator

  write(*,*) ' Seed for random number generator: ', seed
  call init_random_number(seed)

  ! Maximum number of functional evaluations of optimization method

  write(*,*) 
  write(*,*) ' Maximum number of function evaluations of minimization method: ', optpars(1)
  write(*,*) 

  ! Initial guess
  
  if ( .not. read_dihed ) then
    do i = 1, nresidues-1
      call random_number(random)
      phi = (-180.d0 + 360.d0*random)*torad
      call rotate_phi(i+1,x,phi)
      call random_number(random)
      psi = (-180.d0 + 360.d0*random)*torad
      call rotate_psi(i,x,psi)
    end do
  else if ( read_dihed ) then
    write(*,*) ' Reading initial dihedrals from input file: ', trim(read_dihed_file)
    open(10,file=read_dihed_file,status='old',action='read',iostat=ioerr)
    if ( ioerr /= 0 ) then
      write(*,*) ' ERROR: Could not open or read dihedrals input file. ' 
      call stop_all()
    end if
    do i = 1, nresidues-1
      read(10,*,iostat=ioerr) phi, psi
      phi = -1.d0*compute_phi(i+1,x) + phi
      psi = -1.d0*compute_psi(i,x) + psi
      call rotate_phi(i+1,x,phi)
      call rotate_psi(i,x,psi)
    end do
    close(10)
  end if

  !
  ! Start genetic algorithm
  !

  write(*,*) ' Starting genetic algorithm ... '
  write(*,*) ' Number of surviving individuals: ', nsurvive

  ibest = 0
  do i = 1, natoms*3
    xbest(i) = x(i)
  end do
  call computef(xbest,fbest)
  call compute_fconst(xbest,fconstbest1)
  do i = 1, ncgroups
    nctemp(i) = ncconsider(i)
    ncconsider(i) = ncpergroup(i)
  end do
  call compute_fconst(xbest,fconstbest2)
  do i = 1, ncgroups
    ncconsider(i) = nctemp(i)
  end do

  if ( print_pdb .and. .not. restart ) then
    call write_pdb(pdb_file,addbest(pdb_file_out),x)
    call system("cat ./testeBEST.pdb >> ./traj.pdb")
  end if

  ipop = 0
  genetic : do while( ipop < maxpop .or. maxpop == -1 ) 
    ipop = ipop + 1

    write(*,"( 90('-') )")
    write(*,*) ' Structure: ', ipop
    write(*,"( a, tr2, e17.10, a, i4 )") '  Best score up to now: ', &
                                           fbest, ' from individual ', ibest
    write(*,"( a, f17.5 )") '  Constraint violation of best structure (for LOVO) = ', fconstbest1
    write(*,"( a, f17.5 )") '  Constraint violation of best structure (for all constraints) = ', fconstbest2
    write(*,"( 90('-') )")

    !
    ! Rotate one dihedral randomly
    !
    
    do j = 1, natoms*3
      x(j) = xbest(j)
    end do
    call random_number(random)
    if ( random <= probmut ) then
      call random_number(random)
      iresidue = int(random*(nresidues-1))
      write(*,*) ' Rotating dihedral of residue ', iresidue
      if ( muttype == 0 ) then
        call random_number(random)
        ! Phi
        if ( random < 0.5d0 ) then
          phi = compute_phi(iresidue+1,xbest)
          call random_number(random)
          phi = phi + mutbin*(-180. + 360.*random)*torad
          call rotate_phi(iresidue+1,x,phi)
        ! Psi
        else
          psi = compute_psi(iresidue,xbest)
          call random_number(random)
          psi = psi + mutbin*(-180. + 360.*random)*torad
          call rotate_psi(iresidue,x,psi)
        end if
      else if ( muttype == 1 ) then
        call newangles(phi,psi)
        phi = phi - compute_phi(iresidue+1,xbest)
        psi = psi - compute_psi(iresidue,xbest)
        call rotate_phi(iresidue+1,x,phi)
        call rotate_psi(iresidue,x,psi)
      end if
    end if

    ! Minimize the energy of this individual

    if( minimize ) then 

      call computef(x,f)
    
      ! Optimize structure received 

      call optimize(natoms*3,x,f,minoutunit,optpars,minmethod)

      ! Try to escape from trivial local minima by performing small
      ! perturbations (0.5 Angstrom perturbations on the coordinates of each
      ! atom)

      ftemp = f
      do i = 1, nshake
        do j = 1, natoms*3
          xtmp(j) = x(j)
          call random_number(random)
          xtmp(j) = xtmp(j) - 0.5d0 + random
        end do
        call optimize(natoms*3,xtmp,f,minoutunit,optpars,minmethod)
        if ( f < ftemp ) then
          ftemp = f
          do j = 1, natoms*3
            x(j) = xtmp(j)
          end do
        end if
      end do

    end if

    ! Compute the energy of this individual

    call computef(x,f)

    ! Save if best structure so far

    if ( f < fbest ) then
      fbest = f
      do i = 1, natoms*3
        xbest(i) = x(i)
      end do
      if ( print_pdb ) then
        call write_pdb(pdb_file,addbest(pdb_file_out),xbest)
        call system("cat ./testeBEST.pdb >> ./traj.pdb")
      end if
      call compute_fconst(xbest,fconstbest1)
      do i = 1, ncgroups
        nctemp(i) = ncconsider(i)
        ncconsider(i) = ncpergroup(i)
      end do
      call compute_fconst(xbest,fconstbest2)
      do i = 1, ncgroups
        ncconsider(i) = nctemp(i)
      end do
      ibest = ipop

      ! Check if the structure is close enough to the solution

      ica = 0 
      do i = 1, natoms
        ix = (i-1)*3 + 1
        iy = ix + 1
        iz = ix + 2
        if ( name(i) == 'CA' ) then
          ica = ica + 1
          xca(ica,1) = x(ix)
          xca(ica,2) = x(iy)
          xca(ica,3) = x(iz)
        end if
      end do
      call align(nca,xca,xcatrue)
      rmsdmax = 0.d0
      do i = 1, nca
        rmsdmax = max( rmsdmax, ((xca(i,1)-xcatrue(i,1))**2 + &
                                 (xca(i,2)-xcatrue(i,2))**2 + &
                                 (xca(i,3)-xcatrue(i,3))**2) )
      end do
      rmsdmax = dsqrt(rmsdmax)
      write(*,*) ' Maximum error of CA coordinates in new best model: ', rmsdmax
      write(*,*) ' Energy relative to correct model: ', (f/ftrue)*100.d0, '%'
      if ( rmsdmax < 0.5 .and. f < 1.5d0*ftrue ) then
        write(*,*) ' The solution was found. '
        call stop_all()
      end if

    end if
 
  end do genetic
  
  ! Print file containing the final dihedrals

  if ( print_dihed ) then
    write(*,*) ' Writing final dihedrals to file: ', trim(print_dihed_file)
    open(10,file=print_dihed_file)
    do i = 1, nresidues-1
      write(10,*) compute_phi(i+1,x), compute_psi(i,x)
    end do
    close(10)
  end if

  write(*,*) ' Maximum number of populations reached. Stopping. '
  write(*,*) ' Best function value = ', fbest, ' Constraint violation = ', fconstbest1 

end program crfitone
