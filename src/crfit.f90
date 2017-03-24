!
! Program crfit
!
! L. Martinez
!
! Institute of Chemistry - State University of Campinas - Brazil
!

program crfit

  use force_field
  use constraints
  use linkedcells
  use flashsortvars
  use mpi

  implicit none

  double precision, parameter :: torad = 3.141592d0/180.d0 
  integer :: n_parameter_file, narg, i, j, iresidue, iind, nind, irecv, &
             ii, seed, icover, jcover, igen, imut, iparent, nmut, nsurvive, &
             ibest, minoutunit, optpars(10), ioerr, nchange, iprint, iprintlast, &
             minmethod, print, maxmut, nflash, nshake, maxpop, muttype, ix, iy, iz
  integer, allocatable :: nctemp(:)
  character(len=200) :: pdb_file, pdb_file_out, parameter_file(10), psf_file, &
                        addbest, record, record2, value, keyword,&
                        inputfile, minoutputfile
  character(len=200), allocatable :: constraint_file(:)

  double precision :: f, random, phi, psi, fbest, fconstbest1, fconstbest2, frel, ftemp, score,&
                      probmut, mutbin, dist, rmsd, rmsdmax
  double precision, allocatable :: x(:), g(:), xbest(:), xind(:,:), xdihed(:,:,:), xdihedbest(:,:),&
                                   find(:), xref(:), xdihedtemp(:,:), xdihednew(:,:), fsort(:), xtmp(:),&
                                   xdihedinput(:,:), xali(:,:), xaliref(:,:)
  double precision :: compute_psi, compute_phi

  logical :: linearshort, terminate, ffchange, error, useinput,&
             combine, mirror, mutations, minimize, pause, testgradient, restart,&   
             writepop

  ! This are the variables which are specific for the MPI code

  integer :: ierr, rank, myrank, nrank, tag
  integer :: status(mpi_status_size)
  logical, allocatable :: free(:)

  ! These external declarations are required for the gradient checking routine

  external :: computef, computeg 

  ! 
  ! Initializing MPI
  ! 

  call mpi_init(ierr)
  call mpi_comm_rank(mpi_comm_world,myrank,ierr)
  call mpi_comm_size(mpi_comm_world,nrank,ierr)

  !
  ! Activities of the master rank (rank = 0)
  ! 

  if ( myrank == 0 ) then

    write(*,*)
    write(*,"( tr1,80('#') )")
    write(*,*)
    write(*,*) " CRFit - Modelling by LOVO constraint fitting. "
    write(*,*)
    write(*,*) " Institute of Chemistry, State University of Campinas "
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
    nshake = 5
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
    seed = 1431577
    writepop = .false.

    ! 
    ! Read input file 
    ! 

    write(*,*) ' Reading input file name ... '
    error = .false.
    terminate = .false.
    narg = iargc()
    if ( narg == 0 ) then
      write(*,*) ' ERROR: Input file not specified. '
      call stop_all()
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
        case ("restart")
          restart = .true.
        case ("writepop")
          writepop = .true.
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
        case ("terminate")
          write(*,*) ' ERROR: Starting with termination signal. '
          error = .true.
        case default
          write(*,*) ' ERROR: Unrecognized keyword: ', keyword(record)
          error = .true.
       end select
    end do
    close(10)

    ! Check if the number of processes is smaller than the number of individuals

    if ( nrank > nind + 1 .or. nrank < 2 ) then
      write(*,*) ' ERROR: The number of processes (ranks) must be AT MOST the '
      write(*,*) '        the number of individuals plus one, and AT LEAST two. '
      error = .true.
    end if

    ! Stop if errors were found in the input file

    if ( error ) call stop_all()

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

    ! Reading PSF file

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

    ! Compute order of atoms for dihedral rotation

    call nc_dihedrals()

    ! 
    ! Vectors that must be allocated for all ranks
    !

    allocate( x(natoms*3), g(natoms*3), xref(natoms*3) ) 

    ! Reading parameters, initial coordinates, and constraints

    call readprm(parameter_file,n_parameter_file)
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
        call stop_all()
      end if

    else

       ncgroups = 0
       ntotconst = 0

    end if

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
                               
    call build_exclusions()

    ! If the side chains are not polar, update the linear extrapolation of 
    ! short distance interactions
    
    if ( .not. ffcomp(8) ) call updatels()

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

    !
    ! Broadcast force-field parameters to the slave processes
    !

    write(*,*) 
    write(*,*) ' Number of processes: ', nrank
    write(*,*) ' Broadcasting force field data to all slave ranks ... '

    do rank = 1, nrank-1
      call mpi_send( error, 1, mpi_logical, rank, tag, mpi_comm_world, ierr )
    end do
    call mpi_bcast( natoms, 1, mpi_integer, 0, mpi_comm_world, ierr )
    call mpi_bcast( isbackbone, natoms, mpi_logical, 0, mpi_comm_world, ierr )

    call mpi_bcast( n_nonbonded, 1, mpi_integer, 0, mpi_comm_world, ierr )
    call mpi_bcast( ijnonbonded, n_nonbonded*2, mpi_integer, 0, mpi_comm_world, ierr )
    call mpi_bcast( qq, n_nonbonded, mpi_double_precision, 0, mpi_comm_world, ierr )
    call mpi_bcast( epseps, n_nonbonded, mpi_double_precision, 0, mpi_comm_world, ierr )
    call mpi_bcast( ss6, n_nonbonded, mpi_double_precision, 0, mpi_comm_world, ierr )
    call mpi_bcast( cutnb, n_nonbonded*3, mpi_double_precision, 0, mpi_comm_world, ierr )
    call mpi_bcast( nshake, 1, mpi_integer, 0, mpi_comm_world, ierr )
    call mpi_bcast( linearshort, 1, mpi_logical, 0, mpi_comm_world, ierr )
    call mpi_bcast( n_nbbackbone, 1, mpi_integer, 0, mpi_comm_world, ierr ) 
    call mpi_bcast( inbbackbone, n_nbbackbone, mpi_integer, 0, mpi_comm_world, ierr )

    call mpi_bcast( nbonds, 1, mpi_integer, 0, mpi_comm_world, ierr )
    call mpi_bcast( ibond, nbonds*2, mpi_integer, 0, mpi_comm_world, ierr )
    call mpi_bcast( fbond, nbonds*2, mpi_double_precision, 0, mpi_comm_world, ierr )

    call mpi_bcast( nangles, 1, mpi_integer, 0, mpi_comm_world, ierr )
    call mpi_bcast( iangle, nangles*3, mpi_integer , 0, mpi_comm_world, ierr )
    call mpi_bcast( fangle, nangles*2, mpi_double_precision, 0, mpi_comm_world, ierr )

    call mpi_bcast( nureybradley, 1, mpi_integer, 0, mpi_comm_world, ierr )
    call mpi_bcast( iureybradley, nureybradley*2, mpi_integer, 0, mpi_comm_world, ierr )
    call mpi_bcast( fureybradley, nureybradley*2, mpi_double_precision, 0, mpi_comm_world, ierr )

    call mpi_bcast( ndihed, 1, mpi_integer, 0, mpi_comm_world, ierr )
    call mpi_bcast( nfdihed, 1, mpi_integer, 0, mpi_comm_world, ierr )
    call mpi_bcast( idihed, ndihed*6, mpi_integer, 0, mpi_comm_world, ierr )
    call mpi_bcast( fdihed, nfdihed*3, mpi_double_precision, 0, mpi_comm_world, ierr )

    call mpi_bcast( nimpr, 1, mpi_integer, 0, mpi_comm_world, ierr )
    call mpi_bcast( nfimpr, 1, mpi_integer, 0, mpi_comm_world, ierr )
    call mpi_bcast( iimpr, nimpr*6, mpi_integer, 0, mpi_comm_world, ierr )
    call mpi_bcast( fimpr, nfimpr*2, mpi_double_precision, 0, mpi_comm_world, ierr ) 
  
    call mpi_bcast( ncgroups, 1, mpi_integer, 0, mpi_comm_world, ierr )
    call mpi_bcast( ncpergroup, ncgroups, mpi_integer, 0, mpi_comm_world, ierr )
    call mpi_bcast( ncconsider, ncgroups, mpi_integer, 0, mpi_comm_world, ierr )

    call mpi_bcast( ntotconst, 1, mpi_integer, 0, mpi_comm_world, ierr )
    call mpi_bcast( kconstraint, ntotconst, mpi_double_precision, 0, mpi_comm_world, ierr )
    call mpi_bcast( iconstraint, ntotconst*2, mpi_integer, 0, mpi_comm_world, ierr )
    call mpi_bcast( dconstraint, ntotconst*2, mpi_double_precision, 0, mpi_comm_world, ierr )

    call mpi_bcast( minoutunit, 1, mpi_integer, 0, mpi_comm_world, ierr )
    call mpi_bcast( useroutine, 1, mpi_integer, 0, mpi_comm_world, ierr )
    call mpi_bcast( maxboxes, 1, mpi_integer, 0, mpi_comm_world, ierr )
    call mpi_bcast( knonbonded, natoms*natoms, mpi_integer, 0, mpi_comm_world, ierr ) 

    ! Allocate local vectors that are required for the GA evaluation

    allocate( xbest(natoms*3), xind(nind,natoms*3), xdihed(nind,nresidues,2), &
              xdihedbest(nresidues,2), find(nind), xdihedinput(nresidues,2), &
              xdihedtemp(nresidues,2), fsort(nind), xtmp(natoms*3), &
              free(nrank-1), xdihednew(nresidues,2), nctemp(ncgroups), &
              xali(natoms,3), xaliref(natoms,3) )

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

    ! Write the input structure to the the trajectory file (first write)

    if ( .not. restart ) then
      call system("rm -f ./traj.pdb")
      call write_pdb(pdb_file,addbest(pdb_file_out),x)
      call system("cat ./testeBEST.pdb >> ./traj.pdb")
    end if

    ! Perform some local mimimization on the input structure

    write(*,*) 
    write(*,*) ' Minimizing the energy of the input structure ... '
    call optimize(natoms*3,x,f,minoutunit,optpars,minmethod)
    write(*,*) ' Locally minimized input structure: '
    ii = useroutine
    useroutine = 1
    call computef(x,f)
    useroutine = ii
    write(*,*) ' -- Energy: ', f
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

    ! Write locally minimized initial point to trajectory file

    if ( .not. restart ) then
      call write_pdb(pdb_file,addbest(pdb_file_out),x)
      call system("cat ./testeBEST.pdb >> ./traj.pdb")
    end if

    ! Create reference individual with null dihedral angles

    write(*,*) 
    write(*,*) ' Creating null-dihedral reference structure ... '
    do i = 1, nresidues-1
      phi = compute_phi(i+1,x)
      xdihedinput(i+1,1) = phi
      phi = -1.d0*phi
      call rotate_phi(i+1,x,phi)
      psi = compute_psi(i,x)
      xdihedinput(i,2) = psi
      psi = -1.d0*psi
      call rotate_psi(i,x,psi)
    end do
    do i = 1, natoms*3
      xref(i) = x(i)
    end do
    ii = useroutine
    useroutine = 1
    call computef(x,f)
    useroutine = ii
    write(*,*) ' -- Energy: ', f
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

    ! Reconstructing the input structure to test (for testing)

    do i = 1, natoms*3
      x(i) = xref(i)
    end do
    do i = 1, nresidues-1
      phi = xdihedinput(i+1,1)
      call rotate_phi(i+1,x,phi)
      psi = xdihedinput(i,2)
      call rotate_psi(i,x,psi)
    end do
    call computef(x,f)
    write(*,*) ' Reconstructed input structure function: ', f
    call compute_fconst(x,f)
    write(*,*) ' Constraint violation of reconstructed structure: ', f
    
    ! Output files

    if ( .not. restart ) then
      call computef(xref,f)
      call write_pdb(pdb_file,addbest(pdb_file_out),xref)
      call system("cat ./testeBEST.pdb >> ./traj.pdb")
      write(*,*) ' Wrote reference structure to output file. '
    end if

    ! Initialize random number generator

    write(*,*) ' Seed for random number generator: ', seed
    call init_random_number(seed)

    ! Maximum number of functional evaluations of optimization method

    write(*,*) 
    write(*,*) ' Maximum number of function evaluations of minimization method: ', optpars(1)
    write(*,*) 

    ! Initial population
   
    if ( .not. restart ) then

      ! Create random initial individuals
   
      write(*,*) ' Creating random initial individuals ... '
      do i = 1, nind

        do j = 1, natoms*3
          xtmp(j) = xref(j)
        end do
        do j = 1, nresidues - 1

          if ( useinput .and. i == 1 ) then
            phi = xdihedinput(j+1,1)
            psi = xdihedinput(j,2)
          else 
            call newangles(phi,psi)
          end if
          call rotate_phi(j+1,xtmp,phi)
          xdihed(i,j+1,1) = phi
          call rotate_psi(j,xtmp,psi)
          xdihed(i,j,2) = psi

        end do

        call computef(xtmp,f)
        find(i) = f
        do j = 1, natoms*3
          xind(i,j) = xtmp(j)
        end do

      end do
    
    else if ( restart ) then

      ! Restarting from file

      open(13,file="./restart.dat")
        do i = 1, nind
          read(13,*) iind, find(i)
          do j = 1, natoms*3
            read(13,*) xind(i,j)
          end do
        end do
      close(13)

    end if
      
    ! Update the order of the individuals according to their function values

    do i = 1, nind
      fsort(i) = find(i)
    end do
    call flashsort(fsort, nind, lflash, mflash, indflash, indflashback)

    do i = 1, nind

      ftemp = find(i)
      do j = 1, natoms*3
        xtmp(j) = xind(i,j)
      end do
      do j = 1, nresidues-1
        xdihedtemp(j+1,1) = xdihed(i,j+1,1) 
        xdihedtemp(j,2) = xdihed(i,j,2)
      end do

      find(i) = find(indflash(i))
      do j = 1, natoms*3
        xind(i,j) = xind(indflash(i),j)
      end do
      do j = 1, nresidues-1
        xdihed(i,j+1,1) = xdihed(indflash(i),j+1,1) 
        xdihed(i,j,2) = xdihed(indflash(i),j,2)
      end do

      find(indflash(i)) = ftemp
      do j = 1, natoms*3
        xind(indflash(i),j) = xtmp(j)
      end do
      do j = 1, nresidues-1
        xdihed(indflash(i),j+1,1) = xdihedtemp(j+1,1) 
        xdihed(indflash(i),j,2) = xdihedtemp(j,2)
      end do

      do j = i, nind
        if ( indflash(j) == i ) then
          indflash(j) = indflash(i)
          exit
        end if
      end do

    end do

    fbest = find(1)
    do j = 1, natoms*3
      xtmp(j) = xind(1,j)
    end do
    call compute_fconst(xtmp,fconstbest1)
    do i = 1, ncgroups
      nctemp(i) = ncconsider(i)
      ncconsider(i) = ncpergroup(i)
    end do
    call compute_fconst(xtmp,fconstbest2)
    do i = 1, ncgroups
      ncconsider(i) = nctemp(i)
    end do
    write(*,*) ' Function value of best individual: ', fbest

    ! Test analytical gradient against finite difference gradients (testing only)

    if ( testgradient ) then
      do i = 1, natoms*3
        x(i) = xind(1,i)
      end do
      call test_grad(natoms*3,x,g,computef,computeg)
      call stop_all()
    end if

    ! Save the dihedrals of the best individual

    do i = 1, nresidues-1
      xdihedbest(i+1,1) = xdihed(1,i+1,1)
      xdihedbest(i,2) = xdihed(1,i,2)
    end do

    ! Save best coordinates

    do i = 1, natoms*3
      xbest(i) = xind(1,i)
    end do

    !
    ! Start genetic algorithm
    !

    write(*,*) ' Starting genetic algorithm ... '
    write(*,*) ' Number of surviving individuals: ', nsurvive
    open(998,file='./rmsd.dat')
    open(999,file='./generation.dat')

    ibest = 0
    irecv = 0
    iprintlast = -1
    score = 0.d0
    nchange = 0

    ! Set to free the status of all ranks

    do rank = 1, nrank-1 
      free(rank) = .true.
    end do

    genetic : do 

      ! Write the current status of the population

      iprint = mod(irecv,nind)
      if ( iprint ==  0 .and. irecv /= iprintlast ) then

        iprintlast = irecv
        write(*,"( 90('-') )")
        write(*,"( ' Population: ', i5, tr2, 70('-') )") irecv/nind
        write(999,"( i8, 100000(tr2,e17.10) )") irecv, (find(i), i = 1, nind) 
        i = 1
        ii = 0
        do while ( ii < 10 .and. ii < nind ) 
          ii = min0( i + 4, nind )
          write(*,"( 5( i5, tr1, e12.3 ) )") (j, find(j), j = i, ii)
          i = i + 5
        end do
        if ( nind > 10 ) then
          write(*,"( t14,'...  ', 4( i5, tr1, e12.3 ) )") (j, find(j), j=nind-3, nind)
        end if
        write(*,"( a, tr2, e17.10, a, i4 )") ' Best score up to now: ', &
                                               find(1), ' from population ', ibest/nind
        write(*,"( a, tr2, e12.3 )") ' Score of worst surival: ', find(nsurvive)
        write(*,"( a, f17.5 )") ' Constraint violation of best structure (for LOVO) = ', fconstbest1
        write(*,"( a, f17.5 )") ' Constraint violation of best structure (for all constraints) = ', fconstbest2
        score = score / nind
        write(*,"( a, f12.5, a, i5, a )") &
              ' Population improvement score: ', score,' (', nchange, ' individuals changed.)'

        ! Analyze the RMSDs of individuals 

        write(998,"( 90('-') )")
        write(998,*) ' Population: ', irecv/nind
        write(998,*) '        IND                 RMSD                   RMSDMAX'
        do iind = 1, nind
          do i = 1, natoms
            ix = (i-1)*3+1
            iy = ix + 1
            iz = ix + 2
            xali(i,1) = xind(iind,ix)
            xali(i,2) = xind(iind,iy)
            xali(i,3) = xind(iind,iz)
            if ( iind == 1 ) then
              xaliref(i,1) = xind(1,ix)
              xaliref(i,2) = xind(1,iy)
              xaliref(i,3) = xind(1,iz)
            end if
          end do
          call align(natoms,xali,xaliref)
          rmsd = 0.d0
          rmsdmax = 0.d0
          do i = 1, natoms
            dist = (xali(i,1)-xaliref(i,1))**2 + &
                   (xali(i,2)-xaliref(i,2))**2 + &
                   (xali(i,3)-xaliref(i,3))**2
            rmsd = rmsd + dist
            rmsdmax = max( dist, rmsdmax )
          end do
          rmsd = dsqrt(rmsd/natoms)
          rmsdmax = dsqrt(rmsdmax)
          write(998,*) iind, rmsd, rmsdmax
          do i = 1, natoms
            ix = (i-1)*3+1
            iy = ix + 1
            iz = ix + 2
            xind(iind,ix) = xali(i,1)
            xind(iind,iy) = xali(i,2)
            xind(iind,iz) = xali(i,3)
          end do
          if ( writepop ) then
            do i = 1, natoms*3
              xtmp(i) = xind(iind,i)
            end do
            record = './teste_ind.pdb'
            call write_pdb(pdb_file,record,xtmp)
            if ( iind == 1 ) then
              call system("cat ./teste_ind.pdb > ./lastpop.pdb")
            else
              call system("cat ./teste_ind.pdb >> ./lastpop.pdb")
            end if
         end if
        end do

        score = 0.d0
        nchange = 0
        write(*,"( 90('-') )")

        ! Write population to restart file

        open(13,file="./restart.dat")
        do i = 1, nind
          write(13,*) i, find(i)
          do j = 1, natoms*3
            write(13,*) xind(i,j)
          end do
        end do
        close(13)

        ! Stop if maximum number of populations is reached

        if ( irecv/nind == maxpop ) then
          close(998)
          close(999)
          write(*,*) ' Maximum number of populations reached. Stopping. '
          write(*,*) ' Best function value = ', fbest, ' Constraint violation = ', fconstbest1 
          call stop_all()
        end if

        if ( pause ) then
          write(*,*) ' Pause. Type 1 to continue: ' 
          read(*,*) ii 
        end if

      end if

      ! Check if some parameter has changed in the input file

      open(10,file=inputfile,action="read")
      ffchange = .false.
      do
        read(10,"( a200 )",iostat=ioerr) record
        if( ioerr /= 0 ) exit
        if(record(1:1) == "#" .or. len_trim(record) < 1 ) cycle
        select case ( keyword(record) )
          case ("bonds") 
            if ( value(record,1) == "yes" .and. .not. ffcomp(1) ) then
              write(*,*) ' Changed bonds to yes. '
              ffcomp(1) = .true. ; ffchange = .true.
            else if ( value(record,1) == "no" .and. ffcomp(1) ) then
              write(*,*) ' Changed bonds to no. '
              ffcomp(1) = .false. ; ffchange = .true.
            end if
          case ("angles")
            if ( value(record,1) == "yes" .and. .not. ffcomp(2) ) then
              write(*,*) ' Changed angles to yes. '
              ffcomp(2) = .true. ; ffchange = .true.
            else if ( value(record,1) == "no" .and. ffcomp(2) ) then
              write(*,*) ' Changed angles to no. '
              ffcomp(2) = .false. ; ffchange = .true.
            end if
          case ("dihedrals")
            if ( value(record,1) == "yes" .and. .not. ffcomp(3) ) then
              write(*,*) ' Changed dihedrals to yes. '
              ffcomp(3) = .true. ; ffchange = .true.
            else if ( value(record,1) == "no" .and. ffcomp(3) ) then
              write(*,*) ' Changed dihedrals to no. '
              ffcomp(3) = .false. ; ffchange = .true.
            end if
          case ("impropers")
            if ( value(record,1) == "yes" .and. .not. ffcomp(4) ) then
              write(*,*) ' Changed impropers to yes. '
              ffcomp(4) = .true. ; ffchange = .true. 
            else if ( value(record,1) == "no" .and. ffcomp(4) ) then
              write(*,*) ' Changed impropers to no. '
              ffcomp(4) = .false. ; ffchange = .true.
            end if
          case ("nonbonded")
            if ( value(record,1) == "yes" .and. .not. ffcomp(5) ) then
              write(*,*) ' Changed nonbonded to yes. '
              ffcomp(5) = .true. ; ffcomp(7) = .true. ; ffchange = .true.
            else if ( value(record,1) == "yes" .and. ffcomp(5) .and. .not. ffcomp(7) ) then
              write(*,*) ' Changed nonbonded to yes. '
              ffcomp(5) = .true. ; ffcomp(7) = .true. ; ffchange = .true.
            else if ( value(record,1) == "no" .and. ffcomp(5) ) then
              write(*,*) ' Changed nonbonded to no. '
              ffcomp(5) = .false. ; ffchange = .true.
            else if ( value(record,1) == "backbone" .and. .not. ffcomp(5) ) then
              write(*,*) ' Will compute non-bonded backbone interactions now. '
              ffcomp(5) = .true. ;  ffcomp(7) = .false. ; ffchange = .true.
            else if ( value(record,1) == "backbone" .and. ffcomp(5) .and. ffcomp(7) ) then
              write(*,*) ' Will compute ONLY backbone non-bonded interactions now. '
              ffcomp(5) = .true. ;  ffcomp(7) = .false. ; ffchange = .true.
            end if
          case ("sidechains")
            if ( value(record,1) == "polar" .and. .not. ffcomp(8) ) then
              write(*,*) ' Changed sidechains to polar. '
              ffcomp(8) = .true. ; ffchange = .true.
              call updatels()
            else if ( value(record,1) == "neutral" .and. ffcomp(8) ) then
              write(*,*) ' Changed sidechains to neutral. '
              ffcomp(8) = .false. ; ffchange = .true.
              call updatels()
            end if
          case ("useconstraints") 
            if ( value(record,1) == "yes" .and. .not. ffcomp(6) ) then
              write(*,*) ' Changed constraints to yes. '
              ffcomp(6) = .true. ; ffchange = .true. 
            else if ( value(record,1) == "no" .and. ffcomp(6) ) then
              write(*,*) ' Changed constraintsconstraints  to no. '
              ffcomp(6) = .false. ; ffchange = .true.
            end if
          case ("cutoff")
            record2 = value(record,1)
            read(record2,*) ftemp
            if ( abs(ftemp - cutoff) > 1.d-8 ) then
              cutoff = ftemp
              write(*,*) ' Changed cutoff to ', cutoff
              ffchange = .true.
            end if
          case ("combine")
            if ( value(record,1) == "yes" .and. .not. combine ) then
              write(*,*) ' Changed combine to yes. '
              combine = .true.
            else if ( value(record,1) == "no" .and. combine ) then
              write(*,*) ' Changed combine to no. '
              combine = .false. 
            end if
          case ("mirror")
            if ( value(record,1) == "yes" .and. .not. mirror ) then
              write(*,*) ' Changed mirror to yes. '
              mirror = .true.
            else if ( value(record,1) == "no" .and. mirror ) then
              write(*,*) ' Changed mirror to no. '
              mirror = .false. 
            end if
          case ("mutations")
            if ( value(record,1) == "yes" .and. .not. mutations ) then
              write(*,*) ' Changed mutations to yes. '
              mutations = .true.
            else if ( value(record,1) == "no" .and. mutations ) then
              write(*,*) ' Changed mutations to no. '
              mutations = .false. 
            end if
          case ("minimize")
            if ( value(record,1) == "yes" .and. .not. minimize ) then
              write(*,*) ' Changed minimize to yes. '
              minimize = .true.
            else if ( value(record,1) == "no" .and. minimize ) then
              write(*,*) ' Changed minimize to no. '
              minimize = .false. 
            end if
          case ("minmethod")
            if ( value(record,1) == "spg" .and. minmethod == 2 ) then
              write(*,*) ' Changed minmethod to spg. '
              minmethod = 1
            else if ( value(record,1) == "cgnewton" .and. minmethod == 1 ) then
              write(*,*) ' Changed minmethod to cgnewton. '
              minmethod = 2
            end if
          case ("pause")
            if ( value(record,1) == "yes" .and. .not. pause ) then
              write(*,*) ' Changed pause to yes. '
              pause = .true.
            else if ( value(record,1) == "no" .and. pause ) then
              write(*,*) ' Changed pause to no. '
              pause = .false. 
            end if
          case ("print")
            record2 = value(record,1)
            read(record2,*) ii
            if ( ii /= print ) then
              print = ii
              write(*,*) ' Changed print to ', print
            end if
          case ("maxfeval")
            record2 = value(record,1)
            read(record2,*) ii
            if ( ii /= optpars(1) ) then
              optpars(1) = ii
              write(*,*) ' Changed maxopt to ', optpars(1)
            end if
          case ("maxcg")
            record2 = value(record,1)
            read(record2,*) ii
            if ( ii /= optpars(2) ) then
              optpars(2) = ii
              write(*,*) ' Changed maxcg to ', optpars(2)
            end if
          case ("nsurvive")
            record2 = value(record,1)
            read(record2,*) ii
            if ( ii /= nsurvive ) then
              nsurvive = ii
              write(*,*) ' Changed nsurvive to ', nsurvive
            end if
          case ("terminate")
            terminate = .true.
        end select
      end do
      close(10)

      ! If some force-field parameter changed, update the function values 
      
      if ( ffchange ) then 

        write(*,*) ' Updating function values for new force-field. '
        do i = 1, nind
          do j = 1, natoms*3
            x(j) = xind(i,j)
          end do
          call computef(x,f)
          find(i) = f
        end do

        ! Update the order of the individuals according to their function values

        do i = 1, nind
          fsort(i) = find(i)
        end do
        call flashsort(fsort, nind, lflash, mflash, indflash, indflashback)

        do i = 1, nind

          ftemp = find(i)
          do j = 1, natoms*3
            xtmp(j) = xind(i,j)
          end do
          do j = 1, nresidues-1
            xdihedtemp(j+1,1) = xdihed(i,j+1,1) 
            xdihedtemp(j,2) = xdihed(i,j,2)
          end do

          find(i) = find(indflash(i))
          do j = 1, natoms*3
            xind(i,j) = xind(indflash(i),j)
          end do
          do j = 1, nresidues-1
            xdihed(i,j+1,1) = xdihed(indflash(i),j+1,1) 
            xdihed(i,j,2) = xdihed(indflash(i),j,2)
          end do

          find(indflash(i)) = ftemp
          do j = 1, natoms*3
            xind(indflash(i),j) = xtmp(j)
          end do
          do j = 1, nresidues-1
            xdihed(indflash(i),j+1,1) = xdihedtemp(j+1,1) 
            xdihed(indflash(i),j,2) = xdihedtemp(j,2)
          end do

          do j = i, nind
            if ( indflash(j) == i ) then
              indflash(j) = indflash(i)
              exit
            end if
          end do

        end do
        fbest = find(1)
        write(*,*) ' Done. New best function value: ', fbest
        do i = 1, natoms*3
          x(i) = xind(1,i)
        end do
        call compute_fconst(x,fconstbest1)
        do i = 1, ncgroups
          nctemp(i) = ncconsider(i)
          ncconsider(i) = ncpergroup(i)
        end do
        call compute_fconst(x,fconstbest2)
        do i = 1, ncgroups
          ncconsider(i) = nctemp(i)
        end do

      end if

      ! Find the first free rank

      rank = 1
      do i = 1, nrank - 1
        if ( free(rank) ) exit
        rank = rank + 1
      end do

      !
      ! If there is a free rank, create a new individual and send it to the
      ! node for energy minimization
      ! 

      if ( rank <= nrank - 1 ) then

        ! Choose between sexual or assexual reproduction

        if ( combine ) then

          ! Combine best individuals to create a new one

          icover = 1
          jcover = 1 
          do while( jcover < nresidues-1 )

            ! Chose parent

            call random_number(random)
            iparent = int(random*nsurvive)+1
          
            ! Chose up to which point this parent will be copied

            call random_number(random)
            jcover = min0(nresidues-1,icover+int(random*(nresidues-1))+1)
          
            ! Copy parent

            do igen = icover, jcover
              xdihednew(igen+1,1) = xdihed(iparent,igen+1,1)
              xdihednew(igen,2) = xdihed(iparent,igen,2)
            end do

            icover = jcover + 1

          end do

        else if ( .not. combine ) then

          ! Otherwise, simply copy one of the best individuals

          call random_number(random)
          iparent = int(random*nsurvive)+1
          do i = 1, nresidues-1
            xdihednew(i+1,1) = xdihed(iparent,i+1,1)
            xdihednew(i,2) = xdihed(iparent,i,2)
          end do

        end if

        ! Perhaps the individuals should be mirroed

        if ( mirror ) then
          call random_number(random)
          if ( random < 0.1d0 ) then
            do i = 1, nresidues-1
              xdihednew(i+1,1) = -1.d0*xdihednew(i+1,1)
              xdihednew(i,2) = -1.d0*xdihednew(i,2)
            end do
          end if
        end if

        ! Mutate or not mutate?

        if ( mutations ) then
          call random_number(random)
          if ( random <= probmut ) then

            ! How many mutations
       
            call random_number(random)
            nmut = int(random*(maxmut-1)) + 1

            ! Mutate then

            mutate: do imut = 1, nmut

              ! Which residue

              call random_number(random)
              iresidue = int(random*(nresidues-1))+1 ! Between 1 and nresidues-1

              ! Dihedral rotation angle

              if ( muttype == 0 ) then
                call random_number(random)
                ! Phi
                if ( random < 0.5 ) then
                  call random_number(random)
                  phi = xdihednew(iresidue+1,1) + mutbin*(-180. + 360.*random)*torad
                  xdihednew(iresidue+1,1) = phi
                ! Psi
                else
                  call random_number(random)
                  psi = xdihednew(iresidue,2) + mutbin*(-180. + 360.*random)*torad
                  xdihednew(iresidue,2) = psi
                end if
              else if ( muttype == 1 ) then
                call newangles(phi,psi)
                xdihednew(iresidue+1,1) = phi
                xdihednew(iresidue,2) = psi
              end if

            end do mutate
          end if
        end if

        ! Recompute coordinates

        do j = 1, natoms*3
          x(j) = xref(j)
        end do
        do j = 1, nresidues-1
          ! Rotate Phi
          phi = xdihednew(j+1,1)
          call rotate_phi(j+1,x,phi)
          ! Rotate Psi
          psi = xdihednew(j,2)
          call rotate_psi(j,x,psi)
        end do

        ! Send this individual for energy minimization in the free node

        tag = 1
        call mpi_send( minimize, 1, mpi_logical, rank, tag, mpi_comm_world, ierr )
        call mpi_send( minmethod, 1, mpi_integer, rank, tag, mpi_comm_world, ierr )
        call mpi_send( optpars, 10, mpi_integer, rank, tag, mpi_comm_world, ierr )
        call mpi_send( cutoff, 1, mpi_double_precision, rank, tag, mpi_comm_world, ierr )
        call mpi_send( ffcomp, 8, mpi_logical, rank, tag, mpi_comm_world, ierr )
        call mpi_send( qq, n_nonbonded, mpi_double_precision, rank, tag, mpi_comm_world, ierr )
        call mpi_send( cutnb, n_nonbonded*3, mpi_double_precision, rank, tag, mpi_comm_world, ierr )
        call mpi_send( nshake, 1, mpi_integer, rank, tag, mpi_comm_world, ierr )
        call mpi_send( x, natoms*3, mpi_double_precision, rank, tag, mpi_comm_world, ierr )
        free(rank) = .false. 

      end if

      !
      ! Otherwise, if there are no free ranks, just wait for the next answer
      !

      if ( rank > nrank-1 ) then

        call mpi_recv( x, natoms*3, mpi_double_precision, &
                       mpi_any_source, mpi_any_tag, mpi_comm_world, &
                       status, ierr ) 
        irecv = irecv + 1

        ! Set to free the status of this rank

        free(status(mpi_source)) = .true.

        !
        ! Got an answer, therefore, the new individual will be classified
        !
 
        ! Compute the energy of this individual

        call computef(x,f)

        ! Find which is the fit place of this individual in the population

        iind = 1
        do i = 1, nind 
          if ( find(i) > f ) exit
          iind = iind + 1
        end do
        iind = i

        if ( print > 0 ) write(*,*) ' Received: ', irecv,' Position ', iind, ' F = ', f

        ! If the new individual does not have the worst function value,
        ! remove the worst from the population and replace it by this one

        if ( iind <= nind ) then 

          ! Update population improvement score

          score = score + dble(nind-iind+1)/nind
          nchange = nchange + 1

          ! Move all individuals removing the last one

          do i = nind, iind+1, -1
            find(i) = find(i-1)
            do j = 1, natoms*3
              xind(i,j) = xind(i-1,j)
            end do
            do j = 1, nresidues-1
              xdihed(i,j+1,1) = xdihed(i-1,j+1,1)
              xdihed(i,j,2) = xdihed(i-1,j,2)
            end do
          end do

          ! Add new individual to position iind

          find(iind) = f
          do j = 1, natoms*3
            xind(iind,j) = x(j)
          end do
          do i = 1, nresidues-1
            xdihed(iind,i+1,1) = compute_phi(i+1,x)
            xdihed(iind,i,2) = compute_psi(i,x)
          end do

          ! If the best individual is improved, save it 

          if ( find(iind) < fbest ) then
            do j = 1, natoms*3
              xbest(j) = xind(1,j) 
            end do
            do j = 1, nresidues-1
              xdihedbest(j+1,1) = xdihed(1,j+1,1)
              xdihedbest(j,2) = xdihed(1,j,2)
            end do
            frel = (fbest - find(1))/abs(fbest)
            call write_pdb(pdb_file,addbest(pdb_file_out),xbest)
            call system("cat ./testeBEST.pdb >> ./traj.pdb")
            call compute_fconst(xbest,fconstbest1)
            do i = 1, ncgroups
              nctemp(i) = ncconsider(i)
              ncconsider(i) = ncpergroup(i)
            end do
            call compute_fconst(xbest,fconstbest2)
            do i = 1, ncgroups
              ncconsider(i) = nctemp(i)
            end do
            fbest = find(1)
            ibest = irecv
            write(*,"( 90('-') )")
            write(*,"( a, e12.5, a, i6, a, e12.5, a, e12.5 )") &
                    ' New best structure. Improvement: ', frel, &
                    ' I: ', irecv, ' F= ', find(1), ' CV= ', fconstbest1
            write(*,"( 90('-') )")
          end if
        end if

      end if

      ! Terminate now if requested

      if ( terminate ) then

        close(998)
        close(999)
        if ( minoutunit > 0 ) close(minoutunit)

        ! 
        ! Send a termination signal to all processes 
        ! 

        write(*,*) ' Sending termintation signal to all processes...  '
        call stop_all()

        exit genetic
      end if
 
    end do genetic

  end if

  !
  ! Activities of the slave ranks
  !

  if ( myrank > 0 ) then

    ! Receive the signal that all input parameters were read correctly, so
    ! no premature end is needed

    call mpi_recv( error, 1, mpi_logical, 0, mpi_any_tag, mpi_comm_world, status, ierr )

    if ( .not. error ) then 

      !
      ! Receive the force field parameters that were read by the server rank
      !

      call mpi_bcast( natoms, 1, mpi_integer, 0, mpi_comm_world, ierr )
      allocate( x(natoms*3), g(natoms*3), xtmp(natoms*3) )   
      allocate ( isbackbone(natoms) ) 
      call mpi_bcast( isbackbone, natoms, mpi_logical, 0, mpi_comm_world, ierr )

      call mpi_bcast( n_nonbonded, 1, mpi_integer, 0, mpi_comm_world, ierr )
      allocate( ijnonbonded(n_nonbonded,2), qq(n_nonbonded), epseps(n_nonbonded), &
                ss6(n_nonbonded), cutnb(n_nonbonded,3) )
      call mpi_bcast( ijnonbonded, n_nonbonded*2, mpi_integer, 0, mpi_comm_world, ierr )
      call mpi_bcast( qq, n_nonbonded, mpi_double_precision, 0, mpi_comm_world, ierr )
      call mpi_bcast( epseps, n_nonbonded, mpi_double_precision, 0, mpi_comm_world, ierr )
      call mpi_bcast( ss6, n_nonbonded, mpi_double_precision, 0, mpi_comm_world, ierr )
      call mpi_bcast( cutnb, n_nonbonded*3, mpi_double_precision, 0, mpi_comm_world, ierr )
      call mpi_bcast( nshake, 1, mpi_integer, 0, mpi_comm_world, ierr )
      call mpi_bcast( linearshort, 1, mpi_logical, 0, mpi_comm_world, ierr )
      call mpi_bcast( n_nbbackbone, 1, mpi_integer, 0, mpi_comm_world, ierr ) 
      allocate ( inbbackbone(n_nbbackbone) )
      call mpi_bcast( inbbackbone, n_nbbackbone, mpi_integer, 0, mpi_comm_world, ierr )

      call mpi_bcast( nbonds, 1, mpi_integer, 0, mpi_comm_world, ierr )
      allocate( ibond(nbonds,2), fbond(nbonds,2) )
      call mpi_bcast( ibond, nbonds*2, mpi_integer, 0, mpi_comm_world, ierr )
      call mpi_bcast( fbond, nbonds*2, mpi_double_precision, 0, mpi_comm_world, ierr )

      call mpi_bcast( nangles, 1, mpi_integer, 0, mpi_comm_world, ierr )
      allocate( iangle(nangles,3), fangle(nangles,2) )
      call mpi_bcast( iangle, nangles*3, mpi_integer , 0, mpi_comm_world, ierr )
      call mpi_bcast( fangle, nangles*2, mpi_double_precision, 0, mpi_comm_world, ierr )

      call mpi_bcast( nureybradley, 1, mpi_integer, 0, mpi_comm_world, ierr )
      allocate( iureybradley(nureybradley,2), fureybradley(nureybradley,2) )
      call mpi_bcast( iureybradley, nureybradley*2, mpi_integer, 0, mpi_comm_world, ierr )
      call mpi_bcast( fureybradley, nureybradley*2, mpi_double_precision, 0, mpi_comm_world, ierr )

      call mpi_bcast( ndihed, 1, mpi_integer, 0, mpi_comm_world, ierr )
      call mpi_bcast( nfdihed, 1, mpi_integer, 0, mpi_comm_world, ierr )
      allocate( idihed(ndihed,6), fdihed(nfdihed,3) )
      call mpi_bcast( idihed, ndihed*6, mpi_integer, 0, mpi_comm_world, ierr )
      call mpi_bcast( fdihed, nfdihed*3, mpi_double_precision, 0, mpi_comm_world, ierr )

      call mpi_bcast( nimpr, 1, mpi_integer, 0, mpi_comm_world, ierr )
      call mpi_bcast( nfimpr, 1, mpi_integer, 0, mpi_comm_world, ierr )
      allocate( iimpr(nimpr,6), fimpr(nfimpr,2) )
      call mpi_bcast( iimpr, nimpr*6, mpi_integer, 0, mpi_comm_world, ierr )
      call mpi_bcast( fimpr, nfimpr*2, mpi_double_precision, 0, mpi_comm_world, ierr ) 
    
      call mpi_bcast( ncgroups, 1, mpi_integer, 0, mpi_comm_world, ierr )
      allocate( ncpergroup(ncgroups), ncconsider(ncgroups) )
      call mpi_bcast( ncpergroup, ncgroups, mpi_integer, 0, mpi_comm_world, ierr )
      call mpi_bcast( ncconsider, ncgroups, mpi_integer, 0, mpi_comm_world, ierr )

      call mpi_bcast( ntotconst, 1, mpi_integer, 0, mpi_comm_world, ierr )
      allocate( kconstraint(ntotconst), iconstraint(ntotconst,2), dconstraint(ntotconst,2), &
                violation(ntotconst), indflash(ntotconst), lflash(ntotconst), &
                indflashback(ntotconst), viol(ntotconst) )
      call mpi_bcast( kconstraint, ntotconst, mpi_double_precision, 0, mpi_comm_world, ierr )
      call mpi_bcast( iconstraint, ntotconst*2, mpi_integer, 0, mpi_comm_world, ierr )
      call mpi_bcast( dconstraint, ntotconst*2, mpi_double_precision, 0, mpi_comm_world, ierr )
      mflash = 1 + ntotconst / 10

      call mpi_bcast( minoutunit, 1, mpi_integer, 0, mpi_comm_world, ierr )
      call mpi_bcast( useroutine, 1, mpi_integer, 0, mpi_comm_world, ierr )
      call mpi_bcast( maxboxes, 1, mpi_integer, 0, mpi_comm_world, ierr )
      allocate( knonbonded(natoms,natoms), &
                iatomfirst(0:maxboxes+1,0:maxboxes+1,0:maxboxes+1), iatomnext(natoms) )
      call mpi_bcast( knonbonded, natoms*natoms, mpi_integer, 0, mpi_comm_world, ierr ) 

      ! This is an infinite loop that minimies energies on requests, and will be
      ! quited from inside if the termination tag is received

      receive : do
       
        call mpi_recv( minimize, 1, mpi_logical, 0, mpi_any_tag, mpi_comm_world, status, ierr )

        ! If tag == 0, this is a termination call, so exit the infinite loop 

        tag = status(mpi_tag)
        if ( tag == 0 ) then
          exit receive
        end if

        call mpi_recv( minmethod, 1, mpi_integer, 0, mpi_any_tag, mpi_comm_world, status, ierr )
        call mpi_recv( optpars, 10, mpi_integer, 0, mpi_any_tag, mpi_comm_world, status, ierr )
        call mpi_recv( cutoff, 1, mpi_double_precision, 0, mpi_any_tag, mpi_comm_world, status, ierr )
        cutoff2 = cutoff**2
        call mpi_recv( ffcomp, 8, mpi_logical, 0, mpi_any_tag, mpi_comm_world, status, ierr )
        call mpi_recv( qq, n_nonbonded, mpi_double_precision, 0, mpi_any_tag, mpi_comm_world, status, ierr )
        call mpi_recv( cutnb, n_nonbonded*3, mpi_double_precision, 0, mpi_any_tag, mpi_comm_world, status, ierr )
        call mpi_recv( nshake, 1, mpi_integer, 0, mpi_any_tag, mpi_comm_world, status, ierr )

        call mpi_recv( x, natoms*3, mpi_double_precision, 0, mpi_any_tag, mpi_comm_world, status, ierr )  

        if( minimize ) then 

          ! Optimize structure received 

          call optimize(natoms*3,x,f,minoutunit,optpars,minmethod)

          ! Try to escape from trivial local minima by performing small
          ! perturbations (0.5 Angstrom perturbations on the coordinates of each
          ! atom)

          fbest = f
          nshake = 0
          do i = 1, nshake
            do j = 1, natoms*3
              xtmp(j) = x(j)
              call random_number(random)
              xtmp(j) = xtmp(j) - 0.5d0 + random
            end do
            call optimize(natoms*3,xtmp,f,minoutunit,optpars,minmethod)
            if ( f < fbest ) then
              fbest = f
              do j = 1, natoms*3
                x(j) = xtmp(j)
              end do
            end if
          end do

        end if

        ! Return the minimized individual to the master node 
        ! The tag number is the iind number which was sent from the master

        call mpi_send( x, natoms*3, mpi_double_precision, 0, tag, mpi_comm_world, ierr )

      end do receive

    end if

  end if

  call mpi_finalize(ierr)
  if ( myrank == 0 ) write(*,*) ' Exiting gracefully. '

end program crfit
