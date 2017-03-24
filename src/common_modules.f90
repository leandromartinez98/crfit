!
! Module force_field: Contains all arrays that are necessary for energy
!                     calculations
!
! L. Martinez, Nov 11, 2013
! Institute of Chemistry - State University of Campinas - Brazil
!

module force_field

  ! Allocatable arrays
  
  integer :: natoms, nresidues
  character(len=5), allocatable :: segment(:), residue(:), name(:), type(:) 
  integer, allocatable :: residue_number(:), natres(:)
  double precision, allocatable :: charge(:), mass(:), eps(:), sig(:), eps14(:), sig14(:)
  double precision :: scaling14, cutoff, cutoff2, ecomp(8), dielectric
  
  integer :: nbonds
  integer, allocatable :: ibond(:,:), iorder(:)
  double precision, allocatable :: fbond(:,:)
  
  integer :: nangles, nureybradley
  integer, allocatable :: iangle(:,:), iureybradley(:,:)
  double precision, allocatable :: fangle(:,:), fureybradley(:,:)
  
  !
  ! Dihedral angles
  !
  ! idihed, iimpr and icrterm are of dimension (n,6). For each dihedral angle,
  ! the four first components (i,1)...(i,4) are the indexes of the atoms involved
  ! in the dihedral. The fifth (i,5) component is the index in fdihed (etc.) of the first
  ! parameter of the corresponding dihedral, and the sixth component (i,6) is the
  ! number of parameters for this dihedral. Because each dihedral may have any
  ! number of parameters, the fdihed, fimpr, and ficrterm cannot be allocated on
  ! reading the psf file, but only after reading the parameter file.
  
  integer :: ndihed, nfdihed
  integer, allocatable :: idihed(:,:)
  double precision, allocatable :: fdihed(:,:)
  
  integer :: nimpr, nfimpr
  integer, allocatable :: iimpr(:,:)
  double precision, allocatable :: fimpr(:,:)
  
  integer :: ncrterm
  integer, allocatable :: icrterm(:,:)
 
  ! Non-bonded inclusion list (will be allocated when this lists are built)
  ! Since we are doing this with inclusions directly, why not compute the
  ! charge, epsilon and sigma products now only once?

  integer :: n_nonbonded, n_nbbackbone
  integer, allocatable :: ijnonbonded(:,:), inbbackbone(:)
  double precision, allocatable :: qq(:), qqreal(:), epseps(:), ss6(:), cutnb(:,:)
  logical, allocatable :: isbackbone(:)

  ! Variables that define which parts of the force field are used
  
  integer :: useroutine
  logical :: ffcomp(8)

end module force_field

!
! Module constraints: Carries the information of the constraints to the 
!                    optimization routines 
!

module constraints

  integer :: nconstfiles, ntotconst, ncgroups
  integer, allocatable :: ncpergroup(:), iconstraint(:,:)
  double precision, allocatable :: dconstraint(:,:), violation(:), viol(:)

  ! Force constant for the linear potential of the constraints

  double precision, allocatable :: kconstraint(:)

  ! Number of constraints actually considered (LOVO)

  integer, allocatable :: ncconsider(:)

end module constraints

!
! Module dihedrals: Defines the vectors required for the random 
!                   modification of the dihedrals in the genetic algorithm
!                   code
!

module dihedrals

  integer, allocatable :: backbone(:,:)
  
end module dihedrals

!
! Module linkedcells : Defines the vectors and dimensions needed by
!                      the linked cell method of computing the non-bonded
!                      interactions
!

module linkedcells

  integer :: maxboxes, nboxes(3)
  integer, allocatable :: knonbonded(:,:), iatomfirst(:,:,:), iatomnext(:)
  double precision :: boxcut(3)

end module linkedcells

!
! Gets keyword from input file
!

function keyword(string)

  implicit none
  integer :: if, il
  character(len=200) :: keyword, string

  if = 1
  do while(string(if:if) <= ' '.and. if < 200)
    if = if + 1
  end do
  il = if
  do while(string(il:il) > ' '.and.il < 200)
    il = il + 1
  end do
  il = il - 1
  keyword = string(if:il)

return
end function keyword

!
! Gets keyword value from input file
!

function value(string,ivalue)

  implicit none
  integer :: if, il, length, ivalue, i
  character(len=200) :: value, string

  ! Jump keyword 

  if = 1
  do while(string(if:if) <= ' '.and.if < 200)
    if = if + 1
  end do
  il = if
  do while(string(il:il) > ' '.and.il < 200)
    il = il + 1
  end do

  ! The keyword ended, now reading values

  do i = 1, ivalue
    il = il - 1
    if = il + 1
    do while(string(if:if) <= ' '.and.if < 200)
      if = if + 1
    end do
    il = if
    do while(string(il:il) > ' '.and.il < 200)
      il = il + 1
    end do
  end do

  value = string(if:il)
  if(length(value) == 0) then
    write(*,*) ' ERROR: Some keyword without value: '
    write(*,*) string(1:length(string))
    call stop_all()
  end if

return
end function value

!
! Function that sets the length of a string
!

function length(string)

  implicit none
  integer :: length
  character(len=200) :: string

  length = 200
  do while(string(length:length) <= ' ')
    length = length - 1
    if ( length == 0 ) exit
  end do

return
end function length


