! This function adds a "BEST" to the pdb file name to print best structure

character(len=200) function addbest(pdbfile)

  integer :: length, i
  character(len=200) :: pdbfile
  pdbfile = trim(pdbfile)
  length = len(trim(pdbfile))
  addbest = pdbfile(1:length-4)//'BEST'//'.pdb'
  do i = length + 5, 200
    addbest(i:i) = ' '
  end do

end function addbest

