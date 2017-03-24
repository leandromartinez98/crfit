!
! Subroutine write_pdb: Writes solution as pdb file
!
! L. Martinez, Nov 28, 2013
! Institute of Chemistry - State University of Campinas - Brazil
!

subroutine write_pdb(pdb_file,pdb_file_out,x)

  use force_field
  implicit none
  integer :: icount, ioer, ix, iy, iz, i, ii
  character(len=200) :: pdb_file, pdb_file_out, line
  double precision :: x(*), f, energy, xcm, ycm, zcm

  open(10,file=pdb_file,action='read')
  open(20,file=pdb_file_out,action='write')
  icount = 0

  call computef(x,f)
  ii = useroutine
  useroutine = 1
  call computef(x,energy)
  useroutine = ii

  write(20,"( 'REMARK: ENERGY = ',f12.3, ' FVALUE = ', f12.3 )") energy, f
  write(20,"( 'REMARK: BOND: ',f12.3, ' ANG: ',f12.3,' DIHED: ',f12.3,&
              ' IMPR: ',f12.3, ' VDW: ',f12.3,' ELEC: ',f12.3 )") &
              (ecomp(i),i=1,6)
  write(20,"( 'REMARK: Constraint violation: ', f12.3 )") ecomp(7)
  call computef(x,f)
  write(20,"( 'REMARK: Function value: ', f12.3 )") f 

  xcm = 0.d0
  ycm = 0.d0
  zcm = 0.d0
  do i = 1, natoms*3, 3
    xcm = xcm + x(i)
    ycm = ycm + x(i+1)
    zcm = zcm + x(i+2)
  end do
  xcm = xcm / dble(natoms)
  ycm = ycm / dble(natoms)
  zcm = zcm / dble(natoms)

  do
    read(10,"( a200 )",iostat=ioer) line
    if ( ioer /= 0 ) exit
    if ( line(1:4) == "ATOM" .or. line(1:6) == "HETATM" ) then
      ix = icount*3 + 1
      iy = ix + 1
      iz = ix + 2
      write(20,"( a, 3(f8.3), a )") line(1:30),& 
                                    x(ix)-xcm,x(iy)-ycm,x(iz)-zcm,&
                                    line(55:len_trim(line))
      icount = icount + 1
    else
      write(20,"( a )") line(1:len_trim(line))
    end if
  end do
  close(10)
  close(20)

  return
end subroutine write_pdb

