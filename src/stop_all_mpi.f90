!
! subroutine stop_all: Ends all MPI programs
!
! L. Martinez
!
! Institute of Chemistry - State University of Campinas - Brazil
!

subroutine stop_all()

  use mpi
  implicit none
  logical :: error
  integer :: tag, rank, nrank, ierr

  call mpi_comm_size(mpi_comm_world,nrank,ierr)
  error = .true.
  tag = 0
  do rank = 1, nrank-1
    call mpi_send( error, 1, mpi_logical, rank, tag, mpi_comm_world, ierr )
  end do
  call mpi_finalize(ierr)
  stop

end subroutine stop_all

