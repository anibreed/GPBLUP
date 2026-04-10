module M_phimpute_mpi
  implicit none

  logical, save :: mpi_initialized = .false.

  public :: mpi_init_once, mpi_finalize_once, mpi_get_rank_size

contains

  subroutine mpi_init_once()
#ifdef USE_MPI
    use mpi
#endif
    integer :: ierr
    logical :: already_init

    already_init = .false.
#ifdef USE_MPI
    call MPI_Initialized(already_init, ierr)
#endif
    if (.not. already_init .and. .not. mpi_initialized) then
#ifdef USE_MPI
      call MPI_Init(ierr)
#endif
      mpi_initialized = .true.
    end if
  end subroutine mpi_init_once

  subroutine mpi_finalize_once()
#ifdef USE_MPI
    use mpi
#endif
    integer :: ierr
    logical :: already_final

    already_final = .false.
#ifdef USE_MPI
    call MPI_Finalized(already_final, ierr)
#endif
    if (.not. already_final .and. mpi_initialized) then
#ifdef USE_MPI
      call MPI_Finalize(ierr)
#endif
      mpi_initialized = .false.
    end if
  end subroutine mpi_finalize_once

  subroutine mpi_get_rank_size(rank, size)
#ifdef USE_MPI
    use mpi
#endif
    integer, intent(out) :: rank, size
    integer :: ierr

    rank = 0
    size = 1
#ifdef USE_MPI
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)
#endif
  end subroutine mpi_get_rank_size

end module M_phimpute_mpi
