program phimpute_main
  use M_phimpute_service, only: run_phimpute_from_file
  use M_phimpute_mpi, only: mpi_init_once, mpi_finalize_once
  implicit none

  character(len=256) :: par_file
  integer :: arg_count

  arg_count = command_argument_count()
  if (arg_count < 1) then
    write(*,'(A)') 'Usage: phimpute <param_file>'
    stop 1
  end if

  call get_command_argument(1, par_file)
  call mpi_init_once()
  call run_phimpute_from_file(trim(par_file))
  call mpi_finalize_once()
end program phimpute_main
