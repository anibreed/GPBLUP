program relped
  use M_Variables, only: MAX_STR, missC
  use M_relped_service, only: run_relped_from_file
  implicit none

  character(len=MAX_STR) :: PARAM_File
  logical :: use_gpu
  integer :: device_id

  PARAM_File = missC
  use_gpu = .false.
  device_id = 0

  call getarg(1, PARAM_File)
  PARAM_File = trim(adjustl(PARAM_File))
  if (len_trim(PARAM_File) < 1) stop 'Input parameter file'
  print *, "Finish to read parameter file=  ", trim(PARAM_File)

  call run_relped_from_file(PARAM_File, use_gpu, device_id)
end program
