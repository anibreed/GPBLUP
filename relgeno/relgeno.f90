program genoinv
  use M_readpar_relgeno, only: relgeno_params_t, read_relgeno_params
  use M_relgeno_service, only: run_relgeno
  implicit none

  character(len=1), parameter :: missC = achar(32)
  character(len=256) :: param_file
  type(relgeno_params_t) :: params

  param_file = missC
  call getarg(1, param_file)
  param_file = trim(adjustl(param_file))
  if (len_trim(param_file) < 1) stop 'Input parameter file'

  call read_relgeno_params(param_file, params)
  if (.not. params%compute_g) then
    print *, 'INFO: COMPUTE_G=FALSE. Exiting without building G.'
    stop
  end if

  print '(a,1x,a)', 'RELGENO param file:', trim(param_file)
  call run_relgeno(params)
end program genoinv





