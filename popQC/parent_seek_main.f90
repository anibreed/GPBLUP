program parent_seek_main
  implicit none

  character(len=1024) :: param_file, home_dir, popqc_bin, cmd
  character(len=1024) :: param_dir, param_name
  integer :: slash_pos
  integer :: cmdstat, exitstat
  logical :: exe_exists

  param_file = ''
  call getarg(1, param_file)

  if (len_trim(param_file) < 1) then
    print *, 'Usage: parentSeek <parameter_file>'
    print *, 'Example: parentSeek ~/DATA/FinalReport/DD/param'
    stop 1
  end if

  home_dir = ''
  call get_environment_variable('HOME', home_dir)
  if (len_trim(home_dir) > 0) then
    popqc_bin = trim(home_dir) // '/GPBLUP/bin/popQC'
  else
    popqc_bin = 'popQC'
  end if

  inquire(file=trim(popqc_bin), exist=exe_exists)
  if (.not. exe_exists) popqc_bin = 'popQC'

  slash_pos = index(trim(param_file), '/', back=.true.)
  if (slash_pos > 0) then
    param_dir = trim(param_file(1:slash_pos-1))
    param_name = trim(param_file(slash_pos+1:))
  else
    param_dir = '.'
    param_name = trim(param_file)
  end if

  cmd = 'cd "' // trim(param_dir) // '" && ' // trim(popqc_bin) // ' "' // trim(param_name) // '" PARENT_SEEK_ONLY'
  print *, 'Launching:', trim(cmd)

  call execute_command_line(trim(cmd), wait=.true., exitstat=exitstat, cmdstat=cmdstat)
  if (cmdstat /= 0) then
    print *, 'ERROR: failed to execute command. CMDSTAT=', cmdstat
    stop 2
  end if
  if (exitstat /= 0) then
    print *, 'ERROR: popQC exited with status ', exitstat
    stop exitstat
  end if
end program parent_seek_main
