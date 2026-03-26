program test_readpar
  use M_readpar
  implicit none

  character(len=256) :: param_file
  type(FileInfo) :: pedinfo
  integer :: i

  ! 파라미터 파일 경로 지정
  param_file = '~/DATA/FinalReport/param'

  ! 파라미터 파일 읽기
  call read_parameters(param_file)

  ! PEDFile 정보 가져오기
  pedinfo = PEDFile

  print *, 'PEDFile Name:', trim(pedinfo%FileName)
  print *, 'Delimiter:', trim(pedinfo%Delim_char)
  print *, 'Header lines:', pedinfo%Header
  print *, 'Number of Variables:', pedinfo%NVAR
  print *, 'Field Names:'
  do i = 1, pedinfo%NVAR
    print *, '  ', i, ':', trim(pedinfo%FieldName(i)), '(Loc:', pedinfo%FieldLoc(i), ')'
  end do

end program test_readpar