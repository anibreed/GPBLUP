module readpar_mod
  implicit none

  integer, parameter :: RP_STR = 512

  type :: relped_params_t
    character(len=RP_STR) :: pedfile = ''
    character(len=RP_STR) :: reffile = ''
    integer :: pedinfo(5) = 0
    integer :: refinfo(2) = 0
    integer :: inbreed = 0
    integer :: relation = 0
    integer :: relation_inv = 0
    integer :: omp_threads = 32
    logical :: has_reference = .false.
  end type relped_params_t

  public :: relped_params_t, read_relped_params

contains

  subroutine read_relped_params(paramfile, params)
    character(len=*), intent(in) :: paramfile
    type(relped_params_t), intent(inout) :: params

    character(len=RP_STR) :: line, key
    character(len=RP_STR) :: tokens(64)
    character(len=RP_STR) :: base_dir
    integer :: unt, ios, ntok

    base_dir = get_basedir(paramfile)

    open(newunit=unt, file=trim(paramfile), status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'Cannot open parameter file: ', trim(paramfile)
      stop
    end if

    do
      read(unt, '(A)', iostat=ios) line
      if (ios /= 0) exit

      line = strip_comment(line)
      if (len_trim(line) == 0) cycle

      call split_ws(trim(adjustl(line)), tokens, ntok)
      if (ntok <= 0) cycle

      key = to_upper(trim(adjustl(tokens(1))))

      select case (key)
      case ('PEDFILE:')
        if (ntok < 2) then
          write(*,*) 'Pedigree file was not given'
          stop
        end if
        params%pedfile = resolve_path(base_dir, tokens(2))
        write(*,*) 'Pedigree file=  ', trim(params%pedfile)

      case ('PEDINFO:')
        if (ntok < 3) then
          write(*,*) 'Pedigree information was not enough'
          stop
        end if
        params%pedinfo(1:ntok-1) = to_int_slice(tokens(2:ntok))
        write(*,*) 'Column info on pedigree file(ARN, SRN, DRN, BYR, SEX)=', params%pedinfo

      case ('DATAFILE:', 'REFFILE:')
        if (ntok < 2) then
          write(*,*) 'Reference file was not given'
          stop
        end if
        params%reffile = resolve_path(base_dir, tokens(2))
        params%has_reference = .true.
        write(*,*) 'Data file name for renumbering=  ', trim(params%reffile)

      case ('DATAINFO:', 'REFINFO:')
        if (ntok < 3) then
          write(*,*) 'Data file information was not enough'
          stop
        end if
        params%refinfo(1) = to_int(tokens(2))
        params%refinfo(2) = to_int(tokens(3))
        write(*,*) 'Column infos for animal ID and Sex on data file=', params%refinfo

      case ('INBREED:')
        if (ntok < 2) then
          write(*,*) 'Information for calculating inb. was not enough'
          stop
        end if
        params%inbreed = to_int(tokens(2))
        if (params%inbreed < 0 .or. params%inbreed > 1) then
          write(*,*) 'Information for calculating inb. should be 0 or 1', params%inbreed
          stop
        end if
        write(*,*) 'Calculating inb. (No=0, Yes=1)', params%inbreed

      case ('RELATION', 'RELATION:')
        if (ntok < 2) then
          write(*,*) 'Information for calculating A matrix was not enough'
          stop
        end if
        params%relation = to_int(tokens(2))
        if (params%relation < 0 .or. params%relation > 1) then
          write(*,*) 'Information for calculating A matrix should be 0 or 1', params%relation
          stop
        end if
        write(*,*) 'Calculating A matrix (No=0, Yes=1)', params%relation

      case ('RELATION_INV', 'RELATION_INV:', 'AINV', 'AINV:')
        if (ntok < 2) then
          write(*,*) 'Information for calculating A inverse matrix was not enough'
          stop
        end if
        params%relation_inv = to_int(tokens(2))
        if (params%relation_inv < 0 .or. params%relation_inv > 1) then
          write(*,*) 'Information for calculating A inverse matrix should be 0 or 1', params%relation_inv
          stop
        end if
        write(*,*) 'Calculating A inverse matrix (No=0, Yes=1)', params%relation_inv

      case ('OMP_THREADS', 'OMP_THREADS:')
        if (ntok < 2) then
          write(*,*) 'Information for OMP threads was not enough'
          stop
        end if
        params%omp_threads = to_int(tokens(2))
        if (params%omp_threads < 1) then
          write(*,*) 'OMP_THREADS should be >= 1', params%omp_threads
          stop
        end if
        write(*,*) 'OpenMP thread count=', params%omp_threads

      case default
        cycle
      end select
    end do

    close(unt)
  end subroutine read_relped_params

  function get_basedir(path) result(basedir)
    character(len=*), intent(in) :: path
    character(len=RP_STR) :: basedir
    character(len=RP_STR) :: tmp
    integer :: p

    tmp = trim(adjustl(path))
    p = index(tmp, '/', back=.true.)
    if (p > 0) then
      basedir = tmp(1:p)
    else
      basedir = './'
    end if
  end function get_basedir

  function resolve_path(base_dir, path) result(outpath)
    character(len=*), intent(in) :: base_dir
    character(len=*), intent(in) :: path
    character(len=RP_STR) :: outpath
    character(len=RP_STR) :: tmp

    tmp = trim(adjustl(path))
    if (len_trim(tmp) <= 0) then
      outpath = tmp
    else if (tmp(1:1) == '/') then
      outpath = tmp
    else
      outpath = trim(base_dir)//tmp
    end if
  end function resolve_path

  function strip_comment(line) result(out)
    character(len=*), intent(in) :: line
    character(len=RP_STR) :: out
    integer :: p

    out = trim(adjustl(line))
    p = index(out, '#')
    if (p > 0) then
      if (p > 1) then
        out = trim(adjustl(out(1:p-1)))
      else
        out = ''
      end if
    end if
  end function strip_comment

  function to_upper(s) result(u)
    character(len=*), intent(in) :: s
    character(len=len(s)) :: u
    integer :: i, code

    u = s
    do i = 1, len(s)
      code = iachar(u(i:i))
      if (code >= iachar('a') .and. code <= iachar('z')) then
        u(i:i) = achar(code - 32)
      end if
    end do
  end function to_upper

  subroutine split_ws(line, tokens, ntok)
    character(len=*), intent(in) :: line
    character(len=*), intent(out) :: tokens(:)
    integer, intent(out) :: ntok

    integer :: i, n, start

    ntok = 0
    do i = 1, size(tokens)
      tokens(i) = ''
    end do

    n = len_trim(line)
    i = 1
    do while (i <= n)
      do while (i <= n .and. (line(i:i) == ' ' .or. line(i:i) == achar(9)))
        i = i + 1
      end do
      if (i > n) exit

      start = i
      do while (i <= n .and. line(i:i) /= ' ' .and. line(i:i) /= achar(9))
        i = i + 1
      end do

      ntok = ntok + 1
      if (ntok <= size(tokens)) tokens(ntok) = line(start:i-1)
    end do
  end subroutine split_ws

  function to_int(token) result(v)
    character(len=*), intent(in) :: token
    integer :: v
    integer :: ios

    read(token, *, iostat=ios) v
    if (ios /= 0) then
      write(*,*) 'Invalid integer in parameter file: ', trim(token)
      stop
    end if
  end function to_int

  function to_int_slice(tokens) result(vals)
    character(len=*), intent(in) :: tokens(:)
    integer :: vals(size(tokens))
    integer :: i

    do i = 1, size(tokens)
      vals(i) = to_int(tokens(i))
    end do
  end function to_int_slice

end module readpar_mod
