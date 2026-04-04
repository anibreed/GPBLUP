module M_readpar_relgeno
  use M_Kinds
  use M_Variables
  use M_StrEdit
  implicit none

  integer, parameter, public :: FMT_BINARY = 1
  integer, parameter, public :: FMT_TEXT = 2
  integer, parameter, public :: FMT_BOTH = 3

  type :: relgeno_params_t
    character(len=MAX_RECL) :: genofile = ''
    character(len=MAX_RECL) :: outfile_g = ''
    character(len=MAX_RECL) :: outfile_ginv = ''
    character(len=MAX_RECL) :: outfile_gdiag = ''
    character(len=MAX_RECL) :: outfile_gib = ''
    integer :: missing_code = 5
    real(kind=r8) :: weight = 0.95_r8
    real(kind=r8) :: diag_weight = -1.0_r8
    logical :: compute_g = .true.
    logical :: compute_ginv = .true.
    integer :: output_format = FMT_BINARY
    integer :: header_lines = 0
  end type relgeno_params_t

  public :: relgeno_params_t, read_relgeno_params

contains

  function remove_comment(input_line) result(output_line)
    character(len=*), intent(in) :: input_line
    character(len=len(input_line)) :: output_line
    integer :: hash_pos

    hash_pos = index(input_line, '#')
    if (hash_pos > 0) then
      output_line = input_line(1:hash_pos-1)
    else
      output_line = input_line
    end if
  end function remove_comment

  function strip_optional_quotes(input_text) result(output_text)
    character(len=*), intent(in) :: input_text
    character(len=len(input_text)) :: output_text
    integer :: lt

    output_text = trim(adjustl(input_text))
    lt = len_trim(output_text)
    if (lt >= 2) then
      if ((output_text(1:1) == '"' .and. output_text(lt:lt) == '"') .or. &
          (output_text(1:1) == achar(39) .and. output_text(lt:lt) == achar(39))) then
        output_text = output_text(2:lt-1)
      end if
    end if
  end function strip_optional_quotes

  function extract_value_after_keyword(input_line) result(value_text)
    character(len=*), intent(in) :: input_line
    character(len=MAX_RECL) :: value_text
    character(len=MAX_RECL) :: work_line
    integer :: pos_space, pos_tab, pos

    value_text = ''
    work_line = trim(adjustl(input_line))

    pos_space = index(work_line, ' ')
    pos_tab = index(work_line, achar(9))

    if (pos_space == 0 .and. pos_tab == 0) return

    if (pos_space == 0) then
      pos = pos_tab
    else if (pos_tab == 0) then
      pos = pos_space
    else
      pos = min(pos_space, pos_tab)
    end if

    if (pos < len_trim(work_line)) then
      value_text = trim(adjustl(work_line(pos+1:)))
      value_text = trim(strip_optional_quotes(value_text))
    end if
  end function extract_value_after_keyword

  integer function parse_output_format(token) result(fmt)
    character(len=*), intent(in) :: token
    character(len=LEN_STR) :: up

    fmt = FMT_BINARY
    up = to_upper(trim(adjustl(token)))
    select case (trim(up))
    case ('1', 'BINARY', 'BIN')
      fmt = FMT_BINARY
    case ('2', 'TEXT', 'TXT', 'ASCII')
      fmt = FMT_TEXT
    case ('3', 'BOTH', 'ALL')
      fmt = FMT_BOTH
    case default
      fmt = FMT_BINARY
    end select
  end function parse_output_format

  subroutine parse_g_option_line(line, params)
    character(len=*), intent(in) :: line
    type(relgeno_params_t), intent(inout) :: params

    character(len=LEN_STR) :: tokens(MAX_VAR), key
    integer :: ntok

    call split_string(line, tokens, ntok, ' ')
    if (ntok <= 0) return

    key = to_upper(trim(adjustl(tokens(1))))
    select case (key)
    case ('MISSING:', 'MISSING_CODE:')
      if (ntok >= 2) params%missing_code = to_int_relgeno(tokens(2))
    case ('WEIGHT:', 'G_WEIGHT:')
      if (ntok >= 2) params%weight = to_real_relgeno(tokens(2))
    case ('DIAG_WEIGHT:', 'DIAGONAL_WEIGHT:')
      if (ntok >= 2) params%diag_weight = to_real_relgeno(tokens(2))
    case ('COMPUTE_G:')
      if (ntok >= 2) params%compute_g = to_logical(tokens(2))
    case ('COMPUTE_GINV:')
      if (ntok >= 2) params%compute_ginv = to_logical(tokens(2))
    case ('FORMAT:', 'OUTPUT_FORMAT:')
      if (ntok >= 2) params%output_format = parse_output_format(tokens(2))
    case default
      return
    end select
  end subroutine parse_g_option_line

  subroutine read_relgeno_params(paramfile, params)
    character(len=*), intent(in) :: paramfile
    type(relgeno_params_t), intent(inout) :: params

    character(len=MAX_RECL) :: line, raw_line, base_dir
    character(len=LEN_STR) :: tokens(MAX_VAR), key
    character(len=MAX_RECL) :: output_prefix
    integer :: unt, ios, ntok, gopt_remaining
    logical :: in_output_block, has_output_prefix, in_g_options_block, in_input_block

    params = relgeno_params_t()
    base_dir = get_basedir_relgeno(paramfile)
    output_prefix = ''
    in_output_block = .false.
    has_output_prefix = .false.
    in_g_options_block = .false.
    in_input_block = .false.
    gopt_remaining = 0

    open(newunit=unt, file=trim(paramfile), status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'Cannot open parameter file: ', trim(paramfile)
      stop
    end if

    do
      read(unt, '(A)', iostat=ios) raw_line
      if (ios /= 0) exit

      line = trim(adjustl(remove_comment(raw_line)))
      if (len_trim(line) == 0) cycle

      call split_string(line, tokens, ntok, ' ')
      if (ntok <= 0) cycle

      key = to_upper(trim(adjustl(tokens(1))))

      if (trim(key) == 'OUTPUT:') then
        in_output_block = .true.
        cycle
      end if

      if (trim(key) == 'INPUT:') then
        in_input_block = .true.
        cycle
      end if

      if (trim(key) == 'G_OPTIONS:') then
        in_g_options_block = .true.
        gopt_remaining = 0
        cycle
      end if

      if (in_g_options_block) then
        if (trim(key) == 'NO_VARIABLES:' .and. ntok >= 2) then
          gopt_remaining = to_int_relgeno(tokens(2))
        end if

        do while (gopt_remaining > 0)
          read(unt, '(A)', iostat=ios) raw_line
          if (ios /= 0) then
            write(*,*) 'Unexpected EOF while reading G_OPTIONS block'
            stop
          end if

          line = trim(adjustl(remove_comment(raw_line)))
          if (len_trim(line) == 0) cycle
          call parse_g_option_line(line, params)
          gopt_remaining = gopt_remaining - 1
        end do

        in_g_options_block = .false.
        cycle
      end if

      if (in_output_block) then
        select case (trim(key))
        case ('PREFIX:')
          output_prefix = extract_value_after_keyword(line)
          if (len_trim(output_prefix) > 0) has_output_prefix = .true.
          cycle
        case ('FORMAT:', 'OUTPUT_FORMAT:')
          if (ntok >= 2) params%output_format = parse_output_format(tokens(2))
          cycle
        case ('OUT_G:', 'OUT_GFILE:', 'G_OUT:')
          if (ntok >= 2) params%outfile_g = resolve_path_relgeno(base_dir, strip_optional_quotes(tokens(2)))
          cycle
        case ('OUT_GINV:', 'OUT_GI:', 'GINV_OUT:')
          if (ntok >= 2) params%outfile_ginv = resolve_path_relgeno(base_dir, strip_optional_quotes(tokens(2)))
          cycle
        case ('OUT_GDIAG:', 'OUT_G_DIAG:')
          if (ntok >= 2) params%outfile_gdiag = resolve_path_relgeno(base_dir, strip_optional_quotes(tokens(2)))
          cycle
        case ('OUT_GIB:', 'OUT_GIBREED:', 'OUT_GIBREEDING:')
          if (ntok >= 2) params%outfile_gib = resolve_path_relgeno(base_dir, strip_optional_quotes(tokens(2)))
          cycle
        case ('HEADER:', 'DELIM:', 'NO_VARIABLES:')
          cycle
        case default
          in_output_block = .false.
        end select
      end if

      if (in_input_block) then
        select case (trim(key))
        case ('HEADER:')
          if (ntok >= 2) params%header_lines = to_int_relgeno(tokens(2))
          cycle
        case ('DELIM:', 'NO_VARIABLES:')
          cycle
        case default
          in_input_block = .false.
        end select
      end if

      select case (key)
      case ('GENOFILE:')
        params%genofile = resolve_path_relgeno(base_dir, extract_value_after_keyword(line))
        if (len_trim(params%genofile) == 0 .and. ntok >= 2) then
          params%genofile = resolve_path_relgeno(base_dir, strip_optional_quotes(tokens(2)))
        end if
        if (len_trim(params%genofile) == 0) then
          write(*,*) 'GENOFILE was not given'
          stop
        end if

      case ('OUT_G:', 'OUT_GFILE:', 'G_OUT:')
        if (ntok >= 2) params%outfile_g = resolve_path_relgeno(base_dir, strip_optional_quotes(tokens(2)))

      case ('OUT_GINV:', 'OUT_GI:', 'GINV_OUT:')
        if (ntok >= 2) params%outfile_ginv = resolve_path_relgeno(base_dir, strip_optional_quotes(tokens(2)))

      case ('OUT_GDIAG:', 'OUT_G_DIAG:')
        if (ntok >= 2) params%outfile_gdiag = resolve_path_relgeno(base_dir, strip_optional_quotes(tokens(2)))

      case ('OUT_GIB:', 'OUT_GIBREED:', 'OUT_GIBREEDING:')
        if (ntok >= 2) params%outfile_gib = resolve_path_relgeno(base_dir, strip_optional_quotes(tokens(2)))

      case ('FORMAT:', 'OUTPUT_FORMAT:')
        if (ntok >= 2) params%output_format = parse_output_format(tokens(2))

      case ('MISSING:', 'MISSING_CODE:')
        if (ntok < 2) then
          write(*,*) 'MISSING code not given'
          stop
        end if
        params%missing_code = to_int_relgeno(tokens(2))

      case ('WEIGHT:', 'G_WEIGHT:')
        if (ntok < 2) then
          write(*,*) 'WEIGHT not given'
          stop
        end if
        params%weight = to_real_relgeno(tokens(2))

      case ('DIAG_WEIGHT:', 'DIAGONAL_WEIGHT:')
        if (ntok < 2) then
          write(*,*) 'DIAG_WEIGHT not given'
          stop
        end if
        params%diag_weight = to_real_relgeno(tokens(2))

      case ('COMPUTE_G:')
        if (ntok >= 2) params%compute_g = to_logical(tokens(2))

      case ('COMPUTE_GINV:')
        if (ntok >= 2) params%compute_ginv = to_logical(tokens(2))

      case ('HEADER:')
        if (ntok >= 2) params%header_lines = to_int_relgeno(tokens(2))

      case default
        cycle
      end select
    end do
    close(unt)

    if (has_output_prefix .and. len_trim(params%genofile) == 0) then
      params%genofile = resolve_path_relgeno(base_dir, trim(output_prefix) // '_cleaned.txt')
    end if

    if (len_trim(params%genofile) == 0) then
      write(*,*) 'GENOFILE must be provided in parameter file'
      stop
    end if

    if (params%diag_weight < 0.0_r8) then
      params%diag_weight = 1.0_r8 - params%weight
    end if

    if (len_trim(params%outfile_g) == 0) then
      params%outfile_g = trim(params%genofile) // '.G'
    end if
    if (len_trim(params%outfile_ginv) == 0) then
      params%outfile_ginv = trim(params%genofile) // '.Gi'
    end if
    if (len_trim(params%outfile_gdiag) == 0) then
      params%outfile_gdiag = trim(params%genofile) // '.Gdiag'
    end if
    if (len_trim(params%outfile_gib) == 0) then
      params%outfile_gib = trim(params%genofile) // '.Gib'
    end if
  end subroutine read_relgeno_params

  function get_basedir_relgeno(path) result(basedir)
    character(len=*), intent(in) :: path
    character(len=MAX_RECL) :: basedir
    character(len=MAX_RECL) :: tmp
    integer :: p

    tmp = trim(adjustl(path))
    p = index(tmp, '/', back=.true.)
    if (p > 0) then
      basedir = tmp(1:p)
    else
      basedir = './'
    end if
  end function get_basedir_relgeno

  function resolve_path_relgeno(base_dir, path) result(outpath)
    character(len=*), intent(in) :: base_dir
    character(len=*), intent(in) :: path
    character(len=MAX_RECL) :: outpath
    character(len=MAX_RECL) :: tmp

    tmp = trim(adjustl(path))
    if (len_trim(tmp) <= 0) then
      outpath = tmp
    else if (tmp(1:1) == '/') then
      outpath = tmp
    else
      outpath = trim(base_dir)//tmp
    end if
  end function resolve_path_relgeno

  integer function to_int_relgeno(token) result(v)
    character(len=*), intent(in) :: token
    integer :: ios

    read(token, *, iostat=ios) v
    if (ios /= 0) then
      write(*,*) 'Invalid integer in parameter file: ', trim(token)
      stop
    end if
  end function to_int_relgeno

  real(kind=r8) function to_real_relgeno(token) result(v)
    character(len=*), intent(in) :: token
    integer :: ios

    read(token, *, iostat=ios) v
    if (ios /= 0) then
      write(*,*) 'Invalid real in parameter file: ', trim(token)
      stop
    end if
  end function to_real_relgeno

end module M_readpar_relgeno
