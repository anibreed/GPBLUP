module M_readpar_relped
  use M_kinds
  use M_Variables
  use M_renped_config, only: ren_outfile_global
  use M_StrEdit
  implicit none

  integer, parameter, public :: FMT_BINARY = 1
  integer, parameter, public :: FMT_TEXT = 2
  integer, parameter, public :: FMT_BOTH = 3

  type :: relped_params_t
    character(len=MAX_RECL) :: pedfile = ''
    character(len=MAX_RECL) :: reffile = ''
    integer :: pedinfo(5) = 0
    integer :: refinfo(2) = 0
    integer :: inbreed = 0
    integer :: relation = 0
    integer :: relation_inv = 0
    character(len=MAX_RECL) :: relation_file = 'A.mat'
    character(len=MAX_RECL) :: relation_inv_file = 'Ainv.mat'
    character(len=MAX_RECL) :: renumber_file = ''
    integer :: omp_threads = 32
    logical :: has_reference = .false.
  end type relped_params_t

  public :: relped_params_t, read_relped_params

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

  subroutine read_relped_params(paramfile, params)
    character(len=*), intent(in) :: paramfile
    type(relped_params_t), intent(inout) :: params

    character(len=MAX_RECL) :: line, raw_line, base_dir, line_compact
    character(len=LEN_STR) :: tokens(MAX_VAR), key
    integer :: unt, ios, ntok, i
    logical :: in_input_block
    character(len=LEN_STR) :: input_target
    integer :: input_remaining
    integer :: omp_tmp
    logical :: ok_int
    integer :: rel_block

    params = relped_params_t()
    params%inbreed = 1
    ren_outfile_global = ''
    base_dir = get_basedir_relped(paramfile)
    in_input_block = .false.
    input_target = ''
    input_remaining = 0
    rel_block = 0

    open(newunit=unt, file=trim(paramfile), status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(*,*) 'Cannot open parameter file: ', trim(paramfile)
      stop
    end if

    do
      read(unt, '(A)', iostat=ios) raw_line
      if (ios /= 0) exit

      line = trim(adjustl(remove_comment(raw_line)))
      line = strip_cr_relped(line)
      if (len_trim(line) == 0) cycle

      line_compact = line
      if (len_trim(line_compact) > 0) then
        line_compact = adjustl(line_compact)
        line_compact = replace_tabs_relped(line_compact)
        line_compact = strip_cr_relped(line_compact)
      end if
      call remove_duplicate_spaces(line_compact, line)
      call split_string(line, tokens, ntok, ' ')
      if (ntok <= 0) cycle

      if (index(to_upper(line), 'OMP_THREADS') > 0) then
        omp_tmp = 0
        call parse_colon_int_relped(line, omp_tmp)
        if (omp_tmp <= 0 .and. ntok >= 2) then
          do i = 2, ntok
            call try_read_int_relped(tokens(i), omp_tmp, ok_int)
            if (ok_int .and. omp_tmp > 0) exit
          end do
        end if
        if (omp_tmp > 0) params%omp_threads = omp_tmp
        cycle
      end if

      if (input_remaining > 0) then
        if (ntok >= 2) then
          select case (trim(to_upper(tokens(2))))
          case ('ANIMAL_ID', 'ID')
            if (input_target == 'PED') params%pedinfo(1) = to_int_relped(tokens(1))
            if (input_target == 'REF') params%refinfo(1) = to_int_relped(tokens(1))
          case ('SIRE_ID', 'SIRE')
            if (input_target == 'PED') params%pedinfo(2) = to_int_relped(tokens(1))
          case ('DAM_ID', 'DAM')
            if (input_target == 'PED') params%pedinfo(3) = to_int_relped(tokens(1))
          case ('BIRTH', 'BIRTH_YEAR', 'BYEAR', 'BDATE', 'BIRTH_DATE')
            if (input_target == 'PED') params%pedinfo(4) = to_int_relped(tokens(1))
          case ('SEX')
            if (input_target == 'PED') params%pedinfo(5) = to_int_relped(tokens(1))
            if (input_target == 'REF') params%refinfo(2) = to_int_relped(tokens(1))
          case default
            continue
          end select
        end if
        input_remaining = input_remaining - 1
        cycle
      end if

      key = to_upper(trim(adjustl(tokens(1))))

      if (trim(key) == 'INPUT:') then
        in_input_block = .true.
        cycle
      end if
      if (trim(key) == 'NO_VARIABLES:' .and. in_input_block) then
        if (ntok >= 2) input_remaining = to_int_relped(tokens(2))
        in_input_block = .false.
        cycle
      end if
      if (trim(key) == 'HEADER:' .or. trim(key) == 'DELIM:') then
        cycle
      end if
      if (trim(key) == 'OUTPUT:' .or. trim(key) == 'OPTION:' .or. trim(key) == 'FILEFORMAT:') then
        rel_block = 0
        cycle
      end if
      select case (key)
      case ('FILENAME:')
        if (ntok >= 2) then
          params%renumber_file = resolve_path_relped(base_dir, strip_optional_quotes(tokens(2)))
          ren_outfile_global = params%renumber_file
        end if

      case ('PEDFILE:')
        if (ntok < 2) then
          write(*,*) 'Pedigree file was not given'
          stop
        end if
        params%pedfile = resolve_path_relped(base_dir, strip_optional_quotes(tokens(2)))
        input_target = 'PED'

      case ('PEDINFO:')
        if (ntok < 3) then
          write(*,*) 'Pedigree information was not enough'
          stop
        end if
        params%pedinfo = 0
        do i = 2, min(ntok, 6)
          params%pedinfo(i-1) = to_int_relped(tokens(i))
        end do

      case ('DATAFILE:', 'REFFILE:')
        if (ntok < 2) then
          write(*,*) 'Reference file was not given'
          stop
        end if
        params%reffile = resolve_path_relped(base_dir, strip_optional_quotes(tokens(2)))
        params%has_reference = .true.
        input_target = 'REF'

      case ('DATAINFO:', 'REFINFO:')
        if (ntok < 3) then
          write(*,*) 'Data file information was not enough'
          stop
        end if
        params%refinfo(1) = to_int_relped(tokens(2))
        params%refinfo(2) = to_int_relped(tokens(3))

      case ('INBREED:')
        if (ntok >= 2) then
          params%inbreed = 1
        end if

      case ('RELATIONSHIP', 'RELATIONSHIP:')
        rel_block = 1
        if (ntok >= 2) then
          params%relation = to_int_relped(tokens(2))
        end if

      case ('INV_RELATIONSHIP', 'INV_RELATIONSHIP:')
        rel_block = 2
        if (ntok >= 2) then
          params%relation_inv = to_int_relped(tokens(2))
        end if

      case ('RELATION', 'RELATION:')
        if (ntok < 2) then
          write(*,*) 'Information for calculating A matrix was not enough'
          stop
        end if
        params%relation = to_int_relped(tokens(2))
        if (params%relation < 0 .or. params%relation > 1) then
          write(*,*) 'Information for calculating A matrix should be 0 or 1', params%relation
          stop
        end if

      case ('RELATION_FILE', 'RELATION_FILE:')
        if (ntok >= 2) then
          params%relation_file = resolve_path_relped(base_dir, strip_optional_quotes(tokens(2)))
        end if

      case ('FILE', 'FILE:')
        if (ntok >= 2) then
          if (rel_block == 1) then
            params%relation_file = resolve_path_relped(base_dir, strip_optional_quotes(tokens(2)))
          else if (rel_block == 2) then
            params%relation_inv_file = resolve_path_relped(base_dir, strip_optional_quotes(tokens(2)))
          end if
        end if

      case ('RELATION_FORMAT', 'RELATION_FORMAT:', 'A_FORMAT', 'A_FORMAT:', 'OUTPUT_FORMAT', 'OUTPUT_FORMAT:')
        if (ntok >= 2) then
          select case (parse_output_format(tokens(2)))
          case (FMT_TEXT)
            if (index(params%relation_file, '.txt') == 0) then
              params%relation_file = trim(params%relation_file) // '.txt'
            end if
            if (index(params%relation_inv_file, '.txt') == 0) then
              params%relation_inv_file = trim(params%relation_inv_file) // '.txt'
            end if
          case (FMT_BOTH)
            if (index(params%relation_file, '.both') == 0) then
              params%relation_file = trim(params%relation_file) // '.both'
            end if
            if (index(params%relation_inv_file, '.both') == 0) then
              params%relation_inv_file = trim(params%relation_inv_file) // '.both'
            end if
          case default
            continue
          end select
        end if

      case ('FORMAT', 'FORMAT:')
        if (ntok >= 2) then
          select case (parse_output_format(tokens(2)))
          case (FMT_TEXT)
            if (rel_block == 1) then
              if (index(params%relation_file, '.txt') == 0) then
                params%relation_file = trim(params%relation_file) // '.txt'
              end if
            else if (rel_block == 2) then
              if (index(params%relation_inv_file, '.txt') == 0) then
                params%relation_inv_file = trim(params%relation_inv_file) // '.txt'
              end if
            end if
          case (FMT_BOTH)
            if (rel_block == 1) then
              if (index(params%relation_file, '.both') == 0) then
                params%relation_file = trim(params%relation_file) // '.both'
              end if
            else if (rel_block == 2) then
              if (index(params%relation_inv_file, '.both') == 0) then
                params%relation_inv_file = trim(params%relation_inv_file) // '.both'
              end if
            end if
          case default
            continue
          end select
        end if

      case ('RELATION_INV', 'RELATION_INV:', 'AINV', 'AINV:')
        if (ntok < 2) then
          write(*,*) 'Information for calculating A inverse matrix was not enough'
          stop
        end if
        params%relation_inv = to_int_relped(tokens(2))
        if (params%relation_inv < 0 .or. params%relation_inv > 1) then
          write(*,*) 'Information for calculating A inverse matrix should be 0 or 1', params%relation_inv
          stop
        end if

      case ('RELATION_INV_FILE', 'RELATION_INV_FILE:', 'AINV_FILE', 'AINV_FILE:')
        if (ntok >= 2) then
          params%relation_inv_file = resolve_path_relped(base_dir, strip_optional_quotes(tokens(2)))
        end if

      case ('RELATION_INV_FORMAT', 'RELATION_INV_FORMAT:', 'AINV_FORMAT', 'AINV_FORMAT:')
        if (ntok >= 2) then
          select case (parse_output_format(tokens(2)))
          case (FMT_TEXT)
            if (index(params%relation_inv_file, '.txt') == 0) then
              params%relation_inv_file = trim(params%relation_inv_file) // '.txt'
            end if
          case (FMT_BOTH)
            if (index(params%relation_inv_file, '.both') == 0) then
              params%relation_inv_file = trim(params%relation_inv_file) // '.both'
            end if
          case default
            continue
          end select
        end if

      case ('OMP_THREADS', 'OMP_THREADS:')
        if (ntok >= 2) then
          params%omp_threads = to_int_relped(tokens(2))
        else
          call parse_colon_int_relped(line, params%omp_threads)
        end if
        if (params%omp_threads < 1) then
          write(*,*) 'OMP_THREADS should be >= 1', params%omp_threads
          stop
        end if

      case default
        cycle
      end select
    end do

    params%relation_file = resolve_path_relped(base_dir, params%relation_file)
    params%relation_inv_file = resolve_path_relped(base_dir, params%relation_inv_file)

    close(unt)
  end subroutine read_relped_params

  function get_basedir_relped(path) result(basedir)
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
  end function get_basedir_relped

  function resolve_path_relped(base_dir, path) result(outpath)
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
  end function resolve_path_relped

  integer function to_int_relped(token) result(v)
    character(len=*), intent(in) :: token
    integer :: ios

    read(token, *, iostat=ios) v
    if (ios /= 0) then
      write(*,*) 'Invalid integer in parameter file: ', trim(token)
      stop
    end if
  end function to_int_relped

  subroutine parse_colon_int_relped(line, outv)
    character(len=*), intent(in) :: line
    integer, intent(out) :: outv
    integer :: p, ios
    character(len=MAX_RECL) :: rhs

    outv = 0
    p = index(line, ':')
    if (p <= 0) return
    rhs = trim(adjustl(line(p+1:)))
    if (len_trim(rhs) <= 0) return
    read(rhs, *, iostat=ios) outv
    if (ios /= 0) outv = 0
  end subroutine parse_colon_int_relped

  subroutine try_read_int_relped(token, outv, ok)
    character(len=*), intent(in) :: token
    integer, intent(out) :: outv
    logical, intent(out) :: ok
    integer :: ios

    read(token, *, iostat=ios) outv
    ok = (ios == 0)
    if (.not. ok) outv = 0
  end subroutine try_read_int_relped

  function replace_tabs_relped(text) result(out)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: out
    integer :: i

    out = text
    do i = 1, len(out)
      if (iachar(out(i:i)) == 9) out(i:i) = ' '
    end do
  end function replace_tabs_relped

  function strip_cr_relped(text) result(out)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: out
    integer :: i

    out = text
    do i = 1, len(out)
      if (iachar(out(i:i)) == 13) out(i:i) = ' '
    end do
  end function strip_cr_relped

end module M_readpar_relped
