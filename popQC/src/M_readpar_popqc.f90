module M_readpar_popqc
  use M_kinds
  use M_Variables
  use M_StrEdit
  use M_ReadFile
  implicit none

  type :: sentence
    character(len=LEN_STR) :: keyword
    character(len=MAX_RECL), dimension(:), allocatable :: lines
    integer :: n_lines
  end type sentence

  type(sentence), save :: ParSections(5)
  character(len=LEN_STR) :: NAMEKEY(5)
  DATA NAMEKEY /'PEDFILE:', 'DATAFILE:', 'MAPFILE:', 'FRFILE:', 'GENOFILE:'/ 

  public :: read_parameters, read_to_sentences, sentence, ParSections, NAMEKEY, M_readpar_get_popqc_thresholds

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

  subroutine read_to_sentences(Par_file)
    character(len=*), intent(in) :: Par_file
    character(LEN=MAX_RECL) :: line_buffer
    character(LEN=LEN_STR) :: XC(MAX_VAR), keyword_upper
    integer :: n, unitF, i, section_idx, ios
    character(len=MAX_RECL) :: raw_line

    unitF = fopen(Par_file)
    if (unitF < 0) return

    do i = 1, 5
        ParSections(i)%keyword = NAMEKEY(i)
        ParSections(i)%n_lines = 0
        if (allocated(ParSections(i)%lines)) deallocate(ParSections(i)%lines)
        allocate(ParSections(i)%lines(500))
        ParSections(i)%lines = ''
    end do

    section_idx = 0
    rewind(unitF)

    do
      read(unitF, '(A)', iostat=ios) raw_line
      if (ios /= 0) exit

      line_buffer = trim(adjustl(remove_comment(raw_line)))
      if (len_trim(line_buffer) == 0) cycle

      call split_string(line_buffer, XC, n, ' ')
      if (n > 0) then
        keyword_upper = to_upper(trim(adjustl(XC(1))))

        if (trim(keyword_upper) == 'REFFILE:' .or. trim(keyword_upper) == 'DATAFILE:') &
           keyword_upper = 'DATAFILE:'

        do i = 1, 5
           if (trim(keyword_upper) == trim(NAMEKEY(i))) then
               section_idx = i
               exit
           end if
        end do
      end if

      if (section_idx > 0) then
        ParSections(section_idx)%n_lines = ParSections(section_idx)%n_lines + 1
        ParSections(section_idx)%lines(ParSections(section_idx)%n_lines) = line_buffer
      end if
    end do
    close(unitF)
  end subroutine read_to_sentences

  subroutine read_parameters(Par_file)
    character(len=*), intent(in) :: Par_file
    integer :: i

    call init_FRQC(FRQC)
    call init_POPQC(POPQC)
    call init_species_config()

    call read_global_options(Par_file)

    call read_to_sentences(Par_file)

    do i = 1, 5
        if (ParSections(i)%n_lines <= 0) cycle

        select case(trim(ParSections(i)%keyword))
            case('PEDFILE:')
                call parse_ped_section(ParSections(i))
            case('DATAFILE:')
                call parse_data_section(ParSections(i))
            case('MAPFILE:')
                call parse_map_section(ParSections(i))
            case('FRFILE:')
                call parse_fr_section(ParSections(i))
            case('GENOFILE:')
                call parse_geno_section(ParSections(i))
        end select
    end do

    print *, 'INFO: Parameter parsing completed using dedicated section parsers.'
  end subroutine read_parameters

  subroutine read_global_options(par_file)
    character(len=*), intent(in) :: par_file
    character(len=MAX_RECL) :: raw_line, line, value_text
    character(len=LEN_STR) :: XC(MAX_VAR), key
    integer :: n, u_fd, ios
    integer :: autosome_max, chr_x, chr_y, chr_xy, chr_z, chr_w
    logical :: use_zw, any_override
    character(len=LEN_STR) :: species_token

    autosome_max = SPECIES_AUTOSOME_MAX
    chr_x = SPECIES_CHR_X
    chr_y = SPECIES_CHR_Y
    chr_xy = SPECIES_CHR_XY
    chr_z = SPECIES_CHR_Z
    chr_w = SPECIES_CHR_W
    use_zw = SPECIES_USE_ZW
    any_override = .false.

    species_token = ''
    u_fd = fopen(par_file)
    if (u_fd < 0) return

    do
      read(u_fd, '(A)', iostat=ios) raw_line
      if (ios /= 0) exit
      line = trim(adjustl(remove_comment(raw_line)))
      if (len_trim(line) == 0) cycle
      call split_string(line, XC, n, ' ')
      if (n <= 0) cycle
      key = to_upper(trim(adjustl(XC(1))))
      if (trim(key) == 'SPECIES:' .or. trim(key) == 'SPECIES') then
        if (n >= 2) species_token = XC(2)
        exit
      end if
    end do
    close(u_fd)

    if (len_trim(species_token) > 0) then
      call apply_species_defaults(species_token)
      autosome_max = SPECIES_AUTOSOME_MAX
      chr_x = SPECIES_CHR_X
      chr_y = SPECIES_CHR_Y
      chr_xy = SPECIES_CHR_XY
      chr_z = SPECIES_CHR_Z
      chr_w = SPECIES_CHR_W
      use_zw = SPECIES_USE_ZW
    end if

    u_fd = fopen(par_file)
    if (u_fd < 0) return

    do
      read(u_fd, '(A)', iostat=ios) raw_line
      if (ios /= 0) exit
      line = trim(adjustl(remove_comment(raw_line)))
      if (len_trim(line) == 0) cycle
      call split_string(line, XC, n, ' ')
      if (n <= 0) cycle
      key = to_upper(trim(adjustl(XC(1))))

      select case (trim(key))
      case ('AUTOSOME_MAX:', 'AUTOSOME_MAX')
        if (n >= 2) then
          autosome_max = to_int_safe(XC(2))
          any_override = .true.
        end if
      case ('SEXCHR_X:', 'SEXCHR_X', 'X_CHR', 'CHR_X')
        if (n >= 2) then
          chr_x = to_int_safe(XC(2))
          any_override = .true.
        end if
      case ('SEXCHR_Y:', 'SEXCHR_Y', 'Y_CHR', 'CHR_Y')
        if (n >= 2) then
          chr_y = to_int_safe(XC(2))
          any_override = .true.
        end if
      case ('SEXCHR_XY:', 'SEXCHR_XY', 'PAR_CHR', 'CHR_XY')
        if (n >= 2) then
          chr_xy = to_int_safe(XC(2))
          any_override = .true.
        end if
      case ('SEXCHR_Z:', 'SEXCHR_Z', 'Z_CHR', 'CHR_Z')
        if (n >= 2) then
          chr_z = to_int_safe(XC(2))
          any_override = .true.
        end if
      case ('SEXCHR_W:', 'SEXCHR_W', 'W_CHR', 'CHR_W')
        if (n >= 2) then
          chr_w = to_int_safe(XC(2))
          any_override = .true.
        end if
      case ('USE_ZW:', 'USE_ZW')
        if (n >= 2) then
          use_zw = to_logical_safe(XC(2))
          any_override = .true.
        end if
      case ('INCLUDE_SEXCHR:', 'INCLUDE_SEXCHR', 'SEXCHR_INCLUDE', 'USE_SEXCHR')
        if (n >= 2) SPECIES_INCLUDE_SEXCHR = to_logical_safe(XC(2))
      case ('SEXCHR_OVERRIDE:', 'SEXCHR_OVERRIDE')
        value_text = extract_value_after_keyword(line)
        call parse_sexchr_override(value_text, autosome_max, chr_x, chr_y, chr_xy, chr_z, chr_w, use_zw, any_override)
      end select
    end do

    close(u_fd)

    if (any_override) then
      call set_species_config(autosome_max, chr_x, chr_y, chr_xy, chr_z, chr_w, use_zw)
    end if
  end subroutine read_global_options

  subroutine parse_sexchr_override(text, autosome_max, chr_x, chr_y, chr_xy, chr_z, chr_w, use_zw, any_override)
    character(len=*), intent(in) :: text
    integer, intent(inout) :: autosome_max, chr_x, chr_y, chr_xy, chr_z, chr_w
    logical, intent(inout) :: use_zw
    logical, intent(inout) :: any_override
    character(len=MAX_RECL) :: work
    character(len=LEN_STR) :: tokens(MAX_VAR)
    character(len=LEN_STR) :: key, val
    integer :: n, i, eqpos

    work = trim(adjustl(text))
    if (len_trim(work) == 0) return
    do i = 1, len(work)
      if (work(i:i) == ',') work(i:i) = ' '
    end do

    call split_string(work, tokens, n, ' ')
    do i = 1, n
      if (len_trim(tokens(i)) == 0) cycle
      eqpos = index(tokens(i), '=')
      if (eqpos <= 1) cycle
      key = to_upper(trim(adjustl(tokens(i)(1:eqpos-1))))
      val = trim(adjustl(tokens(i)(eqpos+1:)))

      select case (trim(key))
      case ('AUTOSOME_MAX')
        autosome_max = to_int_safe(val)
        any_override = .true.
      case ('X')
        chr_x = to_int_safe(val)
        any_override = .true.
      case ('Y')
        chr_y = to_int_safe(val)
        any_override = .true.
      case ('XY', 'PAR')
        chr_xy = to_int_safe(val)
        any_override = .true.
      case ('Z')
        chr_z = to_int_safe(val)
        any_override = .true.
      case ('W')
        chr_w = to_int_safe(val)
        any_override = .true.
      case ('USE_ZW')
        use_zw = to_logical_safe(val)
        any_override = .true.
      end select
    end do
  end subroutine parse_sexchr_override

  integer function to_int_safe(token) result(val)
    character(len=*), intent(in) :: token
    integer :: ios

    val = 0
    read(token, *, iostat=ios) val
    if (ios /= 0) val = 0
  end function to_int_safe

  logical function to_logical_safe(token) result(val)
    character(len=*), intent(in) :: token
    character(len=LEN_STR) :: up

    up = to_upper(trim(adjustl(token)))
    select case (trim(up))
    case ('1', 'T', 'TRUE', 'YES', 'Y')
      val = .true.
    case default
      val = .false.
    end select
  end function to_logical_safe

  subroutine parse_ped_section(sec)
    type(sentence), intent(in) :: sec
    call parse_section_info(sec, PEDFile)
  end subroutine parse_ped_section

  subroutine parse_data_section(sec)
    type(sentence), intent(in) :: sec
    call parse_section_info(sec, DATAFile)
  end subroutine parse_data_section

  subroutine parse_map_section(sec)
    type(sentence), intent(in) :: sec
    call parse_section_info(sec, MAPFile)
  end subroutine parse_map_section

  subroutine parse_fr_section(sec)
    type(sentence), intent(in) :: sec
    call parse_section_info(sec, FRFile, FR_QC=FRQC, outinfo=FROUTFile)
    call init_pointers(FROUTFile)
  end subroutine parse_fr_section

  subroutine parse_geno_section(sec)
    type(sentence), intent(in) :: sec
    call parse_section_info(sec, GENOFile, POP_QC=POPQC, outinfo=GENOUTFile)
    call init_pointers(GENOUTFile)
  end subroutine parse_geno_section

  subroutine parse_section_info(sec, info, FR_QC, POP_QC, outinfo)
    type(sentence), intent(in) :: sec
    type(FileInfo), intent(inout) :: info
    type(FRQC_matrices), intent(inout), optional :: FR_QC
    type(POPQC_matrices), intent(inout), optional :: POP_QC
    type(FileInfo), intent(inout), optional :: outinfo

    character(LEN=MAX_RECL) :: line
    character(LEN=LEN_STR) :: XC(MAX_VAR), keyword_upper
    integer(kind=ki4) :: XI(MAX_VAR)
    integer :: n, i, L

    info%FileName = trim(extract_value_after_keyword(sec%lines(1)))

    L = 2
    do while (L <= sec%n_lines)
        line = sec%lines(L)
        call split_string(line, XC, n, ' ')
        if (n <= 0) then
            L = L + 1
            cycle
        end if

        keyword_upper = to_upper(trim(adjustl(XC(1))))

        select case(trim(keyword_upper))
          case('INPUT:')
            continue
            case('HEADER:')
                call decode_line(line, XC, XI)
                info%Header = XI(2)
            case('DELIM:')
                info%Delim_char = trim(adjustl(XC(2)))
            case('NO_VARIABLES:')
                call decode_line(line, XC, XI)
                info%NVAR = XI(2)
                do i = 1, info%NVAR
                  L = L + 1
                  call split_string(sec%lines(L), XC, n, ' ')
                  if (n < 2) cycle

                  if (verify(trim(adjustl(XC(1))), '0123456789+-') == 0) then
                    read(XC(1), *) info%FieldLoc(i)
                    info%FieldName(i) = to_upper(trim(adjustl(XC(2))))
                  else if (verify(trim(adjustl(XC(2))), '0123456789+-') == 0) then
                    read(XC(2), *) info%FieldLoc(i)
                    info%FieldName(i) = to_upper(trim(adjustl(XC(1))))
                  else
                    cycle
                  end if

                  if (index(info%FieldName(i), ':', back=.true.) == len_trim(info%FieldName(i))) then
                    info%FieldName(i) = info%FieldName(i)(1:len_trim(info%FieldName(i))-1)
                  end if
                end do
            case('FRQC:')
                call parse_frqc_memory(sec, L, FR_QC)
                L = L - 1
            case('POPQC:')
                call parse_popqc_memory(sec, L, POP_QC)
                L = L - 1
            case('OUTPUT:')
                if (present(outinfo)) then
                    call parse_output_memory(sec, L, outinfo)
                    L = L - 1
                end if
        end select
        L = L + 1
    end do
  end subroutine parse_section_info

  subroutine parse_frqc_memory(sec, current_L, info)
    type(sentence), intent(in) :: sec
    integer, intent(inout) :: current_L
    type(FRQC_matrices), intent(inout) :: info
    character(LEN=LEN_STR) :: XC(MAX_VAR), keyword_upper
    real(kind=r4) :: XR(MAX_VAR)
    integer :: n, i, nvar, ios
    integer :: autosome_max, chr_x, chr_y, chr_xy, chr_z, chr_w
    logical :: use_zw, any_override

    current_L = current_L + 1
    if (current_L > sec%n_lines) return
    call split_string(sec%lines(current_L), XC, n, ' ')
    if (n <= 0) return

    info%gt_score = 0.0
    info%cluster_sep = 0.0
    autosome_max = SPECIES_AUTOSOME_MAX
    chr_x = SPECIES_CHR_X
    chr_y = SPECIES_CHR_Y
    chr_xy = SPECIES_CHR_XY
    chr_z = SPECIES_CHR_Z
    chr_w = SPECIES_CHR_W
    use_zw = SPECIES_USE_ZW
    any_override = .false.

    if (to_upper(trim(XC(1))) == 'NO_VARIABLES:') then
      read(XC(2), *) nvar
    else
      nvar = sec%n_lines - current_L + 1
      current_L = current_L - 1
    end if

    do i = 1, nvar
      current_L = current_L + 1
      if (current_L > sec%n_lines) exit
      call split_string(sec%lines(current_L), XC, n, ' ')
      if (n <= 0) cycle
      keyword_upper = to_upper(trim(XC(1)))
      if (trim(keyword_upper) == 'OUTPUT:' .or. trim(keyword_upper) == 'INPUT:' .or. &
          trim(keyword_upper) == 'FRQC:' .or. trim(keyword_upper) == 'POPQC:' .or. &
          trim(keyword_upper) == 'PEDFILE:' .or. trim(keyword_upper) == 'MAPFILE:' .or. &
          trim(keyword_upper) == 'FRFILE:' .or. trim(keyword_upper) == 'GENOFILE:' .or. &
          trim(keyword_upper) == 'DATAFILE:') then
        current_L = current_L - 1
        exit
      end if
      if (n < 2) cycle
      XR(2) = 0.0
      XR(3) = 0.0
      read(XC(2), *, iostat=ios) XR(2)
      if (n >= 3) read(XC(3), *, iostat=ios) XR(3)

      select case(trim(keyword_upper))
        case('R_INTENSITY:') ; info%r_intensity(1)=XR(2) ; info%r_intensity(2)=XR(3)
        case('GC_SCORE:')    ; info%gc_score=XR(2)
        case('GT_SCORE:')    ; info%gt_score=XR(2)
        case('CLUSTER_SEP:') ; info%cluster_sep=XR(2)
        case('ANIMAL_CR:', 'CALL_RATE:', 'CALL_RATE_THRESHOLD:') ; info%animal_cr=XR(2)
        case('SPECIES:')
          call apply_species_defaults(XC(2))
          autosome_max = SPECIES_AUTOSOME_MAX
          chr_x = SPECIES_CHR_X
          chr_y = SPECIES_CHR_Y
          chr_xy = SPECIES_CHR_XY
          chr_z = SPECIES_CHR_Z
          chr_w = SPECIES_CHR_W
          use_zw = SPECIES_USE_ZW
        case('INCLUDE_SEXCHR:', 'INCLUDE_SEXCHR')
          SPECIES_INCLUDE_SEXCHR = to_logical_safe(XC(2))
        case('SEXCHR_OVERRIDE:', 'SEXCHR_OVERRIDE')
          call parse_sexchr_override(extract_value_after_keyword(sec%lines(current_L)), &
              autosome_max, chr_x, chr_y, chr_xy, chr_z, chr_w, use_zw, any_override)
        case('AUTOSOME_MAX:', 'AUTOSOME_MAX')
          autosome_max = to_int_safe(XC(2))
          any_override = .true.
        case('SEXCHR_X:', 'SEXCHR_X', 'X_CHR', 'CHR_X')
          chr_x = to_int_safe(XC(2))
          any_override = .true.
        case('SEXCHR_Y:', 'SEXCHR_Y', 'Y_CHR', 'CHR_Y')
          chr_y = to_int_safe(XC(2))
          any_override = .true.
        case('SEXCHR_XY:', 'SEXCHR_XY', 'PAR_CHR', 'CHR_XY')
          chr_xy = to_int_safe(XC(2))
          any_override = .true.
        case('SEXCHR_Z:', 'SEXCHR_Z', 'Z_CHR', 'CHR_Z')
          chr_z = to_int_safe(XC(2))
          any_override = .true.
        case('SEXCHR_W:', 'SEXCHR_W', 'W_CHR', 'CHR_W')
          chr_w = to_int_safe(XC(2))
          any_override = .true.
        case('USE_ZW:', 'USE_ZW')
          use_zw = to_logical_safe(XC(2))
          any_override = .true.
      end select
    end do

    if (any_override) then
      call set_species_config(autosome_max, chr_x, chr_y, chr_xy, chr_z, chr_w, use_zw)
    end if
  end subroutine parse_frqc_memory

  subroutine parse_popqc_memory(sec, current_L, info)
    type(sentence), intent(in) :: sec
    integer, intent(inout) :: current_L
    type(POPQC_matrices), intent(inout) :: info
    character(LEN=LEN_STR) :: XC(MAX_VAR), keyword_upper
    integer :: n, i, nvar, ios

    current_L = current_L + 1
    if (current_L > sec%n_lines) return
    call split_string(sec%lines(current_L), XC, n, ' ')
    if (n <= 0) return

    if (to_upper(trim(XC(1))) == 'NO_VARIABLES:') then
      read(XC(2), *) nvar
    else
      nvar = sec%n_lines - current_L + 1
      current_L = current_L - 1
    end if

    info%NVAR = nvar
    do i = 1, nvar
      current_L = current_L + 1
      if (current_L > sec%n_lines) exit
      call split_string(sec%lines(current_L), XC, n, ' ')
      if (n <= 0) cycle
      keyword_upper = to_upper(trim(XC(1)))
      if (trim(keyword_upper) == 'OUTPUT:' .or. trim(keyword_upper) == 'INPUT:' .or. &
          trim(keyword_upper) == 'FRQC:' .or. trim(keyword_upper) == 'POPQC:' .or. &
          trim(keyword_upper) == 'PEDFILE:' .or. trim(keyword_upper) == 'MAPFILE:' .or. &
          trim(keyword_upper) == 'FRFILE:' .or. trim(keyword_upper) == 'GENOFILE:' .or. &
          trim(keyword_upper) == 'DATAFILE:') then
        current_L = current_L - 1
        exit
      end if
      if (n < 2) cycle

      select case(trim(keyword_upper))
                case('SNP_CR:')
                  read(XC(2), *, iostat=ios) info%snp_cr
                case('ANIMAL_CR:')
                  read(XC(2), *, iostat=ios) info%animal_cr
                case('MAF:')
                  read(XC(2), *, iostat=ios) info%maf
                case('MIN_MAF:')
                  read(XC(2), *, iostat=ios) info%maf
                case('HWE:')
                  info%HWE_method=trim(adjustl(XC(2)))
                  if (n >= 3) read(XC(3), *, iostat=ios) info%HWE_threshold
                case('HWE_THRESHOLD:')
                  read(XC(2), *, iostat=ios) info%HWE_threshold
                case('HWE_METHOD:')    ; info%HWE_method=trim(adjustl(XC(2)))
                case('HWE_MIN_N:')
                  read(XC(2), *, iostat=ios) info%HWE_min_n
                  if (ios /= 0 .or. info%HWE_min_n < 1) info%HWE_min_n = 20
                case('HWE_MIN_MAC:')
                  read(XC(2), *, iostat=ios) info%HWE_min_mac
                  if (ios /= 0 .or. info%HWE_min_mac < 0) info%HWE_min_mac = 5
                case('SEX_CONCORDANCE:')
                  read(XC(2), *, iostat=ios) info%sex_concordence_threshold
                  if (ios == 0) then
                    info%sex_check = (info%sex_concordence_threshold > 0.0)
                  end if
                case('MENDELIAN_INCONSISTENCY:')
                  read(XC(2), *, iostat=ios) info%mendel_snp_threshold
                  if (ios == 0) then
                    info%mendel_viol_threshold = info%mendel_snp_threshold
                    if (n >= 3) then
                      read(XC(3), *, iostat=ios) info%mendel_viol_threshold
                      if (ios /= 0) info%mendel_viol_threshold = info%mendel_snp_threshold
                    end if
                  end if
                case('PARENT_SEEK:')
                  read(XC(2), *, iostat=ios) info%parent_seek_threshold
                  if (ios == 0) then
                    info%parent_seek = (info%parent_seek_threshold > 0.0)
                  else
                    info%parent_seek = to_logical(XC(2))
                    if (info%parent_seek) then
                      info%parent_seek_threshold = 0.95
                    else
                      info%parent_seek_threshold = 0.0
                    end if
                  end if
                case('PARENT_SEEK_NSNPS:', 'PARENT_SEEK_NSNP:', 'PARENT_SEEK_SNPS:')
                  read(XC(2), *, iostat=ios) info%parent_seek_nsnps
                  if (ios /= 0 .or. info%parent_seek_nsnps < 1) then
                    info%parent_seek_nsnps = 500
                  end if
                case('PARENT_SEEK_G_TOL:', 'PARENT_SEEK_G_TOLERANCE:')
                  read(XC(2), *, iostat=ios) info%parent_seek_g_tol
                  if (ios /= 0 .or. info%parent_seek_g_tol <= 0.0) then
                    info%parent_seek_g_tol = 0.15
                  end if
                case('ID_MERGE_SAMPLE_FRACTION:')
                  read(XC(2), *, iostat=ios) info%id_merge_sample_fraction
                  if (ios /= 0 .or. info%id_merge_sample_fraction <= 0.0 .or. &
                      info%id_merge_sample_fraction > 1.0) then
                    info%id_merge_sample_fraction = 0.20
                  end if
                case('ID_MERGE_CONFIRM_MATCH_THRESHOLD:')
                  read(XC(2), *, iostat=ios) info%id_merge_confirm_match_threshold
                  if (ios /= 0 .or. info%id_merge_confirm_match_threshold <= 0.0 .or. &
                      info%id_merge_confirm_match_threshold > 1.0) then
                    info%id_merge_confirm_match_threshold = 0.995
                  end if
                case('ID_SWAP_COMPAT_MIN:', 'ID_SWAP_COMPAT_MINIMUM:')
                  read(XC(2), *, iostat=ios) info%id_swap_compat_min
                  if (ios /= 0 .or. info%id_swap_compat_min <= 0.0 .or. &
                      info%id_swap_compat_min > 1.0) then
                    info%id_swap_compat_min = 0.95
                  end if
                case('ID_SWAP_MAX_DEPTH:')
                  read(XC(2), *, iostat=ios) info%id_swap_max_depth
                  if (ios /= 0 .or. info%id_swap_max_depth < 1) then
                    info%id_swap_max_depth = 3
                  end if
        case('MULTI_RESOLVED_ACTION:')
          if (n >= 2) then
            select case (to_upper(trim(adjustl(XC(2)))))
            case ('DELETE', 'DEL', 'REMOVE', '1')
              info%multi_resolved_action = 1
            case default
              info%multi_resolved_action = 0
            end select
          end if
      end select
    end do
  end subroutine parse_popqc_memory

  subroutine parse_output_memory(sec, current_L, info)
    type(sentence), intent(in) :: sec
    integer, intent(inout) :: current_L
    type(FileInfo), intent(inout) :: info
    character(LEN=LEN_STR) :: XC(MAX_VAR), keyword_upper
    integer :: n, i, nvar

    do
        current_L = current_L + 1
        if (current_L > sec%n_lines) exit
        call split_string(sec%lines(current_L), XC, n, ' ')
        if (n <= 0) cycle

        keyword_upper = to_upper(trim(XC(1)))
        select case(trim(keyword_upper))
            case('PREFIX:') ; info%FileName = trim(extract_value_after_keyword(sec%lines(current_L)))
            case('HEADER:') ; read(XC(2), *) info%Header
            case('DELIM:')  ; info%Delim_char = trim(XC(2))
            case('NO_VARIABLES:', 'NO_VAIABLES:')
                if (n >= 2) read(XC(2), *) nvar
                info%NVAR = nvar
                do i = 1, nvar
                  current_L = current_L + 1
                  if (current_L > sec%n_lines) exit
                  call split_string(sec%lines(current_L), XC, n, ' ')
                  if (n < 2) cycle

                  if (verify(trim(adjustl(XC(1))), '0123456789+-') == 0) then
                    read(XC(1), *) info%FieldLoc(i)
                    info%FieldName(i) = to_upper(trim(adjustl(XC(2))))
                  else if (verify(trim(adjustl(XC(2))), '0123456789+-') == 0) then
                    read(XC(2), *) info%FieldLoc(i)
                    info%FieldName(i) = to_upper(trim(adjustl(XC(1))))
                  else
                    cycle
                  end if

                  if (index(info%FieldName(i), ':', back=.true.) == len_trim(info%FieldName(i))) then
                    info%FieldName(i) = info%FieldName(i)(1:len_trim(info%FieldName(i))-1)
                  end if
                end do
                return
            case default
                current_L = current_L - 1
                return
        end select
    end do
  end subroutine parse_output_memory

  subroutine decode_line(line, XC, XI)
    character(len=*), intent(in) :: line
    character(len=LEN_STR), intent(out) :: XC(:)
    integer(kind=ki4), intent(out) :: XI(:)
    integer :: n, i
    call split_string(line, XC, n, ' ')
    XI = 0
    do i = 1, n
        if (len_trim(XC(i)) > 0) then
            if (verify(trim(XC(i)), '0123456789+-') == 0) then
                read(XC(i), *, iostat=n) XI(i)
            end if
        end if
    end do
  end subroutine decode_line

  integer function FindPos(str, Mstr)
    character(len=*) :: str
    character(len=*) :: Mstr(:)
    integer :: k
    do k = 1, size(Mstr)
      if (trim(str) == trim(Mstr(k))) then
        FindPos = k
        return
      endif
    enddo
    FindPos = 0
  end function FindPos

  integer function FieldPos(str, Info)
    character(len=*) :: str
    type(FileInfo), intent(in) :: Info
    integer :: k
    do k = 1, Info%NVAR
      if (trim(str) == trim(Info%FieldName(k))) then
        FieldPos = Info%FieldLoc(k)
        return
      endif
    enddo
    FieldPos = 0
  end function FieldPos

  subroutine print_error_message(err_code)
    integer, intent(in) :: err_code
    select case(err_code)
      case(1)
        print *, '  -> First line read failed (EOF or I/O error)'
      case(2)
        print *, '  -> NO_VARIABLES value is invalid or zero'
      case(3)
        print *, '  -> Memory allocation failed for field arrays'
      case(4)
        print *, '  -> Field definition line read failed'
      case(5)
        print *, '  -> NO_VARIABLES keyword missing from input'
      case default
        print *, '  -> Unknown error'
    end select
  end subroutine print_error_message

  subroutine set_default_parameters()
    write(*,'(A)') 'INFO: Setting default parameter values...'

    PEDFile%FileName = ''
    DATAFile%FileName = ''
    MAPFile%FileName = ''
    FRFile%FileName = ''
    GENOFile%FileName = ''

    call init_FRQC(FRQC)
    call init_POPQC(POPQC)

    write(*,'(A)') '  -> Default parameters loaded. User input may be required.'
  end subroutine set_default_parameters

  subroutine M_readpar_get_popqc_thresholds(out_snp_cr, out_animal_cr, out_maf, out_hwe_method, out_hwe_threshold, &
                                            out_hwe_min_n, out_hwe_min_mac)
    real, intent(out) :: out_snp_cr, out_animal_cr, out_maf, out_hwe_threshold
    character(len=LEN_STR), intent(out) :: out_hwe_method
    integer, intent(out), optional :: out_hwe_min_n, out_hwe_min_mac

    out_snp_cr = POPQC%snp_cr
    out_animal_cr = POPQC%animal_cr
    out_maf = POPQC%maf
    out_hwe_method = POPQC%HWE_method
    out_hwe_threshold = POPQC%HWE_threshold
    if (present(out_hwe_min_n)) out_hwe_min_n = POPQC%HWE_min_n
    if (present(out_hwe_min_mac)) out_hwe_min_mac = POPQC%HWE_min_mac
  end subroutine M_readpar_get_popqc_thresholds

end module M_readpar_popqc
