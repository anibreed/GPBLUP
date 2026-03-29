module M_readpar
  use M_kinds
  use M_Variables
  use M_StrEdit
  use M_ReadFile
  implicit none
  
  type :: sentence
    character(len=LEN_STR) :: keyword ! Keyword
    character(len=MAX_RECL), dimension(:), allocatable :: lines ! text on lines reading
    integer :: n_lines
  end type sentence

  type(sentence), save :: ParSections(5)
  character(len=LEN_STR) :: NAMEKEY(5)
  DATA NAMEKEY /'PEDFILE:', 'DATAFILE:', 'MAPFILE:', 'FRFILE:', 'GENOFILE:'/ 

  ! Public interface  
  public :: read_parameters, read_to_sentences, sentence, ParSections, NAMEKEY, M_readpar_get_popqc_thresholds
CONTAINS

  ! Helper function to remove comments starting with '#'
  function remove_comment(input_line) result(output_line)
    character(len=*), intent(in) :: input_line
    character(len=len(input_line)) :: output_line
    integer :: hash_pos
    
    ! Find position of '#' character
    hash_pos = index(input_line, '#')
    
    if (hash_pos > 0) then
      ! Comment found - keep only text before '#'
      output_line = input_line(1:hash_pos-1)
    else
      ! No comment - return full line
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

    ! Initialize and allocate ParSections
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
      
      ! Remove comments (anything after '#')
      line_buffer = trim(adjustl(remove_comment(raw_line)))
      
      ! Skip empty lines
      if (len_trim(line_buffer) == 0) cycle

      ! Split line into XC to check for keyword
      call split_string(line_buffer, XC, n, ' ')
      if (n > 0) then
        keyword_upper = to_upper(trim(adjustl(XC(1))))

        ! Compatibility alias: REFFILE behaves as DATAFILE.
        if (trim(keyword_upper) == 'REFFILE:') keyword_upper = 'DATAFILE:'
        
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

    ! Start from deterministic defaults before overriding from parameter file.
    call init_FRQC(FRQC)
    call init_POPQC(POPQC)
    
    ! 1. 파일을 읽어 메모리에 섹션별(PEDFILE:, DATAFILE: 등)로 분류 저장
    call read_to_sentences(Par_file)
    
    ! 2. 분류된 각 섹션별로 전용 파서를 호출하여 가독성과 명확성 증대
    do i = 1, 5
        if (ParSections(i)%n_lines <= 0) cycle
        
        select case(trim(ParSections(i)%keyword))
            case("PEDFILE:")
                call parse_ped_section(ParSections(i))
            case("DATAFILE:")
                call parse_data_section(ParSections(i))
            case("MAPFILE:")
                call parse_map_section(ParSections(i))
            case("FRFILE:")
                call parse_fr_section(ParSections(i))
            case("GENOFILE:")
                call parse_geno_section(ParSections(i))
        end select
    end do
    
    print *, "INFO: Parameter parsing completed using dedicated section parsers."
  end subroutine read_parameters

  ! --- 섹션별 전용 파서 정의 ---

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

  ! 메모리에 저장된 sentence 내용을 분석하는 통합 루틴
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
    
    ! 첫번째 라인에서 파일 이름 추출 (예: PEDFILE: test.ped)
    info%FileName = trim(extract_value_after_keyword(sec%lines(1)))

    ! 나머지 라인들 파싱
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
            case("HEADER:")
                call decode_line(line, XC, XI)
                info%Header = XI(2)
            case("DELIM:")
                info%Delim_char = trim(adjustl(XC(2)))
            case("NO_VARIABLES:")
                call decode_line(line, XC, XI)
                info%NVAR = XI(2)
                do i = 1, info%NVAR
                    L = L + 1
                    call split_string(sec%lines(L), XC, n, ' ')
                    info%FieldName(i) = to_upper(trim(adjustl(XC(2))))
                    ! 정수로 변환하여 위치 저장
                    read(XC(1), *) info%FieldLoc(i)
                end do
            case("FRQC:")
                call parse_frqc_memory(sec, L, FR_QC)
                L = L - 1 ! 서브루틴에서 증가된 L을 보정하여 루프 끝의 L = L + 1과 맞춤
            case("POPQC:")
                call parse_popqc_memory(sec, L, POP_QC)
                L = L - 1
            case("OUTPUT:")
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
    
    current_L = current_L + 1
    call split_string(sec%lines(current_L), XC, n, ' ')
    if (to_upper(trim(XC(1))) == "NO_VARIABLES:") then
        read(XC(2), *) nvar
        info%NVAR = nvar
        ! Optional thresholds: if omitted in FRQC section, skip checks by default.
        info%gt_score = 0.0
        info%cluster_sep = 0.0
        do i = 1, nvar
            current_L = current_L + 1
            call split_string(sec%lines(current_L), XC, n, ' ')
            keyword_upper = to_upper(trim(XC(1)))
            XR(2) = 0.0
            XR(3) = 0.0
            read(XC(2), *, iostat=ios) XR(2)
            if (n >= 3) read(XC(3), *, iostat=ios) XR(3)
            
            select case(trim(keyword_upper))
                case("R_INTENSITY:") ; info%r_intensity(1)=XR(2) ; info%r_intensity(2)=XR(3)
                case("GC_SCORE:")    ; info%gc_score=XR(2)
                case("GT_SCORE:")    ; info%gt_score=XR(2)
                case("CLUSTER_SEP:") ; info%cluster_sep=XR(2)
              case("ANIMAL_CR:", "CALL_RATE:", "CALL_RATE_THRESHOLD:") ; info%animal_cr=XR(2)
            end select
        end do
    end if
  end subroutine parse_frqc_memory

  subroutine parse_popqc_memory(sec, current_L, info)
    type(sentence), intent(in) :: sec
    integer, intent(inout) :: current_L
    type(POPQC_matrices), intent(inout) :: info
    character(LEN=LEN_STR) :: XC(MAX_VAR), keyword_upper
    integer :: n, i, nvar, ios
    
    current_L = current_L + 1
    call split_string(sec%lines(current_L), XC, n, ' ')
    if (to_upper(trim(XC(1))) == "NO_VARIABLES:") then
        read(XC(2), *) nvar
        info%NVAR = nvar
        do i = 1, nvar
            current_L = current_L + 1
            call split_string(sec%lines(current_L), XC, n, ' ')
            keyword_upper = to_upper(trim(XC(1)))
            
            select case(trim(keyword_upper))
                case("SNP_CR:")        ; read(XC(2), *) info%snp_cr
                case("ANIMAL_CR:")     ; read(XC(2), *) info%animal_cr
                case("MAF:")           ; read(XC(2), *) info%maf
                case("MIN_MAF:")       ; read(XC(2), *) info%maf
                case("HWE:")           ; info%HWE_method=trim(adjustl(XC(2))) ; read(XC(3), *) info%HWE_threshold
                case("HWE_THRESHOLD:") ; read(XC(2), *) info%HWE_threshold
                case("HWE_METHOD:")    ; info%HWE_method=trim(adjustl(XC(2)))
                case("HWE_MIN_N:")
                  read(XC(2), *, iostat=ios) info%HWE_min_n
                  if (ios /= 0 .or. info%HWE_min_n < 1) info%HWE_min_n = 20
                case("HWE_MIN_MAC:")
                  read(XC(2), *, iostat=ios) info%HWE_min_mac
                  if (ios /= 0 .or. info%HWE_min_mac < 0) info%HWE_min_mac = 5
                case("SEX_CONCORDANCE:")
                  read(XC(2), *, iostat=ios) info%sex_concordence_threshold
                  if (ios == 0) then
                    info%sex_check = (info%sex_concordence_threshold > 0.0)
                  end if
                case("MENDELIAN_INCONSISTENCY:")
                  read(XC(2), *, iostat=ios) info%mendel_snp_threshold
                  if (ios == 0) then
                    info%mendel_viol_threshold = info%mendel_snp_threshold
                    if (n >= 3) then
                      read(XC(3), *, iostat=ios) info%mendel_viol_threshold
                      if (ios /= 0) info%mendel_viol_threshold = info%mendel_snp_threshold
                    end if
                  end if
                case("PARENT_SEEK:")
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
                case("PARENT_SEEK_NSNPS:", "PARENT_SEEK_NSNP:", "PARENT_SEEK_SNPS:")
                  read(XC(2), *, iostat=ios) info%parent_seek_nsnps
                  if (ios /= 0 .or. info%parent_seek_nsnps < 1) then
                    info%parent_seek_nsnps = 500
                  end if
                case("PARENT_SEEK_G_TOL:", "PARENT_SEEK_G_TOLERANCE:")
                  read(XC(2), *, iostat=ios) info%parent_seek_g_tol
                  if (ios /= 0 .or. info%parent_seek_g_tol <= 0.0) then
                    info%parent_seek_g_tol = 0.15
                  end if
                case("ID_MERGE_SAMPLE_FRACTION:")
                  read(XC(2), *, iostat=ios) info%id_merge_sample_fraction
                  if (ios /= 0 .or. info%id_merge_sample_fraction <= 0.0 .or. &
                      info%id_merge_sample_fraction > 1.0) then
                    info%id_merge_sample_fraction = 0.20
                  end if
                case("ID_MERGE_CONFIRM_MATCH_THRESHOLD:")
                  read(XC(2), *, iostat=ios) info%id_merge_confirm_match_threshold
                  if (ios /= 0 .or. info%id_merge_confirm_match_threshold <= 0.0 .or. &
                      info%id_merge_confirm_match_threshold > 1.0) then
                    info%id_merge_confirm_match_threshold = 0.995
                  end if
                case("ID_SWAP_COMPAT_MIN:", "ID_SWAP_COMPAT_MINIMUM:")
                  read(XC(2), *, iostat=ios) info%id_swap_compat_min
                  if (ios /= 0 .or. info%id_swap_compat_min <= 0.0 .or. &
                      info%id_swap_compat_min > 1.0) then
                    info%id_swap_compat_min = 0.95
                  end if
                case("ID_SWAP_MAX_DEPTH:")
                  read(XC(2), *, iostat=ios) info%id_swap_max_depth
                  if (ios /= 0 .or. info%id_swap_max_depth < 1) then
                    info%id_swap_max_depth = 3
                  end if
            end select
        end do
    end if
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
            case("PREFIX:") ; info%FileName = trim(extract_value_after_keyword(sec%lines(current_L)))
            case("HEADER:") ; read(XC(2), *) info%Header
            case("DELIM:")  ; info%Delim_char = trim(XC(2))
            case("NO_VARIABLES:", "NO_VAIABLES:")
                if (n >= 2) read(XC(2), *) nvar
                info%NVAR = nvar
                do i = 1, nvar
                    current_L = current_L + 1
                    if (current_L > sec%n_lines) exit
                    call split_string(sec%lines(current_L), XC, n, ' ')
                    if (n >= 2) then
                        info%FieldName(i) = to_upper(trim(XC(1)))
                        read(XC(2), *) info%FieldLoc(i)
                    end if
                end do
                return
            case default
                current_L = current_L - 1
                return
        end select
    end do
  end subroutine parse_output_memory

  ! Helper routine to decode line into strings and integers (simplified replacement for readline)
  subroutine decode_line(line, XC, XI)
    character(len=*), intent(in) :: line
    character(len=LEN_STR), intent(out) :: XC(:)
    integer(kind=ki4), intent(out) :: XI(:)
    integer :: n, i
    call split_string(line, XC, n, ' ')
    XI = 0
    do i = 1, n
        if (len_trim(XC(i)) > 0) then
            ! Check if it's a number
            if (verify(trim(XC(i)), "0123456789+-") == 0) then
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
        print *, "  -> First line read failed (EOF or I/O error)"
      case(2)
        print *, "  -> NO_VARIABLES value is invalid or zero"
      case(3)
        print *, "  -> Memory allocation failed for field arrays"
      case(4)
        print *, "  -> Field definition line read failed"
      case(5)
        print *, "  -> NO_VARIABLES keyword missing from input"
      case default
        print *, "  -> Unknown error"
    end select
  end subroutine print_error_message

subroutine set_default_parameters()
    ! Set minimal default parameters when file cannot be read
    write(*,'(A)') "INFO: Setting default parameter values..."
    
    ! Initialize file info with empty values
    PEDFile%FileName = ''
    DATAFile%FileName = ''
    MAPFile%FileName = ''
    FRFile%FileName = ''
    GENOFile%FileName = ''
    
    ! Initialize quality control modules with safe defaults
    call init_FRQC(FRQC)
    call init_POPQC(POPQC)
    
    write(*,'(A)') "  -> Default parameters loaded. User input may be required."
  end subroutine set_default_parameters

  subroutine M_readpar_get_popqc_thresholds(out_snp_cr, out_animal_cr, out_maf, out_hwe_method, out_hwe_threshold, &
                                            out_hwe_min_n, out_hwe_min_mac)
    ! Export POPQC thresholds for population QC analysis
    ! These values are populated by read_parameters() from the parameter file
    ! Default values are set by init_POPQC() if parameter file is not provided
    
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
   
end module M_readpar
