module M_IO_Geno
  use M_Kinds
  use M_Variables
  use M_ped, only: CPED, init_cped_record
  use M_ReadFile
  use M_StrEdit
  use M_Hashtable
  use M_IO_Ped, only: PEDHashTable, pht_search
  implicit none

  private
  public :: load_geno_file_mod
  public :: write_valid_geno_for_imputation_mod
  public :: set_qcgenoio_debug_mod

  logical, save :: qcgenoio_debug = .false.

contains

  subroutine set_qcgenoio_debug_mod(in_debug)
    logical, intent(in) :: in_debug
    qcgenoio_debug = in_debug
  end subroutine set_qcgenoio_debug_mod

  subroutine load_geno_file_mod(in_n_animals, in_map_records, in_geno_key_field, out_n_snps, &
                                io_geno_matrix, io_geno_animal_ids)
    integer, intent(in) :: in_n_animals
    type(SNPInfo), intent(in) :: in_map_records(:)
    character(len=*), intent(in) :: in_geno_key_field
    integer, intent(out) :: out_n_snps
    integer(kind=ki1), allocatable, intent(inout) :: io_geno_matrix(:,:)
    character(len=LEN_STR), allocatable, intent(inout) :: io_geno_animal_ids(:)

    integer :: i, j, u_fd, k
    character(len=:), allocatable :: a
    integer :: io_stat, geno_len
    integer :: geno_col_idx, geno_col_loc
    integer :: geno_key_col_idx, geno_key_col_loc
    character(len=1) :: dlm
    logical :: found_field
    integer :: cv
    integer :: a_len, geno_s, geno_e, key_s, key_e
    integer :: log_interval
    logical :: positions_cached

    if (len_trim(GENOFile%FileName) == 0) then
       print *, 'ERROR: GENOFile not specified.'
       out_n_snps = 0
       return
    end if

    out_n_snps = 0
    do i = 1, size(in_map_records)
      if (in_map_records(i)%Array_All > out_n_snps) out_n_snps = in_map_records(i)%Array_All
    end do

    geno_col_idx = 0
    do j = 1, GENOFile%NVAR
      if (trim(to_upper(GENOFile%FieldName(j))) == 'GENO') then
        geno_col_idx = j
        if (qcgenoio_debug) print *, '  Found GENO field in metadata index: ', geno_col_idx
        exit
      end if
    end do
    if (geno_col_idx == 0) then
      print *, "ERROR: 'GENO' field not found in GENOFILE definition."
      return
    end if
    geno_col_loc = GENOFile%FieldLoc(geno_col_idx)

    geno_key_col_idx = 0
    do j = 1, GENOFile%NVAR
      if (trim(to_upper(GENOFile%FieldName(j))) == trim(to_upper(in_geno_key_field))) then
        geno_key_col_idx = j
        exit
      end if
    end do
    if (geno_key_col_idx == 0) then
      print *, 'ERROR: GENO key field not found: ', trim(in_geno_key_field)
      return
    end if
    geno_key_col_loc = GENOFile%FieldLoc(geno_key_col_idx)

    if (allocated(io_geno_matrix)) deallocate(io_geno_matrix)
    allocate(io_geno_matrix(in_n_animals, out_n_snps))
    io_geno_matrix = 9_ki1
    if (allocated(io_geno_animal_ids)) deallocate(io_geno_animal_ids)
    allocate(io_geno_animal_ids(in_n_animals))
    io_geno_animal_ids = ''

    print *, '  Loading genotypes from: ', trim(GENOFile%FileName)
    dlm = FindDelim(trim(GENOFile%Delim_char))
    if (qcgenoio_debug) then
      print *, '  Delimiter: [', dlm, ']'
      print *, '  GENO data column position: ', geno_col_loc
      print *, '  GENO animal key field/position: ', trim(in_geno_key_field), ' / ', geno_key_col_loc
    end if
    u_fd = fopen(trim(GENOFile%FileName), recl=80000)

    allocate(character(len=80000) :: a)

    do i = 1, GENOFile%Header
      read(u_fd, '(a)', iostat=io_stat) a
    end do

    log_interval = max(100, in_n_animals / 10)
    positions_cached = .false.

    do i = 1, in_n_animals
      read(u_fd, '(a)', iostat=io_stat) a
      if (io_stat /= 0) then
        print *, '  ERROR: Failed to read line ', i, ' iostat=', io_stat
        exit
      end if

      if (mod(i, log_interval) == 0 .or. i == in_n_animals) then
        write(*, '(a,i0,a,i0,a,i0,a)') &
          '  GENO reading: ', i, ' / ', in_n_animals, &
          '  (', i * 100 / in_n_animals, '% done)'
      end if

      a_len = len_trim(a)

      if (.not. positions_cached) then
        call get_field_pos(a, a_len, dlm, geno_col_loc, geno_s, geno_e, found_field)
        if (.not. found_field) then
          print *, '  WARNING: Could not read GENO column at line ', i
          cycle
        end if
        call get_field_pos(a, a_len, dlm, geno_key_col_loc, key_s, key_e, found_field)
        positions_cached = .true.
      end if

      if (.not. found_field) cycle
      if (found_field .and. key_e >= key_s) io_geno_animal_ids(i) = a(key_s:key_e)

      if (qcgenoio_debug .and. i == 1) then
        print *, '  --- DEBUG: First Animal Key/GENO Extraction ---'
        if (found_field .and. key_e >= key_s) print *, '  Animal key: [', a(key_s:key_e), ']'
        geno_len = geno_e - geno_s + 1
        print *, '  GENO string length: ', geno_len
        if (geno_len > 0) print *, '  GENO first 50 chars: [', a(geno_s:min(geno_s+49, geno_e)), ']'
        print *, '  -----------------------------------------------'
      end if

      geno_len = min(geno_e - geno_s + 1, out_n_snps)
      do k = 1, geno_len
        cv = ichar(a(geno_s + k - 1 : geno_s + k - 1)) - 48
        if (cv >= 0 .and. cv <= 2) then
          io_geno_matrix(i, k) = int(cv, ki1)
        else
          io_geno_matrix(i, k) = 9_ki1
        end if
      end do
    end do

    close(u_fd)
    if (allocated(a)) deallocate(a)

    print *, '  Finished loading GENO matrix: ', in_n_animals, ' x ', out_n_snps
  end subroutine load_geno_file_mod

  subroutine write_valid_geno_for_imputation_mod(in_n_animals, in_n_snps, in_geno_matrix, in_geno_animal_ids, &
                                                 in_pht, in_n_snps_pass_qc, in_snp_qc_status, in_snp_is_autosome, &
                                                 in_action)
    integer, intent(in) :: in_n_animals, in_n_snps
    integer(kind=ki1), intent(in) :: in_geno_matrix(in_n_animals, in_n_snps)
    character(len=*), intent(in) :: in_geno_animal_ids(in_n_animals)
    type(PEDHashTable), intent(in) :: in_pht
    integer, intent(in) :: in_n_snps_pass_qc
    integer, intent(in), optional :: in_snp_qc_status(:)
    logical, intent(in), optional :: in_snp_is_autosome(:)
    integer, intent(in), optional :: in_action(:)

    integer :: i, j, s, out_unit
    integer :: nvar, n_written, n_excluded_delete, n_excluded_empty
    integer :: width, n_non_missing
    integer :: header_flag
    character(len=MAX_RECL) :: out_file
    character(len=1) :: out_delim
    character(len=LEN_STR) :: field_names(MAX_VAR)
    integer :: field_widths(MAX_VAR)
    character(len=LEN_STR) :: field_up
    character(len=MAX_STR) :: out_value
    logical :: use_output_block, has_geno_field, keep_row
    logical :: found
    type(CPED) :: ped_rec
    integer, allocatable :: auto_col(:)
    integer :: n_auto, buf_pos
    character(len=in_n_snps) :: geno_buf

    n_auto = in_n_snps
    if (present(in_snp_is_autosome) .and. present(in_snp_qc_status)) then
      if (size(in_snp_is_autosome) >= in_n_snps .and. size(in_snp_qc_status) >= in_n_snps) then
        n_auto = count(in_snp_is_autosome(1:in_n_snps) .and. in_snp_qc_status(1:in_n_snps) == 0)
      end if
    else if (present(in_snp_qc_status)) then
      if (size(in_snp_qc_status) >= in_n_snps) n_auto = count(in_snp_qc_status(1:in_n_snps) == 0)
    end if

    allocate(auto_col(n_auto))
    buf_pos = 0
    if (present(in_snp_is_autosome) .and. present(in_snp_qc_status)) then
      if (size(in_snp_is_autosome) >= in_n_snps .and. size(in_snp_qc_status) >= in_n_snps) then
        do s = 1, in_n_snps
          if (in_snp_is_autosome(s) .and. in_snp_qc_status(s) == 0) then
            buf_pos = buf_pos + 1
            auto_col(buf_pos) = s
          end if
        end do
      end if
    else if (present(in_snp_qc_status)) then
      if (size(in_snp_qc_status) >= in_n_snps) then
        do s = 1, in_n_snps
          if (in_snp_qc_status(s) == 0) then
            buf_pos = buf_pos + 1
            auto_col(buf_pos) = s
          end if
        end do
      end if
    else
      do s = 1, in_n_snps
        buf_pos = buf_pos + 1
        auto_col(buf_pos) = s
      end do
    end if

    field_names = ''
    field_widths = 0
    nvar = 0

    use_output_block = (len_trim(GENOUTFile%FileName) > 0 .and. GENOUTFile%NVAR > 0)
    if (len_trim(GENOUTFile%FileName) > 0) then
      out_file = trim(GENOUTFile%FileName) // '_cleaned.txt'
    else
      out_file = 'popQC_cleaned.txt'
    end if

    if (use_output_block) then
      out_delim = FindDelim(trim(GENOUTFile%Delim_char))
      header_flag = GENOUTFile%Header
      nvar = GENOUTFile%NVAR
      if (nvar > MAX_VAR) nvar = MAX_VAR
      do j = 1, nvar
        field_names(j) = trim(GENOUTFile%FieldName(j))
        field_widths(j) = max(0, GENOUTFile%FieldLoc(j))
      end do
    else
      out_delim = FindDelim(trim(GENOFile%Delim_char))
      header_flag = GENOFile%Header
      if (header_flag <= 0) header_flag = 1
      nvar = GENOFile%NVAR
      if (nvar > MAX_VAR) nvar = MAX_VAR
      if (nvar <= 0) then
        nvar = 2
        field_names(1) = 'ANIMAL_ID'
        field_names(2) = 'GENO'
      else
        do j = 1, nvar
          field_names(j) = trim(GENOFile%FieldName(j))
          field_widths(j) = 0
        end do
      end if
    end if

    has_geno_field = .false.
    do j = 1, nvar
      field_up = to_upper(trim(field_names(j)))
      if (trim(field_up) == 'GENO') then
        has_geno_field = .true.
        exit
      end if
    end do
    if (.not. has_geno_field .and. nvar < MAX_VAR) then
      nvar = nvar + 1
      field_names(nvar) = 'GENO'
      field_widths(nvar) = 0
    end if

    open(newunit=out_unit, file=trim(out_file), status='replace', action='write')

    if (header_flag > 0) then
      do j = 1, nvar
        field_up = to_upper(trim(field_names(j)))
        select case (trim(field_up))
        case ('ID', 'ANIMAL_ID', 'ANIMAL_ARN', 'ARN')
          out_value = 'ANIMAL_ID'
          call write_fixed_popqc(out_unit, out_value, field_widths(j))
        case ('GENO')
          write(out_unit, '(A)', advance='no') 'GENO'
        case default
          out_value = trim(field_names(j))
          call write_fixed_popqc(out_unit, out_value, field_widths(j))
        end select
        if (j < nvar) write(out_unit, '(A1)', advance='no') out_delim
      end do
      write(out_unit, '()')
    end if

    n_written = 0
    n_excluded_delete = 0
    n_excluded_empty = 0

    do i = 1, in_n_animals
      keep_row = .true.
      if (present(in_action)) then
        if (i <= size(in_action)) then
          if (in_action(i) == 2) then
            keep_row = .false.
            n_excluded_delete = n_excluded_delete + 1
          end if
        end if
      end if
      if (keep_row) then
        if (len_trim(in_geno_animal_ids(i)) == 0) then
          keep_row = .false.
          n_excluded_empty = n_excluded_empty + 1
        end if
      end if
      if (keep_row) then
        n_non_missing = 0
        do s = 1, n_auto
          if (int(in_geno_matrix(i, auto_col(s))) /= 9) n_non_missing = n_non_missing + 1
        end do
        if (n_non_missing <= 0) then
          keep_row = .false.
          n_excluded_empty = n_excluded_empty + 1
        end if
      end if
      if (.not. keep_row) cycle

      found = pht_search(in_pht, trim(in_geno_animal_ids(i)), ped_rec)
      if (.not. found) call init_cped_record(ped_rec)

      do j = 1, nvar
        field_up = to_upper(trim(field_names(j)))
        select case (trim(field_up))
        case ('ANIMAL_ID', 'ID')
          out_value = trim(in_geno_animal_ids(i))
          call write_fixed_popqc(out_unit, out_value, field_widths(j))
        case ('ANIMAL_ARN', 'ARN')
          out_value = trim(ped_rec%ARN)
          call write_fixed_popqc(out_unit, out_value, field_widths(j))
        case ('BREED')
          out_value = trim(ped_rec%BREED)
          call write_fixed_popqc(out_unit, out_value, field_widths(j))
        case ('SIRE', 'SIRE_ID')
          out_value = trim(ped_rec%SIRE)
          call write_fixed_popqc(out_unit, out_value, field_widths(j))
        case ('DAM', 'DAM_ID')
          out_value = trim(ped_rec%DAM)
          call write_fixed_popqc(out_unit, out_value, field_widths(j))
        case ('SEX')
          write(out_value, '(I0)') ped_rec%SEX
          out_value = trim(adjustl(out_value))
          call write_fixed_popqc(out_unit, out_value, field_widths(j))
        case ('BDATE', 'BIRTH', 'BIRTH_DATE', 'BIRTHDATE')
          write(out_value, '(I0)') ped_rec%BDate
          out_value = trim(adjustl(out_value))
          call write_fixed_popqc(out_unit, out_value, field_widths(j))
        case ('LOC')
          out_value = trim(ped_rec%LOC)
          call write_fixed_popqc(out_unit, out_value, field_widths(j))
        case ('GENO')
          do s = 1, n_auto
            geno_buf(s:s) = achar(iachar('0') + int(in_geno_matrix(i, auto_col(s))))
          end do
          write(out_unit, '(A)', advance='no') geno_buf(1:n_auto)
        case default
          out_value = ''
          width = field_widths(j)
          call write_fixed_popqc(out_unit, out_value, width)
        end select
        if (j < nvar) write(out_unit, '(A1)', advance='no') out_delim
      end do
      write(out_unit, '()')
      n_written = n_written + 1
    end do

    close(out_unit)
    if (allocated(auto_col)) deallocate(auto_col)

    print '(a,1x,a)',  '  Imputation GENO output file:', trim(out_file)
    print '(a, i8)',   '  Valid animals written:      ', n_written
    print '(a, i8)',   '  Excluded (action=2):        ', n_excluded_delete
    print '(a, i8)',   '  Excluded (empty/dropped):   ', n_excluded_empty
    print '(a, i8)',   '  Valid SNPs (QC passed):     ', in_n_snps_pass_qc
    print '(a, i8)',   '  Valid SNPs written (Chr1-18, QC pass):', n_auto
    print '(a, i8)',   '  Total SNPs in GENO:         ', in_n_snps
  end subroutine write_valid_geno_for_imputation_mod

  subroutine write_fixed_popqc(unit_no, value, width)
    integer, intent(in) :: unit_no, width
    character(len=*), intent(in) :: value
    integer :: w, n
    character(len=MAX_STR) :: tmp

    tmp = adjustl(value)
    w = max(0, width)
    n = len_trim(tmp)

    if (w == 0) then
      write(unit_no, '(A)', advance='no') trim(tmp)
    else if (n >= w) then
      write(unit_no, '(A)', advance='no') tmp(1:w)
    else
      write(unit_no, '(A)', advance='no') trim(tmp) // repeat(' ', w - n)
    end if
  end subroutine write_fixed_popqc

end module M_IO_Geno
