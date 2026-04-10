module M_phimpute_io_geno
  use M_kinds
  use M_Variables, only: FileInfo, LEN_STR, MAX_VAR
  use M_ReadFile, only: fopen, N_recf
  use M_StrEdit, only: to_upper, FindDelim
  use M_phimpute_data, only: geno_matrix_t, init_geno_matrix
  use M_phimpute_bitpack, only: pack_geno_string
  implicit none

  public :: load_candidate_geno, load_reference_geno

contains

  subroutine load_candidate_geno(info, n_snps, geno_mat, animal_ids, stop_on_mismatch)
    type(FileInfo), intent(in) :: info
    integer, intent(in) :: n_snps
    type(geno_matrix_t), intent(inout) :: geno_mat
    character(len=LEN_STR), allocatable, intent(out) :: animal_ids(:)
    logical, intent(in), optional :: stop_on_mismatch

    integer :: n_animals, u_fd, i, ios
    integer :: idx_animal, idx_geno
    character(len=:), allocatable :: raw_line
    character(len=1) :: delim
    integer :: line_len
    integer :: a_start, a_end, g_start, g_end
    logical :: found_field
    integer :: first_pos
    integer :: short_cnt, long_cnt
    integer :: geno_len
    logical :: strict

    n_animals = N_recf(trim(info%FileName)) - info%Header
    strict = .false.
    if (present(stop_on_mismatch)) strict = stop_on_mismatch
    if (n_animals <= 0) then
      call init_geno_matrix(geno_mat, 0, n_snps)
      allocate(animal_ids(0))
      return
    end if

    call init_geno_matrix(geno_mat, n_animals, n_snps)
    allocate(animal_ids(n_animals))
    animal_ids = ''

    idx_animal = find_field_index(info, 'ANIMAL_ID')
    idx_geno = find_field_index(info, 'GENO')
    if (idx_animal <= 0 .or. idx_geno <= 0) then
      write(*,*) 'CANDIDATE_GENOFILE missing ANIMAL_ID or GENO field'
      stop 1
    end if

    u_fd = fopen(trim(info%FileName), recl=80000)
    allocate(character(len=80000) :: raw_line)
    do i = 1, info%Header
      read(u_fd, '(A)', iostat=ios) raw_line
      if (ios /= 0) exit
    end do

    delim = FindDelim(trim(info%Delim_char))

    i = 0
    short_cnt = 0
    long_cnt = 0
    do
      read(u_fd, '(A)', iostat=ios) raw_line
      if (ios /= 0) exit
      line_len = len_trim(raw_line)
      if (line_len == 0) cycle
      first_pos = 1
      do while (first_pos <= line_len .and. raw_line(first_pos:first_pos) == ' ')
        first_pos = first_pos + 1
      end do
      if (first_pos > line_len) cycle
      if (raw_line(first_pos:first_pos) == '#') cycle

      call get_field_pos_line(raw_line, line_len, delim, info%FieldLoc(idx_animal), a_start, a_end, found_field)
      if (.not. found_field) cycle
      call get_field_pos_line(raw_line, line_len, delim, info%FieldLoc(idx_geno), g_start, g_end, found_field)
      if (.not. found_field) cycle

      i = i + 1
      if (i > n_animals) exit

      animal_ids(i) = trim(raw_line(a_start:a_end))
      geno_len = g_end - g_start + 1
      if (geno_len < n_snps) short_cnt = short_cnt + 1
      if (geno_len > n_snps) long_cnt = long_cnt + 1
      if (i <= 2) then
        write(*,'(A,I0,A,A,A,I0,A,I0,A,I0)') 'GENO PREVIEW [', i, '] ID=', trim(animal_ids(i)), &
          ' LEN=', geno_len, ' START=', g_start, ' END=', g_end
      end if
      call pack_geno_string(raw_line(g_start:g_end), n_snps, geno_mat%packed(:, i))
    end do

    close(u_fd)
    if (allocated(raw_line)) deallocate(raw_line)

    write(*,'(A,I0,A,I0)') 'CANDIDATE_GENOFILE records read: ', i, ' of ', n_animals

    if (short_cnt > 0 .or. long_cnt > 0) then
      write(*,'(A,I0,A,I0)') 'GENO length check: short=', short_cnt, ' long=', long_cnt
      if (strict) then
        write(*,'(A)') 'STOP_ON_MISMATCH=1: abort due to GENO length mismatch.'
        stop 3
      end if
    else
      write(*,'(A,I0)') 'GENO length check: all records match n_snps=', n_snps
    end if
  end subroutine load_candidate_geno

  subroutine load_reference_geno(info, n_snps, geno_mat, animal_ids, stop_on_mismatch)
    type(FileInfo), intent(in) :: info
    integer, intent(in) :: n_snps
    type(geno_matrix_t), intent(inout) :: geno_mat
    character(len=LEN_STR), allocatable, intent(out) :: animal_ids(:)
    logical, intent(in), optional :: stop_on_mismatch

    call load_candidate_geno(info, n_snps, geno_mat, animal_ids, stop_on_mismatch)
  end subroutine load_reference_geno

  integer function find_field_index(info, name) result(idx)
    type(FileInfo), intent(in) :: info
    character(len=*), intent(in) :: name
    character(len=LEN_STR) :: up
    integer :: i

    idx = 0
    up = to_upper(trim(adjustl(name)))
    do i = 1, info%NVAR
      if (to_upper(trim(adjustl(info%FieldName(i)))) == up) then
        idx = i
        return
      end if
    end do
  end function find_field_index

  subroutine get_field_pos_line(line, line_len, delim, target_col, start_pos, end_pos, found)
    character(len=*), intent(in) :: line
    integer, intent(in) :: line_len
    character(len=1), intent(in) :: delim
    integer, intent(in) :: target_col
    integer, intent(out) :: start_pos, end_pos
    logical, intent(out) :: found

    integer :: i, col

    start_pos = 0
    end_pos = 0
    found = .false.
    if (target_col <= 0) return

    i = 1
    col = 0
    do while (i <= line_len)
      if (delim == ' ') then
        do while (i <= line_len .and. line(i:i) == ' ')
          i = i + 1
        end do
        if (i > line_len) exit
      end if

      col = col + 1
      start_pos = i

      do while (i <= line_len .and. line(i:i) /= delim)
        i = i + 1
      end do
      end_pos = i - 1

      if (col == target_col) then
        found = .true.
        return
      end if

      i = i + 1
    end do
  end subroutine get_field_pos_line

end module M_phimpute_io_geno
