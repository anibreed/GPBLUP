module M_phimpute_map
  use M_kinds
  use M_Variables, only: FileInfo, LEN_STR, MAX_VAR
  use M_ReadFile, only: fopen, readline, N_recf
  use M_phimpute_data, only: map_info_t, init_map_info
  use M_StrEdit, only: to_upper
  implicit none

  public :: load_map_sorted, lookup_array_all

contains

  subroutine load_map_sorted(map_info, map, stop_on_mismatch, max_autosome)
    type(FileInfo), intent(in) :: map_info
    type(map_info_t), intent(inout) :: map
    logical, intent(in), optional :: stop_on_mismatch
    integer, intent(in), optional :: max_autosome

    integer :: n_snp, u_fd, n, i, j
    integer(kind=ki4) :: XI(MAX_VAR)
    character(len=LEN_STR) :: XC(MAX_VAR)
    integer, allocatable :: order(:)
    integer(kind=ki8), allocatable :: keys(:)
    integer :: current_chr, chr_rank
    logical :: strict
    integer :: c, k, max_chr, dup_cnt, miss_cnt, bad_cnt
    integer, allocatable :: chr_ids(:), chr_max(:), chr_counts(:), tmp_counts(:)
    integer :: max_chr_filter, valid_cnt
    character(len=LEN_STR), allocatable :: tmp_snp_id(:)
    integer, allocatable :: tmp_chr(:), tmp_pos(:)
    integer, allocatable :: tmp_array_all(:), tmp_array_chr(:)
    character(len=1), allocatable :: tmp_a1(:), tmp_a2(:)
    integer :: idx_snp, idx_chr, idx_pos, idx_all, idx_chrarr, idx_a1, idx_a2
    logical :: has_array_all, has_array_chr, has_alleles
    logical :: ok_all, ok_chr
    logical :: use_input_arrays, need_chrpos_sort

    strict = .false.
    if (present(stop_on_mismatch)) strict = stop_on_mismatch
    max_chr_filter = 0
    if (present(max_autosome)) max_chr_filter = max_autosome

    n_snp = N_recf(trim(map_info%FileName)) - map_info%Header
    if (n_snp <= 0) then
      call init_map_info(map, 0)
      return
    end if

    allocate(tmp_snp_id(n_snp), tmp_chr(n_snp), tmp_pos(n_snp))
    allocate(tmp_array_all(n_snp), tmp_array_chr(n_snp))
    allocate(tmp_a1(n_snp), tmp_a2(n_snp))
    tmp_snp_id = ''
    tmp_chr = 0
    tmp_pos = 0
    tmp_array_all = 0
    tmp_array_chr = 0
    tmp_a1 = ''
    tmp_a2 = ''

    idx_snp = find_field_index(map_info, 'SNP_ID')
    idx_chr = find_field_index(map_info, 'CHR')
    idx_pos = find_field_index(map_info, 'POS')
    idx_all = find_field_index(map_info, 'ARRAY_ALL')
    idx_chrarr = find_field_index(map_info, 'ARRAY_CHR')
    idx_a1 = find_field_index(map_info, 'ALLELE1')
    if (idx_a1 <= 0) idx_a1 = find_field_index(map_info, 'A1')
    idx_a2 = find_field_index(map_info, 'ALLELE2')
    if (idx_a2 <= 0) idx_a2 = find_field_index(map_info, 'A2')
    has_array_all = (idx_all > 0)
    has_array_chr = (idx_chrarr > 0)
    has_alleles = (idx_a1 > 0 .and. idx_a2 > 0)

    u_fd = fopen(trim(map_info%FileName))
    do i = 1, map_info%Header
      call readline(u_fd, n, XC, XI, dlm_str=trim(map_info%Delim_char))
    end do

    do i = 1, n_snp
      call readline(u_fd, n, XC, XI, dlm_str=trim(map_info%Delim_char))
      if (n < 0) exit

      if (idx_snp > 0) tmp_snp_id(i) = trim(XC(map_info%FieldLoc(idx_snp)))
      if (idx_chr > 0) tmp_chr(i) = XI(map_info%FieldLoc(idx_chr))
      if (idx_pos > 0) tmp_pos(i) = XI(map_info%FieldLoc(idx_pos))
      if (has_array_all) tmp_array_all(i) = XI(map_info%FieldLoc(idx_all))
      if (has_array_chr) tmp_array_chr(i) = XI(map_info%FieldLoc(idx_chrarr))
      if (has_alleles) then
        tmp_a1(i) = XC(map_info%FieldLoc(idx_a1))(1:1)
        tmp_a2(i) = XC(map_info%FieldLoc(idx_a2))(1:1)
      end if
    end do
    close(u_fd)

    valid_cnt = 0
    do i = 1, n_snp
      if (tmp_chr(i) <= 0) cycle
      if (max_chr_filter > 0 .and. tmp_chr(i) > max_chr_filter) cycle
      valid_cnt = valid_cnt + 1
    end do

    write(*,'(A,I0,A,I0)') 'MAPFILE records read: ', n_snp, ' valid: ', valid_cnt

    if (valid_cnt <= 0) then
      call init_map_info(map, 0)
      deallocate(tmp_snp_id, tmp_chr, tmp_pos)
      return
    end if

    call init_map_info(map, valid_cnt)
    k = 0
    do i = 1, n_snp
      if (tmp_chr(i) <= 0) cycle
      if (max_chr_filter > 0 .and. tmp_chr(i) > max_chr_filter) cycle
      k = k + 1
      map%snp_id(k) = tmp_snp_id(i)
      map%chr(k) = tmp_chr(i)
      map%pos(k) = tmp_pos(i)
      if (has_alleles) then
        map%allele1(k) = tmp_a1(i)
        map%allele2(k) = tmp_a2(i)
      end if
      if (has_array_all) map%array_all(k) = tmp_array_all(i)
      if (has_array_chr) map%array_chr(k) = tmp_array_chr(i)
    end do

    deallocate(tmp_snp_id, tmp_chr, tmp_pos, tmp_array_all, tmp_array_chr, tmp_a1, tmp_a2)

    use_input_arrays = has_array_all .and. has_array_chr
    ok_all = .false.
    ok_chr = .false.
    if (use_input_arrays) then
      call validate_array_all(map, strict, ok_all)
      call validate_array_chr(map, strict, ok_chr)
    end if

    need_chrpos_sort = .not. (use_input_arrays .and. ok_all .and. ok_chr)
    if (use_input_arrays) then
      if (need_chrpos_sort) then
        write(*,'(A)') 'MAP arrays invalid; rebuilding ARRAY_ALL/ARRAY_CHR by CHR/POS sort.'
      else
        write(*,'(A)') 'MAP arrays valid; using input ARRAY_ALL/ARRAY_CHR without sorting.'
      end if
    else
      write(*,'(A)') 'MAP arrays not provided; building ARRAY_ALL/ARRAY_CHR by CHR/POS sort.'
    end if
    if (need_chrpos_sort) then
      allocate(order(map%n_snps), keys(map%n_snps))
      do i = 1, map%n_snps
        order(i) = i
        keys(i) = chrpos_key(map%chr(i), map%pos(i))
      end do

      call sort_keys(keys, order, 1, map%n_snps)

      do i = 1, map%n_snps
        map%array_all(order(i)) = i
      end do

      current_chr = -1
      chr_rank = 0
      do i = 1, map%n_snps
        if (map%chr(order(i)) /= current_chr) then
          current_chr = map%chr(order(i))
          chr_rank = 0
        end if
        chr_rank = chr_rank + 1
        map%array_chr(order(i)) = chr_rank
      end do
    end if

    do i = 1, min(2, map%n_snps)
      write(*,'(A,I0,A,A)') 'MAP PREVIEW [', i, '] SNP_ID=', trim(map%snp_id(i))
      write(*,'(A,I0,A,I0,A,I0,A,I0)') '  CHR=', map%chr(i), ' POS=', map%pos(i), &
        ' ARRAY_ALL=', map%array_all(i), ' ARRAY_CHR=', map%array_chr(i)
    end do

    if (need_chrpos_sort) then
      if (allocated(map%keys_sorted)) deallocate(map%keys_sorted)
      if (allocated(map%array_all_sorted)) deallocate(map%array_all_sorted)
      allocate(map%keys_sorted(map%n_snps), map%array_all_sorted(map%n_snps))
      do i = 1, map%n_snps
        map%keys_sorted(i) = keys(i)
        map%array_all_sorted(i) = map%array_all(order(i))
      end do

      deallocate(order, keys)

      if (.not. use_input_arrays) then
        call validate_array_chr(map, strict, ok_chr)
      end if
    end if
  end subroutine load_map_sorted

  subroutine validate_array_all(map, strict, ok)
    type(map_info_t), intent(in) :: map
    logical, intent(in) :: strict
    logical, intent(out) :: ok
    integer :: n_snp, i
    integer :: dup_cnt, miss_cnt, bad_cnt
    integer(kind=ki8), allocatable :: keys(:)
    integer, allocatable :: order(:)

    n_snp = map%n_snps
    ok = .true.
    if (n_snp <= 0) return

    allocate(keys(n_snp), order(n_snp))
    do i = 1, n_snp
      keys(i) = int(map%array_all(i), kind=ki8)
      order(i) = i
    end do

    call sort_keys(keys, order, 1, n_snp)

    dup_cnt = 0
    miss_cnt = 0
    bad_cnt = 0

    if (keys(1) /= 1_ki8) miss_cnt = miss_cnt + 1
    do i = 2, n_snp
      if (keys(i) == keys(i-1)) dup_cnt = dup_cnt + 1
      if (keys(i) /= keys(i-1) + 1_ki8) miss_cnt = miss_cnt + 1
    end do
    do i = 1, n_snp
      if (map%array_all(i) <= 0) bad_cnt = bad_cnt + 1
    end do

    if (dup_cnt > 0 .or. miss_cnt > 0 .or. bad_cnt > 0) then
      ok = .false.
      write(*,'(A,I0,A,I0,A,I0)') 'MAP ARRAY_ALL check: dup=', dup_cnt, ' missing=', miss_cnt, ' invalid=', bad_cnt
      if (strict) then
        write(*,'(A)') 'STOP_ON_MISMATCH=1: abort due to MAP ARRAY_ALL inconsistencies.'
        stop 2
      end if
    end if

    deallocate(keys, order)
  end subroutine validate_array_all

  subroutine validate_array_chr(map, strict, ok)
    type(map_info_t), intent(in) :: map
    logical, intent(in) :: strict
    logical, intent(out) :: ok
    integer :: n_snp, i, c, idx, n_chr
    integer, allocatable :: chr_ids(:), chr_max(:), chr_counts(:), tmp_counts(:)
    integer :: dup_cnt, miss_cnt, bad_cnt, k

    n_snp = map%n_snps
    ok = .true.
    if (n_snp <= 0) return

    allocate(chr_ids(n_snp))
    n_chr = 0
    do i = 1, n_snp
      if (map%chr(i) <= 0) cycle
      do c = 1, n_chr
        if (chr_ids(c) == map%chr(i)) exit
      end do
      if (c > n_chr) then
        n_chr = n_chr + 1
        chr_ids(n_chr) = map%chr(i)
      end if
    end do

    allocate(chr_max(n_chr), chr_counts(n_chr))
    chr_max = 0
    chr_counts = 0
    bad_cnt = 0

    do i = 1, n_snp
      if (map%chr(i) <= 0) cycle
      c = find_chr_index(chr_ids, n_chr, map%chr(i))
      if (c <= 0) cycle
      if (map%array_chr(i) <= 0) then
        bad_cnt = bad_cnt + 1
        cycle
      end if
      if (map%array_chr(i) > chr_max(c)) chr_max(c) = map%array_chr(i)
      chr_counts(c) = chr_counts(c) + 1
    end do

    dup_cnt = 0
    miss_cnt = 0
    do c = 1, n_chr
      if (chr_max(c) <= 0) cycle
      allocate(tmp_counts(chr_max(c)))
      tmp_counts = 0
      do i = 1, n_snp
        if (map%chr(i) /= chr_ids(c)) cycle
        if (map%array_chr(i) <= 0) cycle
        if (map%array_chr(i) <= chr_max(c)) then
          tmp_counts(map%array_chr(i)) = tmp_counts(map%array_chr(i)) + 1
        end if
      end do
      do k = 1, chr_max(c)
        if (tmp_counts(k) == 0) miss_cnt = miss_cnt + 1
        if (tmp_counts(k) > 1) dup_cnt = dup_cnt + (tmp_counts(k) - 1)
      end do
      deallocate(tmp_counts)
    end do

    if (dup_cnt > 0 .or. miss_cnt > 0 .or. bad_cnt > 0) then
      ok = .false.
      write(*,'(A,I0,A,I0,A,I0)') 'MAP ARRAY_CHR check: dup=', dup_cnt, ' missing=', miss_cnt, ' invalid=', bad_cnt
      if (strict) then
        write(*,'(A)') 'STOP_ON_MISMATCH=1: abort due to MAP ARRAY_CHR inconsistencies.'
        stop 2
      end if
    end if

    deallocate(chr_ids, chr_max, chr_counts)
  end subroutine validate_array_chr

  integer function find_field_index(info, name) result(idx)
    type(FileInfo), intent(in) :: info
    character(len=*), intent(in) :: name
    character(len=LEN_STR) :: up
    integer :: i

    idx = 0
    up = to_upper(trim(adjustl(name)))
    do i = 1, info%NVAR
      if (to_upper(trim(adjustl(info%FieldName(i)))) == trim(adjustl(up))) then
        idx = i
        return
      end if
    end do
  end function find_field_index

  integer function find_chr_index(chr_ids, n_chr, chr_value) result(idx)
    integer, intent(in) :: chr_ids(:)
    integer, intent(in) :: n_chr
    integer, intent(in) :: chr_value
    integer :: i

    idx = 0
    do i = 1, n_chr
      if (chr_ids(i) == chr_value) then
        idx = i
        return
      end if
    end do
  end function find_chr_index

  function lookup_array_all(map, chr, pos) result(array_all)
    type(map_info_t), intent(in) :: map
    integer, intent(in) :: chr, pos
    integer :: array_all
    integer :: lo, hi, mid
    integer(kind=ki8) :: key

    array_all = 0
    if (.not. allocated(map%keys_sorted)) return

    key = chrpos_key(chr, pos)
    lo = 1
    hi = size(map%keys_sorted)

    do while (lo <= hi)
      mid = (lo + hi) / 2
      if (map%keys_sorted(mid) == key) then
        array_all = map%array_all_sorted(mid)
        return
      else if (map%keys_sorted(mid) < key) then
        lo = mid + 1
      else
        hi = mid - 1
      end if
    end do
  end function lookup_array_all

  integer(kind=ki8) function chrpos_key(chr, pos) result(key)
    integer, intent(in) :: chr, pos

    key = int(chr, kind=ki8) * 1000000000_ki8 + int(pos, kind=ki8)
  end function chrpos_key

  recursive subroutine sort_keys(keys, order, left, right)
    integer(kind=ki8), intent(inout) :: keys(:)
    integer, intent(inout) :: order(:)
    integer, intent(in) :: left, right
    integer :: i, j, tmp
    integer(kind=ki8) :: pivot, tkey

    if (left >= right) return

    pivot = keys((left + right) / 2)
    i = left
    j = right

    do
      do while (keys(i) < pivot)
        i = i + 1
      end do
      do while (keys(j) > pivot)
        j = j - 1
      end do
      if (i <= j) then
        tkey = keys(i)
        keys(i) = keys(j)
        keys(j) = tkey

        tmp = order(i)
        order(i) = order(j)
        order(j) = tmp

        i = i + 1
        j = j - 1
      end if
      if (i > j) exit
    end do

    if (left < j) call sort_keys(keys, order, left, j)
    if (i < right) call sort_keys(keys, order, i, right)
  end subroutine sort_keys

end module M_phimpute_map
