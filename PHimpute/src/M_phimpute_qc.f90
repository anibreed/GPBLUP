module M_phimpute_qc
  use M_kinds
  use M_Variables, only: LEN_STR
  use M_phimpute_data, only: phimpute_state_t, geno_matrix_t, map_info_t, init_geno_matrix, init_map_info, &
    copy_geno_matrix
  use M_phimpute_bitpack, only: get_geno_value, set_geno_value, GENO_MISSING
  use M_phimpute_haplotype_lib, only: count_missing
  use M_phimpute_params, only: phimpute_params_t
  implicit none

  public :: apply_qc_filters

contains

  subroutine apply_qc_filters(params, state)
    type(phimpute_params_t), intent(in) :: params
    type(phimpute_state_t), intent(inout) :: state
    logical, allocatable :: keep_snps(:)
    integer :: kept_snps

    if (.not. params%qc_enable) then
      write(*,'(A)') 'QC disabled: skip sample/SNP filtering.'
      return
    end if

    if (state%cand%n_snps <= 0) return

    call filter_samples(state%cand, state%cand_ids, params%qc_sample_miss_max, 'cand', params%qc_log_prefix)
    call filter_samples(state%ref, state%ref_ids, params%qc_sample_miss_max, 'ref', params%qc_log_prefix)
    if (state%has_reference .and. state%ref%n_animals <= 0) then
      write(*,'(A)') 'REFERENCE panel empty after QC; fallback to candidate GENO.'
      call copy_geno_matrix(state%cand, state%ref)
      if (allocated(state%ref_ids)) deallocate(state%ref_ids)
      if (allocated(state%cand_ids)) then
        allocate(state%ref_ids(size(state%cand_ids)))
        state%ref_ids = state%cand_ids
      end if
      state%has_reference = .false.
    end if

    call filter_snps(state, params%qc_snp_miss_max, params%qc_min_maf, keep_snps, kept_snps)
    call write_removed_snps(state%map%snp_id, keep_snps, params%qc_log_prefix, 'removed_snps.txt')
    if (params%qc_remove_ambiguous) then
      call apply_ambiguous_filter(state, keep_snps, kept_snps, params%qc_log_prefix)
    end if
    if (kept_snps == 0) then
      call clear_all_snps(state)
    else if (kept_snps < state%map%n_snps) then
      call apply_snp_filter(state, keep_snps, kept_snps)
    end if

    if (allocated(keep_snps)) deallocate(keep_snps)
    call write_reference_map(state%map, params%qc_log_prefix)
  end subroutine apply_qc_filters

  subroutine clear_all_snps(state)
    type(phimpute_state_t), intent(inout) :: state

    call init_geno_matrix(state%cand, state%cand%n_animals, 0)
    call init_geno_matrix(state%ref, state%ref%n_animals, 0)
    call init_map_info(state%map, 0)
    state%ref_panel%built = .false.
    write(*,'(A)') 'QC SNP filter: removed all SNPs.'
  end subroutine clear_all_snps

  subroutine filter_samples(geno, ids, max_missing, label, log_prefix)
    type(geno_matrix_t), intent(inout) :: geno
    character(len=*), allocatable, intent(inout) :: ids(:)
    real, intent(in) :: max_missing
    character(len=*), intent(in) :: label
    character(len=*), intent(in) :: log_prefix
    integer, allocatable :: keep_idx(:)
    integer :: a, kept, new_idx
    integer(kind=ki8) :: miss_count
    real :: miss_rate
    type(geno_matrix_t) :: tmp
    character(len=:), allocatable :: tmp_ids(:)
    character(len=:), allocatable :: orig_ids(:)
    integer :: id_len

    if (geno%n_animals <= 0 .or. geno%n_snps <= 0) return

    if (max_missing < 0.0 .or. max_missing >= 1.0) return

    allocate(keep_idx(geno%n_animals))
    keep_idx = 0
    id_len = 1
    if (allocated(ids)) id_len = len(ids)
    allocate(character(len=id_len) :: orig_ids(geno%n_animals))
    orig_ids = ''
    if (allocated(ids)) orig_ids = ids
    kept = 0
    do a = 1, geno%n_animals
      miss_count = count_missing(geno%packed(:, a), 1, geno%n_snps)
      miss_rate = real(miss_count) / real(geno%n_snps)
      if (miss_rate <= max_missing) then
        kept = kept + 1
        keep_idx(kept) = a
      end if
    end do

    if (kept == geno%n_animals) then
      write(*,'(A,I0,A,A,A)') 'QC sample filter: kept all ', kept, ' (', trim(label), ')'
      return
    end if

    if (kept <= 0) then
      call init_geno_matrix(geno, 0, geno%n_snps)
      if (allocated(ids)) deallocate(ids)
      allocate(ids(0))
      write(*,'(A,A,A)') 'QC sample filter: removed all samples (', trim(label), ')'
      return
    end if

    call init_geno_matrix(tmp, kept, geno%n_snps)
    allocate(character(len=id_len) :: tmp_ids(kept))
    do new_idx = 1, kept
      a = keep_idx(new_idx)
      tmp%packed(:, new_idx) = geno%packed(:, a)
      if (allocated(ids)) tmp_ids(new_idx) = ids(a)
    end do

    geno = tmp
    if (allocated(ids)) then
      deallocate(ids)
      allocate(character(len=id_len) :: ids(kept))
      ids = tmp_ids
    end if
    call write_removed_samples(orig_ids, keep_idx, size(keep_idx), trim(label), trim(log_prefix))
    write(*,'(A,I0,A,I0,A)') 'QC sample filter: kept ', kept, ' of ', size(keep_idx), &
      ' (' // trim(label) // ')'
    if (allocated(orig_ids)) deallocate(orig_ids)
    if (allocated(keep_idx)) deallocate(keep_idx)
  end subroutine filter_samples

  subroutine write_removed_samples(orig_ids, keep_idx, total, label, log_prefix)
    character(len=*), intent(in) :: orig_ids(:)
    integer, intent(in) :: keep_idx(:)
    integer, intent(in) :: total
    character(len=*), intent(in) :: label, log_prefix
    logical, allocatable :: is_kept(:)
    integer :: i, u_fd
    character(len=LEN_STR) :: file_name

    if (total <= 0) return

    allocate(is_kept(total))
    is_kept = .false.
    do i = 1, size(keep_idx)
      if (keep_idx(i) > 0 .and. keep_idx(i) <= total) is_kept(keep_idx(i)) = .true.
    end do

    file_name = trim(log_prefix) // '_removed_samples_' // trim(label) // '.txt'
    u_fd = 98
    open(unit=u_fd, file=trim(file_name), status='replace', action='write')
    do i = 1, total
      if (.not. is_kept(i)) write(u_fd,'(A)') trim(orig_ids(i))
    end do
    close(u_fd)
    deallocate(is_kept)
  end subroutine write_removed_samples

  subroutine filter_snps(state, max_missing, min_maf, keep_snps, kept_snps)
    type(phimpute_state_t), intent(in) :: state
    real, intent(in) :: max_missing, min_maf
    logical, allocatable, intent(out) :: keep_snps(:)
    integer, intent(out) :: kept_snps
    integer :: s, n_snps
    integer(kind=ki8), allocatable :: miss_counts(:), alt_counts(:), nonmiss_counts(:)
    integer(kind=ki8) :: total_geno, total_alleles
    real :: miss_rate, maf, p

    n_snps = state%map%n_snps
    kept_snps = 0
    if (n_snps <= 0) then
      allocate(keep_snps(0))
      return
    end if

    allocate(miss_counts(n_snps), alt_counts(n_snps), nonmiss_counts(n_snps))
    miss_counts = 0
    alt_counts = 0
    nonmiss_counts = 0

    call update_snp_stats(state%cand, miss_counts, alt_counts, nonmiss_counts)
    if (state%has_reference .and. state%ref%n_animals > 0) then
      call update_snp_stats(state%ref, miss_counts, alt_counts, nonmiss_counts)
    end if

    allocate(keep_snps(n_snps))
    keep_snps = .false.
    do s = 1, n_snps
      total_geno = miss_counts(s) + nonmiss_counts(s)
      if (total_geno <= 0) cycle
      miss_rate = real(miss_counts(s)) / real(total_geno)
      total_alleles = 2_ki8 * nonmiss_counts(s)
      if (total_alleles <= 0) cycle
      p = real(alt_counts(s)) / real(total_alleles)
      maf = p
      if (maf > 0.5) maf = 1.0 - maf

      if (miss_rate <= max_missing .and. maf >= min_maf) then
        keep_snps(s) = .true.
        kept_snps = kept_snps + 1
      end if
    end do

    write(*,'(A,I0,A,I0,A,F6.3,A,F6.3,A)') 'QC SNP filter: kept ', kept_snps, ' of ', n_snps, &
      ' (miss<=', max_missing, ', maf>=', min_maf, ')'

    deallocate(miss_counts, alt_counts, nonmiss_counts)
  end subroutine filter_snps

  subroutine write_removed_snps(snp_ids, keep_snps, log_prefix, file_tail)
    character(len=*), intent(in) :: snp_ids(:)
    logical, intent(in) :: keep_snps(:)
    character(len=*), intent(in) :: log_prefix, file_tail
    integer :: s, u_fd
    character(len=LEN_STR) :: file_name

    if (size(keep_snps) <= 0) return
    file_name = trim(log_prefix) // '_' // trim(file_tail)
    u_fd = 97
    open(unit=u_fd, file=trim(file_name), status='replace', action='write')
    do s = 1, min(size(keep_snps), size(snp_ids))
      if (.not. keep_snps(s)) write(u_fd,'(A)') trim(snp_ids(s))
    end do
    close(u_fd)
  end subroutine write_removed_snps

  subroutine apply_ambiguous_filter(state, keep_snps, kept_snps, log_prefix)
    type(phimpute_state_t), intent(in) :: state
    logical, intent(inout) :: keep_snps(:)
    integer, intent(inout) :: kept_snps
    character(len=*), intent(in) :: log_prefix
    integer :: s, removed
    character(len=1) :: a1, a2
    logical :: any_alleles

    if (size(keep_snps) <= 0) return

    any_alleles = .false.
    do s = 1, min(size(keep_snps), state%map%n_snps)
      if (state%map%allele1(s) /= '' .or. state%map%allele2(s) /= '') then
        any_alleles = .true.
        exit
      end if
    end do
    if (.not. any_alleles) then
      write(*,'(A)') 'QC ambiguous SNP filter: skipped (allele info missing).'
      return
    end if

    removed = 0
    do s = 1, min(size(keep_snps), state%map%n_snps)
      if (.not. keep_snps(s)) cycle
      a1 = state%map%allele1(s)
      a2 = state%map%allele2(s)
      if (is_ambiguous(a1, a2)) then
        keep_snps(s) = .false.
        removed = removed + 1
        kept_snps = kept_snps - 1
      end if
    end do

    if (removed > 0) then
      write(*,'(A,I0)') 'QC ambiguous SNP filter: removed ', removed
      call write_removed_snps(state%map%snp_id, keep_snps, log_prefix, 'removed_snps_ambiguous.txt')
    else
      write(*,'(A)') 'QC ambiguous SNP filter: none removed.'
    end if
  end subroutine apply_ambiguous_filter

  logical function is_ambiguous(a1, a2) result(flag)
    character(len=1), intent(in) :: a1, a2
    character(len=1) :: x1, x2

    x1 = to_upper_char(a1)
    x2 = to_upper_char(a2)
    flag = .false.
    if ((x1 == 'A' .and. x2 == 'T') .or. (x1 == 'T' .and. x2 == 'A')) flag = .true.
    if ((x1 == 'C' .and. x2 == 'G') .or. (x1 == 'G' .and. x2 == 'C')) flag = .true.
  end function is_ambiguous

  character(len=1) function to_upper_char(ch) result(out)
    character(len=1), intent(in) :: ch

    out = ch
    select case (ch)
    case ('a')
      out = 'A'
    case ('c')
      out = 'C'
    case ('g')
      out = 'G'
    case ('t')
      out = 'T'
    case default
      continue
    end select
  end function to_upper_char

  subroutine update_snp_stats(geno, miss_counts, alt_counts, nonmiss_counts)
    type(geno_matrix_t), intent(in) :: geno
    integer(kind=ki8), intent(inout) :: miss_counts(:)
    integer(kind=ki8), intent(inout) :: alt_counts(:)
    integer(kind=ki8), intent(inout) :: nonmiss_counts(:)
    integer :: a, s
    integer(kind=ki1) :: g

    if (geno%n_animals <= 0 .or. geno%n_snps <= 0) return

    do a = 1, geno%n_animals
      do s = 1, geno%n_snps
        g = get_geno_value(geno%packed(:, a), s)
        if (g == GENO_MISSING) then
          miss_counts(s) = miss_counts(s) + 1
        else
          nonmiss_counts(s) = nonmiss_counts(s) + 1
          alt_counts(s) = alt_counts(s) + g
        end if
      end do
    end do
  end subroutine update_snp_stats

  subroutine apply_snp_filter(state, keep_snps, kept_snps)
    type(phimpute_state_t), intent(inout) :: state
    logical, intent(in) :: keep_snps(:)
    integer, intent(in) :: kept_snps
    integer :: s, new_s, a
    type(geno_matrix_t) :: new_cand, new_ref

    call init_geno_matrix(new_cand, state%cand%n_animals, kept_snps)
    call init_geno_matrix(new_ref, state%ref%n_animals, kept_snps)

    new_s = 0
    do s = 1, state%map%n_snps
      if (.not. keep_snps(s)) cycle
      new_s = new_s + 1

      do a = 1, state%cand%n_animals
        call set_geno_value(new_cand%packed(:, a), new_s, get_geno_value(state%cand%packed(:, a), s))
      end do
      do a = 1, state%ref%n_animals
        call set_geno_value(new_ref%packed(:, a), new_s, get_geno_value(state%ref%packed(:, a), s))
      end do
    end do

    state%cand = new_cand
    state%ref = new_ref

    call filter_map(state, keep_snps, kept_snps)
    state%ref_panel%built = .false.
  end subroutine apply_snp_filter

  subroutine filter_map(state, keep_snps, kept_snps)
    type(phimpute_state_t), intent(inout) :: state
    logical, intent(in) :: keep_snps(:)
    integer, intent(in) :: kept_snps
    integer :: s, new_s
    type(map_info_t) :: old_map

    old_map = state%map
    call init_map_info(state%map, kept_snps)
    new_s = 0
    do s = 1, size(keep_snps)
      if (.not. keep_snps(s)) cycle
      new_s = new_s + 1
      state%map%snp_id(new_s) = old_map%snp_id(s)
      state%map%chr(new_s) = old_map%chr(s)
      state%map%pos(new_s) = old_map%pos(s)
      state%map%allele1(new_s) = old_map%allele1(s)
      state%map%allele2(new_s) = old_map%allele2(s)
    end do

    call rebuild_map_arrays(state%map)

    if (allocated(state%map%keys_sorted)) deallocate(state%map%keys_sorted)
    if (allocated(state%map%array_all_sorted)) deallocate(state%map%array_all_sorted)
  end subroutine filter_map

  subroutine rebuild_map_arrays(map)
    type(map_info_t), intent(inout) :: map
    integer :: n, i, current_chr, chr_rank
    integer(kind=ki8), allocatable :: keys(:)
    integer, allocatable :: order(:)

    n = map%n_snps
    if (n <= 0) return

    allocate(keys(n), order(n))
    do i = 1, n
      keys(i) = chrpos_key(map%chr(i), map%pos(i))
      order(i) = i
    end do

    call sort_keys(keys, order, 1, n)

    do i = 1, n
      map%array_all(order(i)) = i
    end do

    current_chr = -1
    chr_rank = 0
    do i = 1, n
      if (map%chr(order(i)) /= current_chr) then
        current_chr = map%chr(order(i))
        chr_rank = 0
      end if
      chr_rank = chr_rank + 1
      map%array_chr(order(i)) = chr_rank
    end do

    deallocate(keys, order)
  end subroutine rebuild_map_arrays

  subroutine write_reference_map(map, log_prefix)
    type(map_info_t), intent(in) :: map
    character(len=*), intent(in) :: log_prefix
    integer :: u_fd, s
    character(len=LEN_STR) :: file_name

    if (map%n_snps <= 0) return

    file_name = trim(log_prefix) // '_reference_map.txt'
    u_fd = 96
    open(unit=u_fd, file=trim(file_name), status='replace', action='write')
    write(u_fd,'(A)') 'SNP_ID	CHR	POS	ARRAY_ALL	ARRAY_CHR	ALLELE1	ALLELE2'
    do s = 1, map%n_snps
      write(u_fd,'(A,1X,I0,1X,I0,1X,I0,1X,I0,1X,A,1X,A)') trim(map%snp_id(s)), map%chr(s), map%pos(s), &
        map%array_all(s), map%array_chr(s), map%allele1(s), map%allele2(s)
    end do
    close(u_fd)
  end subroutine write_reference_map

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

end module M_phimpute_qc
