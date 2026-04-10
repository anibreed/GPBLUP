module M_phimpute_family_phasing
  use M_kinds
  use M_Variables, only: PEDFile, LEN_STR, FileInfo
  use M_IO_Ped, only: PEDHashTable, pht_load_ped_file, pht_free, pht_search
  use M_ped, only: CPED, str_hash_t
  use M_phimpute_params, only: phimpute_params_t
  use M_phimpute_data, only: phimpute_state_t, phasing_stats_t
  use M_phimpute_bitpack, only: GENO_MISSING
  use M_phimpute_haplotype_lib, only: build_haplotype_panel, summarize_haplotype_panel
  use omp_lib
  implicit none

  public :: run_family_phasing

contains

  subroutine run_family_phasing(params, state, stats, verbose)
    type(phimpute_params_t), intent(in) :: params
    type(phimpute_state_t), intent(inout) :: state
    type(phasing_stats_t), intent(out), optional :: stats
    logical, intent(in), optional :: verbose
    integer :: i, n_animals, n_trios
    integer, allocatable :: child_idx(:), sire_idx(:), dam_idx(:)
    integer(kind=ki8) :: filled_total, checked_total
    real :: t0, t1
    logical :: do_log
    integer, allocatable :: chr_by_snp(:), chr_index_by_snp(:), chr_ids(:)
    integer(kind=ki8), allocatable :: checked_by_chr(:), filled_by_chr(:)
    integer(kind=ki8), allocatable :: checked_by_snp(:), filled_by_snp(:)
    integer :: n_chr

    do_log = .true.
    if (present(verbose)) do_log = verbose

    n_animals = state%cand%n_animals
    if (n_animals <= 0 .or. state%cand%n_snps <= 0) then
      return
    end if

    call build_trio_index(params%pedfile, state, child_idx, sire_idx, dam_idx, n_trios)
    if (n_trios <= 0) then
      state%family_phased = .true.
      return
    end if

    if (do_log) write(*,'(A,I0)') 'FAMILY TRIOS: ', n_trios

    if (.not. state%ref_panel%built) then
      call cpu_time(t0)
      call build_haplotype_panel(state%ref_panel, state%ref, pick_hap_window(params), pick_hap_shift(params))
      call cpu_time(t1)
      if (do_log) write(*,'(A,F8.3,A)') 'TIME haplotype panel build: ', t1 - t0, ' s'
    end if

    if (do_log) call summarize_haplotype_panel(state%ref_panel)

    call omp_set_num_threads(params%omp_threads)

    filled_total = 0
    checked_total = 0

    call build_chr_index_map(state, chr_by_snp, chr_index_by_snp, chr_ids, n_chr)
    allocate(checked_by_chr(n_chr), filled_by_chr(n_chr))
    allocate(checked_by_snp(state%cand%n_snps), filled_by_snp(state%cand%n_snps))
    checked_by_chr = 0
    filled_by_chr = 0
    checked_by_snp = 0
    filled_by_snp = 0

    !$omp parallel do default(shared) private(i) reduction(+:filled_total,checked_total) schedule(static)
    do i = 1, n_trios
      call phase_one_trio(state, child_idx(i), sire_idx(i), dam_idx(i), chr_index_by_snp, &
        checked_by_chr, filled_by_chr, checked_by_snp, filled_by_snp, filled_total, checked_total)
    end do
    !$omp end parallel do

    if (present(stats)) then
      stats%n_trios = n_trios
      stats%n_snps = state%cand%n_snps
      stats%n_animals = state%cand%n_animals
      stats%checked_total = checked_total
      stats%filled_total = filled_total
      stats%n_chr = n_chr
      if (allocated(stats%chr_ids)) deallocate(stats%chr_ids)
      if (allocated(stats%checked_by_chr)) deallocate(stats%checked_by_chr)
      if (allocated(stats%filled_by_chr)) deallocate(stats%filled_by_chr)
      if (allocated(stats%checked_by_snp)) deallocate(stats%checked_by_snp)
      if (allocated(stats%filled_by_snp)) deallocate(stats%filled_by_snp)
      allocate(stats%chr_ids(n_chr), stats%checked_by_chr(n_chr), stats%filled_by_chr(n_chr))
      allocate(stats%checked_by_snp(state%cand%n_snps), stats%filled_by_snp(state%cand%n_snps))
      stats%chr_ids = chr_ids
      stats%checked_by_chr = checked_by_chr
      stats%filled_by_chr = filled_by_chr
      stats%checked_by_snp = checked_by_snp
      stats%filled_by_snp = filled_by_snp
    end if

    if (do_log) then
      call log_phasing_stats(state, chr_ids, checked_by_chr, filled_by_chr, checked_by_snp, filled_by_snp, &
        checked_total, filled_total)
    end if

    if (allocated(child_idx)) deallocate(child_idx)
    if (allocated(sire_idx)) deallocate(sire_idx)
    if (allocated(dam_idx)) deallocate(dam_idx)
    if (allocated(chr_by_snp)) deallocate(chr_by_snp)
    if (allocated(chr_index_by_snp)) deallocate(chr_index_by_snp)
    if (allocated(chr_ids)) deallocate(chr_ids)
    if (allocated(checked_by_chr)) deallocate(checked_by_chr)
    if (allocated(filled_by_chr)) deallocate(filled_by_chr)
    if (allocated(checked_by_snp)) deallocate(checked_by_snp)
    if (allocated(filled_by_snp)) deallocate(filled_by_snp)

    state%family_phased = .true.
  end subroutine run_family_phasing

  integer function pick_hap_window(params) result(win)
    type(phimpute_params_t), intent(in) :: params

    win = params%window_sizes(1)
    if (params%hap_window_count > 0) then
      if (params%hap_window_candidates(1) > 0) win = params%hap_window_candidates(1)
    end if
  end function pick_hap_window

  integer function pick_hap_shift(params) result(sh)
    type(phimpute_params_t), intent(in) :: params

    sh = 0
    if (params%hap_shift_count > 0) then
      if (params%hap_shift_candidates(1) >= 0) sh = params%hap_shift_candidates(1)
    end if
  end function pick_hap_shift

  subroutine phase_one_trio(state, child, sire, dam, chr_index_by_snp, checked_by_chr, filled_by_chr, &
      checked_by_snp, filled_by_snp, filled_total, checked_total)
    type(phimpute_state_t), intent(inout) :: state
    integer, intent(in) :: child, sire, dam
    integer, intent(in) :: chr_index_by_snp(:)
    integer(kind=ki8), intent(inout) :: checked_by_chr(:), filled_by_chr(:)
    integer(kind=ki8), intent(inout) :: checked_by_snp(:), filled_by_snp(:)
    integer(kind=ki8), intent(inout) :: filled_total, checked_total
    integer :: w, offset, snp_idx
    integer(kind=ki4) :: child_word, sire_word, dam_word
    integer(kind=ki1) :: g_child, g_sire, g_dam, g_out
    logical :: ok
    integer :: chr_idx

    do w = 1, state%cand%n_words
      child_word = state%cand%packed(w, child)
      sire_word = state%cand%packed(w, sire)
      dam_word = state%cand%packed(w, dam)

      do offset = 0, 15
        snp_idx = (w - 1) * 16 + offset + 1
        if (snp_idx > state%cand%n_snps) exit

        g_child = int(iand(ishft(child_word, -2 * offset), 3_ki4), kind=ki1)
        if (g_child /= GENO_MISSING) cycle

        checked_total = checked_total + 1
        !$omp atomic
        checked_by_snp(snp_idx) = checked_by_snp(snp_idx) + 1
        chr_idx = chr_index_by_snp(snp_idx)
        if (chr_idx > 0) then
          !$omp atomic
          checked_by_chr(chr_idx) = checked_by_chr(chr_idx) + 1
        end if
        g_sire = int(iand(ishft(sire_word, -2 * offset), 3_ki4), kind=ki1)
        g_dam = int(iand(ishft(dam_word, -2 * offset), 3_ki4), kind=ki1)

        call mendel_phase_rules(g_sire, g_dam, g_out, ok)
        if (ok) then
          call set_word_geno(child_word, offset, g_out)
          filled_total = filled_total + 1
          !$omp atomic
          filled_by_snp(snp_idx) = filled_by_snp(snp_idx) + 1
          if (chr_idx > 0) then
            !$omp atomic
            filled_by_chr(chr_idx) = filled_by_chr(chr_idx) + 1
          end if
        end if
      end do

      state%cand%packed(w, child) = child_word
    end do
  end subroutine phase_one_trio

  subroutine build_chr_index_map(state, chr_by_snp, chr_index_by_snp, chr_ids, n_chr)
    type(phimpute_state_t), intent(in) :: state
    integer, allocatable, intent(out) :: chr_by_snp(:)
    integer, allocatable, intent(out) :: chr_index_by_snp(:)
    integer, allocatable, intent(out) :: chr_ids(:)
    integer, intent(out) :: n_chr

    integer :: n_snps, i, idx, c
    logical :: seen

    n_snps = state%cand%n_snps
    allocate(chr_by_snp(n_snps), chr_index_by_snp(n_snps))
    chr_by_snp = 0
    chr_index_by_snp = 0

    do i = 1, state%map%n_snps
      idx = state%map%array_all(i)
      if (idx > 0 .and. idx <= n_snps) then
        chr_by_snp(idx) = state%map%chr(i)
      end if
    end do

    allocate(chr_ids(n_snps))
    n_chr = 0
    do i = 1, n_snps
      if (chr_by_snp(i) <= 0) cycle
      seen = .false.
      do c = 1, n_chr
        if (chr_ids(c) == chr_by_snp(i)) then
          seen = .true.
          exit
        end if
      end do
      if (.not. seen) then
        n_chr = n_chr + 1
        chr_ids(n_chr) = chr_by_snp(i)
      end if
    end do

    do i = 1, n_snps
      if (chr_by_snp(i) <= 0) cycle
      chr_index_by_snp(i) = find_chr_index(chr_ids, n_chr, chr_by_snp(i))
    end do

    if (n_chr > 0 .and. n_chr < n_snps) then
      call shrink_chr_ids(chr_ids, n_chr)
    end if
  end subroutine build_chr_index_map

  subroutine shrink_chr_ids(chr_ids, n_chr)
    integer, allocatable, intent(inout) :: chr_ids(:)
    integer, intent(in) :: n_chr
    integer, allocatable :: tmp(:)

    allocate(tmp(n_chr))
    tmp = chr_ids(1:n_chr)
    deallocate(chr_ids)
    allocate(chr_ids(n_chr))
    chr_ids = tmp
    deallocate(tmp)
  end subroutine shrink_chr_ids

  integer function find_chr_index(chr_ids, n_chr, chr_val) result(idx)
    integer, intent(in) :: chr_ids(:)
    integer, intent(in) :: n_chr, chr_val
    integer :: i

    idx = 0
    do i = 1, n_chr
      if (chr_ids(i) == chr_val) then
        idx = i
        return
      end if
    end do
  end function find_chr_index

  subroutine log_phasing_stats(state, chr_ids, checked_by_chr, filled_by_chr, checked_by_snp, filled_by_snp, &
      checked_total, filled_total)
    type(phimpute_state_t), intent(in) :: state
    integer, intent(in) :: chr_ids(:)
    integer(kind=ki8), intent(in) :: checked_by_chr(:), filled_by_chr(:)
    integer(kind=ki8), intent(in) :: checked_by_snp(:), filled_by_snp(:)
    integer(kind=ki8), intent(in) :: checked_total, filled_total
    integer :: i, n_show
    integer, allocatable :: order(:)
    integer(kind=ki8), allocatable :: counts(:)

    write(*,'(A,I0)') 'MENDEL checked missing genotypes: ', checked_total
    write(*,'(A,I0)') 'MENDEL filled genotypes: ', filled_total
    if (checked_total > 0) then
      write(*,'(A,F8.3)') 'MENDEL fill rate: ', real(filled_total) / real(checked_total)
    end if

    if (size(chr_ids) > 0) then
      write(*,'(A)') 'MENDEL per-CHR summary:'
      do i = 1, size(chr_ids)
        if (checked_by_chr(i) > 0) then
          write(*,'(A,I0,A,I0,A,I0,A,F6.3)') '  CHR=', chr_ids(i), ' checked=', checked_by_chr(i), &
            ' filled=', filled_by_chr(i), ' rate=', real(filled_by_chr(i)) / real(checked_by_chr(i))
        end if
      end do
    end if

    allocate(order(state%cand%n_snps), counts(state%cand%n_snps))
    do i = 1, state%cand%n_snps
      order(i) = i
      counts(i) = filled_by_snp(i)
    end do

    call sort_counts(counts, order, 1, state%cand%n_snps)
    n_show = min(5, state%cand%n_snps)
    write(*,'(A)') 'MENDEL top SNP fills:'
    do i = 1, n_show
      if (counts(i) <= 0) exit
      write(*,'(A,I0,A,I0)') '  SNP_IDX=', order(i), ' filled=', counts(i)
    end do

    deallocate(order, counts)
  end subroutine log_phasing_stats

  recursive subroutine sort_counts(counts, order, left, right)
    integer(kind=ki8), intent(inout) :: counts(:)
    integer, intent(inout) :: order(:)
    integer, intent(in) :: left, right
    integer :: i, j, tmp
    integer(kind=ki8) :: pivot, tcount

    if (left >= right) return

    pivot = counts((left + right) / 2)
    i = left
    j = right

    do
      do while (counts(i) > pivot)
        i = i + 1
      end do
      do while (counts(j) < pivot)
        j = j - 1
      end do
      if (i <= j) then
        tcount = counts(i)
        counts(i) = counts(j)
        counts(j) = tcount

        tmp = order(i)
        order(i) = order(j)
        order(j) = tmp

        i = i + 1
        j = j - 1
      end if
      if (i > j) exit
    end do

    if (left < j) call sort_counts(counts, order, left, j)
    if (i < right) call sort_counts(counts, order, i, right)
  end subroutine sort_counts

  subroutine set_word_geno(word_val, offset, geno)
    integer(kind=ki4), intent(inout) :: word_val
    integer, intent(in) :: offset
    integer(kind=ki1), intent(in) :: geno
    integer(kind=ki4) :: mask

    mask = ishft(3_ki4, 2 * offset)
    word_val = iand(word_val, not(mask))
    word_val = ior(word_val, ishft(int(geno, ki4), 2 * offset))
  end subroutine set_word_geno

  subroutine mendel_phase_rules(g_sire, g_dam, g_out, ok)
    integer(kind=ki1), intent(in) :: g_sire, g_dam
    integer(kind=ki1), intent(out) :: g_out
    logical, intent(out) :: ok

    ok = .false.
    g_out = GENO_MISSING

    if (g_sire == GENO_MISSING .or. g_dam == GENO_MISSING) return

    if (g_sire == 0 .and. g_dam == 0) then
      g_out = 0
      ok = .true.
    else if (g_sire == 2 .and. g_dam == 2) then
      g_out = 2
      ok = .true.
    else if ((g_sire == 0 .and. g_dam == 2) .or. (g_sire == 2 .and. g_dam == 0)) then
      g_out = 1
      ok = .true.
    end if
  end subroutine mendel_phase_rules

  subroutine build_trio_index(pedinfo, state, child_idx, sire_idx, dam_idx, n_trios)
    type(FileInfo), intent(in) :: pedinfo
    type(phimpute_state_t), intent(in) :: state
    integer, allocatable, intent(out) :: child_idx(:), sire_idx(:), dam_idx(:)
    integer, intent(out) :: n_trios

    type(PEDHashTable) :: pht
    type(CPED) :: ped_rec
    type(str_hash_t) :: id_to_idx
    integer :: n_rec, unit_fd
    integer :: i, idx_out
    character(len=LEN_STR) :: sire_id, dam_id
    logical :: found

    n_trios = 0
    if (.not. allocated(state%cand_ids)) return
    if (state%cand%n_animals <= 0) return

    PEDFile = pedinfo
    call pht_load_ped_file(pht, n_rec, unit_fd, 'ID', quiet=.true.)

    call id_to_idx%reset()
    do i = 1, state%cand%n_animals
      call id_to_idx%put(trim(state%cand_ids(i)), i)
    end do

    allocate(child_idx(state%cand%n_animals))
    allocate(sire_idx(state%cand%n_animals))
    allocate(dam_idx(state%cand%n_animals))
    child_idx = 0
    sire_idx = 0
    dam_idx = 0

    do i = 1, state%cand%n_animals
      found = pht_search(pht, trim(state%cand_ids(i)), ped_rec)
      if (.not. found) cycle
      sire_id = trim(ped_rec%SIRE)
      dam_id = trim(ped_rec%DAM)
      if (len_trim(sire_id) == 0 .or. sire_id == '0') cycle
      if (len_trim(dam_id) == 0 .or. dam_id == '0') cycle

      call id_to_idx%get(trim(sire_id), idx_out, found)
      if (.not. found) cycle
      sire_idx(i) = idx_out

      call id_to_idx%get(trim(dam_id), idx_out, found)
      if (.not. found) cycle
      dam_idx(i) = idx_out

      n_trios = n_trios + 1
      child_idx(n_trios) = i
    end do

    if (n_trios > 0 .and. n_trios < state%cand%n_animals) then
      call shrink_trios(child_idx, sire_idx, dam_idx, n_trios)
    end if

    call pht_free(pht)
    call id_to_idx%reset()
  end subroutine build_trio_index

  subroutine shrink_trios(child_idx, sire_idx, dam_idx, n_trios)
    integer, intent(inout), allocatable :: child_idx(:), sire_idx(:), dam_idx(:)
    integer, intent(in) :: n_trios
    integer, allocatable :: tmp_child(:), tmp_sire(:), tmp_dam(:)

    allocate(tmp_child(n_trios), tmp_sire(n_trios), tmp_dam(n_trios))
    tmp_child = child_idx(1:n_trios)
    tmp_sire = sire_idx(1:n_trios)
    tmp_dam = dam_idx(1:n_trios)

    deallocate(child_idx, sire_idx, dam_idx)
    allocate(child_idx(n_trios), sire_idx(n_trios), dam_idx(n_trios))
    child_idx = tmp_child
    sire_idx = tmp_sire
    dam_idx = tmp_dam

    deallocate(tmp_child, tmp_sire, tmp_dam)
  end subroutine shrink_trios

end module M_phimpute_family_phasing
