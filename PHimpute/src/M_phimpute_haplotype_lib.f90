module M_phimpute_haplotype_lib
  use M_kinds
  use M_phimpute_data, only: geno_matrix_t, haplotype_panel_t, hap_window_t
  implicit none

  public :: build_haplotype_panel
  public :: summarize_haplotype_panel
  public :: hash_window
  public :: find_key_index
  public :: count_missing

  integer(kind=ki8), parameter :: FNV_OFFSET = 1469598103934665603_ki8
  integer(kind=ki8), parameter :: FNV_PRIME  = 1099511628211_ki8

contains

  subroutine build_haplotype_panel(panel, geno_mat, window_size, shift_in)
    type(haplotype_panel_t), intent(inout) :: panel
    type(geno_matrix_t), intent(in) :: geno_mat
    integer, intent(in) :: window_size
    integer, intent(in) :: shift_in

    integer :: n_snps, n_animals
    integer :: n_windows, n_shift, w, idx
    integer :: start_snp, end_snp
    integer :: shift
    integer(kind=ki8), allocatable :: keys(:)
    integer, allocatable :: order(:)
    integer, allocatable :: miss_counts(:)

    n_snps = geno_mat%n_snps
    n_animals = geno_mat%n_animals
    if (n_snps <= 0 .or. n_animals <= 0) return
    if (window_size <= 0) return

    n_windows = (n_snps + window_size - 1) / window_size
    shift = max(0, min(shift_in, window_size / 2))
    n_shift = 0
    if (shift > 0) then
      do w = 1, n_windows
        end_snp = min(w * window_size, n_snps)
        if (end_snp + shift <= n_snps) n_shift = n_shift + 1
      end do
    end if

    panel%window_size = window_size
    panel%shift_size = shift
    panel%n_windows = n_windows + n_shift
    panel%built = .false.
    if (allocated(panel%windows)) deallocate(panel%windows)
    allocate(panel%windows(panel%n_windows))

    idx = 0
    do w = 1, n_windows
      start_snp = (w - 1) * window_size + 1
      end_snp = min(w * window_size, n_snps)

      idx = idx + 1
      call build_one_window(panel%windows(idx), geno_mat, start_snp, end_snp, n_animals)
    end do

    if (shift > 0 .and. n_shift > 0) then
      do w = 1, n_windows
        start_snp = (w - 1) * window_size + 1
        end_snp = min(w * window_size, n_snps)
        if (end_snp + shift > n_snps) cycle

        idx = idx + 1
        call build_one_window(panel%windows(idx), geno_mat, start_snp + shift, end_snp + shift, n_animals)
      end do
    end if

    panel%built = .true.
  end subroutine build_haplotype_panel

  subroutine build_one_window(window, geno_mat, start_snp, end_snp, n_animals)
    type(hap_window_t), intent(inout) :: window
    type(geno_matrix_t), intent(in) :: geno_mat
    integer, intent(in) :: start_snp, end_snp
    integer, intent(in) :: n_animals
    integer(kind=ki8), allocatable :: keys(:)
    integer, allocatable :: order(:)
    integer, allocatable :: miss_counts(:)

    window%start_snp = start_snp
    window%end_snp = end_snp

    allocate(keys(n_animals), miss_counts(n_animals))
    call build_window_hashes(geno_mat, start_snp, end_snp, keys, miss_counts)

    allocate(order(n_animals))
    call init_order(order)
    call sort_keys(keys, order, 1, n_animals)
    call compress_keys(keys, order, miss_counts, window)

    deallocate(order, keys, miss_counts)
  end subroutine build_one_window

  subroutine summarize_haplotype_panel(panel)
    type(haplotype_panel_t), intent(in) :: panel
    integer :: w, n_show
    integer, allocatable :: order(:)
    integer, allocatable :: counts(:)

    if (.not. panel%built) return
    if (panel%n_windows <= 0) return

    allocate(order(panel%n_windows), counts(panel%n_windows))
    do w = 1, panel%n_windows
      order(w) = w
      counts(w) = panel%windows(w)%n_keys
    end do

    call sort_counts(counts, order, 1, panel%n_windows)

    n_show = min(5, panel%n_windows)
    write(*,'(A,I0,A,I0,A,I0)') 'HAP PANEL windows=', panel%n_windows, ' window_size=', panel%window_size, &
      ' shift=', panel%shift_size
    do w = 1, n_show
      write(*,'(A,I0,A,I0,A,I0)') '  window=', order(w), ' keys=', panel%windows(order(w))%n_keys, &
        ' snps=', panel%windows(order(w))%end_snp - panel%windows(order(w))%start_snp + 1
    end do

    deallocate(order, counts)
  end subroutine summarize_haplotype_panel

  subroutine build_window_hashes(geno_mat, start_snp, end_snp, keys, miss_counts)
    type(geno_matrix_t), intent(in) :: geno_mat
    integer, intent(in) :: start_snp, end_snp
    integer(kind=ki8), intent(out) :: keys(:)
    integer, intent(out) :: miss_counts(:)

    integer :: a
    integer(kind=ki8) :: hash

    do a = 1, geno_mat%n_animals
      call hash_window(geno_mat%packed(:, a), start_snp, end_snp, hash)
      keys(a) = hash
      miss_counts(a) = count_missing(geno_mat%packed(:, a), start_snp, end_snp)
    end do
  end subroutine build_window_hashes

  subroutine hash_window(packed_row, start_snp, end_snp, hash)
    integer(kind=ki4), intent(in) :: packed_row(:)
    integer, intent(in) :: start_snp, end_snp
    integer(kind=ki8), intent(out) :: hash

    integer :: start_word, end_word, w
    integer :: start_bit, end_bit
    integer :: offset, snp_idx
    integer(kind=ki4) :: word_val, clean_val
    integer(kind=ki1) :: g

    hash = FNV_OFFSET

    start_word = (start_snp - 1) / 16 + 1
    end_word = (end_snp - 1) / 16 + 1
    start_bit = 2 * mod(start_snp - 1, 16)
    end_bit = 2 * mod(end_snp - 1, 16) + 1

    do w = start_word, end_word
      word_val = packed_row(w)
      clean_val = 0_ki4
      do offset = 0, 15
        snp_idx = (w - 1) * 16 + offset + 1
        if (snp_idx < start_snp .or. snp_idx > end_snp) cycle
        g = int(iand(ishft(word_val, -2 * offset), 3_ki4), kind=ki1)
        if (g == 3_ki1) g = 0_ki1
        clean_val = ior(clean_val, ishft(int(g, ki4), 2 * offset))
      end do

      if (w == start_word .and. w == end_word) then
        clean_val = iand(clean_val, mask_bits(start_bit, end_bit))
      else if (w == start_word) then
        clean_val = iand(clean_val, mask_bits(start_bit, 31))
      else if (w == end_word) then
        clean_val = iand(clean_val, mask_bits(0, end_bit))
      end if

      hash = ieor(hash, int(clean_val, ki8))
      hash = hash * FNV_PRIME
    end do
  end subroutine hash_window

  integer(kind=ki4) function mask_bits(start_bit, end_bit) result(mask)
    integer, intent(in) :: start_bit, end_bit
    integer(kind=ki4) :: left_mask, right_mask

    if (start_bit < 0 .or. end_bit > 31 .or. start_bit > end_bit) then
      mask = -1_ki4
      return
    end if

    left_mask = ishft(-1_ki4, start_bit)
    right_mask = not(ishft(-1_ki4, end_bit + 1))
    mask = iand(left_mask, right_mask)
  end function mask_bits

  subroutine init_order(order)
    integer, intent(out) :: order(:)
    integer :: i

    do i = 1, size(order)
      order(i) = i
    end do
  end subroutine init_order

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

  recursive subroutine sort_counts(counts, order, left, right)
    integer, intent(inout) :: counts(:)
    integer, intent(inout) :: order(:)
    integer, intent(in) :: left, right
    integer :: i, j, tmp
    integer :: pivot, tcount

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

  subroutine compress_keys(keys, order, miss_counts, window)
    integer(kind=ki8), intent(in) :: keys(:)
    integer, intent(in) :: order(:)
    integer, intent(in) :: miss_counts(:)
    type(hap_window_t), intent(inout) :: window

    integer :: n, i, k
    integer(kind=ki8) :: prev
    integer :: rep_miss
    integer, parameter :: K_TOP = 5
    integer :: kslot

    n = size(keys)
    if (n <= 0) return

    allocate(window%keys(n), window%counts(n), window%rep_idx(n))
    allocate(window%topk_idx(n, K_TOP), window%topk_miss(n, K_TOP), window%topk_count(n))
    window%topk_idx = 0
    window%topk_miss = 0
    window%topk_count = 0

    k = 0
    prev = keys(1)
    k = k + 1
    window%keys(k) = prev
    window%counts(k) = 1
    window%rep_idx(k) = order(1)
    rep_miss = miss_counts(order(1))
    window%topk_count(k) = 1
    window%topk_idx(k, 1) = order(1)
    window%topk_miss(k, 1) = miss_counts(order(1))

    do i = 2, n
      if (keys(i) == prev) then
        window%counts(k) = window%counts(k) + 1
        if (miss_counts(order(i)) < rep_miss) then
          window%rep_idx(k) = order(i)
          rep_miss = miss_counts(order(i))
        end if
        call update_topk(window%topk_idx(k, :), window%topk_miss(k, :), window%topk_count(k), order(i), miss_counts, K_TOP)
      else
        prev = keys(i)
        k = k + 1
        window%keys(k) = prev
        window%counts(k) = 1
        window%rep_idx(k) = order(i)
        rep_miss = miss_counts(order(i))
        window%topk_count(k) = 1
        window%topk_idx(k, 1) = order(i)
        window%topk_miss(k, 1) = miss_counts(order(i))
      end if
    end do

    window%n_keys = k
    if (k < n) then
      call shrink_arrays(window, k, K_TOP)
    end if
  end subroutine compress_keys

  subroutine shrink_arrays(window, n_keep, k_top)
    type(hap_window_t), intent(inout) :: window
    integer, intent(in) :: n_keep
    integer, intent(in) :: k_top
    integer(kind=ki8), allocatable :: tmp_keys(:)
    integer, allocatable :: tmp_counts(:)
    integer, allocatable :: tmp_rep(:)
    integer, allocatable :: tmp_topk(:,:)
    integer, allocatable :: tmp_topk_miss(:,:)
    integer, allocatable :: tmp_topk_count(:)

    allocate(tmp_keys(n_keep), tmp_counts(n_keep), tmp_rep(n_keep))
    allocate(tmp_topk(n_keep, k_top), tmp_topk_miss(n_keep, k_top), tmp_topk_count(n_keep))
    tmp_keys = window%keys(1:n_keep)
    tmp_counts = window%counts(1:n_keep)
    tmp_rep = window%rep_idx(1:n_keep)
    tmp_topk = window%topk_idx(1:n_keep, :)
    tmp_topk_miss = window%topk_miss(1:n_keep, :)
    tmp_topk_count = window%topk_count(1:n_keep)

    deallocate(window%keys, window%counts, window%rep_idx, window%topk_idx, window%topk_miss, window%topk_count)
    allocate(window%keys(n_keep), window%counts(n_keep), window%rep_idx(n_keep))
    allocate(window%topk_idx(n_keep, k_top), window%topk_miss(n_keep, k_top), window%topk_count(n_keep))
    window%keys = tmp_keys
    window%counts = tmp_counts
    window%rep_idx = tmp_rep
    window%topk_idx = tmp_topk
    window%topk_miss = tmp_topk_miss
    window%topk_count = tmp_topk_count

    deallocate(tmp_keys, tmp_counts, tmp_rep, tmp_topk, tmp_topk_miss, tmp_topk_count)
  end subroutine shrink_arrays

  subroutine update_topk(topk_idx, topk_miss, topk_count, cand_idx, miss_counts, k_top)
    integer, intent(inout) :: topk_idx(:)
    integer, intent(inout) :: topk_miss(:)
    integer, intent(inout) :: topk_count
    integer, intent(in) :: cand_idx
    integer, intent(in) :: miss_counts(:)
    integer, intent(in) :: k_top
    integer :: i, worst_pos, worst_miss

    if (cand_idx <= 0) return

    if (topk_count < k_top) then
      topk_count = topk_count + 1
      topk_idx(topk_count) = cand_idx
      topk_miss(topk_count) = miss_counts(cand_idx)
      return
    end if

    worst_pos = 1
    worst_miss = miss_counts(topk_idx(1))
    do i = 2, k_top
      if (miss_counts(topk_idx(i)) > worst_miss) then
        worst_miss = miss_counts(topk_idx(i))
        worst_pos = i
      end if
    end do

    if (miss_counts(cand_idx) < worst_miss) then
      topk_idx(worst_pos) = cand_idx
      topk_miss(worst_pos) = miss_counts(cand_idx)
    end if
  end subroutine update_topk

  integer function count_missing(packed_row, start_snp, end_snp) result(cnt)
    integer(kind=ki4), intent(in) :: packed_row(:)
    integer, intent(in) :: start_snp, end_snp
    integer :: w, offset, snp_idx
    integer :: start_word, end_word
    integer(kind=ki4) :: word_val
    integer(kind=ki1) :: g

    cnt = 0
    start_word = (start_snp - 1) / 16 + 1
    end_word = (end_snp - 1) / 16 + 1

    do w = start_word, end_word
      word_val = packed_row(w)
      do offset = 0, 15
        snp_idx = (w - 1) * 16 + offset + 1
        if (snp_idx < start_snp .or. snp_idx > end_snp) cycle
        g = int(iand(ishft(word_val, -2 * offset), 3_ki4), kind=ki1)
        if (g == 3_ki1) cnt = cnt + 1
      end do
    end do
  end function count_missing

  integer function find_key_index(window, key) result(idx)
    type(hap_window_t), intent(in) :: window
    integer(kind=ki8), intent(in) :: key
    integer :: lo, hi, mid

    idx = 0
    if (window%n_keys <= 0) return
    lo = 1
    hi = window%n_keys
    do while (lo <= hi)
      mid = (lo + hi) / 2
      if (window%keys(mid) == key) then
        idx = mid
        return
      else if (window%keys(mid) < key) then
        lo = mid + 1
      else
        hi = mid - 1
      end if
    end do
  end function find_key_index

end module M_phimpute_haplotype_lib
