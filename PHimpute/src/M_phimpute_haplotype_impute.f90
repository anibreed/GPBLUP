module M_phimpute_haplotype_impute
  use M_kinds
  use M_phimpute_data, only: phimpute_state_t
  use M_phimpute_haplotype_lib, only: hash_window, find_key_index, count_missing
  use M_phimpute_bitpack, only: GENO_MISSING
  use omp_lib
  implicit none

  public :: run_haplotype_impute

contains

  subroutine run_haplotype_impute(state, checked_total, filled_total)
    type(phimpute_state_t), intent(inout) :: state
    integer(kind=ki8), intent(out), optional :: checked_total, filled_total

    integer :: a, w
    integer :: start_snp, end_snp
    integer :: window_len, min_non_missing
    integer :: miss_cand, miss_ref
    integer :: rep_idx, key_idx
    integer(kind=ki8) :: key
    integer(kind=ki8) :: checked_local, filled_local
    integer(kind=ki8) :: matched_windows

    if (.not. state%ref_panel%built) return
    if (state%cand%n_animals <= 0 .or. state%cand%n_snps <= 0) return

    checked_local = 0
    filled_local = 0
    matched_windows = 0

    !$omp parallel do default(shared) private(a, w, start_snp, end_snp, rep_idx, key_idx, key) &
    !$omp reduction(+:checked_local, filled_local, matched_windows) schedule(static)
    do a = 1, state%cand%n_animals
      do w = 1, state%ref_panel%n_windows
        start_snp = state%ref_panel%windows(w)%start_snp
        end_snp = state%ref_panel%windows(w)%end_snp

        window_len = end_snp - start_snp + 1
        min_non_missing = (window_len * 7) / 10
        if (min_non_missing < 10) min_non_missing = 10

        miss_cand = count_missing(state%cand%packed(:, a), start_snp, end_snp)
        if (window_len - miss_cand < min_non_missing) cycle

        call hash_window(state%cand%packed(:, a), start_snp, end_snp, key)
        key_idx = find_key_index(state%ref_panel%windows(w), key)
        if (key_idx <= 0) cycle

        matched_windows = matched_windows + 1

        rep_idx = state%ref_panel%windows(w)%rep_idx(key_idx)
        if (rep_idx <= 0) cycle

        miss_ref = count_missing(state%ref%packed(:, rep_idx), start_snp, end_snp)
        if (window_len - miss_ref < min_non_missing) cycle

        call fill_window_topk(state, a, state%ref_panel%windows(w)%topk_idx(key_idx, :), &
          state%ref_panel%windows(w)%topk_miss(key_idx, :), state%ref_panel%windows(w)%topk_count(key_idx), &
          start_snp, end_snp, checked_local, filled_local)
      end do
    end do
    !$omp end parallel do

    write(*,'(A,I0)') 'HAP IMPUTE checked missing genotypes: ', checked_local
    write(*,'(A,I0)') 'HAP IMPUTE filled genotypes: ', filled_local
    if (checked_local > 0) then
      write(*,'(A,F8.3)') 'HAP IMPUTE fill rate: ', real(filled_local) / real(checked_local)
    end if
    write(*,'(A,I0)') 'HAP IMPUTE matched windows: ', matched_windows

    if (present(checked_total)) checked_total = checked_local
    if (present(filled_total)) filled_total = filled_local
  end subroutine run_haplotype_impute

  subroutine fill_window_topk(state, cand_idx, topk_idx, topk_miss, topk_count, start_snp, end_snp, checked_total, filled_total)
    type(phimpute_state_t), intent(inout) :: state
    integer, intent(in) :: cand_idx
    integer, intent(in) :: topk_idx(:)
    integer, intent(in) :: topk_miss(:)
    integer, intent(in) :: topk_count
    integer, intent(in) :: start_snp, end_snp
    integer(kind=ki8), intent(inout) :: checked_total, filled_total

    integer :: w, offset, snp_idx
    integer :: start_word, end_word, window_len
    integer(kind=ki4) :: cand_word
    integer(kind=ki1) :: g_cand, g_vote

    start_word = (start_snp - 1) / 16 + 1
    end_word = (end_snp - 1) / 16 + 1
    window_len = end_snp - start_snp + 1

    do w = start_word, end_word
      cand_word = state%cand%packed(w, cand_idx)
      do offset = 0, 15
        snp_idx = (w - 1) * 16 + offset + 1
        if (snp_idx < start_snp .or. snp_idx > end_snp) cycle
        if (snp_idx > state%cand%n_snps) exit

        g_cand = int(iand(ishft(cand_word, -2 * offset), 3_ki4), kind=ki1)
        if (g_cand /= GENO_MISSING) cycle

        checked_total = checked_total + 1
        call majority_vote(state, topk_idx, topk_miss, topk_count, snp_idx, window_len, g_vote)
        if (g_vote /= GENO_MISSING) then
          call set_word_geno(cand_word, offset, g_vote)
          filled_total = filled_total + 1
        end if
      end do

      state%cand%packed(w, cand_idx) = cand_word
    end do
  end subroutine fill_window_topk

  subroutine majority_vote(state, topk_idx, topk_miss, topk_count, snp_idx, window_len, g_out)
    type(phimpute_state_t), intent(in) :: state
    integer, intent(in) :: topk_idx(:)
    integer, intent(in) :: topk_miss(:)
    integer, intent(in) :: topk_count
    integer, intent(in) :: snp_idx
    integer, intent(in) :: window_len
    integer(kind=ki1), intent(out) :: g_out
    integer :: i
    integer(kind=ki1) :: g
    integer(kind=ki8) :: cnt0, cnt1, cnt2
    integer :: weight

    cnt0 = 0
    cnt1 = 0
    cnt2 = 0
    do i = 1, topk_count
      if (topk_idx(i) <= 0) cycle
      g = int(iand(ishft(state%ref%packed((snp_idx - 1) / 16 + 1, topk_idx(i)), -2 * mod(snp_idx - 1, 16)), 3_ki4), kind=ki1)
      weight = window_len - topk_miss(i)
      if (weight < 1) weight = 1
      select case (g)
      case (0)
        cnt0 = cnt0 + weight
      case (1)
        cnt1 = cnt1 + weight
      case (2)
        cnt2 = cnt2 + weight
      case default
        continue
      end select
    end do

    g_out = GENO_MISSING
    if (cnt0 == 0 .and. cnt1 == 0 .and. cnt2 == 0) return

    if (cnt0 >= cnt1 .and. cnt0 >= cnt2) then
      g_out = 0
    else if (cnt1 >= cnt0 .and. cnt1 >= cnt2) then
      g_out = 1
    else
      g_out = 2
    end if
  end subroutine majority_vote

  subroutine set_word_geno(word_val, offset, geno)
    integer(kind=ki4), intent(inout) :: word_val
    integer, intent(in) :: offset
    integer(kind=ki1), intent(in) :: geno
    integer(kind=ki4) :: mask

    mask = ishft(3_ki4, 2 * offset)
    word_val = iand(word_val, not(mask))
    word_val = ior(word_val, ishft(int(geno, ki4), 2 * offset))
  end subroutine set_word_geno

end module M_phimpute_haplotype_impute
