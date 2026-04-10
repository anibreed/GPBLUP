module M_phimpute_population_impute
  use M_kinds
  use M_phimpute_data, only: phimpute_state_t, geno_matrix_t
  use M_phimpute_bitpack, only: get_geno_value, set_geno_value, GENO_MISSING
  implicit none

  public :: run_population_impute
  integer, parameter :: MAX_LD_WINDOW = 25

contains

  subroutine run_population_impute(state, ld_window, ld_window_max, checked_total, filled_total)
    type(phimpute_state_t), intent(inout) :: state
    integer, intent(in) :: ld_window
    integer, intent(in) :: ld_window_max
    integer(kind=ki8), intent(out), optional :: checked_total, filled_total

    integer :: n_snps, a, s
    integer(kind=ki1) :: g, g_left, g_right
    integer(kind=ki8), allocatable :: freq_counts(:,:)
    integer(kind=ki8), allocatable :: left_counts(:,:,:,:)
    integer(kind=ki8), allocatable :: right_counts(:,:,:,:)
    integer(kind=ki8) :: checked_local, filled_local
    integer(kind=ki8) :: counts(3)
    integer :: pick_idx
    integer :: d, window

    if (state%cand%n_animals <= 0 .or. state%cand%n_snps <= 0) return

    n_snps = state%cand%n_snps
    window = max(1, min(ld_window, min(ld_window_max, MAX_LD_WINDOW)))
    call build_ld_counts(state%ref, window, freq_counts, left_counts, right_counts)

    checked_local = 0_ki8
    filled_local = 0_ki8

    do a = 1, state%cand%n_animals
      do s = 1, n_snps
        g = get_geno_value(state%cand%packed(:, a), s)
        if (g /= GENO_MISSING) cycle

        checked_local = checked_local + 1_ki8
        counts = freq_counts(s, :)

        do d = 1, window
          if (s > d) then
            g_left = get_geno_value(state%cand%packed(:, a), s - d)
            if (g_left /= GENO_MISSING) then
              counts = counts + left_counts(s, d, g_left + 1, :)
            end if
          end if

          if (s + d <= n_snps) then
            g_right = get_geno_value(state%cand%packed(:, a), s + d)
            if (g_right /= GENO_MISSING) then
              counts = counts + right_counts(s, d, g_right + 1, :)
            end if
          end if
        end do

        if (counts(1) + counts(2) + counts(3) <= 0_ki8) cycle

        pick_idx = maxloc(counts, dim=1)
        call set_geno_value(state%cand%packed(:, a), s, int(pick_idx - 1, kind=ki1))
        filled_local = filled_local + 1_ki8
      end do
    end do

    write(*,'(A,I0,A,I0)') 'LD IMPUTE checked missing genotypes: ', checked_local, ' window=', window
    write(*,'(A,I0)') 'LD IMPUTE filled genotypes: ', filled_local
    if (checked_local > 0) then
      write(*,'(A,F8.3)') 'LD IMPUTE fill rate: ', real(filled_local) / real(checked_local)
    end if

    if (present(checked_total)) checked_total = checked_local
    if (present(filled_total)) filled_total = filled_local

    if (allocated(freq_counts)) deallocate(freq_counts)
    if (allocated(left_counts)) deallocate(left_counts)
    if (allocated(right_counts)) deallocate(right_counts)
  end subroutine run_population_impute

  subroutine build_ld_counts(geno, window, freq_counts, left_counts, right_counts)
    type(geno_matrix_t), intent(in) :: geno
    integer, intent(in) :: window
    integer(kind=ki8), allocatable, intent(out) :: freq_counts(:,:)
    integer(kind=ki8), allocatable, intent(out) :: left_counts(:,:,:,:)
    integer(kind=ki8), allocatable, intent(out) :: right_counts(:,:,:,:)

    integer :: n_snps, a, s, d, left_idx, right_idx
    integer(kind=ki1) :: g, g_left, g_right

    n_snps = geno%n_snps
    if (n_snps <= 0 .or. geno%n_animals <= 0) then
      allocate(freq_counts(0, 0), left_counts(0, 0, 0, 0), right_counts(0, 0, 0, 0))
      return
    end if

    allocate(freq_counts(n_snps, 3))
    allocate(left_counts(n_snps, window, 3, 3))
    allocate(right_counts(n_snps, window, 3, 3))

    freq_counts = 0_ki8
    left_counts = 0_ki8
    right_counts = 0_ki8

    do a = 1, geno%n_animals
      do s = 1, n_snps
        g = get_geno_value(geno%packed(:, a), s)
        if (g == GENO_MISSING) cycle
        freq_counts(s, g + 1) = freq_counts(s, g + 1) + 1_ki8

        do d = 1, window
          left_idx = s - d
          if (left_idx >= 1) then
            g_left = get_geno_value(geno%packed(:, a), left_idx)
            if (g_left /= GENO_MISSING) then
              left_counts(s, d, g_left + 1, g + 1) = left_counts(s, d, g_left + 1, g + 1) + 1_ki8
            end if
          end if

          right_idx = s + d
          if (right_idx <= n_snps) then
            g_right = get_geno_value(geno%packed(:, a), right_idx)
            if (g_right /= GENO_MISSING) then
              right_counts(s, d, g_right + 1, g + 1) = right_counts(s, d, g_right + 1, g + 1) + 1_ki8
            end if
          end if
        end do
      end do
    end do
  end subroutine build_ld_counts

end module M_phimpute_population_impute
