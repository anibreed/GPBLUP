module M_phimpute_steps
  use M_kinds
  use M_phimpute_params, only: phimpute_params_t
  use M_phimpute_data, only: phimpute_state_t, geno_matrix_t, copy_geno_matrix, init_geno_matrix, phasing_stats_t
  use M_phimpute_map, only: load_map_sorted
  use M_phimpute_io_geno, only: load_candidate_geno, load_reference_geno
  use M_phimpute_family_phasing, only: run_family_phasing
  use M_phimpute_haplotype_lib, only: build_haplotype_panel, count_missing
  use M_phimpute_qc, only: apply_qc_filters
  use M_phimpute_haplotype_impute, only: run_haplotype_impute
  use M_phimpute_population_impute, only: run_population_impute
  use M_phimpute_pedigree_deterministic, only: run_pedigree_deterministic
  use M_ReadFile, only: fopen
  use M_Variables, only: FileInfo, MAX_RECL, LEN_STR
  use M_phimpute_output, only: write_geno_file
  implicit none

contains

  subroutine run_phimpute_pipeline(params)
    type(phimpute_params_t), intent(in) :: params
    type(phimpute_state_t) :: state
    real :: t0, t1, t2, t3
    type(geno_matrix_t) :: cand_backup, ref_backup
    type(phasing_stats_t) :: stats_ref, stats_noref
    integer(kind=ki8) :: base_missing, mendel_missing, hap_missing

    call cpu_time(t0)
    call load_input_files(params, state)
    call apply_qc_filters(params, state)
    call cpu_time(t1)
    write(*,'(A,F8.3,A)') 'TIME load_input_files: ', t1 - t0, ' s'
    base_missing = count_total_missing(state%cand)
    mendel_missing = base_missing

    if (params%enable_pedigree) then
      call step_pedigree_deterministic(params, state)
    end if
    if (params%enable_family_phasing) then
      call copy_geno_matrix(state%cand, cand_backup)
      call copy_geno_matrix(state%ref, ref_backup)

      if (params%has_reference .and. len_trim(params%ref_genofile%FileName) > 0) then
        if (state%ref%n_animals <= 0 .or. state%ref%n_snps /= state%cand%n_snps) then
          write(*,'(A)') 'REFERENCE comparison: skipped (reference size mismatch).'
        else
          write(*,'(A)') 'REFERENCE comparison: run with reference panel'
          call run_family_phasing(params, state, stats_ref, .false.)
          call copy_geno_matrix(cand_backup, state%cand)
          call copy_geno_matrix(cand_backup, state%ref)
          state%ref_panel%built = .false.
          state%has_reference = .false.

          write(*,'(A)') 'REFERENCE comparison: run without reference panel'
          call run_family_phasing(params, state, stats_noref, .false.)

          write(*,'(A,F8.3,A,F8.3)') 'REFERENCE comparison fill rate (ref/no-ref): ', &
            real(stats_ref%filled_total) / max(1.0, real(stats_ref%checked_total)), ' / ', &
            real(stats_noref%filled_total) / max(1.0, real(stats_noref%checked_total))

          call copy_geno_matrix(cand_backup, state%cand)
          call copy_geno_matrix(ref_backup, state%ref)
          state%ref_panel%built = .false.
          state%has_reference = .true.
        end if
      end if

      call cpu_time(t2)
      call run_family_phasing(params, state, stats_ref, .true.)
      call cpu_time(t3)
      write(*,'(A,F8.3,A)') 'TIME family_phasing: ', t3 - t2, ' s'

      mendel_missing = count_total_missing(state%cand)

    end if
    if (params%enable_haplotype) then
      call step_haplotype_matching(params, state)
      hap_missing = count_total_missing(state%cand)
      call log_cumulative_fill(base_missing, mendel_missing, hap_missing)
    end if
    if (params%enable_population) then
      call step_population_impute(params, state)
      call write_geno_file(trim(params%output_genofile), state%cand, state%cand_ids)
    else if (params%enable_haplotype) then
      call write_geno_file(trim(params%output_genofile), state%cand, state%cand_ids)
    end if
  end subroutine run_phimpute_pipeline

  subroutine load_input_files(params, state)
    type(phimpute_params_t), intent(in) :: params
    type(phimpute_state_t), intent(inout) :: state
    integer :: n_snps_map, n_snps_geno

    state%has_reference = params%has_reference

    call load_map_sorted(params%mapfile, state%map, params%stop_on_mismatch, params%max_autosome)
    n_snps_map = state%map%n_snps

    call load_candidate_geno(params%cand_genofile, state%map%n_snps, state%cand, state%cand_ids, &
      params%stop_on_mismatch)
    n_snps_geno = state%cand%n_snps

    if (params%has_reference .and. len_trim(params%ref_genofile%FileName) > 0) then
      call load_reference_geno(params%ref_genofile, state%map%n_snps, state%ref, state%ref_ids, &
        params%stop_on_mismatch)
      state%has_reference = .true.
      write(*,'(A,I0,A,I0)') 'REFERENCE panel: animals=', state%ref%n_animals, ' snps=', state%ref%n_snps
    else
      call copy_geno_matrix(state%cand, state%ref)
      if (allocated(state%ref_ids)) deallocate(state%ref_ids)
      if (allocated(state%cand_ids)) then
        allocate(state%ref_ids(size(state%cand_ids)))
        state%ref_ids = state%cand_ids
      end if
      state%has_reference = .false.
      write(*,'(A)') 'REFERENCE panel: using candidate GENO'
      write(*,'(A,I0,A,I0)') 'REFERENCE panel: animals=', state%ref%n_animals, ' snps=', state%ref%n_snps
    end if

    call preview_ped_records(params%pedfile)

    if (n_snps_map /= n_snps_geno) then
      write(*,'(A,I0,A,I0)') 'MAP/GENO SNP count mismatch: MAP=', n_snps_map, ' GENO=', n_snps_geno
    else
      write(*,'(A,I0)') 'MAP/GENO SNP count match: ', n_snps_map
    end if

    ! TODO: Load PED and optional reference GENO.
    !       Candidate GENO follows MAP ARRAY_ALL order.
  end subroutine load_input_files

  subroutine preview_ped_records(info)
    type(FileInfo), intent(in) :: info
    integer :: u_fd, ios, i, shown
    integer :: total
    character(len=MAX_RECL) :: raw_line
    integer :: line_len, first_pos

    if (len_trim(info%FileName) == 0) return

    u_fd = fopen(trim(info%FileName))
    shown = 0
    total = 0
    do i = 1, info%Header
      read(u_fd, '(A)', iostat=ios) raw_line
      if (ios /= 0) exit
    end do

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

      total = total + 1
      if (shown < 2) then
        shown = shown + 1
        write(*,'(A,I0,A)') 'PED PREVIEW [', shown, ']: '
        write(*,'(A)') trim(raw_line)
      end if
    end do

    write(*,'(A,I0)') 'PEDFILE records read: ', total

    close(u_fd)
  end subroutine preview_ped_records

  subroutine step_pedigree_deterministic(params, state)
    type(phimpute_params_t), intent(in) :: params
    type(phimpute_state_t), intent(inout) :: state

    call run_pedigree_deterministic(params, state)
  end subroutine step_pedigree_deterministic

  subroutine step_family_phasing(params, state)
    type(phimpute_params_t), intent(in) :: params
    type(phimpute_state_t), intent(inout) :: state

    call run_family_phasing(params, state)
  end subroutine step_family_phasing

  subroutine step_haplotype_matching(params, state)
    type(phimpute_params_t), intent(in) :: params
    type(phimpute_state_t), intent(inout) :: state
    integer(kind=ki8) :: hap_checked, hap_filled
    type(geno_matrix_t) :: ref_backup
    character(len=LEN_STR), allocatable :: ref_ids_backup(:)
    logical :: filtered
    integer :: best_window, best_shift
    real, parameter :: REF_MAX_MISSING = 0.10

    filtered = .false.
    call filter_reference_panel(state, REF_MAX_MISSING, ref_backup, ref_ids_backup, filtered)

    call get_haplotype_defaults(params, best_window, best_shift)
    if (params%hap_tune) then
      call tune_haplotype_params(params, state, best_window, best_shift)
    end if

    if (.not. state%ref_panel%built .or. state%ref_panel%window_size /= best_window .or. &
        state%ref_panel%shift_size /= best_shift) then
      call build_haplotype_panel(state%ref_panel, state%ref, best_window, best_shift)
    end if
    call run_haplotype_impute(state, hap_checked, hap_filled)
    if (filtered) call restore_reference_panel(state, ref_backup, ref_ids_backup)
  end subroutine step_haplotype_matching

  subroutine step_population_impute(params, state)
    type(phimpute_params_t), intent(in) :: params
    type(phimpute_state_t), intent(inout) :: state
    integer(kind=ki8) :: ld_checked, ld_filled

    call run_population_impute(state, params%ld_window, params%ld_window_max, ld_checked, ld_filled)
  end subroutine step_population_impute

  subroutine log_cumulative_fill(base_missing, mendel_missing, hap_missing)
    integer(kind=ki8), intent(in) :: base_missing, mendel_missing, hap_missing
    integer(kind=ki8) :: mendel_filled, hap_filled, total_filled
    real :: mendel_rate, hap_rate, total_rate

    mendel_filled = max(0_ki8, base_missing - mendel_missing)
    hap_filled = max(0_ki8, mendel_missing - hap_missing)
    total_filled = max(0_ki8, base_missing - hap_missing)
    mendel_rate = 0.0
    hap_rate = 0.0
    total_rate = 0.0
    if (base_missing > 0) mendel_rate = real(mendel_filled) / real(base_missing)
    if (mendel_missing > 0) hap_rate = real(hap_filled) / real(mendel_missing)
    if (base_missing > 0) total_rate = real(total_filled) / real(base_missing)

    write(*,'(A,I0)') 'CUMULATIVE total missing (base): ', base_missing
    write(*,'(A,I0)') 'CUMULATIVE filled by Mendel: ', mendel_filled
    write(*,'(A,I0)') 'CUMULATIVE filled by HAP: ', hap_filled
    write(*,'(A,I0)') 'CUMULATIVE total filled: ', total_filled
    write(*,'(A,I0)') 'CUMULATIVE remaining missing: ', hap_missing
    if (base_missing > 0) write(*,'(A,F8.3)') 'CUMULATIVE Mendel fill rate: ', mendel_rate
    if (mendel_missing > 0) write(*,'(A,F8.3)') 'CUMULATIVE HAP fill rate (on remaining): ', hap_rate
    if (base_missing > 0) write(*,'(A,F8.3)') 'CUMULATIVE fill rate (Mendel+HAP): ', total_rate
  end subroutine log_cumulative_fill

  subroutine get_haplotype_defaults(params, window_size, shift_size)
    type(phimpute_params_t), intent(in) :: params
    integer, intent(out) :: window_size, shift_size

    window_size = 20
    shift_size = 8
    if (params%hap_window_count > 0) then
      if (params%hap_window_candidates(1) > 0) window_size = params%hap_window_candidates(1)
    end if
    if (params%hap_shift_count > 0) then
      if (params%hap_shift_candidates(1) >= 0) shift_size = params%hap_shift_candidates(1)
    end if
  end subroutine get_haplotype_defaults

  subroutine tune_haplotype_params(params, state, best_window, best_shift)
    type(phimpute_params_t), intent(in) :: params
    type(phimpute_state_t), intent(inout) :: state
    integer, intent(inout) :: best_window, best_shift
    type(geno_matrix_t) :: cand_backup
    integer :: i, j, win, sh, sh_eff, tried
    integer(kind=ki8) :: checked, filled, best_filled
    real :: rate, best_rate

    if (state%cand%n_animals <= 0 .or. state%cand%n_snps <= 0) return
    if (params%hap_window_count <= 0 .or. params%hap_shift_count <= 0) return

    call copy_geno_matrix(state%cand, cand_backup)
    best_filled = -1_ki8
    best_rate = -1.0
    tried = 0

    do i = 1, params%hap_window_count
      win = params%hap_window_candidates(i)
      if (win <= 0) cycle
      do j = 1, params%hap_shift_count
        sh = params%hap_shift_candidates(j)
        if (sh < 0) cycle
        sh_eff = max(0, min(sh, win / 2))
        tried = tried + 1

        call copy_geno_matrix(cand_backup, state%cand)
        call build_haplotype_panel(state%ref_panel, state%ref, win, sh_eff)
        call run_haplotype_impute(state, checked, filled)

        rate = 0.0
        if (checked > 0) rate = real(filled) / real(checked)

        if (filled > best_filled .or. (filled == best_filled .and. rate > best_rate)) then
          best_filled = filled
          best_rate = rate
          best_window = win
          best_shift = sh_eff
        end if
      end do
    end do

    call copy_geno_matrix(cand_backup, state%cand)
    write(*,'(A,I0,A,I0,A,I0,A,F8.3)') 'HAP TUNE best window=', best_window, ' shift=', best_shift, &
      ' filled=', best_filled, ' rate=', best_rate
  end subroutine tune_haplotype_params

  subroutine filter_reference_panel(state, max_missing_rate, ref_backup, ref_ids_backup, filtered)
    type(phimpute_state_t), intent(inout) :: state
    real, intent(in) :: max_missing_rate
    type(geno_matrix_t), intent(out) :: ref_backup
    character(len=LEN_STR), allocatable, intent(out) :: ref_ids_backup(:)
    logical, intent(out) :: filtered
    integer :: a, kept, new_idx
    integer, allocatable :: keep_idx(:)
    integer(kind=ki8) :: miss_count
    real :: miss_rate

    filtered = .false.
    if (state%ref%n_animals <= 0 .or. state%ref%n_snps <= 0) return
    if (max_missing_rate <= 0.0 .or. max_missing_rate >= 1.0) return

    allocate(keep_idx(state%ref%n_animals))
    kept = 0
    do a = 1, state%ref%n_animals
      miss_count = count_missing(state%ref%packed(:, a), 1, state%ref%n_snps)
      miss_rate = real(miss_count) / real(state%ref%n_snps)
      if (miss_rate <= max_missing_rate) then
        kept = kept + 1
        keep_idx(kept) = a
      end if
    end do

    if (kept <= 0 .or. kept == state%ref%n_animals) then
      deallocate(keep_idx)
      return
    end if

    call copy_geno_matrix(state%ref, ref_backup)
    if (allocated(ref_ids_backup)) deallocate(ref_ids_backup)
    if (allocated(state%ref_ids)) then
      allocate(ref_ids_backup(size(state%ref_ids)))
      ref_ids_backup = state%ref_ids
    end if

    call init_geno_matrix(state%ref, kept, ref_backup%n_snps)
    if (allocated(state%ref_ids)) deallocate(state%ref_ids)
    allocate(state%ref_ids(kept))

    do new_idx = 1, kept
      a = keep_idx(new_idx)
      state%ref%packed(:, new_idx) = ref_backup%packed(:, a)
      if (allocated(ref_ids_backup)) state%ref_ids(new_idx) = ref_ids_backup(a)
    end do

    state%ref_panel%built = .false.
    filtered = .true.
    write(*,'(A,I0,A,I0,A,F5.3,A)') 'REFERENCE filter: kept ', kept, ' of ', ref_backup%n_animals, &
      ' (max_miss=', max_missing_rate, ')'

    deallocate(keep_idx)
  end subroutine filter_reference_panel

  subroutine restore_reference_panel(state, ref_backup, ref_ids_backup)
    type(phimpute_state_t), intent(inout) :: state
    type(geno_matrix_t), intent(in) :: ref_backup
    character(len=LEN_STR), allocatable, intent(in) :: ref_ids_backup(:)

    call copy_geno_matrix(ref_backup, state%ref)
    if (allocated(state%ref_ids)) deallocate(state%ref_ids)
    if (allocated(ref_ids_backup)) then
      allocate(state%ref_ids(size(ref_ids_backup)))
      state%ref_ids = ref_ids_backup
    end if
    state%ref_panel%built = .false.
  end subroutine restore_reference_panel

  function count_total_missing(geno) result(total_missing)
    type(geno_matrix_t), intent(in) :: geno
    integer(kind=ki8) :: total_missing
    integer :: a

    total_missing = 0_ki8
    if (geno%n_snps <= 0 .or. geno%n_animals <= 0) return
    do a = 1, geno%n_animals
      total_missing = total_missing + count_missing(geno%packed(:, a), 1, geno%n_snps)
    end do
  end function count_total_missing

end module M_phimpute_steps
