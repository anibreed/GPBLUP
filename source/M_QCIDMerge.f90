module M_QCIDMerge
  use M_Kinds
  use M_Variables
  use M_Hashtable
  use H_PedRead
  use M_QCMendel, only: is_duo_mendel_error, is_trio_mendel_error
  use M_QCSNPSelect, only: select_parent_seek_snps
  use M_QCIDMergePairEval, only: compute_pair_similarity_sampled_mod, count_non_missing_sampled_mod
  implicit none

  private
  public :: set_id_merge_context
  public :: merge_highly_identical_animals
  public :: finalize_deferred_identical_pairs

  type(PEDHashTable), pointer, save :: PHT => null()
  type(SNPInfo), pointer, save :: map_records(:) => null()
  integer(kind=ki1), pointer, save :: geno_matrix(:,:) => null()
  character(len=LEN_STR), pointer, save :: geno_animal_ids(:) => null()
  character(len=LEN_STR), save :: ped_animal_key_field = 'ANIMAL_ID'
  real(r8), pointer, save :: snp_call_rate(:) => null()
  real(r8), pointer, save :: snp_maf(:) => null()
  real(r8), pointer, save :: snp_hwe_pvalue(:) => null()
  real(r8), save :: hwe_threshold = 0.0_r8
  integer, pointer, save :: id_merge_seek_flag(:) => null()

contains

  subroutine set_id_merge_context(in_pht, in_map_records, in_geno_matrix, in_geno_animal_ids, &
                                  in_ped_key_field, in_snp_call_rate, in_snp_maf, &
                                  in_snp_hwe_pvalue, in_hwe_threshold, in_id_merge_seek_flag)
    type(PEDHashTable), target, intent(inout) :: in_pht
    type(SNPInfo), target, intent(in) :: in_map_records(:)
    integer(kind=ki1), target, intent(inout) :: in_geno_matrix(:,:)
    character(len=LEN_STR), target, intent(inout) :: in_geno_animal_ids(:)
    character(len=*), intent(in) :: in_ped_key_field
    real(r8), target, intent(in), optional :: in_snp_call_rate(:), in_snp_maf(:), in_snp_hwe_pvalue(:)
    real(r8), intent(in), optional :: in_hwe_threshold
    integer, target, intent(inout), optional :: in_id_merge_seek_flag(:)

    PHT => in_pht
    map_records => in_map_records
    geno_matrix => in_geno_matrix
    geno_animal_ids => in_geno_animal_ids
    ped_animal_key_field = trim(in_ped_key_field)

    if (present(in_snp_call_rate)) then
      snp_call_rate => in_snp_call_rate
    else
      nullify(snp_call_rate)
    end if
    if (present(in_snp_maf)) then
      snp_maf => in_snp_maf
    else
      nullify(snp_maf)
    end if
    if (present(in_snp_hwe_pvalue)) then
      snp_hwe_pvalue => in_snp_hwe_pvalue
    else
      nullify(snp_hwe_pvalue)
    end if
    if (present(in_hwe_threshold)) hwe_threshold = in_hwe_threshold
    if (present(in_id_merge_seek_flag)) then
      id_merge_seek_flag => in_id_merge_seek_flag
    else
      nullify(id_merge_seek_flag)
    end if
  end subroutine set_id_merge_context

  subroutine merge_highly_identical_animals(in_n_animals, in_n_snps, in_snp_qc_status, in_match_threshold)
    integer, intent(in) :: in_n_animals, in_n_snps
    integer, intent(in) :: in_snp_qc_status(*)
    real(r8), intent(in) :: in_match_threshold

    integer,  parameter :: MIN_COMPARE_SNPS  = 500
    integer,  parameter :: FINGERPRINT_SNPS_MIN  = 16
    integer,  parameter :: FINGERPRINT_SNPS_MAX  = 512
    integer,  parameter :: CONFIRM_MIN_SNPS  = 2000

    integer :: i, j, k, s, p, n_comp
    integer :: keep_row, drop_row
    integer :: n_merged, n_conflict_to_missing, n_ref_updates
    integer :: n_sample, n_fp, n_valid_rows, n_key_slots
    integer :: idx, n_bucket, pos, g, call_bin
    integer(kind=8) :: n_pairs_before
    integer(kind=8) :: keyv
    integer :: n_candidate_pairs, n_pass_pairs
    integer :: n_confirm_pass, n_confirm_fail, n_confirm_valid
    real(r8) :: match_rate, rate_local, confirm_rate, confirm_match_threshold
    real(r8) :: t_start, t_bucket_done, t_eval_start, t_eval_done, t_end
    integer(kind=8), allocatable :: fp_key(:), key_slot(:)
    logical, allocatable :: key_used(:), row_dropped(:), pair_pass(:)
    integer, allocatable :: key_head(:), row_next(:), bucket_rows(:)
    integer, allocatable :: pair_i(:), pair_j(:), sample_idx(:), non_missing_qc(:)
    real(r8), allocatable :: pair_rate(:), sample_maf(:)
    logical :: resolved_by_ped
    character(len=LEN_STR) :: keep_id, drop_id

    if (.false.) then
      if (in_snp_qc_status(1) /= 0) then
      end if
    end if

    t_start = 0.0_r8
    t_bucket_done = 0.0_r8
    t_eval_start = 0.0_r8
    t_eval_done = 0.0_r8
    t_end = 0.0_r8
    call cpu_time(t_start)
    confirm_match_threshold = real(POPQC%id_merge_confirm_match_threshold, r8)
    if (confirm_match_threshold <= 0.0_r8 .or. confirm_match_threshold > 1.0_r8) then
      confirm_match_threshold = 0.995_r8
    end if
    if (.not. associated(map_records)) then
      print *, '  No MAP records in ID merge context. Step 6 skipped.'
      return
    end if
    if (.not. associated(snp_maf)) then
      print *, '  No SNP MAF profile in ID merge context. Step 6 skipped.'
      return
    end if
    if (.not. associated(snp_call_rate)) then
      print *, '  No SNP call-rate profile in ID merge context. Step 6 skipped.'
      return
    end if
    if (associated(snp_hwe_pvalue)) then
      call select_parent_seek_snps(in_n_snps, map_records, snp_call_rate, snp_maf, sample_idx, sample_maf, &
                                   n_sample, hwe_threshold, snp_hwe_pvalue, 'ID_MERGE')
    else
      call select_parent_seek_snps(in_n_snps, map_records, snp_call_rate, snp_maf, sample_idx, sample_maf, &
                                   n_sample, hwe_threshold, in_context='ID_MERGE')
    end if

    allocate(row_dropped(in_n_animals), non_missing_qc(in_n_animals), row_next(in_n_animals))
    allocate(fp_key(in_n_animals), bucket_rows(in_n_animals))
    row_dropped = .false.
    non_missing_qc = 0
    row_next = 0
    fp_key = -1_8

    if (n_sample <= 0) then
      print *, '  No SNPs sampled for ID merge. Step 6 skipped.'
      if (allocated(sample_idx)) deallocate(sample_idx)
      if (allocated(sample_maf)) deallocate(sample_maf)
      if (allocated(row_dropped)) deallocate(row_dropped)
      if (allocated(non_missing_qc)) deallocate(non_missing_qc)
      if (allocated(row_next)) deallocate(row_next)
      if (allocated(fp_key)) deallocate(fp_key)
      if (allocated(bucket_rows)) deallocate(bucket_rows)
      return
    end if

    print '(a, i8)', '  SNPs sampled for ID-merge concordance: ', n_sample

    do i = 1, in_n_animals
      if (len_trim(geno_animal_ids(i)) == 0) cycle
      do k = 1, n_sample
        s = sample_idx(k)
        if (int(geno_matrix(i, s)) /= 9) non_missing_qc(i) = non_missing_qc(i) + 1
      end do
    end do

    n_valid_rows = count([(len_trim(geno_animal_ids(i)) > 0, i=1,in_n_animals)])
    n_pairs_before = int(n_valid_rows, 8) * int(max(0, n_valid_rows - 1), 8) / 2_8

    ! Adaptive fingerprint size based on number of sampled SNPs
    ! For large SNP sets (>2400), use more SNPs for better discrimination power
    ! Otherwise use minimum 16 SNPs
    n_fp = min(FINGERPRINT_SNPS_MAX, max(FINGERPRINT_SNPS_MIN, n_sample / 150))
    do i = 1, in_n_animals
      if (len_trim(geno_animal_ids(i)) == 0) cycle
      keyv = 0_8
      do k = 1, n_fp
        s = sample_idx(k)
        g = int(geno_matrix(i, s))
        if (g < 0 .or. g > 2) g = 3
        keyv = ishft(keyv, 2) + int(g, 8)
      end do
      call_bin = min(255, non_missing_qc(i) / max(1, n_fp))
      keyv = ishft(keyv, 8) + int(call_bin, 8)
      fp_key(i) = keyv
    end do

    n_key_slots = htab_get_prime_size(max(101, 2 * max(1, n_valid_rows) + 1))
    allocate(key_used(n_key_slots), key_slot(n_key_slots), key_head(n_key_slots))
    key_used = .false.
    key_slot = 0_8
    key_head = 0

    do i = 1, in_n_animals
      if (len_trim(geno_animal_ids(i)) == 0) cycle
      if (fp_key(i) < 0_8) cycle
      idx = int(modulo(fp_key(i), int(n_key_slots, 8))) + 1
      do
        if (.not. key_used(idx)) then
          key_used(idx) = .true.
          key_slot(idx) = fp_key(i)
          key_head(idx) = i
          row_next(i) = 0
          exit
        else if (key_slot(idx) == fp_key(i)) then
          row_next(i) = key_head(idx)
          key_head(idx) = i
          exit
        else
          idx = idx + 1
          if (idx > n_key_slots) idx = 1
        end if
      end do
    end do

    n_candidate_pairs = 0
    do idx = 1, n_key_slots
      if (.not. key_used(idx)) cycle
      n_bucket = 0
      i = key_head(idx)
      do while (i > 0)
        n_bucket = n_bucket + 1
        bucket_rows(n_bucket) = i
        i = row_next(i)
      end do
      if (n_bucket >= 2) n_candidate_pairs = n_candidate_pairs + (n_bucket * (n_bucket - 1)) / 2
    end do

    if (n_candidate_pairs > 0) then
      allocate(pair_i(n_candidate_pairs), pair_j(n_candidate_pairs), pair_rate(n_candidate_pairs), pair_pass(n_candidate_pairs))
      pair_rate = 0.0_r8
      pair_pass = .false.

      pos = 0
      do idx = 1, n_key_slots
        if (.not. key_used(idx)) cycle
        n_bucket = 0
        i = key_head(idx)
        do while (i > 0)
          n_bucket = n_bucket + 1
          bucket_rows(n_bucket) = i
          i = row_next(i)
        end do
        if (n_bucket < 2) cycle
        do i = 1, n_bucket - 1
          do j = i + 1, n_bucket
            pos = pos + 1
            pair_i(pos) = min(bucket_rows(i), bucket_rows(j))
            pair_j(pos) = max(bucket_rows(i), bucket_rows(j))
          end do
        end do
      end do
    else
      allocate(pair_i(0), pair_j(0), pair_rate(0), pair_pass(0))
    end if

    call cpu_time(t_bucket_done)
    print '(a, i0)', '  ID concordance pair count (full):    ', n_pairs_before
    print '(a, i0)', '  ID concordance pair count (bucket):  ', n_candidate_pairs

    call cpu_time(t_eval_start)
!$omp parallel do default(shared) private(p,i,j,n_comp,rate_local) schedule(dynamic)
    do p = 1, n_candidate_pairs
      i = pair_i(p)
      j = pair_j(p)
      call compute_pair_similarity_sampled_mod(geno_matrix, i, j, n_sample, sample_idx, n_comp, rate_local)
      if (n_comp < MIN_COMPARE_SNPS) cycle
      if (rate_local >= in_match_threshold) then
        pair_pass(p) = .true.
        pair_rate(p) = rate_local
      end if
    end do
!$omp end parallel do
    call cpu_time(t_eval_done)

    n_pass_pairs = count(pair_pass)
    n_confirm_pass = 0
    n_confirm_fail = 0
    print '(a, i0)', '  ID concordance pairs passing cutoff: ', n_pass_pairs

    n_merged = 0
    n_conflict_to_missing = 0
    n_ref_updates = 0

    do p = 1, n_candidate_pairs
      if (.not. pair_pass(p)) cycle
      i = pair_i(p)
      j = pair_j(p)
      if (len_trim(geno_animal_ids(i)) == 0) cycle
      if (len_trim(geno_animal_ids(j)) == 0) cycle
      if (row_dropped(i) .or. row_dropped(j)) cycle
      if (trim(geno_animal_ids(i)) == trim(geno_animal_ids(j))) cycle

      match_rate = pair_rate(p)
      call compute_pair_similarity_confirm(i, j, in_n_snps, n_confirm_valid, confirm_rate)
      if (n_confirm_valid >= CONFIRM_MIN_SNPS) then
        if (confirm_rate < confirm_match_threshold) then
          n_confirm_fail = n_confirm_fail + 1
          cycle
        else
          n_confirm_pass = n_confirm_pass + 1
        end if
      end if

      call resolve_identical_pair_id(i, j, in_n_animals, n_sample, sample_idx, &
                                     keep_row, drop_row, resolved_by_ped)
      if (.not. resolved_by_ped) then
        if (associated(id_merge_seek_flag)) then
          id_merge_seek_flag(i) = 1
          id_merge_seek_flag(j) = 1
        end if
        print '(a,1x,a,1x,a,1x,f7.4)', &
          '  ID merge deferred (send to parent-seek):', trim(geno_animal_ids(i)), &
          trim(geno_animal_ids(j)), match_rate
        cycle
      end if

      keep_id = trim(geno_animal_ids(keep_row))
      drop_id = trim(geno_animal_ids(drop_row))

      call merge_genotype_rows(keep_row, drop_row, in_n_snps, n_conflict_to_missing)
      call rewrite_parent_references(drop_id, keep_id, in_n_animals, n_ref_updates)

      geno_matrix(drop_row, :) = 9_ki1
      geno_animal_ids(drop_row) = ''
      row_dropped(drop_row) = .true.

      call count_non_missing_sampled_mod(geno_matrix, keep_row, n_sample, sample_idx, non_missing_qc(keep_row))

      n_merged = n_merged + 1
      print '(a,1x,a,1x,a,1x,f7.4,1x,a,1x,a)', &
        '  ID merge:', trim(drop_id), '->', match_rate, 'kept as', trim(keep_id)
      print '(a,1x,a,1x,a,1x,a)', '    ID remap:', trim(drop_id), '=>', trim(keep_id)
    end do

    call cpu_time(t_end)

    print '(a, i8)', '  Merged duplicate-ID records:        ', n_merged
    print '(a, i8)', '  Confirm pass (full-autosome):       ', n_confirm_pass
    print '(a, i8)', '  Confirm fail (full-autosome):       ', n_confirm_fail
    print '(a, f6.3)', '  Confirm threshold (match):           ', confirm_match_threshold
    print '(a, i8)', '  Conflicting SNPs set to missing:    ', n_conflict_to_missing
    print '(a, i8)', '  Pedigree parent references rewired: ', n_ref_updates
    print '(a, f10.3)', '  ID concordance bucketization time (s): ', t_bucket_done - t_start
    print '(a, f10.3)', '  ID concordance pair-eval time (s):     ', t_eval_done - t_eval_start
    print '(a, f10.3)', '  ID concordance total time (s):         ', t_end - t_start

    if (allocated(sample_idx)) deallocate(sample_idx)
    if (allocated(sample_maf)) deallocate(sample_maf)
    if (allocated(row_dropped)) deallocate(row_dropped)
    if (allocated(non_missing_qc)) deallocate(non_missing_qc)
    if (allocated(row_next)) deallocate(row_next)
    if (allocated(fp_key)) deallocate(fp_key)
    if (allocated(bucket_rows)) deallocate(bucket_rows)
    if (allocated(key_used)) deallocate(key_used)
    if (allocated(key_slot)) deallocate(key_slot)
    if (allocated(key_head)) deallocate(key_head)
    if (allocated(pair_i)) deallocate(pair_i)
    if (allocated(pair_j)) deallocate(pair_j)
    if (allocated(pair_rate)) deallocate(pair_rate)
    if (allocated(pair_pass)) deallocate(pair_pass)
  end subroutine merge_highly_identical_animals

  subroutine compute_pair_similarity_confirm(in_row_a, in_row_b, in_n_snps, out_valid, out_match_rate)
    integer, intent(in) :: in_row_a, in_row_b, in_n_snps
    integer, intent(out) :: out_valid
    real(r8), intent(out) :: out_match_rate

    integer, parameter :: AUTOSOME_MIN = 1, AUTOSOME_MAX = 18
    integer :: k, j, g1, g2, n_match

    out_valid = 0
    out_match_rate = 0.0_r8
    n_match = 0

    if (.not. associated(map_records)) return

    do k = 1, size(map_records)
      if (map_records(k)%Chr < AUTOSOME_MIN .or. map_records(k)%Chr > AUTOSOME_MAX) cycle
      j = map_records(k)%Array_All
      if (j < 1 .or. j > in_n_snps) cycle

      g1 = int(geno_matrix(in_row_a, j))
      if (g1 == 9) cycle
      g2 = int(geno_matrix(in_row_b, j))
      if (g2 == 9) cycle
      out_valid = out_valid + 1
      if (g1 == g2) n_match = n_match + 1
    end do

    if (out_valid > 0) out_match_rate = real(n_match, r8) / real(out_valid, r8)
  end subroutine compute_pair_similarity_confirm

  subroutine resolve_identical_pair_id(in_row_a, in_row_b, in_n_animals, in_n_snps, in_snp_idx, &
                                       out_keep_row, out_drop_row, out_resolved)
    integer, intent(in)  :: in_row_a, in_row_b, in_n_animals
    integer, intent(in)  :: in_n_snps
    integer, intent(in)  :: in_snp_idx(in_n_snps)
    integer, intent(out) :: out_keep_row, out_drop_row
    logical, intent(out) :: out_resolved

    integer,  parameter :: MIN_EVAL_SNPS = 100
    real(r8), parameter :: TIE_TOL = 0.005_r8

    integer :: c, r_other, row_ref
    integer :: sire_row, dam_row, k, snp, g_c, g_s, g_d
    integer :: n_valid, n_err
    real(r8) :: score(2), err_rate
    integer :: rows(2)
    character(len=LEN_STR) :: id_c
    type(PEDInfo) :: ped_c
    logical :: found
    integer :: nm_a, nm_b
    integer :: par_row

    rows(1) = in_row_a
    rows(2) = in_row_b
    row_ref = in_row_a
    score = -1.0_r8

    do c = 1, 2
      id_c = trim(geno_animal_ids(rows(c)))
      if (len_trim(id_c) == 0) cycle

      found = pht_search(PHT, id_c, ped_c)
      if (.not. found) cycle

      sire_row = 0
      dam_row = 0
      if (len_trim(ped_c%SIRE) > 0 .and. trim(ped_c%SIRE) /= '0') then
        do r_other = 1, in_n_animals
          if (len_trim(geno_animal_ids(r_other)) == 0) cycle
          if (trim(geno_animal_ids(r_other)) == trim(ped_c%SIRE)) then
            sire_row = r_other
            exit
          end if
        end do
      end if
      if (len_trim(ped_c%DAM) > 0 .and. trim(ped_c%DAM) /= '0') then
        do r_other = 1, in_n_animals
          if (len_trim(geno_animal_ids(r_other)) == 0) cycle
          if (trim(geno_animal_ids(r_other)) == trim(ped_c%DAM)) then
            dam_row = r_other
            exit
          end if
        end do
      end if

      if (sire_row == 0 .and. dam_row == 0) cycle

      n_valid = 0
      n_err = 0
      if (sire_row > 0 .and. dam_row > 0) then
        do k = 1, in_n_snps
          snp = in_snp_idx(k)
          g_c = int(geno_matrix(row_ref, snp))
          if (g_c == 9) cycle
          g_s = int(geno_matrix(sire_row, snp))
          if (g_s == 9) cycle
          g_d = int(geno_matrix(dam_row, snp))
          if (g_d == 9) cycle
          n_valid = n_valid + 1
          if (is_trio_mendel_error(g_c, g_s, g_d)) n_err = n_err + 1
        end do
      else
        par_row = sire_row
        if (par_row == 0) par_row = dam_row
        do k = 1, in_n_snps
          snp = in_snp_idx(k)
          g_c = int(geno_matrix(row_ref, snp))
          if (g_c == 9) cycle
          g_s = int(geno_matrix(par_row, snp))
          if (g_s == 9) cycle
          n_valid = n_valid + 1
          if (is_duo_mendel_error(g_c, g_s)) n_err = n_err + 1
        end do
      end if

      if (n_valid >= MIN_EVAL_SNPS) then
        err_rate = real(n_err, r8) / real(n_valid, r8)
        score(c) = 1.0_r8 - err_rate
      end if
    end do

    out_resolved = .false.
    if (score(2) > score(1) + TIE_TOL) then
      out_keep_row = in_row_b
      out_drop_row = in_row_a
      out_resolved = .true.
    else if (score(1) > score(2) + TIE_TOL) then
      out_keep_row = in_row_a
      out_drop_row = in_row_b
      out_resolved = .true.
    else
      nm_a = 0
      nm_b = 0
      do k = 1, in_n_snps
        snp = in_snp_idx(k)
        if (int(geno_matrix(in_row_a, snp)) /= 9) nm_a = nm_a + 1
        if (int(geno_matrix(in_row_b, snp)) /= 9) nm_b = nm_b + 1
      end do
      if (nm_b > nm_a) then
        out_keep_row = in_row_b
        out_drop_row = in_row_a
      else
        out_keep_row = in_row_a
        out_drop_row = in_row_b
      end if
      out_resolved = .false.
    end if
  end subroutine resolve_identical_pair_id

  subroutine merge_genotype_rows(in_keep_row, in_drop_row, in_n_snps, io_conflict_count)
    integer, intent(in) :: in_keep_row, in_drop_row, in_n_snps
    integer, intent(inout) :: io_conflict_count
    integer :: s, gk, gd

    do s = 1, in_n_snps
      gk = int(geno_matrix(in_keep_row, s))
      gd = int(geno_matrix(in_drop_row, s))
      if (gk == 9 .and. gd /= 9) then
        geno_matrix(in_keep_row, s) = int(gd, ki1)
      else if (gk /= 9 .and. gd /= 9 .and. gk /= gd) then
        geno_matrix(in_keep_row, s) = 9_ki1
        io_conflict_count = io_conflict_count + 1
      end if
    end do
  end subroutine merge_genotype_rows

  subroutine rewrite_parent_references(in_old_id, in_new_id, in_n_animals, io_updates)
    character(len=*), intent(in) :: in_old_id, in_new_id
    integer, intent(in) :: in_n_animals
    integer, intent(inout) :: io_updates

    integer :: i
    logical :: found, changed
    type(PEDInfo) :: ped_rec

    if (len_trim(in_old_id) == 0) return
    if (len_trim(in_new_id) == 0) return
    if (trim(in_old_id) == trim(in_new_id)) return

    do i = 1, in_n_animals
      if (len_trim(geno_animal_ids(i)) == 0) cycle
      found = pht_search(PHT, trim(geno_animal_ids(i)), ped_rec)
      if (.not. found) cycle

      changed = .false.
      if (trim(ped_rec%SIRE) == trim(in_old_id)) then
        ped_rec%SIRE = trim(in_new_id)
        changed = .true.
      end if
      if (trim(ped_rec%DAM) == trim(in_old_id)) then
        ped_rec%DAM = trim(in_new_id)
        changed = .true.
      end if

      if (changed) then
        call pht_insert_with_key(PHT, ped_rec, trim(ped_animal_key_field))
        io_updates = io_updates + 1
      end if
    end do
  end subroutine rewrite_parent_references

  subroutine finalize_deferred_identical_pairs(in_n_animals, in_n_snps, in_match_threshold)
    integer, intent(in) :: in_n_animals, in_n_snps
    real(r8), intent(in) :: in_match_threshold

    integer, parameter :: MIN_COMPARE_SNPS = 500
    integer :: i, j, n_comp
    integer :: keep_row, drop_row
    integer :: n_final_merged, n_ref_updates, n_conflict_to_missing
    integer :: n_final_ped_resolved, n_final_fallback
    integer :: n_sample, nm_keep, nm_drop
    real(r8) :: match_rate
    logical :: resolved_by_ped
    logical, allocatable :: row_dropped(:)
    integer, allocatable :: sample_idx(:)
    real(r8), allocatable :: sample_maf(:)
    character(len=LEN_STR) :: keep_id, drop_id

    if (.not. associated(id_merge_seek_flag)) then
      print *, '  INFO: No deferred identical-ID state. Step 9 skipped.'
      return
    end if
    if (count(id_merge_seek_flag(1:in_n_animals) == 1) < 2) then
      print *, '  INFO: No deferred identical-ID pairs to finalize.'
      return
    end if
    if (.not. associated(map_records)) then
      print *, '  WARNING: Step 9 skipped (no MAP records in ID merge context).'
      return
    end if
    if (.not. associated(snp_maf)) then
      print *, '  WARNING: Step 9 skipped (no SNP MAF profile in ID merge context).'
      return
    end if
    if (.not. associated(snp_call_rate)) then
      print *, '  WARNING: Step 9 skipped (no SNP call-rate profile in ID merge context).'
      return
    end if

    if (associated(snp_hwe_pvalue)) then
      call select_parent_seek_snps(in_n_snps, map_records, snp_call_rate, snp_maf, sample_idx, sample_maf, &
                                   n_sample, hwe_threshold, snp_hwe_pvalue, 'ID_MERGE_FINAL')
    else
      call select_parent_seek_snps(in_n_snps, map_records, snp_call_rate, snp_maf, sample_idx, sample_maf, &
                                   n_sample, hwe_threshold, in_context='ID_MERGE_FINAL')
    end if
    if (n_sample <= 0) then
      print *, '  WARNING: Step 9 skipped (no SNP panel).'
      if (allocated(sample_idx)) deallocate(sample_idx)
      if (allocated(sample_maf)) deallocate(sample_maf)
      return
    end if

    allocate(row_dropped(in_n_animals))
    row_dropped = .false.
    n_final_merged = 0
    n_ref_updates = 0
    n_conflict_to_missing = 0
    n_final_ped_resolved = 0
    n_final_fallback = 0

    do i = 1, in_n_animals - 1
      if (id_merge_seek_flag(i) /= 1) cycle
      if (len_trim(geno_animal_ids(i)) == 0) cycle
      if (row_dropped(i)) cycle

      do j = i + 1, in_n_animals
        if (id_merge_seek_flag(j) /= 1) cycle
        if (len_trim(geno_animal_ids(j)) == 0) cycle
        if (row_dropped(j)) cycle
        if (trim(geno_animal_ids(i)) == trim(geno_animal_ids(j))) cycle

        call compute_pair_similarity_sampled_mod(geno_matrix, i, j, n_sample, sample_idx, n_comp, match_rate)

        if (n_comp < MIN_COMPARE_SNPS) cycle
        if (match_rate < in_match_threshold) cycle

        call resolve_identical_pair_id(i, j, in_n_animals, n_sample, sample_idx, &
                                       keep_row, drop_row, resolved_by_ped)
        if (resolved_by_ped) then
          n_final_ped_resolved = n_final_ped_resolved + 1
        else
          n_final_fallback = n_final_fallback + 1
          call count_non_missing_sampled_mod(geno_matrix, i, n_sample, sample_idx, nm_keep)
          call count_non_missing_sampled_mod(geno_matrix, j, n_sample, sample_idx, nm_drop)
          if (nm_drop > nm_keep) then
            keep_row = j
            drop_row = i
          else
            keep_row = i
            drop_row = j
          end if
        end if

        keep_id = trim(geno_animal_ids(keep_row))
        drop_id = trim(geno_animal_ids(drop_row))

        call merge_genotype_rows(keep_row, drop_row, in_n_snps, n_conflict_to_missing)
        call rewrite_parent_references(drop_id, keep_id, in_n_animals, n_ref_updates)

        geno_matrix(drop_row, :) = 9_ki1
        geno_animal_ids(drop_row) = ''
        row_dropped(drop_row) = .true.
        id_merge_seek_flag(drop_row) = 0
        id_merge_seek_flag(keep_row) = 1

        n_final_merged = n_final_merged + 1
        print '(a,1x,a,1x,a,1x,f7.4,1x,a,1x,a)', &
          '  Final ID collapse:', trim(drop_id), '->', match_rate, 'kept as', trim(keep_id)
      end do
    end do

    print '(a, i8)', '  Step 9 collapsed records:           ', n_final_merged
    print '(a, i8)', '    - PED-resolved merges:           ', n_final_ped_resolved
    print '(a, i8)', '    - fallback merges:               ', n_final_fallback
    print '(a, i8)', '  Step 9 parent refs rewired:         ', n_ref_updates
    print '(a, i8)', '  Step 9 SNP conflicts -> missing:    ', n_conflict_to_missing

    if (allocated(row_dropped)) deallocate(row_dropped)
    if (allocated(sample_idx)) deallocate(sample_idx)
    if (allocated(sample_maf)) deallocate(sample_maf)
  end subroutine finalize_deferred_identical_pairs

end module M_QCIDMerge
