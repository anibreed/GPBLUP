module M_QCParentSeek
  use M_Kinds
  use M_Variables
  use M_ped, only: CPED
  use M_StrEdit, only: bdate_to_jdn
  use M_Hashtable
  use H_PedRead
  use M_QCMendel, only: build_geno_row_map, find_geno_row_map, is_duo_mendel_error, is_trio_mendel_error
  use M_QCSNPSelect, only: select_parent_seek_snps_shared => select_parent_seek_snps
  implicit none
  private
  public :: set_parent_seek_context
  public :: identify_correct_parentage_mod

  type(PEDHashTable), pointer, save :: PHT => null()
  type(SNPInfo), pointer, save :: map_records(:) => null()
  integer(kind=ki1), pointer, save :: geno_matrix(:,:) => null()
  character(len=LEN_STR), pointer, save :: geno_animal_ids(:) => null()
  character(len=LEN_STR), save :: ped_animal_key_field = "ANIMAL_ID"
  real(r8), pointer, save :: snp_call_rate(:) => null()
  real(r8), pointer, save :: snp_maf(:) => null()
  real(r8), pointer, save :: snp_hwe_pvalue(:) => null()
  real(r8), save :: hwe_threshold = 0.0_r8

contains

  subroutine set_parent_seek_context(in_pht, in_map_records, in_geno_matrix, in_geno_animal_ids, &
                                     in_ped_key_field, in_snp_call_rate, in_snp_maf, &
                                     in_snp_hwe_pvalue, in_hwe_threshold)
    type(PEDHashTable), target, intent(inout) :: in_pht
    type(SNPInfo), target, intent(in) :: in_map_records(:)
    integer(kind=ki1), target, intent(inout) :: in_geno_matrix(:,:)
    character(len=LEN_STR), target, intent(inout) :: in_geno_animal_ids(:)
    character(len=*), intent(in) :: in_ped_key_field
    real(r8), target, intent(in), optional :: in_snp_call_rate(:), in_snp_maf(:), &
                         in_snp_hwe_pvalue(:)
    real(r8), intent(in), optional :: in_hwe_threshold
    PHT => in_pht
    map_records => in_map_records
    geno_matrix => in_geno_matrix
    geno_animal_ids => in_geno_animal_ids
    ped_animal_key_field = trim(in_ped_key_field)
    if (present(in_snp_call_rate)) snp_call_rate => in_snp_call_rate
    if (present(in_snp_maf)) snp_maf => in_snp_maf
    if (present(in_snp_hwe_pvalue)) snp_hwe_pvalue => in_snp_hwe_pvalue
    if (present(in_hwe_threshold)) hwe_threshold = in_hwe_threshold
  end subroutine set_parent_seek_context


  subroutine identify_correct_parentage_mod(in_n_animals, in_n_snps, in_action, in_candidate_threshold, &
      in_popqc, in_sex_qc_status, in_corrected_sex)
    !-------------------------------------------------------------------
    ! For each animal flagged as action=1 (CLEAR_PARENT):
    !
    ! Candidate evaluation follows a strict logical order:
    !
    !   [1] Sex filter
    !       SIRE candidates: MALE or UNKNOWN sex only
    !       DAM  candidates: FEMALE or UNKNOWN sex only
    !
    !   [2] Age filter (BDate)
    !       Candidate must be born >= MIN_PARENT_AGE_DAYS before the child.
    !       Skipped when either BDate is unknown (0).
    !
    !   [3] Mendelian compatibility (Mendel score)
    !       1 - duo_error_rate over all autosome SNPs.
    !       Candidates below in_candidate_threshold are rejected.
    !       When multiple G≈0.5 candidates exist, lowest error count wins
    !       (full siblings share G≈0.5 but produce Mendel conflicts).
    !
    !   [4] Genomic relationship (G filter)
    !       Up to MAX_G_SNPS_PER_CHR high-quality autosomal SNPs per
    !       chromosome (Chr 1-18; PAR=19, X=20, Y=21, unknown=23 excluded).
    !       G(child,cand) computed via VanRaden normalisation.
    !       Only candidates with |G - 0.5| <= G_TOLERANCE pass.
    !
    !   [5] Outcome classification (per animal):
    !       CASE A – ID_SWAP
    !         No candidate matches when child is treated as itself,
    !         but ANOTHER existing animal ID matches both a sire AND dam.
    !         => Report swapped ID; update PHT with genetic parents.
    !       CASE B – PARENT_UPDATE
    !         Normal case: true parents found in the population.
    !         => Update PHT SIRE/DAM with found candidates.
    !       CASE C – MULTI_CANDIDATE
    !         Multiple G≈0.5 candidates after sex/age/G filters.
    !         => Select the one with minimum Mendel error (siblings have errors).
    !       CASE D – NOT_FOUND
    !         No G≈0.5 candidate exists in the population.
    !         => SIRE/DAM remain '0' (unknown) for AlphaImpute2 LD resolution.
    !
    ! Animals with action=2 (DELETE) are listed but not updated.
    !-------------------------------------------------------------------
    integer,   intent(in) :: in_n_animals
    integer,   intent(in) :: in_n_snps
    integer,   intent(inout) :: in_action(:)
    real(r8),  intent(in) :: in_candidate_threshold
    type(POPQC_matrices), intent(in) :: in_popqc
    integer,   intent(in), optional :: in_sex_qc_status(:)  ! 0=pass, 1=fail, 2=ambiguous, -1=no_data
    integer,   intent(in), optional :: in_corrected_sex(:)  ! post-sex-switch corrected sex per geno row

    ! ── Thresholds ────────────────────────────────────────────────────
    integer,  parameter :: MIN_VALID_SNPS      = 100    ! min overlapping SNPs for Mendel
    integer,  parameter :: MIN_PARENT_AGE_DAYS = 365    ! min age gap child-parent (days)
    real(r8), parameter :: G_TARGET            = 0.5_r8
    real(r8), parameter :: POST_RECHECK_FACTOR = 5.0_r8

    ! ── Local variables ───────────────────────────────────────────────
    integer :: i, c, j, k
    integer :: n_seek, n_assigned_bc, n_no_match, n_id_swap
    integer :: n_parent_update, n_multi_resolved, n_multi_resolved_deleted, n_not_found
    integer :: n_id_swap_direct, n_id_swap_recursive, swap_depth_found
    integer :: n_id_swap_changed, n_id_swap_nochange
    integer :: n_parent_update_changed, n_parent_update_nochange
    integer :: n_multi_resolved_changed, n_multi_resolved_nochange
    integer :: n_not_found_changed, n_not_found_nochange
    integer :: g_child, g_cand, n_valid, n_error
    integer :: child_jdn, cand_jdn, cand_sex_val
    real(r8) :: compat, g_rel, g_tolerance
    ! Best-candidate tracking (lowest error count among G≈0.5 candidates)
    integer  :: best_m_errors, best_f_errors
    real(r8) :: best_m_compat, best_f_compat
    real(r8) :: best_m_g,      best_f_g
    integer  :: best_m_cand,   best_f_cand
    integer  :: n_m_cands,     n_f_cands     ! count of G-passing candidates per slot
    type(CPED) :: ped_target, ped_cand, ped_new
    logical  :: found, age_ok, g_ok
    logical  :: upd_sire, upd_dam
    character(len=LEN_STR) :: sire_id, dam_id
    character(len=16)      :: outcome_tag
    character(len=4)       :: outcome_code
    character(len=LEN_STR) :: orig_id, final_id, swap_target_id
    character(len=LEN_STR) :: animal_disp
    character(len=LEN_STR) :: norm_id
    character(len=LEN_STR) :: sire_old_disp, sire_new_disp, dam_old_disp, dam_new_disp
    character(len=10)      :: change_class
    character(len=3)       :: change_code
    character(len=MAX_RECL) :: parent_change_disp
    character(len=MAX_RECL) :: change_field_disp
    logical                :: id_changed, sire_changed, dam_changed, any_changed
    logical                :: src_mendel, src_sex, post_eval, post_pass
    character(len=3)       :: source_code
    character(len=4)       :: mendel_code
    integer                :: src_idx, out_idx, mdl_idx
    integer                :: sire_row_post, dam_row_post, parent_row_post
    integer                :: n_valid_post, n_err_post, gp_child, gp_sire, gp_dam
    integer                :: os, oo, om
    integer                :: summary_src_out_mdl(3,4,3)
    real(r8)               :: post_err_rate, mendel_thr, post_keep_thr
    character(len=16)      :: src_label
    character(len=20)      :: out_label
    character(len=22)      :: mdl_label
    character(len=LEN_STR), allocatable :: map_keys(:)
    integer,                allocatable :: map_rows(:)
    logical,                allocatable :: map_used(:)
    integer                :: n_log_rows, l, ord_out, ord_chg, sec_count
    character(len=MAX_RECL), allocatable :: log_lines(:)
    integer,                 allocatable :: log_row(:)
    character(len=2),       allocatable :: log_out(:)
    character(len=3),       allocatable :: log_chg(:)
    character(len=LEN_STR), allocatable :: log_animal(:)
    character(len=3),       allocatable :: log_src(:)
    character(len=4),       allocatable :: log_mdl(:)
    character(len=LEN_STR), allocatable :: log_sire_old(:)
    character(len=LEN_STR), allocatable :: log_dam_old(:)
    real(r8) :: num_g, denom_g, p_k, xi_k, xc_k

    ! Parent-seek SNP panel selected by dedicated criteria
    integer,  allocatable :: seek_snp_idx(:)
    real(r8), allocatable :: seek_snp_maf(:)
    integer               :: n_seek_snps

    ! ID-swap scan state (recursive + bidirectional identity check)
    integer  :: swap_m_cand, swap_f_cand
    integer  :: swap_curr, swap_depth, swap_best_q, swap_back_best
    integer  :: q, r, nv_q, ne_q, gq, gi
    real(r8) :: swap_m_compat, swap_f_compat
    real(r8) :: compat_q, compat_q_adj, swap_best_compat, swap_back_compat
    
    ! Top-2 margin tracking: best vs second-best compatibility score
    real(r8) :: best_m_compat_2nd, best_f_compat_2nd, best_m_margin, best_f_margin
    
    ! BDATE/SEX validation for ID swap
    real(r8) :: m_top_compat, f_top_compat
    logical  :: swap_bdate_ok, swap_sex_ok, q_bdate_ok, q_sex_ok
    logical  :: is_id_swap, swap_found, fq
    logical, allocatable :: swap_visited(:)
    type(CPED) :: ped_q

    ! Pre-initialise allocatable arrays so gfortran array-descriptor bounds are
    ! always valid even before select_parent_seek_snps is called.
    allocate(seek_snp_idx(0))
    allocate(seek_snp_maf(0))
    n_seek_snps = 0

    if (.not. associated(geno_matrix) .or. .not. associated(geno_animal_ids)) then
      print *, "  WARNING: GENO matrix/ID array not ready. Step 7 skipped."
      return
    end if

    ! Count animals needing parent seek: Mendelian error (action=1) OR sex mismatch (status=1)
    n_seek = 0
    do i = 1, in_n_animals
      if (in_action(i) == 1) then
        n_seek = n_seek + 1
      else if (present(in_sex_qc_status)) then
        if (i <= size(in_sex_qc_status)) then
          if (in_sex_qc_status(i) == 1) n_seek = n_seek + 1
        end if
      end if
    end do
    n_assigned_bc = 0
    n_no_match = 0
    n_id_swap  = 0
    n_parent_update = 0
    n_multi_resolved = 0
    n_not_found = 0
    n_id_swap_direct = 0
    n_id_swap_recursive = 0
    n_id_swap_changed = 0
    n_id_swap_nochange = 0
    n_parent_update_changed = 0
    n_parent_update_nochange = 0
    n_multi_resolved_changed = 0
    n_multi_resolved_nochange = 0
    n_not_found_changed = 0
    n_not_found_nochange = 0
    summary_src_out_mdl = 0
    n_multi_resolved_deleted = 0
    g_tolerance = real(POPQC%parent_seek_g_tol, r8)
    if (g_tolerance <= 0.0_r8) g_tolerance = 0.15_r8
    mendel_thr = real(POPQC%mendel_viol_threshold, r8)
    post_keep_thr = POST_RECHECK_FACTOR * mendel_thr

    print '(a, i8)',   "  Animals seeking parent (Mendel/SexFail): ", n_seek
    print '(a, i8)',   "  Animals marked for deletion  (=2): ", count(in_action(1:in_n_animals) == 2)
    print '(a, f8.4)', "  Candidate compatibility threshold: ", in_candidate_threshold
    print '(a, f8.4, a, f8.4)', "  G relationship cutoff:          |G - ", G_TARGET, "| <= ", g_tolerance

    if (count(in_action(1:in_n_animals) == 2) > 0) then
      print *, ""
      print *, "  Animals marked for deletion:"
      do i = 1, in_n_animals
        if (in_action(i) == 2) then
          if (len_trim(geno_animal_ids(i)) > 0) then
            print '(4x, a32)', adjustl(geno_animal_ids(i))
          end if
        end if
      end do
    end if

    if (n_seek == 0) then
      print *, "  No animals require parent-seeking."
      return
    end if

    ! Select MAP-based SNP panel: autosomal 10% bp-uniform sampling.
    call select_parent_seek_snps_local(in_n_snps, seek_snp_idx, seek_snp_maf, n_seek_snps, 'PARENT_SEEK')
    print '(a, i6)', "  Parent-seek SNPs selected: ", n_seek_snps
    if (n_seek_snps < MIN_VALID_SNPS) then
      print *, "  WARNING: Too few valid SNPs for parent-seek. Step 7 cannot assign parents."
      if (allocated(seek_snp_idx)) deallocate(seek_snp_idx)
      if (allocated(seek_snp_maf)) deallocate(seek_snp_maf)
      return
    end if

    allocate(log_lines(n_seek))
    allocate(log_row(n_seek))
    allocate(log_out(n_seek))
    allocate(log_chg(n_seek))
    allocate(log_animal(n_seek))
    allocate(log_src(n_seek))
    allocate(log_mdl(n_seek))
    allocate(log_sire_old(n_seek))
    allocate(log_dam_old(n_seek))
    call build_geno_row_map(in_n_animals, geno_animal_ids, map_keys, map_rows, map_used)
    n_log_rows = 0

    ! ──────────────────────────────────────────────────────────────────
    print *, ""
    print *, "  Parent seek results (grouped):"
    print *, "  Row  Animal ID         Src  Mdl   Parent changes"
    print *, "  " // repeat('-', 140)

    do i = 1, in_n_animals
      ! Process animals with Mendelian parent error (action=1) OR sex check mismatch (status=1)
      if (in_action(i) /= 1) then
        if (.not. present(in_sex_qc_status)) cycle
        if (i > size(in_sex_qc_status)) cycle
        if (in_sex_qc_status(i) /= 1) cycle
      end if

      best_m_compat = 0.0_r8;  best_f_compat = 0.0_r8
      best_m_compat_2nd = 0.0_r8; best_f_compat_2nd = 0.0_r8
      best_m_margin = 0.0_r8; best_f_margin = 0.0_r8
      m_top_compat = 0.0_r8; f_top_compat = 0.0_r8
      best_m_g      = -9.9_r8; best_f_g      = -9.9_r8
      best_m_errors = huge(0); best_f_errors = huge(0)
      best_m_cand   = 0;       best_f_cand   = 0
      n_m_cands     = 0;       n_f_cands     = 0
      orig_id       = trim(geno_animal_ids(i))
      final_id      = trim(geno_animal_ids(i))
      swap_target_id = '0'
      swap_bdate_ok = .false.
      swap_sex_ok = .false.
      src_mendel    = (in_action(i) == 1)
      src_sex       = .false.
      if (present(in_sex_qc_status)) then
        if (i <= size(in_sex_qc_status)) src_sex = (in_sex_qc_status(i) == 1)
      end if
      if (src_mendel .and. src_sex) then
        source_code = 'BTH'
        src_idx = 3
      else if (src_mendel) then
        source_code = 'MEN'
        src_idx = 1
      else
        source_code = 'SEX'
        src_idx = 2
      end if

      ! PED lookup for child
      found = pht_search(PHT, trim(geno_animal_ids(i)), ped_target)
      if (.not. found) cycle
      child_jdn = bdate_to_jdn(ped_target%BDate)

      ! ── [PRE] ID Swap Pre-detection (선탐지) ──────────────────────────────────
      !   Scan for high-confidence ID swaps BEFORE normal parent-seeking.
      !   If i's genotype matches another q's recorded SIRE/DAM pair with high
      !   Mendel compatibility, assume swap and proceed with q's parents directly.
      is_id_swap  = .false.
      swap_m_cand = 0;  swap_f_cand = 0
      swap_m_compat = 0.0_r8;  swap_f_compat = 0.0_r8

      allocate(swap_visited(in_n_animals))
      swap_visited = .false.
      swap_curr = i
      swap_visited(swap_curr) = .true.
      swap_found = .false.

      do swap_depth = 1, in_popqc%id_swap_max_depth
        swap_best_q = 0
        swap_best_compat = 0.0_r8

        do q = 1, in_n_animals
          if (q == swap_curr) cycle
          if (in_action(q) == 2) cycle
          if (len_trim(geno_animal_ids(q)) == 0) cycle

          nv_q = 0
          ne_q = 0
          do k = 1, n_seek_snps
            j = seek_snp_idx(k)
            gi = int(geno_matrix(swap_curr, j));  if (gi == 9) cycle
            gq = int(geno_matrix(q, j));          if (gq == 9) cycle
            nv_q = nv_q + 1
            if (gi /= gq) ne_q = ne_q + 1
          end do
          if (nv_q < MIN_VALID_SNPS) cycle
          compat_q = 1.0_r8 - real(ne_q, r8) / real(nv_q, r8)
          if (compat_q < in_popqc%id_swap_compat_min) cycle

          call evaluate_swap_profile(q, child_jdn, q_bdate_ok, q_sex_ok)
          compat_q_adj = compat_q
          if (q_bdate_ok) compat_q_adj = compat_q_adj + 0.01_r8
          if (q_sex_ok)   compat_q_adj = compat_q_adj + 0.01_r8

          if (compat_q_adj > swap_best_compat) then
            swap_best_compat = compat_q_adj
            swap_best_q = q
            swap_bdate_ok = q_bdate_ok
            swap_sex_ok = q_sex_ok
          end if
        end do

        if (swap_best_q == 0) exit

        ! Bidirectional check: best_q must point back to swap_curr.
        swap_back_best = 0
        swap_back_compat = 0.0_r8
        do r = 1, in_n_animals
          if (r == swap_best_q) cycle
          if (in_action(r) == 2) cycle
          if (len_trim(geno_animal_ids(r)) == 0) cycle

          nv_q = 0
          ne_q = 0
          do k = 1, n_seek_snps
            j = seek_snp_idx(k)
            gi = int(geno_matrix(swap_best_q, j));  if (gi == 9) cycle
            gq = int(geno_matrix(r, j));            if (gq == 9) cycle
            nv_q = nv_q + 1
            if (gi /= gq) ne_q = ne_q + 1
          end do
          if (nv_q < MIN_VALID_SNPS) cycle
          compat_q = 1.0_r8 - real(ne_q, r8) / real(nv_q, r8)
          if (compat_q < in_popqc%id_swap_compat_min) cycle

          if (compat_q > swap_back_compat) then
            swap_back_compat = compat_q
            swap_back_best = r
          end if
        end do
        if (swap_back_best /= swap_curr) exit

        fq = pht_search(PHT, trim(geno_animal_ids(swap_best_q)), ped_q)
        if (fq) then
          if (len_trim(ped_q%SIRE) > 0 .and. trim(ped_q%SIRE) /= '0' .and. &
              len_trim(ped_q%DAM)  > 0 .and. trim(ped_q%DAM)  /= '0') then
            is_id_swap = .true.
            swap_m_cand = swap_best_q
            swap_m_compat = swap_best_compat
            swap_depth_found = swap_depth
            swap_found = .true.
            exit
          end if
        end if

        if (swap_visited(swap_best_q)) exit
        swap_visited(swap_best_q) = .true.
        swap_curr = swap_best_q
      end do

      deallocate(swap_visited)

      ! ── Normal parent-seeking only if ID swap NOT detected (선탐지 미해결)
      if (.not. is_id_swap) then

      do c = 1, in_n_animals
        if (c == i) cycle
        if (len_trim(geno_animal_ids(c)) == 0) cycle

        if (present(in_sex_qc_status)) then
          if (c <= size(in_sex_qc_status)) then
            if (in_sex_qc_status(c) /= 0) cycle
          end if
        end if

        ! ── [2] Age filter ──────────────────────────────────────────
        found = pht_search(PHT, trim(geno_animal_ids(c)), ped_cand)
        age_ok = .true.
        if (found) then
          cand_jdn = bdate_to_jdn(ped_cand%BDate)
          if (child_jdn > 0 .and. cand_jdn > 0) &
            age_ok = (child_jdn - cand_jdn) >= MIN_PARENT_AGE_DAYS
        else
          cand_jdn = 0
        end if
        if (.not. age_ok) cycle

        ! ── [3] Mendelian compatibility ──────────────────────────────
        n_valid = 0;  n_error = 0
        do k = 1, n_seek_snps
          j = seek_snp_idx(k)
          g_child = int(geno_matrix(i, j));  if (g_child == 9) cycle
          g_cand  = int(geno_matrix(c, j));  if (g_cand  == 9) cycle
          n_valid = n_valid + 1
          if (is_duo_mendel_error(g_child, g_cand)) n_error = n_error + 1
        end do
        if (n_valid < MIN_VALID_SNPS) cycle
        compat = 1.0_r8 - real(n_error, r8) / real(n_valid, r8)
        if (compat < in_candidate_threshold) cycle

        ! Genomic relationship filter on selected SNP panel
        g_ok  = .false.
        g_rel = -9.9_r8
        num_g = 0.0_r8;  denom_g = 0.0_r8
        do k = 1, n_seek_snps
          j    = seek_snp_idx(k)
          p_k  = seek_snp_maf(k)
          xi_k = real(int(geno_matrix(i, j)), r8)
          xc_k = real(int(geno_matrix(c, j)), r8)
          if (xi_k > 8.0_r8 .or. xc_k > 8.0_r8) cycle
          num_g   = num_g   + (xi_k - 2.0_r8*p_k) * (xc_k - 2.0_r8*p_k)
          denom_g = denom_g + 2.0_r8 * p_k * (1.0_r8 - p_k)
        end do
        if (denom_g > 0.0_r8) then
          g_rel = num_g / denom_g
          g_ok = (abs(g_rel - G_TARGET) <= g_tolerance)
        end if
        if (.not. g_ok) cycle

        ! ── [1] Sex routing + CASE C: keep minimum-error among multi-candidates
        !   When G≈0.5 candidates are multiple (e.g., full siblings also have
        !   G≈0.5 but produce Mendel conflicts), prefer lowest n_error.
        ! Use corrected sex (post-sex-switch) if available, otherwise PHT sex
        if (present(in_corrected_sex) .and. c <= size(in_corrected_sex)) then
          cand_sex_val = in_corrected_sex(c)
        else if (found) then
          cand_sex_val = ped_cand%SEX
        else
          cand_sex_val = UNKNOWN
        end if

        if (found) then
          select case (cand_sex_val)
          case (MALE)
            n_m_cands = n_m_cands + 1
            if (compat > m_top_compat) then
              best_m_compat_2nd = m_top_compat
              m_top_compat = compat
            else if (compat > best_m_compat_2nd) then
              best_m_compat_2nd = compat
            end if
            if (n_error < best_m_errors .or. &
                (n_error == best_m_errors .and. compat > best_m_compat)) then
              best_m_errors = n_error;  best_m_compat = compat
              best_m_cand   = c;        best_m_g      = g_rel
            end if
          case (FEMALE)
            n_f_cands = n_f_cands + 1
            if (compat > f_top_compat) then
              best_f_compat_2nd = f_top_compat
              f_top_compat = compat
            else if (compat > best_f_compat_2nd) then
              best_f_compat_2nd = compat
            end if
            if (n_error < best_f_errors .or. &
                (n_error == best_f_errors .and. compat > best_f_compat)) then
              best_f_errors = n_error;  best_f_compat = compat
              best_f_cand   = c;        best_f_g      = g_rel
            end if
          case default  ! UNKNOWN sex – compete for both slots
            n_m_cands = n_m_cands + 1
            n_f_cands = n_f_cands + 1
            if (compat > m_top_compat) then
              best_m_compat_2nd = m_top_compat
              m_top_compat = compat
            else if (compat > best_m_compat_2nd) then
              best_m_compat_2nd = compat
            end if
            if (compat > f_top_compat) then
              best_f_compat_2nd = f_top_compat
              f_top_compat = compat
            else if (compat > best_f_compat_2nd) then
              best_f_compat_2nd = compat
            end if
            if (n_error < best_m_errors .or. &
                (n_error == best_m_errors .and. compat > best_m_compat)) then
              best_m_errors = n_error;  best_m_compat = compat
              best_m_cand   = c;        best_m_g      = g_rel
            end if
            if (n_error < best_f_errors .or. &
                (n_error == best_f_errors .and. compat > best_f_compat)) then
              best_f_errors = n_error;  best_f_compat = compat
              best_f_cand   = c;        best_f_g      = g_rel
            end if
          end select
        else  ! not in PED: UNKNOWN sex
          n_m_cands = n_m_cands + 1
          n_f_cands = n_f_cands + 1
          if (compat > m_top_compat) then
            best_m_compat_2nd = m_top_compat
            m_top_compat = compat
          else if (compat > best_m_compat_2nd) then
            best_m_compat_2nd = compat
          end if
          if (compat > f_top_compat) then
            best_f_compat_2nd = f_top_compat
            f_top_compat = compat
          else if (compat > best_f_compat_2nd) then
            best_f_compat_2nd = compat
          end if
          if (n_error < best_m_errors .or. &
              (n_error == best_m_errors .and. compat > best_m_compat)) then
            best_m_errors = n_error;  best_m_compat = compat
            best_m_cand   = c;        best_m_g      = g_rel
          end if
          if (n_error < best_f_errors .or. &
              (n_error == best_f_errors .and. compat > best_f_compat)) then
            best_f_errors = n_error;  best_f_compat = compat
            best_f_cand   = c;        best_f_g      = g_rel
          end if
        end if
      end do  ! candidate loop

      if (n_m_cands > 1 .and. m_top_compat > 0.0_r8) best_m_margin = m_top_compat - best_m_compat_2nd
      if (n_f_cands > 1 .and. f_top_compat > 0.0_r8) best_f_margin = f_top_compat - best_f_compat_2nd

      end if  ! if (.not. is_id_swap) - end of normal parent-seek block

      ! ── [5] Outcome classification ───────────────────────────────
      !
      ! CASE A – ID_SWAP detection
      !   ID swap may be detected via prescanning (선탐지) or when normal
      !   parent-seeking finds no candidates (best_m_cand==0 && best_f_cand==0).
      !   Confirm swap by checking if genotype has high Mendel compatibility
      !   with another individual's recorded parents.
      !   'q' (q≠i) has compat(i,q) >= ID_SWAP_COMPAT_MIN for BOTH its
      !   recorded sire AND dam, the sample is likely a swap with 'q'.
      ! 
      ! If prescanning (선탐지) already found ID swap, is_id_swap=true and
      ! swap_m_cand will be set. (normal parent-seeking loop was skipped)
      ! 
      ! If prescanning did not find swap (is_id_swap=false) but parent-seeking
      ! also failed (best_m/f_cand=0), attempt fallback ID swap detection.

      if (.not. is_id_swap .and. best_m_cand == 0 .and. best_f_cand == 0) then
        ! Fallback: No parent found by normal seeking AND no swap by prescanning
        ! -> perform a second ID-swap scan as fallback.
        !
        ! 1) For current row, find the best identity mate (>= ID_SWAP_COMPAT_MIN).
        ! 2) Require bidirectional consistency: mate's best identity mate must point back.
        ! 3) If mate has valid parents, use that row as swap target.
        ! 4) Otherwise recurse (follow identity chain) up to ID_SWAP_MAX_DEPTH.
        allocate(swap_visited(in_n_animals))
        swap_visited = .false.
        swap_curr = i
        swap_visited(swap_curr) = .true.
        swap_found = .false.

        do swap_depth = 1, in_popqc%id_swap_max_depth
          swap_best_q = 0
          swap_best_compat = 0.0_r8

          do q = 1, in_n_animals
            if (q == swap_curr) cycle
            if (in_action(q) == 2) cycle
            if (len_trim(geno_animal_ids(q)) == 0) cycle

            nv_q = 0
            ne_q = 0
            do k = 1, n_seek_snps
              j = seek_snp_idx(k)
              gi = int(geno_matrix(swap_curr, j));  if (gi == 9) cycle
              gq = int(geno_matrix(q, j));          if (gq == 9) cycle
              nv_q = nv_q + 1
              if (gi /= gq) ne_q = ne_q + 1
            end do
            if (nv_q < MIN_VALID_SNPS) cycle
            compat_q = 1.0_r8 - real(ne_q, r8) / real(nv_q, r8)
            if (compat_q < in_popqc%id_swap_compat_min) cycle

            call evaluate_swap_profile(q, child_jdn, q_bdate_ok, q_sex_ok)
            compat_q_adj = compat_q
            if (q_bdate_ok) compat_q_adj = compat_q_adj + 0.01_r8
            if (q_sex_ok)   compat_q_adj = compat_q_adj + 0.01_r8

            if (compat_q_adj > swap_best_compat) then
              swap_best_compat = compat_q_adj
              swap_best_q = q
              swap_bdate_ok = q_bdate_ok
              swap_sex_ok = q_sex_ok
            end if
          end do

          if (swap_best_q == 0) exit

          ! Bidirectional check: best_q must point back to swap_curr.
          swap_back_best = 0
          swap_back_compat = 0.0_r8
          do r = 1, in_n_animals
            if (r == swap_best_q) cycle
            if (in_action(r) == 2) cycle
            if (len_trim(geno_animal_ids(r)) == 0) cycle

            nv_q = 0
            ne_q = 0
            do k = 1, n_seek_snps
              j = seek_snp_idx(k)
              gi = int(geno_matrix(swap_best_q, j));  if (gi == 9) cycle
              gq = int(geno_matrix(r, j));            if (gq == 9) cycle
              nv_q = nv_q + 1
              if (gi /= gq) ne_q = ne_q + 1
            end do
            if (nv_q < MIN_VALID_SNPS) cycle
            compat_q = 1.0_r8 - real(ne_q, r8) / real(nv_q, r8)
            if (compat_q < in_popqc%id_swap_compat_min) cycle

            if (compat_q > swap_back_compat) then
              swap_back_compat = compat_q
              swap_back_best = r
            end if
          end do
          if (swap_back_best /= swap_curr) exit

          fq = pht_search(PHT, trim(geno_animal_ids(swap_best_q)), ped_q)
          if (fq) then
            if (len_trim(ped_q%SIRE) > 0 .and. trim(ped_q%SIRE) /= '0' .and. &
                len_trim(ped_q%DAM)  > 0 .and. trim(ped_q%DAM)  /= '0') then
              is_id_swap = .true.
              swap_m_cand = swap_best_q
              swap_m_compat = swap_best_compat
              swap_depth_found = swap_depth
              swap_found = .true.
              exit
            end if
          end if

          if (swap_visited(swap_best_q)) exit
          swap_visited(swap_best_q) = .true.
          swap_curr = swap_best_q
        end do

        deallocate(swap_visited)
      end if

      ! ── Build PHT update ─────────────────────────────────────────
      ped_new  = ped_target
      upd_sire = .false.;  upd_dam = .false.
      sire_id  = '0';  dam_id = '0'   ! default: unknown

      if (is_id_swap) then
        ! CASE A: ID_SWAP – take the parents of the swapped individual
        outcome_tag = 'ID_SWAP'
        n_id_swap = n_id_swap + 1
        if (swap_m_cand > 0) swap_target_id = trim(geno_animal_ids(swap_m_cand))
        block
          type(CPED) :: ped_swap
          logical :: fs
          fs = pht_search(PHT, trim(geno_animal_ids(swap_m_cand)), ped_swap)
          if (fs) then
            sire_id = trim(ped_swap%SIRE)
            dam_id  = trim(ped_swap%DAM)
            ped_new%SIRE = sire_id
            ped_new%DAM  = dam_id
            upd_sire = (len_trim(sire_id) > 0 .and. trim(sire_id) /= '0')
            upd_dam  = (len_trim(dam_id)  > 0 .and. trim(dam_id)  /= '0')
          end if
        end block
        if (swap_depth_found <= 1) then
          n_id_swap_direct = n_id_swap_direct + 1
        else
          n_id_swap_recursive = n_id_swap_recursive + 1
        end if

      else if (best_m_cand > 0 .or. best_f_cand > 0) then
        ! CASE B/C: PARENT_UPDATE (single or multi-candidate resolved)
        if (n_m_cands > 1 .or. n_f_cands > 1) then
          outcome_tag = 'MULTI_RESOLVED'
          n_multi_resolved = n_multi_resolved + 1
        else
          outcome_tag = 'PARENT_UPDATE'
          n_parent_update = n_parent_update + 1
        end if

        if (best_m_cand > 0) then
          sire_id      = trim(geno_animal_ids(best_m_cand))
          ped_new%SIRE = sire_id
          upd_sire     = .true.
        end if
        if (best_f_cand > 0) then
          dam_id      = trim(geno_animal_ids(best_f_cand))
          ped_new%DAM = dam_id
          upd_dam     = .true.
        end if

      else
        ! CASE D: NOT_FOUND – leave SIRE/DAM = '0' for AlphaImpute2
        outcome_tag = 'NOT_FOUND'
        n_not_found = n_not_found + 1
        ped_new%SIRE = '0'
        ped_new%DAM  = '0'
        n_no_match   = n_no_match + 1
        sire_id = '0';  dam_id = '0'
      end if

      ! ── GENO ID normalization (PED immutable) ────────────────────
      norm_id = ''
      if (trim(outcome_tag) == 'ID_SWAP' .and. swap_m_cand > 0) then
        block
          type(CPED) :: ped_swap
          logical :: fs
          fs = pht_search(PHT, trim(geno_animal_ids(swap_m_cand)), ped_swap)
          if (fs .and. len_trim(ped_swap%ID) > 0) then
            norm_id = trim(ped_swap%ID)
          else
            norm_id = trim(geno_animal_ids(swap_m_cand))
          end if
        end block
      else
        if (len_trim(ped_target%ID) > 0) norm_id = trim(ped_target%ID)
      end if
      if (len_trim(norm_id) > 0) geno_animal_ids(i) = trim(norm_id)

      ! ── Multi-resolved handling option (update vs delete) ─────────
      if (trim(outcome_tag) == 'MULTI_RESOLVED' .and. in_popqc%multi_resolved_action == 1) then
        in_action(i) = 2
        geno_matrix(i, :) = 9_ki1
        geno_animal_ids(i) = ''
        n_multi_resolved_deleted = n_multi_resolved_deleted + 1
        cycle
      end if

      ! ── Geno-data corrections and pedigree updates ────────────────
      select case (trim(outcome_tag))
      case ('ID_SWAP')
        ! Correct animal ID in geno data; normalized to PED ID
        final_id = trim(geno_animal_ids(i))
      case ('PARENT_UPDATE', 'MULTI_RESOLVED')
        ! True parents identified; keep PED hash unchanged (input PED immutable).
        n_assigned_bc = n_assigned_bc + 1
      case ('NOT_FOUND')
        ! Keep original ID/genotype to maximize retention.
        ! Parent IDs are set to unknown so downstream imputation can still use the animal.
        final_id = trim(geno_animal_ids(i))
      end select

      if (len_trim(orig_id) == 0) orig_id = '-'
      if (len_trim(final_id) == 0) final_id = '-'
      animal_disp = trim(final_id)

      sire_old_disp = trim(ped_target%SIRE)
      dam_old_disp  = trim(ped_target%DAM)
      sire_new_disp = trim(ped_new%SIRE)
      dam_new_disp  = trim(ped_new%DAM)
      if (len_trim(sire_old_disp) == 0 .or. trim(sire_old_disp) == '0') sire_old_disp = '-'
      if (len_trim(dam_old_disp)  == 0 .or. trim(dam_old_disp)  == '0') dam_old_disp  = '-'
      if (len_trim(sire_new_disp) == 0 .or. trim(sire_new_disp) == '0') sire_new_disp = '-'
      if (len_trim(dam_new_disp)  == 0 .or. trim(dam_new_disp)  == '0') dam_new_disp  = '-'

      id_changed = (trim(orig_id) /= trim(final_id))
      sire_changed = (trim(sire_old_disp) /= trim(sire_new_disp))
      dam_changed = (trim(dam_old_disp) /= trim(dam_new_disp))
      any_changed = id_changed .or. sire_changed .or. dam_changed
      if (any_changed) then
        change_class = 'CHANGED'
        change_code = 'CHG'
      else
        change_class = 'NOCHANGE'
        change_code = 'NCH'
      end if

      select case (trim(outcome_tag))
      case ('PARENT_UPDATE')
        outcome_code = 'PU'
        out_idx = 1
      case ('MULTI_RESOLVED')
        outcome_code = 'MR'
        out_idx = 2
      case ('ID_SWAP')
        outcome_code = 'IS'
        out_idx = 3
      case ('NOT_FOUND')
        outcome_code = 'NF'
        out_idx = 4
      case default
        outcome_code = 'NA'
        out_idx = 4
      end select

      post_eval = .false.
      post_pass = .false.
      mendel_code = 'NA'
      mdl_idx = 3
      if (src_mendel) then
        mendel_code = 'NA'
        if (mendel_thr > 0.0_r8) then
          sire_row_post = 0
          dam_row_post = 0
          if (len_trim(ped_new%SIRE) > 0 .and. trim(ped_new%SIRE) /= '0') then
            sire_row_post = find_geno_row_map(trim(ped_new%SIRE), map_keys, map_rows, map_used)
          end if
          if (len_trim(ped_new%DAM) > 0 .and. trim(ped_new%DAM) /= '0') then
            dam_row_post = find_geno_row_map(trim(ped_new%DAM), map_keys, map_rows, map_used)
          end if
          if (sire_row_post > 0 .or. dam_row_post > 0) then
            parent_row_post = sire_row_post
            if (parent_row_post == 0) parent_row_post = dam_row_post
            n_valid_post = 0
            n_err_post = 0
            do j = 1, in_n_snps
              gp_child = int(geno_matrix(i, j)); if (gp_child == 9) cycle
              if (sire_row_post > 0 .and. dam_row_post > 0) then
                gp_sire = int(geno_matrix(sire_row_post, j))
                gp_dam  = int(geno_matrix(dam_row_post, j))
                if (gp_sire == 9 .or. gp_dam == 9) cycle
                n_valid_post = n_valid_post + 1
                if (is_trio_mendel_error(gp_child, gp_sire, gp_dam)) n_err_post = n_err_post + 1
              else
                gp_sire = int(geno_matrix(parent_row_post, j)); if (gp_sire == 9) cycle
                n_valid_post = n_valid_post + 1
                if (is_duo_mendel_error(gp_child, gp_sire)) n_err_post = n_err_post + 1
              end if
            end do
            if (n_valid_post >= MIN_VALID_SNPS) then
              post_eval = .true.
              post_err_rate = real(n_err_post, r8) / real(n_valid_post, r8)
              post_pass = (post_err_rate <= post_keep_thr)
            end if
          end if
        end if
        if (post_eval) then
          if (post_pass) then
            mendel_code = 'PASS'
            mdl_idx = 1
          else
            mendel_code = 'FAIL'
            mdl_idx = 2
          end if
        end if
      end if
      summary_src_out_mdl(src_idx, out_idx, mdl_idx) = summary_src_out_mdl(src_idx, out_idx, mdl_idx) + 1

      select case (trim(outcome_tag))
      case ('ID_SWAP')
        if (any_changed) then
          n_id_swap_changed = n_id_swap_changed + 1
        else
          n_id_swap_nochange = n_id_swap_nochange + 1
        end if
      case ('PARENT_UPDATE')
        if (any_changed) then
          n_parent_update_changed = n_parent_update_changed + 1
        else
          n_parent_update_nochange = n_parent_update_nochange + 1
        end if
      case ('MULTI_RESOLVED')
        if (any_changed) then
          n_multi_resolved_changed = n_multi_resolved_changed + 1
        else
          n_multi_resolved_nochange = n_multi_resolved_nochange + 1
        end if
      case ('NOT_FOUND')
        if (any_changed) then
          n_not_found_changed = n_not_found_changed + 1
        else
          n_not_found_nochange = n_not_found_nochange + 1
        end if
      end select

      parent_change_disp = ''
      if (id_changed) then
        parent_change_disp = trim(parent_change_disp) // 'ID[' // trim(orig_id) // '->' // trim(final_id) // ']'
      end if
      if (sire_changed) then
        change_field_disp = 'SIRE[' // trim(sire_old_disp) // '->' // trim(sire_new_disp) // ']'
        if (len_trim(parent_change_disp) > 0) then
          parent_change_disp = trim(parent_change_disp) // ' ' // trim(change_field_disp)
        else
          parent_change_disp = trim(change_field_disp)
        end if
      end if
      if (dam_changed) then
        change_field_disp = 'DAM[' // trim(dam_old_disp) // '->' // trim(dam_new_disp) // ']'
        if (len_trim(parent_change_disp) > 0) then
          parent_change_disp = trim(parent_change_disp) // ' ' // trim(change_field_disp)
        else
          parent_change_disp = trim(change_field_disp)
        end if
      end if
      if (trim(outcome_tag) == 'MULTI_RESOLVED') then
        if (m_top_compat > 0.0_r8) then
          write(change_field_disp, '(a, f7.4)') 'TOP1_SIRE_COMPAT=', m_top_compat
          if (len_trim(parent_change_disp) > 0) then
            parent_change_disp = trim(parent_change_disp) // ' ' // trim(change_field_disp)
          else
            parent_change_disp = trim(change_field_disp)
          end if
        end if
        if (best_m_margin > 0.0_r8) then
          write(change_field_disp, '(a, f7.4)') 'TOP2_SIRE_MARGIN=', best_m_margin
          if (len_trim(parent_change_disp) > 0) then
            parent_change_disp = trim(parent_change_disp) // ' ' // trim(change_field_disp)
          else
            parent_change_disp = trim(change_field_disp)
          end if
        end if
        if (f_top_compat > 0.0_r8) then
          write(change_field_disp, '(a, f7.4)') 'TOP1_DAM_COMPAT=', f_top_compat
          if (len_trim(parent_change_disp) > 0) then
            parent_change_disp = trim(parent_change_disp) // ' ' // trim(change_field_disp)
          else
            parent_change_disp = trim(change_field_disp)
          end if
        end if
        if (best_f_margin > 0.0_r8) then
          write(change_field_disp, '(a, f7.4)') 'TOP2_DAM_MARGIN=', best_f_margin
          if (len_trim(parent_change_disp) > 0) then
            parent_change_disp = trim(parent_change_disp) // ' ' // trim(change_field_disp)
          else
            parent_change_disp = trim(change_field_disp)
          end if
        end if
      end if
      if (trim(outcome_tag) == 'ID_SWAP') then
        write(change_field_disp, '(a, l1, a, l1, a)') 'SWAPCHK[BDATE=', swap_bdate_ok, ' SEX=', swap_sex_ok, ']'
        if (len_trim(parent_change_disp) > 0) then
          parent_change_disp = trim(parent_change_disp) // ' ' // trim(change_field_disp)
        else
          parent_change_disp = trim(change_field_disp)
        end if
      end if
      if (len_trim(parent_change_disp) == 0) parent_change_disp = '-'
      if (trim(outcome_code) == 'IS' .and. len_trim(swap_target_id) > 0 .and. trim(swap_target_id) /= '0') then
        parent_change_disp = trim(parent_change_disp) // ' T[' // trim(swap_target_id) // ']'
      end if

      n_log_rows = n_log_rows + 1
      write(log_lines(n_log_rows), '(2x, i4, 2x, a, 2x, a3, 2x, a4, 2x, a)') &
        i, trim(animal_disp), trim(source_code), trim(mendel_code), trim(parent_change_disp)
      log_row(n_log_rows) = i
      log_out(n_log_rows) = trim(outcome_code)
      log_chg(n_log_rows) = trim(change_code)
      log_animal(n_log_rows) = trim(animal_disp)
      log_src(n_log_rows) = trim(source_code)
      log_mdl(n_log_rows) = trim(mendel_code)
      log_sire_old(n_log_rows) = trim(sire_old_disp)
      log_dam_old(n_log_rows) = trim(dam_old_disp)
    end do  ! animal loop

    do ord_out = 1, 4
      do ord_chg = 1, 2
        sec_count = 0
        do l = 1, n_log_rows
          select case (ord_out)
          case (1)
            if (trim(log_out(l)) /= 'PU') cycle
          case (2)
            if (trim(log_out(l)) /= 'MR') cycle
          case (3)
            if (trim(log_out(l)) /= 'IS') cycle
          case (4)
            if (trim(log_out(l)) /= 'NF') cycle
          end select
          select case (ord_chg)
          case (1)
            if (trim(log_chg(l)) /= 'CHG') cycle
          case (2)
            if (trim(log_chg(l)) /= 'NCH') cycle
          end select
          sec_count = sec_count + 1
        end do

        if (sec_count > 0) then
          select case (ord_out)
          case (1)
            if (ord_chg == 1) print '(a)', '  [Parent update - Changed]'
            if (ord_chg == 2) print '(a)', '  [Parent update - No change]'
          case (2)
            if (ord_chg == 1) print '(a)', '  [Multi resolved - Changed]'
            if (ord_chg == 2) print '(a)', '  [Multi resolved - No change]'
          case (3)
            if (ord_chg == 1) print '(a)', '  [ID swap - Changed]'
            if (ord_chg == 2) print '(a)', '  [ID swap - No change]'
          case (4)
            if (ord_chg == 1) print '(a)', '  [Not found - Changed]'
            if (ord_chg == 2) print '(a)', '  [Not found - No change]'
          end select
          do l = 1, n_log_rows
            select case (ord_out)
            case (1)
              if (trim(log_out(l)) /= 'PU') cycle
            case (2)
              if (trim(log_out(l)) /= 'MR') cycle
            case (3)
              if (trim(log_out(l)) /= 'IS') cycle
            case (4)
              if (trim(log_out(l)) /= 'NF') cycle
            end select
            select case (ord_chg)
            case (1)
              if (trim(log_chg(l)) /= 'CHG') cycle
            case (2)
              if (trim(log_chg(l)) /= 'NCH') cycle
            end select
            if (ord_out == 4) then
              print '(2x, i4, 2x, a, 2x, a, a, a, a)', log_row(l), trim(log_animal(l)), &
                'PED[SIRE=', trim(log_sire_old(l)), ' DAM=', trim(log_dam_old(l)) // ']'
            else if (ord_chg == 2) then
              cycle
            else
              print '(a)', trim(log_lines(l))
            end if
          end do
          if (ord_chg == 2 .and. ord_out /= 4) then
            print '(a)', '    (No-change rows omitted; see CHG/NCH summary counts below)'
          end if
        end if
      end do
    end do

    if (allocated(seek_snp_idx)) deallocate(seek_snp_idx)
    if (allocated(seek_snp_maf)) deallocate(seek_snp_maf)
    if (allocated(log_lines)) deallocate(log_lines)
    if (allocated(log_row)) deallocate(log_row)
    if (allocated(log_out)) deallocate(log_out)
    if (allocated(log_chg)) deallocate(log_chg)
    if (allocated(log_animal)) deallocate(log_animal)
    if (allocated(log_src)) deallocate(log_src)
    if (allocated(log_mdl)) deallocate(log_mdl)
    if (allocated(log_sire_old)) deallocate(log_sire_old)
    if (allocated(log_dam_old)) deallocate(log_dam_old)
    if (allocated(map_keys)) deallocate(map_keys)
    if (allocated(map_rows)) deallocate(map_rows)
    if (allocated(map_used)) deallocate(map_used)

    print *, ""
    print '(a, i8)', "  Parent update total:             ", n_parent_update
    print '(a, i8, a, i8)', "    - Changed / No change:          ", n_parent_update_changed, " /", n_parent_update_nochange
    print '(a, i8)', "  Multi resolved total:            ", n_multi_resolved
    if (n_multi_resolved_deleted > 0) then
      print '(a, i8)', "    - Deleted by parameter:         ", n_multi_resolved_deleted
    end if
    print '(a, i8, a, i8)', "    - Changed / No change:          ", n_multi_resolved_changed, " /", n_multi_resolved_nochange
    print '(a, i8)', "  ID swap total:                   ", n_id_swap
    print '(a, i8, a, i8)', "    - Changed / No change:          ", n_id_swap_changed, " /", n_id_swap_nochange
    print '(a, i8)', "    - Direct match (depth=1):        ", n_id_swap_direct
    print '(a, i8)', "    - Recursive match (depth>=2):    ", n_id_swap_recursive
    print '(a, i8)', "  Not found total:                 ", n_not_found
    print '(a, i8, a, i8)', "    - Changed / No change:          ", n_not_found_changed, " /", n_not_found_nochange
    print *, ""
    print *, "  Source x Outcome x Post-Mendel status:"
    do os = 1, 3
      select case (os)
      case (1)
        src_label = 'MENDEL'
      case (2)
        src_label = 'SEX'
      case default
        src_label = 'BOTH'
      end select
      do oo = 1, 4
        select case (oo)
        case (1)
          out_label = 'Parent update'
        case (2)
          out_label = 'Multi resolved'
        case (3)
          out_label = 'ID swap'
        case default
          out_label = 'Not found'
        end select
        do om = 1, 3
          if (summary_src_out_mdl(os, oo, om) <= 0) cycle
          select case (om)
          case (1)
            mdl_label = 'PRE_FAIL->POST_PASS'
          case (2)
            mdl_label = 'PRE_FAIL->POST_FAIL'
          case default
            mdl_label = 'PRE_NA'
          end select
          print '(4x,"SRC=",a,"  OUT=",a,"  MDL=",a,"  N=",i8)', &
            trim(src_label), trim(out_label), trim(mdl_label), summary_src_out_mdl(os, oo, om)
        end do
      end do
    end do
    print '(a, i8)', "  Animals with parent assigned (B/C): ", n_assigned_bc
    print '(a, i8)', "  No candidate found          (D):   ", n_no_match
    if (n_no_match > 0) then
      print '(a)', "    (D animals are retained with original ID; PED parents set to 0/0)"
    end if

  end subroutine identify_correct_parentage_mod

  subroutine select_parent_seek_snps_local(in_n_snps, out_idx, out_maf, out_n, in_context)
    integer, intent(in) :: in_n_snps
    integer, allocatable, intent(out) :: out_idx(:)
    real(r8), allocatable, intent(out) :: out_maf(:)
    integer, intent(out) :: out_n
    character(len=*), intent(in), optional :: in_context

    if (.not. associated(map_records)) return
    if (.not. associated(snp_maf)) return
    if (.not. associated(snp_call_rate)) return

    if (associated(snp_hwe_pvalue)) then
      call select_parent_seek_snps_shared(in_n_snps, map_records, snp_call_rate, snp_maf, out_idx, out_maf, &
                                          out_n, hwe_threshold, snp_hwe_pvalue, in_context)
    else
      call select_parent_seek_snps_shared(in_n_snps, map_records, snp_call_rate, snp_maf, out_idx, out_maf, &
                                          out_n, hwe_threshold, in_context=in_context)
    end if
  end subroutine select_parent_seek_snps_local

  subroutine evaluate_swap_profile(in_row_q, in_child_jdn, out_bdate_ok, out_sex_ok)
    integer, intent(in) :: in_row_q
    integer, intent(in) :: in_child_jdn
    logical, intent(out) :: out_bdate_ok, out_sex_ok

    type(CPED) :: ped_q_local, ped_sire, ped_dam
    logical :: fq_local, fs_local, fd_local
    integer :: sire_jdn, dam_jdn

    out_bdate_ok = .true.
    out_sex_ok = .true.

    if (in_row_q < 1 .or. in_row_q > size(geno_animal_ids)) return
    if (len_trim(geno_animal_ids(in_row_q)) == 0) return

    fq_local = pht_search(PHT, trim(geno_animal_ids(in_row_q)), ped_q_local)
    if (.not. fq_local) return

    if (len_trim(ped_q_local%SIRE) > 0 .and. trim(ped_q_local%SIRE) /= '0') then
      fs_local = pht_search(PHT, trim(ped_q_local%SIRE), ped_sire)
      if (fs_local) then
        out_sex_ok = out_sex_ok .and. (ped_sire%SEX == MALE .or. ped_sire%SEX == UNKNOWN)
        if (in_child_jdn > 0) then
          sire_jdn = bdate_to_jdn(ped_sire%BDate)
          if (sire_jdn > 0) out_bdate_ok = out_bdate_ok .and. (in_child_jdn - sire_jdn >= 365)
        end if
      end if
    end if

    if (len_trim(ped_q_local%DAM) > 0 .and. trim(ped_q_local%DAM) /= '0') then
      fd_local = pht_search(PHT, trim(ped_q_local%DAM), ped_dam)
      if (fd_local) then
        out_sex_ok = out_sex_ok .and. (ped_dam%SEX == FEMALE .or. ped_dam%SEX == UNKNOWN)
        if (in_child_jdn > 0) then
          dam_jdn = bdate_to_jdn(ped_dam%BDate)
          if (dam_jdn > 0) out_bdate_ok = out_bdate_ok .and. (in_child_jdn - dam_jdn >= 365)
        end if
      end if
    end if

  end subroutine evaluate_swap_profile

end module M_QCParentSeek
