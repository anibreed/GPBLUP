module M_QCSNP
  ! SNP-level Quality Control module
  !
  ! Provides:
  !   calculate_snp_statistics    - call rate, allele freq, MAF, HWE p-value
  !   calculate_hwe_test          - exact HWE test (Wigginton)
  !   report_snp_qc               - print QC report (autosome-filtered)
  !   identify_monomorphic_snps   - flag MAF=0 SNPs
  !   determine_snp_qc_status     - assign per-SNP QC status (0-5)
  !   apply_bh_fdr                - Benjamini-Hochberg FDR correction
  !
  ! Note: calculate_snp_statistics requires the genotype matrix (in_geno)
  !       passed explicitly.  report_snp_qc and determine_snp_qc_status
  !       accept an optional logical mask (in_snp_is_autosome) to restrict
  !       statistics to Chr 1-18.
  !
  use M_Kinds
  use M_Variables
  use M_StrEdit
  use M_readpar_popqc, only: M_readpar_get_popqc_thresholds
  implicit none
  private

  public :: calculate_snp_statistics
  public :: calculate_hwe_test
  public :: hwe_exact_pvalue
  public :: report_snp_qc
  public :: identify_monomorphic_snps
  public :: determine_snp_qc_status
  public :: apply_bh_fdr

contains

  ! ============================================================
  ! calculate_snp_statistics
  !   Computes per-SNP call rate, allele frequencies, MAF, and
  !   per-animal call rate.  Optionally runs HWE exact test and
  !   BH-FDR correction.
  !
  !   in_geno(in_n_animals, in_n_snps) : genotype matrix (ki1)
  !     0=AA, 1=AB, 2=BB, 9=missing
  ! ============================================================
  subroutine calculate_snp_statistics(param_file, in_n_snps, in_n_animals, in_geno, &
                                       call_rate, allele_freq, maf,            &
                                       out_animal_call_rate, hwe_pvalue, out_hwe_threshold)
    character(len=*), intent(in)  :: param_file
    integer,          intent(in)  :: in_n_snps, in_n_animals
    integer(kind=ki1),intent(in)  :: in_geno(in_n_animals, in_n_snps)
    real(r8),         intent(out) :: call_rate(in_n_snps)
    real(r8),         intent(out) :: allele_freq(in_n_snps, 2)
    real(r8),         intent(out) :: maf(in_n_snps)
    real(r8),         intent(out) :: out_animal_call_rate(in_n_animals)
    real(r8),         intent(out) :: hwe_pvalue(in_n_snps)
    real(r8),         intent(out) :: out_hwe_threshold

    integer :: i_snp, i_animal, n_valid, n_calls
    integer, allocatable :: n_valid_per_snp(:), n_0_per_snp(:), &
                            n_1_per_snp(:), n_2_per_snp(:)
    real(r8) :: freq_A, freq_B
    real(r8) :: maf_threshold, snp_cr_threshold, animal_cr_threshold
    real(r8) :: hwe_alpha
    logical  :: hwe_enabled
    real     :: popqc_snp_cr, popqc_animal_cr, popqc_maf, popqc_hwe_threshold
    integer  :: popqc_hwe_min_n, popqc_hwe_min_mac
    character(len=LEN_STR) :: popqc_hwe_method

    call M_readpar_get_popqc_thresholds(popqc_snp_cr, popqc_animal_cr, popqc_maf, &
                        popqc_hwe_method, popqc_hwe_threshold,                   &
                        popqc_hwe_min_n, popqc_hwe_min_mac)

    maf_threshold        = real(popqc_maf,            r8)
    snp_cr_threshold     = real(popqc_snp_cr,         r8)
    animal_cr_threshold  = real(popqc_animal_cr,      r8)

    hwe_alpha = real(popqc_hwe_threshold, r8)
    if (hwe_alpha > 0.0_r8) then
      out_hwe_threshold = hwe_alpha / real(max(1, in_n_snps), r8)
      hwe_enabled = .true.
    else
      out_hwe_threshold = 0.0_r8
      hwe_enabled = .false.
    end if

    print *, "  ======================================"
    print *, "  Source parameter file:     ", trim(param_file)
    print *, "  QC Thresholds Applied:"
    print *, "    MAF threshold:         ", maf_threshold
    print *, "    SNP call rate:         ", snp_cr_threshold
    print *, "    Animal call rate:      ", animal_cr_threshold
    print *, "    HWE min N:             ", popqc_hwe_min_n
    print *, "    HWE min MAC:           ", popqc_hwe_min_mac
    if (.not. hwe_enabled) then
      print *, "    HWE test:              Disabled (HWE_THRESHOLD=0.0)"
      print *, "    HWE threshold:         N/A"
    else
      write(*, '(a, es12.4, a, i0, a)') &
        "    HWE FWER (alpha):      ", hwe_alpha, "  (÷", in_n_snps, " SNPs)"
      write(*, '(a, es12.4)') &
        "    HWE per-SNP threshold: ", out_hwe_threshold
      print *, "    HWE correction:        Bonferroni"
    end if
    print *, "  ======================================"
    print *, ""

    ! ---------- per-SNP statistics ----------
    ! Optimized: Initialize SNP arrays first, then row-major loop
    print *, "  Calculating SNP-level statistics..."
    call_rate = 0.0_r8
    allele_freq = 0.0_r8
    maf = 0.0_r8
    
    ! Temporary working arrays (allocated once)
    allocate(n_valid_per_snp(in_n_snps), n_0_per_snp(in_n_snps), &
             n_1_per_snp(in_n_snps), n_2_per_snp(in_n_snps))
    n_valid_per_snp = 0
    n_0_per_snp = 0
    n_1_per_snp = 0
    n_2_per_snp = 0

    ! Row-major access (animals outer loop for better cache locality)
    do i_animal = 1, in_n_animals
      do i_snp = 1, in_n_snps
        if (in_geno(i_animal, i_snp) /= 9) then
          n_valid_per_snp(i_snp) = n_valid_per_snp(i_snp) + 1
          select case (int(in_geno(i_animal, i_snp)))
          case (0)
            n_0_per_snp(i_snp) = n_0_per_snp(i_snp) + 1
          case (1)
            n_1_per_snp(i_snp) = n_1_per_snp(i_snp) + 1
          case (2)
            n_2_per_snp(i_snp) = n_2_per_snp(i_snp) + 1
          end select
        end if
      end do
    end do

    ! Finalize SNP statistics
    do i_snp = 1, in_n_snps
      n_valid = n_valid_per_snp(i_snp)
      call_rate(i_snp) = merge(real(n_valid, r8) / real(in_n_animals, r8), &
                               0.0_r8, in_n_animals > 0)
      if (n_valid > 0) then
        freq_A = real(2 * n_0_per_snp(i_snp) + n_1_per_snp(i_snp), r8) / real(2 * n_valid, r8)
        freq_B = real(2 * n_2_per_snp(i_snp) + n_1_per_snp(i_snp), r8) / real(2 * n_valid, r8)
        allele_freq(i_snp, 1) = freq_A
        allele_freq(i_snp, 2) = freq_B
        maf(i_snp) = min(freq_A, freq_B)
      else
        allele_freq(i_snp, 1) = 0.0_r8
        allele_freq(i_snp, 2) = 0.0_r8
        maf(i_snp)            = 0.0_r8
      end if
    end do
    deallocate(n_valid_per_snp, n_0_per_snp, n_1_per_snp, n_2_per_snp)

    ! ---------- per-animal call rates ----------
    print *, "  Calculating animal-level call rates..."
    do i_animal = 1, in_n_animals
      n_calls = 0
      do i_snp = 1, in_n_snps
        if (in_geno(i_animal, i_snp) /= 9) n_calls = n_calls + 1
      end do
      out_animal_call_rate(i_animal) = merge(real(n_calls, r8) / real(in_n_snps, r8), &
                                             0.0_r8, in_n_snps > 0)
    end do

    ! ---------- HWE ----------
    hwe_pvalue = 1.0_r8
    if (hwe_enabled) then
      call calculate_hwe_test(in_n_snps, in_geno, allele_freq, hwe_pvalue, &
                              popqc_hwe_min_n, popqc_hwe_min_mac)
      if (index(to_upper(trim(popqc_hwe_method)), 'FDR') > 0 .or. &
          index(to_upper(trim(popqc_hwe_method)), 'BH')  > 0) then
        print *, "  Applying Benjamini-Hochberg FDR correction to HWE p-values..."
        block
          real(r8), allocatable :: tmp(:)
          allocate(tmp(in_n_snps))
          call apply_bh_fdr(hwe_pvalue, tmp)
          hwe_pvalue = tmp
          deallocate(tmp)
        end block
        print *, "    FDR-adjusted p-values calculated."
      end if
    else
      print *, "  HWE test skipped because HWE_THRESHOLD is 0.0"
    end if

    ! ---------- summary ----------
    print *, ""
    print *, "  ======================================"
    print *, "  SNP Statistics Summary:"
    print *, "    Total SNPs analyzed:              ", in_n_snps
    print *, "    SNPs with MAF < threshold:        ", count(maf < maf_threshold)
    print *, "    SNPs with call rate < threshold:  ", count(call_rate < snp_cr_threshold)
    print *, "    Monomorphic SNPs (MAF=0):         ", count(maf == 0.0_r8)
    if (hwe_enabled) then
      print *, "    SNPs with HWE p < threshold:      ", count(hwe_pvalue < out_hwe_threshold)
    else
      print *, "    SNPs with HWE p < threshold:      ", 0
    end if
    print *, ""
    print *, "  Animal Statistics Summary:"
    print *, "    Total animals analyzed:           ", in_n_animals
    print *, "    Animals with call rate < threshold:", &
             count(out_animal_call_rate < animal_cr_threshold)
    print *, "    Average animal call rate:         ", &
             sum(out_animal_call_rate) / real(in_n_animals, r8)
    print *, "  ======================================"
    print *, ""

  end subroutine calculate_snp_statistics

  ! ============================================================
  ! calculate_hwe_test
  !   Exact HWE test (Wigginton et al.).  Called internally by
  !   calculate_snp_statistics; also public for standalone use.
  ! ============================================================
  subroutine calculate_hwe_test(n_snps, in_geno, allele_freq, pvalues, &
                                 in_hwe_min_n, in_hwe_min_mac)
    integer,          intent(in)  :: n_snps
    integer(kind=ki1),intent(in)  :: in_geno(:,:)
    real(r8),         intent(in)  :: allele_freq(n_snps, 2)
    real(r8),         intent(out) :: pvalues(n_snps)
    integer,          intent(in)  :: in_hwe_min_n, in_hwe_min_mac

    integer :: j, i, n_animals_local
    integer :: n_AA, n_AB, n_BB, n_valid
    integer :: minor_ac, n_skipped_low_n, n_skipped_mono, n_skipped_low_mac
    integer :: hwe_min_n, hwe_min_mac

    print *, "  Calculating HWE p-values using exact test (allele-count conditioned)..."

    hwe_min_n   = max(1, in_hwe_min_n)
    hwe_min_mac = max(0, in_hwe_min_mac)
    n_animals_local   = size(in_geno, 1)
    n_skipped_low_n   = 0
    n_skipped_mono    = 0
    n_skipped_low_mac = 0

    do j = 1, n_snps
      n_AA = 0; n_AB = 0; n_BB = 0; n_valid = 0
      do i = 1, n_animals_local
        if (in_geno(i, j) /= 9) then
          n_valid = n_valid + 1
          if (in_geno(i, j) == 0) n_AA = n_AA + 1
          if (in_geno(i, j) == 1) n_AB = n_AB + 1
          if (in_geno(i, j) == 2) n_BB = n_BB + 1
        end if
      end do

      if (n_valid < hwe_min_n) then
        pvalues(j) = 1.0_r8; n_skipped_low_n = n_skipped_low_n + 1; cycle
      end if
      if (allele_freq(j, 1) <= 0.0_r8 .or. allele_freq(j, 2) <= 0.0_r8) then
        pvalues(j) = 1.0_r8; n_skipped_mono = n_skipped_mono + 1; cycle
      end if
      minor_ac = min(2 * n_AA + n_AB, 2 * n_BB + n_AB)
      if (minor_ac < hwe_min_mac) then
        pvalues(j) = 1.0_r8; n_skipped_low_mac = n_skipped_low_mac + 1; cycle
      end if

      pvalues(j) = hwe_exact_pvalue(n_AB, n_AA, n_BB)
    end do

    print *, "  HWE p-value calculation complete."
    print *, "    Skipped (N <", hwe_min_n, "): ", n_skipped_low_n
    print *, "    Skipped (monomorphic): ", n_skipped_mono
    print *, "    Skipped (minor AC <", hwe_min_mac, "): ", n_skipped_low_mac
    print *, "    SNPs with p < 0.001: ", count(pvalues < 0.001_r8)
    print *, "    SNPs with p < 0.01:  ", count(pvalues < 0.01_r8)
    print *, "    SNPs with p < 0.05:  ", count(pvalues < 0.05_r8)

  end subroutine calculate_hwe_test

  ! hwe_exact_pvalue  (Wigginton two-sided exact test)
  ! ============================================================
  function hwe_exact_pvalue(obs_hets, obs_hom1, obs_hom2) result(pval)
    integer, intent(in) :: obs_hets, obs_hom1, obs_hom2
    real(r8) :: pval

    integer :: n, rare_copies, h_obs
    integer :: mid, h, hr, hc
    real(r8), allocatable :: probs(:)
    real(r8) :: sum_prob, p_obs

    n = obs_hets + obs_hom1 + obs_hom2
    if (n <= 0) then; pval = 1.0_r8; return; end if

    rare_copies = 2 * min(obs_hom1, obs_hom2) + obs_hets
    h_obs = obs_hets
    if (h_obs < 0 .or. h_obs > rare_copies) then; pval = 0.0_r8; return; end if

    allocate(probs(0:rare_copies))
    probs = 0.0_r8

    mid = int(real(rare_copies, r8) * real(2 * n - rare_copies, r8) / real(2 * n, r8))
    if (mod(mid, 2) /= mod(rare_copies, 2)) mid = mid + 1

    probs(mid) = 1.0_r8
    sum_prob   = probs(mid)

    h = mid
    hr = (rare_copies - h) / 2
    hc = n - h - hr
    do while (h >= 2)
      probs(h - 2) = probs(h) * real(h * (h - 1), r8) / real(4 * (hr + 1) * (hc + 1), r8)
      sum_prob = sum_prob + probs(h - 2)
      h = h - 2; hr = hr + 1; hc = hc + 1
    end do

    h = mid
    hr = (rare_copies - h) / 2
    hc = n - h - hr
    do while (h <= rare_copies - 2)
      probs(h + 2) = probs(h) * real(4 * hr * hc, r8) / real((h + 2) * (h + 1), r8)
      sum_prob = sum_prob + probs(h + 2)
      h = h + 2; hr = hr - 1; hc = hc - 1
    end do

    if (sum_prob <= 0.0_r8) then
      pval = 1.0_r8; deallocate(probs); return
    end if

    probs  = probs / sum_prob
    p_obs  = probs(h_obs)
    pval   = sum(probs, mask=probs <= p_obs + 1.0e-14_r8)
    pval   = max(0.0_r8, min(1.0_r8, pval))
    deallocate(probs)
  end function hwe_exact_pvalue

  ! ============================================================
  ! report_snp_qc
  !   Prints QC summary report.
  !   in_snp_is_autosome (optional): if present, statistics are
  !   restricted to Chr 1-18 SNPs; otherwise all SNPs are used.
  ! ============================================================
  subroutine report_snp_qc(n_snps, in_n_animals, call_rate, maf, monomorphic, hwe_p, &
                            status, a_call_rate, n_pass, n_low_maf, n_low_cr,         &
                            n_maf, n_mono, n_hwe, in_snp_is_autosome)
    integer,  intent(in)           :: n_snps, in_n_animals
    integer,  intent(in)           :: n_pass, n_low_maf, n_low_cr, n_maf, n_mono, n_hwe
    real(r8), intent(in)           :: call_rate(n_snps), maf(n_snps)
    real(r8), intent(in)           :: hwe_p(n_snps), a_call_rate(in_n_animals)
    integer,  intent(in)           :: monomorphic(n_snps), status(n_snps)
    logical,  intent(in), optional :: in_snp_is_autosome(n_snps)

    real(r8) :: avg_cr, avg_maf_v, min_cr, max_cr
    real(r8) :: avg_acr, min_acr, max_acr
    integer  :: n_auto, n_total_filtered

    ! Suppress unused-argument warnings
    if (n_snps > 0) then
      if (monomorphic(1) == -999) print *, hwe_p(1)
    end if

    ! Stats restricted to autosome SNPs when mask is provided
    if (present(in_snp_is_autosome)) then
      n_auto = count(in_snp_is_autosome)
      if (n_auto > 0) then
        avg_cr   = sum(call_rate, mask=in_snp_is_autosome) / real(n_auto, r8)
        min_cr   = minval(call_rate, mask=in_snp_is_autosome)
        max_cr   = maxval(call_rate, mask=in_snp_is_autosome)
        avg_maf_v = sum(maf, mask=in_snp_is_autosome) / real(n_auto, r8)
      else
        avg_cr = 0.0_r8; min_cr = 0.0_r8; max_cr = 0.0_r8; avg_maf_v = 0.0_r8
      end if
    else
      n_auto   = n_snps
      avg_cr   = sum(call_rate) / real(n_snps, r8)
      min_cr   = minval(call_rate)
      max_cr   = maxval(call_rate)
      avg_maf_v = sum(maf) / real(n_snps, r8)
    end if

    avg_acr = sum(a_call_rate) / real(in_n_animals, r8)
    min_acr = minval(a_call_rate)
    max_acr = maxval(a_call_rate)
    n_total_filtered = n_auto - n_pass

    print *, ""
    print *, "========================================="
    print *, "SNP QUALITY CONTROL REPORT"
    print *, "========================================="
    print *, ""
    print *, "--- SNP QC Summary ---"
    print *, "  Total SNPs in GENO:            ", n_snps
    print *, "  Autosome SNPs (Chr 1-18):      ", n_auto
    print *, "  Sex/Unknown SNPs (Chr >18):    ", n_snps - n_auto
    print *, "  ** VALID SNPs (PASS QC):       ", n_pass
    print *, "  Total filtered autosome SNPs:  ", n_total_filtered
    print *, ""
    print *, "--- Filtering Details (autosome only) ---"
    print *, "  Monomorphic SNPs:              ", n_mono
    print *, "  Low call rate SNPs:            ", n_low_cr
    print *, "  Low MAF SNPs:                  ", n_low_maf
    print *, "  MAF-fail SNPs (all):           ", n_maf
    print *, "  HWE test failed SNPs:          ", n_hwe
    print *, ""
    print *, "--- SNP Statistics (autosome only) ---"
    print *, "  Average SNP call rate:         ", avg_cr
    print *, "  Min SNP call rate:             ", min_cr
    print *, "  Max SNP call rate:             ", max_cr
    print *, "  Average MAF:                   ", avg_maf_v
    print *, ""
    print *, "--- Animal Statistics ---"
    print *, "  Total animals analyzed:        ", in_n_animals
    print *, "  Average animal call rate:      ", avg_acr
    print *, "  Min animal call rate:          ", min_acr
    print *, "  Max animal call rate:          ", max_acr
    print *, ""
    print *, "--- QC Status Distribution ---"
    print *, "  Status 0 (Pass):               ", count(status == 0)
    print *, "  Status 1 (Low MAF):            ", count(status == 1)
    print *, "  Status 2 (Monomorphic):        ", count(status == 2)
    print *, "  Status 3 (HWE fail):           ", count(status == 3)
    print *, "  Status 4 (Low call rate):      ", count(status == 4)
    print *, "  Status 5 (Sex/Unk chr, skip):  ", count(status == 5)
    print *, ""
    print *, "========================================="
    print *, ""
  end subroutine report_snp_qc

  ! ============================================================
  ! identify_monomorphic_snps
  ! ============================================================
  subroutine identify_monomorphic_snps(n_snps, maf, monomorphic)
    integer,  intent(in)  :: n_snps
    real(r8), intent(in)  :: maf(n_snps)
    integer,  intent(out) :: monomorphic(n_snps)

    integer :: j, n_monomorphic

    n_monomorphic = 0
    do j = 1, n_snps
      if (maf(j) == 0.0_r8) then
        monomorphic(j) = 1
        n_monomorphic  = n_monomorphic + 1
      else
        monomorphic(j) = 0
      end if
    end do
    print *, "  Monomorphic SNP detection complete."
    print *, "    Total monomorphic SNPs: ", n_monomorphic
  end subroutine identify_monomorphic_snps

  ! ============================================================
  ! determine_snp_qc_status
  !   Status codes:
  !     0 = Pass
  !     1 = Low MAF
  !     2 = Monomorphic
  !     3 = HWE fail
  !     4 = Low call rate
  !     5 = Non-autosome (Chr >18); only set when in_snp_is_autosome provided
  !
  !   in_snp_is_autosome (optional): if present, non-autosome SNPs
  !     receive status 5 and are excluded from all other counts.
  ! ============================================================
  subroutine determine_snp_qc_status(n_snps, call_rate, maf, monomorphic, hwe_p, &
                                      in_hwe_threshold, status, n_pass,           &
                                      n_low_maf, n_low_cr, n_maf_fail, n_mono, n_hwe, &
                                      in_snp_is_autosome)
    integer,  intent(in)           :: n_snps
    real(r8), intent(in)           :: call_rate(n_snps), maf(n_snps), hwe_p(n_snps)
    real(r8), intent(in)           :: in_hwe_threshold
    integer,  intent(in)           :: monomorphic(n_snps)
    integer,  intent(out)          :: status(n_snps)
    integer,  intent(out)          :: n_pass, n_low_maf, n_low_cr, n_maf_fail, n_mono, n_hwe
    logical,  intent(in), optional :: in_snp_is_autosome(n_snps)

    integer :: j
    real(r8) :: maf_threshold, snp_cr_threshold
    logical  :: hwe_enabled
    real :: popqc_snp_cr, popqc_animal_cr, popqc_maf, popqc_hwe_threshold
    character(len=LEN_STR) :: popqc_hwe_method

    call M_readpar_get_popqc_thresholds(popqc_snp_cr, popqc_animal_cr, popqc_maf, &
                                        popqc_hwe_method, popqc_hwe_threshold)
    snp_cr_threshold = real(popqc_snp_cr, r8)
    maf_threshold    = real(popqc_maf,    r8)
    hwe_enabled      = (in_hwe_threshold > 0.0_r8)

    n_pass = 0; n_low_maf = 0; n_low_cr = 0
    n_maf_fail = 0; n_mono = 0; n_hwe = 0

    do j = 1, n_snps
      ! Non-autosome: skip QC counts
      if (present(in_snp_is_autosome)) then
        if (.not. in_snp_is_autosome(j)) then
          status(j) = 5; cycle
        end if
      end if

      if (monomorphic(j) == 1) then
        status(j) = 2; n_mono = n_mono + 1
      else if (call_rate(j) < snp_cr_threshold) then
        status(j) = 4; n_low_cr = n_low_cr + 1
      else if (maf(j) < maf_threshold) then
        status(j) = 1; n_low_maf = n_low_maf + 1; n_maf_fail = n_maf_fail + 1
      else if (hwe_enabled .and. hwe_p(j) < in_hwe_threshold) then
        status(j) = 3; n_hwe = n_hwe + 1
      else
        status(j) = 0; n_pass = n_pass + 1
      end if
    end do
  end subroutine determine_snp_qc_status

  ! ============================================================
  ! apply_bh_fdr  — Benjamini-Hochberg FDR correction
  ! ============================================================
  subroutine apply_bh_fdr(inp_p, out_p)
    real(r8), intent(in)  :: inp_p(:)
    real(r8), intent(out) :: out_p(size(inp_p))

    integer :: n, i, j
    integer, allocatable  :: order(:), rank(:)
    real(r8), allocatable :: sorted_p(:), adjusted_p(:)
    real(r8) :: cummin

    n = size(inp_p)
    if (n == 0) return

    allocate(order(n), rank(n), sorted_p(n), adjusted_p(n))

    do i = 1, n; order(i) = i; end do

    do i = 1, n - 1
      do j = i + 1, n
        if (inp_p(order(j)) < inp_p(order(i))) call swap_int(order(i), order(j))
      end do
    end do

    do i = 1, n
      sorted_p(i)   = inp_p(order(i))
      rank(order(i)) = i
    end do

    do i = 1, n
      adjusted_p(i) = min(sorted_p(i) * real(n, r8) / real(i, r8), 1.0_r8)
    end do

    cummin = adjusted_p(n)
    do i = n - 1, 1, -1
      cummin = min(cummin, adjusted_p(i))
      adjusted_p(i) = cummin
    end do

    do i = 1, n; out_p(order(i)) = adjusted_p(i); end do

    deallocate(order, rank, sorted_p, adjusted_p)
  end subroutine apply_bh_fdr

  ! ============================================================
  ! swap_int  — private integer swap helper
  ! ============================================================
  subroutine swap_int(a, b)
    integer, intent(inout) :: a, b
    integer :: tmp
    tmp = a; a = b; b = tmp
  end subroutine swap_int

end module M_QCSNP
