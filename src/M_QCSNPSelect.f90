module M_QCSNPSelect
  use M_Kinds
  use M_Variables
  implicit none

  private
  public :: select_parent_seek_snps

contains

  subroutine select_parent_seek_snps(in_n_snps, in_map_records, in_snp_call_rate, in_snp_maf, &
                                     out_idx, out_maf, out_n, in_hwe_threshold, in_snp_hwe_pvalue, in_context)
    integer, intent(in) :: in_n_snps
    type(SNPInfo), intent(in) :: in_map_records(:)
    real(r8), intent(in) :: in_snp_call_rate(:)
    real(r8), intent(in) :: in_snp_maf(:)
    integer, allocatable, intent(out) :: out_idx(:)
    real(r8), allocatable, intent(out) :: out_maf(:)
    integer, intent(out) :: out_n
    real(r8), intent(in) :: in_hwe_threshold
    real(r8), intent(in), optional :: in_snp_hwe_pvalue(:)
    character(len=*), intent(in), optional :: in_context

    integer,  parameter :: AUTOSOME_MIN = 1, AUTOSOME_MAX = 18
    real(r8), parameter :: MIN_MAF = 0.30_r8
    real(r8), parameter :: MIN_CALL_RATE = 0.99_r8
    real(r8), parameter :: MIN_HWE_P = 0.001_r8
    real(r8), parameter :: SAMPLE_FRACTION = 0.10_r8

    integer :: i, j, k, chr
    integer :: n_candidates, cnt
    integer :: nc, qi, quota, b
    integer :: pos_min, pos_max, bin_lo, bin_hi, best_i
    real(r8) :: best_maf, bin_width_r8, sample_fraction_local
    real(r8) :: bin_bp_avg, bin_bp_min, bin_bp_max
    integer :: bin_stat_n
    logical :: hwe_enabled
    logical, allocatable :: map_pos_used(:)
    integer, allocatable :: chr_candidates(:), chr_selected(:)
    integer, allocatable :: cand_idx(:), cand_chr(:), cand_pos(:)
    real(r8), allocatable :: cand_maf(:)
    integer, allocatable :: c_idx(:), c_pos(:)
    real(r8), allocatable :: c_maf(:)
    character(len=24) :: context_label

    out_n = 0
    if (allocated(out_idx)) deallocate(out_idx)
    if (allocated(out_maf)) deallocate(out_maf)

    context_label = 'UNKNOWN'
    if (present(in_context)) context_label = trim(in_context)

    sample_fraction_local = SAMPLE_FRACTION
    if (trim(context_label) == 'ID_MERGE' .or. trim(context_label) == 'ID_MERGE_FINAL') then
      sample_fraction_local = real(POPQC%id_merge_sample_fraction, r8)
      if (sample_fraction_local <= 0.0_r8 .or. sample_fraction_local > 1.0_r8) then
        sample_fraction_local = 0.20_r8
      end if
    end if

    if (in_n_snps <= 0) return

    allocate(cand_idx(in_n_snps), cand_chr(in_n_snps), cand_pos(in_n_snps), cand_maf(in_n_snps))
    allocate(map_pos_used(in_n_snps))
    allocate(chr_candidates(AUTOSOME_MAX), chr_selected(AUTOSOME_MAX))
    map_pos_used = .false.
    chr_candidates = 0
    chr_selected = 0
    bin_bp_avg = 0.0_r8
    bin_bp_min = huge(0.0_r8)
    bin_bp_max = 0.0_r8
    bin_stat_n = 0
    n_candidates = 0

    hwe_enabled = (in_hwe_threshold > 0.0_r8)

    do k = 1, size(in_map_records)
      chr = in_map_records(k)%Chr
      if (chr < AUTOSOME_MIN .or. chr > AUTOSOME_MAX) cycle
      j = in_map_records(k)%Array_All
      if (j < 1 .or. j > in_n_snps) cycle
      if (map_pos_used(j)) cycle
      map_pos_used(j) = .true.
      if (in_snp_call_rate(j) < MIN_CALL_RATE) cycle
      if (in_snp_maf(j) < MIN_MAF) cycle
      if (hwe_enabled .and. present(in_snp_hwe_pvalue)) then
        if (in_snp_hwe_pvalue(j) < MIN_HWE_P) cycle
      end if

      n_candidates = n_candidates + 1
      cand_idx(n_candidates) = j
      cand_chr(n_candidates) = chr
      cand_pos(n_candidates) = in_map_records(k)%Pos
      cand_maf(n_candidates) = in_snp_maf(j)
      chr_candidates(chr) = chr_candidates(chr) + 1
    end do

    if (n_candidates <= 0) then
      print '(a,1x,a,1x,a)', '  SNP sampling profile [', trim(context_label), ']'
      print '(a)',           '    Method: autosome bp-bin sampling (highest MAF per bin)'
      print '(a, f6.3)',     '    Fraction per chromosome: ', sample_fraction_local
      print '(a, f8.4)',     '    Threshold CallRate >= ', MIN_CALL_RATE
      print '(a, f8.4)',     '    Threshold MAF      >= ', MIN_MAF
      if (hwe_enabled .and. present(in_snp_hwe_pvalue)) then
        print '(a, f8.4)',   '    Threshold HWE p    >= ', MIN_HWE_P
      else
        print '(a)',         '    Threshold HWE p    : disabled'
      end if
      print '(a, i8)',       '    Candidate SNPs: ', n_candidates
      print '(a, i8)',       '    Selected SNPs : ', 0
      deallocate(cand_idx, cand_chr, cand_pos, cand_maf, map_pos_used, chr_candidates, chr_selected)
      return
    end if

    allocate(out_idx(n_candidates), out_maf(n_candidates))
    out_idx = 0
    out_maf = 0.0_r8
    cnt = 0

    do chr = AUTOSOME_MIN, AUTOSOME_MAX
      nc = 0
      do i = 1, n_candidates
        if (cand_chr(i) == chr) nc = nc + 1
      end do
      if (nc == 0) cycle

      allocate(c_idx(nc), c_pos(nc), c_maf(nc))
      nc = 0
      pos_min = huge(0)
      pos_max = -huge(0)

      do i = 1, n_candidates
        if (cand_chr(i) /= chr) cycle
        nc = nc + 1
        c_idx(nc) = cand_idx(i)
        c_pos(nc) = cand_pos(i)
        c_maf(nc) = cand_maf(i)
        if (c_pos(nc) < pos_min) pos_min = c_pos(nc)
        if (c_pos(nc) > pos_max) pos_max = c_pos(nc)
      end do

      quota = max(1, nint(real(nc, r8) * sample_fraction_local))
      if (pos_max <= pos_min) then
        best_i = 0
        best_maf = -1.0_r8
        do b = 1, nc
          if (c_maf(b) > best_maf) then
            best_maf = c_maf(b)
            best_i = b
          end if
        end do
        if (best_i > 0 .and. cnt < size(out_idx)) then
          cnt = cnt + 1
          out_idx(cnt) = c_idx(best_i)
          out_maf(cnt) = c_maf(best_i)
          chr_selected(chr) = chr_selected(chr) + 1
        end if
        deallocate(c_idx, c_pos, c_maf)
        cycle
      end if

      bin_width_r8 = real(pos_max - pos_min, r8) / real(quota, r8)
      if (bin_width_r8 > 0.0_r8) then
        bin_bp_avg = bin_bp_avg + bin_width_r8
        bin_bp_min = min(bin_bp_min, bin_width_r8)
        bin_bp_max = max(bin_bp_max, bin_width_r8)
        bin_stat_n = bin_stat_n + 1
      end if
      do qi = 1, quota
        bin_lo = pos_min + nint(real(qi - 1, r8) * bin_width_r8)
        bin_hi = pos_min + nint(real(qi, r8) * bin_width_r8) - 1
        if (qi == quota) bin_hi = pos_max

        best_i = 0
        best_maf = -1.0_r8
        do b = 1, nc
          if (c_pos(b) < bin_lo .or. c_pos(b) > bin_hi) cycle
          if (c_maf(b) > best_maf) then
            best_maf = c_maf(b)
            best_i = b
          end if
        end do

        if (best_i > 0 .and. cnt < size(out_idx)) then
          cnt = cnt + 1
          out_idx(cnt) = c_idx(best_i)
          out_maf(cnt) = c_maf(best_i)
          chr_selected(chr) = chr_selected(chr) + 1
        end if
      end do

      deallocate(c_idx, c_pos, c_maf)
    end do

    out_n = cnt
    if (out_n <= 0) then
      if (allocated(out_idx)) deallocate(out_idx)
      if (allocated(out_maf)) deallocate(out_maf)
    else if (out_n < size(out_idx)) then
      block
        integer, allocatable :: tmp_idx(:)
        real(r8), allocatable :: tmp_maf(:)
        allocate(tmp_idx(out_n), tmp_maf(out_n))
        tmp_idx = out_idx(1:out_n)
        tmp_maf = out_maf(1:out_n)
        call move_alloc(tmp_idx, out_idx)
        call move_alloc(tmp_maf, out_maf)
      end block
    end if

    print '(a,1x,a,1x,a)', '  SNP sampling profile [', trim(context_label), ']'
    print '(a)',           '    Method: autosome bp-bin sampling (highest MAF per bin)'
    print '(a, f6.3)',     '    Fraction per chromosome: ', sample_fraction_local
    print '(a, f6.3)',     '    Threshold CallRate >= ', MIN_CALL_RATE
    print '(a, f6.3)',     '    Threshold MAF      >= ', MIN_MAF
    if (hwe_enabled .and. present(in_snp_hwe_pvalue)) then
      print '(a, f8.4)',   '    Threshold HWE p    >= ', MIN_HWE_P
    else
      print '(a)',         '    Threshold HWE p    : disabled'
    end if
    print '(a, i8)',       '    Candidate SNPs: ', n_candidates
    print '(a, i8)',       '    Selected SNPs : ', out_n
    if (bin_stat_n > 0) then
      print '(a, f12.1)',  '    Window(bin) bp avg: ', bin_bp_avg / real(bin_stat_n, r8)
      print '(a, f12.1)',  '    Window(bin) bp min: ', bin_bp_min
      print '(a, f12.1)',  '    Window(bin) bp max: ', bin_bp_max
    end if
    print *, '    Chr   Candidates   Selected'
    do chr = AUTOSOME_MIN, AUTOSOME_MAX
      if (chr_candidates(chr) <= 0) cycle
      print '(4x, i2, 6x, i8, 4x, i8)', chr, chr_candidates(chr), chr_selected(chr)
    end do

    deallocate(cand_idx, cand_chr, cand_pos, cand_maf, map_pos_used, chr_candidates, chr_selected)
  end subroutine select_parent_seek_snps

end module M_QCSNPSelect
