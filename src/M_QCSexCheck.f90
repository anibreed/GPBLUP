module M_QCSexCheck
  use M_Kinds
  use M_Variables
  use M_ped, only: CPED
  use M_Hashtable
  use H_PedRead
  implicit none

  private
  public :: detect_x_chromosome_snps
  public :: fill_ped_sex_from_hashtable
  public :: calculate_sex_concordance
  public :: report_sex_concordance

contains

  subroutine detect_x_chromosome_snps(in_map, in_n_map, out_indices, out_n)
    ! Scan map records for Chr==20 (X chromosome) and collect Array_All indices.
    type(SNPInfo), intent(in) :: in_map(:)
    integer, intent(in) :: in_n_map
    integer, allocatable, intent(out) :: out_indices(:)
    integer, intent(out) :: out_n

    integer, parameter :: X_CHR = 20
    integer :: i, cnt

    cnt = 0
    do i = 1, in_n_map
      if (in_map(i)%Chr == X_CHR) cnt = cnt + 1
    end do

    out_n = cnt
    if (cnt == 0) return

    allocate(out_indices(cnt))
    cnt = 0
    do i = 1, in_n_map
      if (in_map(i)%Chr == X_CHR) then
        cnt = cnt + 1
        out_indices(cnt) = in_map(i)%Array_All
      end if
    end do
  end subroutine detect_x_chromosome_snps

  subroutine fill_ped_sex_from_hashtable(in_pht, in_n, in_ids, out_sex)
    type(PEDHashTable), intent(in) :: in_pht
    integer, intent(in) :: in_n
    character(len=*), intent(in) :: in_ids(:)
    integer, intent(out) :: out_sex(:)

    integer :: i
    type(CPED) :: ped_rec
    logical :: found

    out_sex = UNKNOWN
    do i = 1, in_n
      if (len_trim(in_ids(i)) > 0) then
        found = pht_search(in_pht, trim(in_ids(i)), ped_rec)
        if (found) out_sex(i) = ped_rec%SEX
      end if
    end do
  end subroutine fill_ped_sex_from_hashtable

  subroutine calculate_sex_concordance(in_n_animals, in_n_x_snps, in_x_indices, in_geno_matrix, &
                                       in_ped_sex, in_female_f_threshold, out_f_values, out_qc_status)
    ! Compute X-chromosome F-statistic per animal.
    integer, intent(in) :: in_n_animals, in_n_x_snps
    integer, intent(in) :: in_x_indices(:)
    integer(kind=ki1), intent(in) :: in_geno_matrix(:,:)
    integer, intent(in) :: in_ped_sex(:)
    real(r8), intent(in) :: in_female_f_threshold
    real(r8), intent(out) :: out_f_values(:)
    integer, intent(out) :: out_qc_status(:)

    real(r8), allocatable :: x_allele_freq(:)
    integer :: i, k, col
    integer :: n_valid, n_het
    integer :: n_all, n_AA, n_AB, n_BB
    real(r8) :: p_k, q_k, sum_exp_het, f_val
    integer :: predicted_sex
    integer :: n_x_valid  ! Count of valid X SNP indices

    allocate(x_allele_freq(in_n_x_snps))
    x_allele_freq = 0.5_r8

    ! Pre-validate X SNP indices
    n_x_valid = 0
    do k = 1, in_n_x_snps
      if (in_x_indices(k) >= 1 .and. in_x_indices(k) <= size(in_geno_matrix, 2)) then
        n_x_valid = n_x_valid + 1
      end if
    end do

    ! Calculate X allele frequencies (loop 1: SNP-outer for column computation)
    do k = 1, in_n_x_snps
      col = in_x_indices(k)
      if (col < 1 .or. col > size(in_geno_matrix, 2)) cycle
      n_all = 0
      n_AA = 0
      n_AB = 0
      n_BB = 0
      do i = 1, in_n_animals
        if (in_geno_matrix(i, col) /= 9) then
          select case (in_geno_matrix(i, col))
            case (0); n_AA = n_AA + 1
            case (1); n_AB = n_AB + 1
            case (2); n_BB = n_BB + 1
          end select
          n_all = n_all + 1
        end if
      end do
      if (n_all > 0) x_allele_freq(k) = real(2 * n_AA + n_AB, r8) / real(2 * n_all, r8)
    end do

    ! Calculate F-statistic per individual (loop 2: individual-outer for row efficiency)
    do i = 1, in_n_animals
      n_valid = 0
      n_het = 0
      sum_exp_het = 0.0_r8

      do k = 1, in_n_x_snps
        col = in_x_indices(k)
        if (col < 1 .or. col > size(in_geno_matrix, 2)) cycle
        if (in_geno_matrix(i, col) == 9) cycle

        n_valid = n_valid + 1
        if (in_geno_matrix(i, col) == 1) n_het = n_het + 1
        p_k = x_allele_freq(k)
        q_k = 1.0_r8 - p_k
        sum_exp_het = sum_exp_het + 2.0_r8 * p_k * q_k
      end do

      if (n_valid >= 10 .and. sum_exp_het > 0.0_r8) then
        f_val = 1.0_r8 - real(n_het, r8) / sum_exp_het
      else
        f_val = -9.0_r8
      end if
      out_f_values(i) = f_val

      if (f_val < -8.0_r8) then
        predicted_sex = UNKNOWN
      else if (f_val > 0.8_r8) then
        predicted_sex = MALE
      else if (f_val < in_female_f_threshold) then
        predicted_sex = FEMALE
      else
        predicted_sex = UNKNOWN
      end if

      if (f_val < -8.0_r8) then
        out_qc_status(i) = -1
      else if (predicted_sex == UNKNOWN) then
        out_qc_status(i) = 2
      else if (in_ped_sex(i) /= UNKNOWN .and. in_ped_sex(i) /= predicted_sex) then
        out_qc_status(i) = 1
      else
        out_qc_status(i) = 0
      end if
    end do

    deallocate(x_allele_freq)
  end subroutine calculate_sex_concordance

  subroutine report_sex_concordance(in_n_animals, in_ids, in_ped_sex, in_f_values, in_qc_status, in_female_f_threshold)
    integer, intent(in) :: in_n_animals
    character(len=*), intent(in) :: in_ids(:)
    integer, intent(in) :: in_ped_sex(:), in_qc_status(:)
    real(r8), intent(in) :: in_f_values(:)
    real(r8), intent(in) :: in_female_f_threshold

    integer :: i, n_ok, n_mismatch, n_ambiguous, n_nodata, n_issues, predicted_sex
    character(len=8) :: ped_sex_str, pred_sex_str, status_str

    n_ok = 0
    n_mismatch = 0
    n_ambiguous = 0
    n_nodata = 0
    n_issues = 0

    print *, "  Sex concordance mismatch details (MISMATCH only):"
    print *, "  Animal                   Ped_Sex  Pred_Sex  F_value  Status"
    print *, "  ---------------------------------------------------------------"

    do i = 1, in_n_animals
      if (in_f_values(i) < -8.0_r8) then
        predicted_sex = UNKNOWN
      else if (in_f_values(i) > 0.8_r8) then
        predicted_sex = MALE
      else if (in_f_values(i) < in_female_f_threshold) then
        predicted_sex = FEMALE
      else
        predicted_sex = UNKNOWN
      end if

      if (in_ped_sex(i) == FEMALE) then
        ped_sex_str = 'FEMALE'
      else if (in_ped_sex(i) == MALE) then
        ped_sex_str = 'MALE'
      else
        ped_sex_str = 'UNKNOWN'
      end if

      if (predicted_sex == FEMALE) then
        pred_sex_str = 'FEMALE'
      else if (predicted_sex == MALE) then
        pred_sex_str = 'MALE'
      else
        pred_sex_str = 'AMBIG'
      end if

      select case (in_qc_status(i))
      case (0)
        status_str = 'OK'
        n_ok = n_ok + 1
      case (1)
        status_str = 'MISMATCH'
        n_mismatch = n_mismatch + 1
      case (2)
        status_str = 'AMBIG'
        n_ambiguous = n_ambiguous + 1
      case default
        status_str = 'NO_DATA'
        n_nodata = n_nodata + 1
      end select

      if (in_qc_status(i) == 1) then
        n_issues = n_issues + 1
        print '(2x, a24, 2x, a7, 2x, a6, 2x, f7.4, 2x, a8)', &
          trim(in_ids(i)), trim(ped_sex_str), trim(pred_sex_str), in_f_values(i), trim(status_str)
      end if
    end do

    if (n_issues == 0) print *, "  (no MISMATCH animals found)"

    print *, "  ---------------------------------------------------------------"
    print *, "  Sex concordance summary:"
    print '(a, f8.4)', "    Female F threshold: ", in_female_f_threshold
    print '(a, i6)', "    OK:          ", n_ok
    print '(a, i6)', "    Mismatch:    ", n_mismatch
    print '(a, i6)', "    Ambiguous:   ", n_ambiguous
    print '(a, i6)', "    No data:     ", n_nodata
    print '(a, i6)', "    Total:       ", n_ok + n_mismatch + n_ambiguous + n_nodata
  end subroutine report_sex_concordance

end module M_QCSexCheck
