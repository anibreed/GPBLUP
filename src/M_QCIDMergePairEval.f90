module M_QCIDMergePairEval
  use M_Kinds
  implicit none

  private
  public :: compute_pair_similarity_sampled_mod
  public :: count_non_missing_sampled_mod

contains

  subroutine compute_pair_similarity_sampled_mod(in_geno_matrix, in_row_a, in_row_b, in_n_sample, in_sample_idx, &
                                                 out_n_comp, out_match_rate)
    integer(kind=ki1), intent(in) :: in_geno_matrix(:,:)
    integer, intent(in) :: in_row_a, in_row_b, in_n_sample
    integer, intent(in) :: in_sample_idx(in_n_sample)
    integer, intent(out) :: out_n_comp
    real(r8), intent(out) :: out_match_rate

    integer :: k, s, g1, g2, n_match

    out_n_comp = 0
    out_match_rate = 0.0_r8
    n_match = 0

    do k = 1, in_n_sample
      s = in_sample_idx(k)
      g1 = int(in_geno_matrix(in_row_a, s))
      if (g1 == 9) cycle
      g2 = int(in_geno_matrix(in_row_b, s))
      if (g2 == 9) cycle
      out_n_comp = out_n_comp + 1
      if (g1 == g2) n_match = n_match + 1
    end do

    if (out_n_comp > 0) out_match_rate = real(n_match, r8) / real(out_n_comp, r8)
  end subroutine compute_pair_similarity_sampled_mod

  subroutine count_non_missing_sampled_mod(in_geno_matrix, in_row, in_n_sample, in_sample_idx, out_count)
    integer(kind=ki1), intent(in) :: in_geno_matrix(:,:)
    integer, intent(in) :: in_row, in_n_sample
    integer, intent(in) :: in_sample_idx(in_n_sample)
    integer, intent(out) :: out_count

    integer :: k, s

    out_count = 0
    do k = 1, in_n_sample
      s = in_sample_idx(k)
      if (int(in_geno_matrix(in_row, s)) /= 9) out_count = out_count + 1
    end do
  end subroutine count_non_missing_sampled_mod

end module M_QCIDMergePairEval
