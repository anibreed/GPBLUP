module M_QCStepRunner
  use M_Kinds
  use M_Variables
  implicit none

  private
  public :: run_mendelian_step
  public :: run_parent_seek_step

  abstract interface
    subroutine mendel_proc_iface(in_n_animals, in_n_snps, in_snp_error_threshold, in_error_threshold, &
                                 in_sex_qc_status, out_action)
      import :: r8
      integer, intent(in) :: in_n_animals
      integer, intent(in) :: in_n_snps
      real, intent(in) :: in_snp_error_threshold
      real, intent(in) :: in_error_threshold
      integer, intent(in), optional :: in_sex_qc_status(:)
      integer, intent(out) :: out_action(in_n_animals)
    end subroutine mendel_proc_iface
  end interface

  abstract interface
    subroutine parent_seek_proc_iface(in_n_animals, in_n_snps, in_action, in_candidate_threshold, &
                                      in_popqc, in_sex_qc_status, in_corrected_sex)
      import :: r8, POPQC_matrices
      integer, intent(in) :: in_n_animals
      integer, intent(in) :: in_n_snps
      integer, intent(in) :: in_action(:)
      real(r8), intent(in) :: in_candidate_threshold
      type(POPQC_matrices), intent(in) :: in_popqc
      integer, intent(in), optional :: in_sex_qc_status(:)
      integer, intent(in), optional :: in_corrected_sex(:)
    end subroutine parent_seek_proc_iface
  end interface

contains

  subroutine run_mendelian_step(in_n_animals, in_n_snps, in_snp_error_threshold, in_error_threshold, &
                                out_action, in_mendel_proc, in_sex_qc_status)
    integer, intent(in) :: in_n_animals
    integer, intent(in) :: in_n_snps
    real, intent(in) :: in_snp_error_threshold
    real, intent(in) :: in_error_threshold
    integer, intent(out) :: out_action(:)
    procedure(mendel_proc_iface) :: in_mendel_proc
    integer, intent(in), optional :: in_sex_qc_status(:)

    if (present(in_sex_qc_status)) then
      call in_mendel_proc(in_n_animals, in_n_snps, in_snp_error_threshold, in_error_threshold, &
                          in_sex_qc_status=in_sex_qc_status, out_action=out_action)
    else
      call in_mendel_proc(in_n_animals, in_n_snps, in_snp_error_threshold, in_error_threshold, &
                          out_action=out_action)
    end if
  end subroutine run_mendelian_step

  subroutine run_parent_seek_step(in_n_animals, in_n_snps, in_action, in_candidate_threshold, in_popqc, &
                                  in_parent_seek_proc, in_sex_qc_status, in_corrected_sex)
    integer, intent(in) :: in_n_animals
    integer, intent(in) :: in_n_snps
    integer, intent(in) :: in_action(:)
    real(r8), intent(in) :: in_candidate_threshold
    type(POPQC_matrices), intent(in) :: in_popqc
    procedure(parent_seek_proc_iface) :: in_parent_seek_proc
    integer, intent(in), optional :: in_sex_qc_status(:)
    integer, intent(in), optional :: in_corrected_sex(:)

    if (present(in_sex_qc_status) .and. present(in_corrected_sex)) then
      call in_parent_seek_proc(in_n_animals, in_n_snps, in_action, in_candidate_threshold, in_popqc, &
                               in_sex_qc_status=in_sex_qc_status, in_corrected_sex=in_corrected_sex)
    else if (present(in_sex_qc_status)) then
      call in_parent_seek_proc(in_n_animals, in_n_snps, in_action, in_candidate_threshold, in_popqc, &
                               in_sex_qc_status=in_sex_qc_status)
    else
      call in_parent_seek_proc(in_n_animals, in_n_snps, in_action, in_candidate_threshold, in_popqc)
    end if
  end subroutine run_parent_seek_step

end module M_QCStepRunner
