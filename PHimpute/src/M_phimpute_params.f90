module M_phimpute_params
  use M_kinds
  use M_Variables, only: FileInfo, init_File, LEN_STR
  implicit none

  type :: phimpute_params_t
    type(FileInfo) :: pedfile
    type(FileInfo) :: mapfile
    type(FileInfo) :: cand_genofile
    type(FileInfo) :: ref_genofile
    logical :: has_reference = .false.
    integer :: omp_threads = 32
    integer :: window_sizes(3) = (/1000, 200, 50/)
    integer :: min_match_snps = 12
    integer :: max_iter_phasing = 3
    logical :: enable_pedigree = .true.
    logical :: enable_family_phasing = .true.
    logical :: enable_haplotype = .true.
    logical :: enable_population = .true.
    integer :: ld_window = 2
    integer :: ld_window_max = 5
    logical :: stop_on_mismatch = .false.
    integer :: max_autosome = 0
    logical :: hap_tune = .true.
    integer :: hap_window_candidates(7) = (/24, 28, 32, 0, 0, 0, 0/)
    integer :: hap_shift_candidates(7) = (/8, 10, 12, 0, 0, 0, 0/)
    integer :: hap_window_count = 3
    integer :: hap_shift_count = 3
    logical :: qc_enable = .false.
    real :: qc_sample_miss_max = 0.10
    real :: qc_snp_miss_max = 0.02
    real :: qc_min_maf = 0.01
    logical :: qc_remove_ambiguous = .false.
    character(len=LEN_STR) :: qc_log_prefix = 'qc'
    character(len=LEN_STR) :: output_genofile = 'phimpute_imputed.geno'
  end type phimpute_params_t

contains

  subroutine init_phimpute_params(params)
    type(phimpute_params_t), intent(out) :: params

    call init_File(params%pedfile)
    call init_File(params%mapfile)
    call init_File(params%cand_genofile)
    call init_File(params%ref_genofile)

    params%has_reference = .false.
    params%omp_threads = 32
    params%window_sizes = (/1000, 200, 50/)
    params%min_match_snps = 12
    params%max_iter_phasing = 3
    params%enable_pedigree = .true.
    params%enable_family_phasing = .true.
    params%enable_haplotype = .true.
    params%enable_population = .true.
    params%ld_window = 2
    params%ld_window_max = 5
    params%stop_on_mismatch = .false.
    params%max_autosome = 0
    params%hap_tune = .true.
    params%hap_window_candidates = (/24, 28, 32, 0, 0, 0, 0/)
    params%hap_shift_candidates = (/8, 10, 12, 0, 0, 0, 0/)
    params%hap_window_count = 3
    params%hap_shift_count = 3
    params%qc_enable = .false.
    params%qc_sample_miss_max = 0.10
    params%qc_snp_miss_max = 0.02
    params%qc_min_maf = 0.01
    params%qc_remove_ambiguous = .false.
    params%qc_log_prefix = 'qc'
    params%output_genofile = 'phimpute_imputed.geno'
  end subroutine init_phimpute_params

end module M_phimpute_params
