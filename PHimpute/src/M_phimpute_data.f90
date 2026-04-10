module M_phimpute_data
  use M_kinds
  use M_Variables, only: LEN_STR
  implicit none

  type :: geno_matrix_t
    integer :: n_animals = 0
    integer :: n_snps = 0
    integer :: n_words = 0
    integer(kind=ki4), allocatable :: packed(:,:)
  end type geno_matrix_t

  type :: map_info_t
    integer :: n_snps = 0
    character(len=LEN_STR), allocatable :: snp_id(:)
    integer, allocatable :: chr(:)
    integer, allocatable :: pos(:)
    character(len=1), allocatable :: allele1(:)
    character(len=1), allocatable :: allele2(:)
    integer, allocatable :: array_all(:)
    integer, allocatable :: array_chr(:)
    integer(kind=ki8), allocatable :: keys_sorted(:)
    integer, allocatable :: array_all_sorted(:)
  end type map_info_t

  type :: hap_window_t
    integer :: start_snp = 0
    integer :: end_snp = 0
    integer :: n_keys = 0
    integer(kind=ki8), allocatable :: keys(:)
    integer, allocatable :: counts(:)
    integer, allocatable :: rep_idx(:)
    integer, allocatable :: topk_idx(:,:)
    integer, allocatable :: topk_miss(:,:)
    integer, allocatable :: topk_count(:)
  end type hap_window_t

  type :: haplotype_panel_t
    integer :: window_size = 0
    integer :: shift_size = 0
    integer :: n_windows = 0
    logical :: built = .false.
    type(hap_window_t), allocatable :: windows(:)
  end type haplotype_panel_t

  type :: phasing_stats_t
    integer :: n_trios = 0
    integer :: n_snps = 0
    integer :: n_animals = 0
    integer(kind=ki8) :: checked_total = 0
    integer(kind=ki8) :: filled_total = 0
    integer :: n_chr = 0
    integer, allocatable :: chr_ids(:)
    integer(kind=ki8), allocatable :: checked_by_chr(:)
    integer(kind=ki8), allocatable :: filled_by_chr(:)
    integer(kind=ki8), allocatable :: checked_by_snp(:)
    integer(kind=ki8), allocatable :: filled_by_snp(:)
  end type phasing_stats_t

  type :: phimpute_state_t
    type(geno_matrix_t) :: cand
    type(geno_matrix_t) :: ref
    type(map_info_t) :: map
    logical :: has_reference = .false.
    character(len=LEN_STR), allocatable :: cand_ids(:)
    character(len=LEN_STR), allocatable :: ref_ids(:)
    logical :: family_phased = .false.
    type(haplotype_panel_t) :: ref_panel
  end type phimpute_state_t

contains

  subroutine init_geno_matrix(mat, n_animals, n_snps)
    type(geno_matrix_t), intent(inout) :: mat
    integer, intent(in) :: n_animals, n_snps

    mat%n_animals = n_animals
    mat%n_snps = n_snps
    mat%n_words = (n_snps + 15) / 16
    if (allocated(mat%packed)) deallocate(mat%packed)
    allocate(mat%packed(mat%n_words, n_animals))
    mat%packed = 0
  end subroutine init_geno_matrix

  subroutine init_map_info(map, n_snps)
    type(map_info_t), intent(inout) :: map
    integer, intent(in) :: n_snps

    map%n_snps = n_snps
    if (allocated(map%snp_id)) deallocate(map%snp_id)
    if (allocated(map%chr)) deallocate(map%chr)
    if (allocated(map%pos)) deallocate(map%pos)
    if (allocated(map%allele1)) deallocate(map%allele1)
    if (allocated(map%allele2)) deallocate(map%allele2)
    if (allocated(map%array_all)) deallocate(map%array_all)
    if (allocated(map%array_chr)) deallocate(map%array_chr)
    if (allocated(map%keys_sorted)) deallocate(map%keys_sorted)
    if (allocated(map%array_all_sorted)) deallocate(map%array_all_sorted)
    allocate(map%snp_id(n_snps), map%chr(n_snps), map%pos(n_snps))
    allocate(map%allele1(n_snps), map%allele2(n_snps))
    allocate(map%array_all(n_snps), map%array_chr(n_snps))
    map%chr = 0
    map%pos = 0
    map%allele1 = ''
    map%allele2 = ''
    map%array_all = 0
    map%array_chr = 0
  end subroutine init_map_info

  subroutine copy_geno_matrix(src, dst)
    type(geno_matrix_t), intent(in) :: src
    type(geno_matrix_t), intent(inout) :: dst

    call init_geno_matrix(dst, src%n_animals, src%n_snps)
    if (dst%n_words > 0 .and. dst%n_animals > 0) then
      dst%packed = src%packed
    end if
  end subroutine copy_geno_matrix

end module M_phimpute_data
