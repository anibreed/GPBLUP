module M_phimpute_pedigree_deterministic
  use M_kinds
  use M_Variables, only: PEDFile, LEN_STR
  use M_IO_Ped, only: PEDHashTable, pht_load_ped_file, pht_free, pht_search
  use M_ped, only: CPED, str_hash_t
  use M_phimpute_params, only: phimpute_params_t
  use M_phimpute_data, only: phimpute_state_t
  use M_phimpute_bitpack, only: get_geno_value, set_geno_value, GENO_MISSING
  use M_phimpute_mpi, only: mpi_get_rank_size
  use omp_lib
  implicit none

  public :: run_pedigree_deterministic

contains

  subroutine run_pedigree_deterministic(params, state)
    type(phimpute_params_t), intent(in) :: params
    type(phimpute_state_t), intent(inout) :: state
    type(PEDHashTable) :: pht
    type(CPED) :: ped_rec
    type(str_hash_t) :: id_to_idx
    integer :: n_rec, unit_fd
    integer :: i, n_animals, n_snps
    integer, allocatable :: sire_idx(:), dam_idx(:)
    character(len=LEN_STR) :: sire_id, dam_id
    integer :: idx_out
    logical :: found
    integer, allocatable :: chr_by_idx(:)
    integer, allocatable :: chr_ids(:), chr_counts(:), chr_offsets(:), chr_snps(:)
    integer :: n_chr, c, pos, total_snps
    integer :: rank, n_ranks, c_start, c_end

    n_animals = state%cand%n_animals
    n_snps = state%cand%n_snps
    if (n_animals <= 0 .or. n_snps <= 0) return
    if (.not. allocated(state%cand_ids)) return

    PEDFile = params%pedfile

    call pht_load_ped_file(pht, n_rec, unit_fd, 'ID', quiet=.true.)

    call id_to_idx%reset()
    do i = 1, n_animals
      call id_to_idx%put(trim(state%cand_ids(i)), i)
    end do

    allocate(sire_idx(n_animals), dam_idx(n_animals))
    sire_idx = 0
    dam_idx = 0

    do i = 1, n_animals
      found = pht_search(pht, trim(state%cand_ids(i)), ped_rec)
      if (.not. found) cycle
      sire_id = trim(ped_rec%SIRE)
      dam_id = trim(ped_rec%DAM)
      if (len_trim(sire_id) > 0 .and. sire_id /= '0') then
        call id_to_idx%get(trim(sire_id), idx_out, found)
        if (found) sire_idx(i) = idx_out
      end if
      if (len_trim(dam_id) > 0 .and. dam_id /= '0') then
        call id_to_idx%get(trim(dam_id), idx_out, found)
        if (found) dam_idx(i) = idx_out
      end if
    end do

    call build_chr_index(state, chr_by_idx, chr_ids, chr_counts, chr_offsets, chr_snps, n_chr)
    call mpi_get_rank_size(rank, n_ranks)
    call split_chr_range(n_chr, rank, n_ranks, c_start, c_end)

    if (c_start <= c_end) then
      !$omp parallel do collapse(2) default(shared) private(i, c, pos) schedule(static)
      do c = c_start, c_end
        do i = 1, n_animals
          if (sire_idx(i) <= 0 .or. dam_idx(i) <= 0) cycle
          do pos = chr_offsets(c), chr_offsets(c+1) - 1
            call impute_one_snp(state, i, sire_idx(i), dam_idx(i), chr_snps(pos))
          end do
        end do
      end do
      !$omp end parallel do
    end if

    deallocate(sire_idx, dam_idx, chr_by_idx, chr_ids, chr_counts, chr_offsets, chr_snps)
    call pht_free(pht)
    call id_to_idx%reset()
  end subroutine run_pedigree_deterministic

  subroutine impute_one_snp(state, child_idx, sire_idx, dam_idx, snp_idx)
    type(phimpute_state_t), intent(inout) :: state
    integer, intent(in) :: child_idx, sire_idx, dam_idx, snp_idx
    integer(kind=ki1) :: g_child, g_sire, g_dam
    integer(kind=ki1) :: g_out
    logical :: ok

    g_child = get_geno_value(state%cand%packed(:, child_idx), snp_idx)
    if (g_child /= GENO_MISSING) return
    g_sire = get_geno_value(state%cand%packed(:, sire_idx), snp_idx)
    g_dam = get_geno_value(state%cand%packed(:, dam_idx), snp_idx)
    call mendel_deterministic(g_sire, g_dam, g_out, ok)
    if (ok) call set_geno_value(state%cand%packed(:, child_idx), snp_idx, g_out)
  end subroutine impute_one_snp

  subroutine build_chr_index(state, chr_by_idx, chr_ids, chr_counts, chr_offsets, chr_snps, n_chr)
    type(phimpute_state_t), intent(in) :: state
    integer, allocatable, intent(out) :: chr_by_idx(:)
    integer, allocatable, intent(out) :: chr_ids(:)
    integer, allocatable, intent(out) :: chr_counts(:)
    integer, allocatable, intent(out) :: chr_offsets(:)
    integer, allocatable, intent(out) :: chr_snps(:)
    integer, intent(out) :: n_chr

    integer :: n_snps, i, idx, c, pos, total_snps
    logical :: seen

    n_snps = state%cand%n_snps
    allocate(chr_by_idx(n_snps))
    chr_by_idx = 0

    do i = 1, state%map%n_snps
      idx = state%map%array_all(i)
      if (idx > 0 .and. idx <= n_snps) then
        chr_by_idx(idx) = state%map%chr(i)
      end if
    end do

    allocate(chr_ids(n_snps))
    n_chr = 0
    do i = 1, n_snps
      if (chr_by_idx(i) <= 0) cycle
      seen = .false.
      do c = 1, n_chr
        if (chr_ids(c) == chr_by_idx(i)) then
          seen = .true.
          exit
        end if
      end do
      if (.not. seen) then
        n_chr = n_chr + 1
        chr_ids(n_chr) = chr_by_idx(i)
      end if
    end do
    if (n_chr <= 0) then
      allocate(chr_counts(1), chr_offsets(2), chr_snps(1))
      chr_counts = 0
      chr_offsets = 1
      chr_snps = 0
      return
    end if

    allocate(chr_counts(n_chr))
    chr_counts = 0
    do i = 1, state%map%n_snps
      idx = state%map%array_all(i)
      if (idx <= 0 .or. idx > n_snps) cycle
      c = find_chr_index(chr_ids, n_chr, state%map%chr(i))
      if (c <= 0) cycle
      if (state%map%array_chr(i) > chr_counts(c)) chr_counts(c) = state%map%array_chr(i)
    end do

    allocate(chr_offsets(n_chr + 1))
    chr_offsets(1) = 1
    do c = 1, n_chr
      chr_offsets(c + 1) = chr_offsets(c) + chr_counts(c)
    end do
    total_snps = chr_offsets(n_chr + 1) - 1
    allocate(chr_snps(total_snps))
    chr_counts = 0

    do i = 1, state%map%n_snps
      idx = state%map%array_all(i)
      if (idx <= 0 .or. idx > n_snps) cycle
      c = find_chr_index(chr_ids, n_chr, state%map%chr(i))
      if (c <= 0) cycle
      if (state%map%array_chr(i) > 0) then
        pos = chr_offsets(c) + state%map%array_chr(i) - 1
      else
        chr_counts(c) = chr_counts(c) + 1
        pos = chr_offsets(c) + chr_counts(c) - 1
      end if
      if (pos >= chr_offsets(c) .and. pos < chr_offsets(c + 1)) then
        chr_snps(pos) = idx
      end if
    end do
  end subroutine build_chr_index

  integer function find_chr_index(chr_ids, n_chr, chr_val) result(idx)
    integer, intent(in) :: chr_ids(:)
    integer, intent(in) :: n_chr, chr_val
    integer :: i

    idx = 0
    do i = 1, n_chr
      if (chr_ids(i) == chr_val) then
        idx = i
        return
      end if
    end do
  end function find_chr_index

  subroutine split_chr_range(n_chr, rank, n_ranks, c_start, c_end)
    integer, intent(in) :: n_chr, rank, n_ranks
    integer, intent(out) :: c_start, c_end
    integer :: base, rem

    if (n_ranks <= 1) then
      c_start = 1
      c_end = n_chr
      return
    end if

    base = n_chr / n_ranks
    rem = mod(n_chr, n_ranks)
    if (rank < rem) then
      c_start = rank * (base + 1) + 1
      c_end = c_start + base
    else
      c_start = rank * base + rem + 1
      c_end = c_start + base - 1
    end if
  end subroutine split_chr_range

  subroutine mendel_deterministic(g_sire, g_dam, g_out, ok)
    integer(kind=ki1), intent(in) :: g_sire, g_dam
    integer(kind=ki1), intent(out) :: g_out
    logical, intent(out) :: ok

    ok = .false.
    g_out = GENO_MISSING

    if (g_sire == GENO_MISSING .or. g_dam == GENO_MISSING) return

    if (g_sire == 0 .and. g_dam == 0) then
      g_out = 0
      ok = .true.
    else if (g_sire == 2 .and. g_dam == 2) then
      g_out = 2
      ok = .true.
    else if ((g_sire == 0 .and. g_dam == 2) .or. (g_sire == 2 .and. g_dam == 0)) then
      g_out = 1
      ok = .true.
    end if
  end subroutine mendel_deterministic

end module M_phimpute_pedigree_deterministic
